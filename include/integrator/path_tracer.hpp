#pragma once

#include "geometry/scene.hpp"
#include "core/camera.hpp"
#include "integrator/sampling.hpp"
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <cmath>
#include <algorithm>

namespace kajiya {

/**
 * @brief Acumulador estadístico en línea por píxel mediante el algoritmo de Welford (1962).
 *        Calcula simultáneamente el valor esperado <L_o> (media muestral) y la varianza
 *        muestral insesgada s^2_N de la radiancia conforme N (SPP) se incrementa.
 */
struct PixelAccumulator {
    Color mean = Color(0.0, 0.0, 0.0);
    Color m2   = Color(0.0, 0.0, 0.0); // Suma de cuadrados de diferencias respecto a la media
    double lum_mean = 0.0;
    double lum_m2   = 0.0;
    int samples = 0;

    void reset() {
        mean = Color(0.0, 0.0, 0.0);
        m2   = Color(0.0, 0.0, 0.0);
        lum_mean = 0.0;
        lum_m2   = 0.0;
        samples  = 0;
    }

    inline void add_sample(const Color& sample_radiance) {
        samples += 1;
        double inv_n = 1.0 / static_cast<double>(samples);

        // Actualización de Welford por canal RGB
        Color delta = sample_radiance - mean;
        mean += delta * inv_n;
        Color delta2 = sample_radiance - mean;
        m2 += delta * delta2;

        // Actualización de Welford en luminancia escalar Y = 0.2126 R + 0.7152 G + 0.0722 B
        double y = sample_radiance.luminance();
        double dy = y - lum_mean;
        lum_mean += dy * inv_n;
        double dy2 = y - lum_mean;
        lum_m2 += dy * dy2;
    }

    /**
     * @brief Varianza muestral insesgada de una muestra individual X_k: Var(X) = M2 / (N - 1).
     */
    double sample_variance() const {
        if (samples < 2) return 0.0;
        return (m2.x + m2.y + m2.z) / (3.0 * static_cast<double>(samples - 1));
    }

    /**
     * @brief Varianza del estimador de Monte Carlo <L_o>_N: Var(<L_o>_N) = Var(X) / N.
     */
    double estimator_variance() const {
        if (samples < 2) return 0.0;
        return sample_variance() / static_cast<double>(samples);
    }
};

/**
 * @brief Parámetros de configuración del integrador de Monte Carlo.
 */
struct RenderConfig {
    int width = 256;
    int height = 256;
    int spp = 64;
    int max_bounces = 5;
    bool enable_russian_roulette = true;
    bool jitter_primary_rays = true;
    SamplingStrategy strategy = SamplingStrategy::CosineWeightedHemisphere;
    int num_threads = 0; // 0 = detección automática de hilos de hardware
};

/**
 * @brief Buffer de imagen con acumuladores estadísticos por píxel y métricas globales.
 */
struct RenderBuffer {
    int width;
    int height;
    std::vector<PixelAccumulator> pixels;

    RenderBuffer(int w, int h) : width(w), height(h), pixels(static_cast<size_t>(w * h)) {}

    void clear() {
        for (auto& p : pixels) p.reset();
    }

    PixelAccumulator& at(int x, int y) {
        return pixels[static_cast<size_t>(y * width + x)];
    }

    const PixelAccumulator& at(int x, int y) const {
        return pixels[static_cast<size_t>(y * width + x)];
    }

    /**
     * @brief Varianza muestral promedio Var(X) en toda la imagen.
     */
    double average_sample_variance() const {
        double sum = 0.0;
        for (const auto& p : pixels) {
            sum += p.sample_variance();
        }
        return sum / static_cast<double>(pixels.size());
    }

    /**
     * @brief Varianza promedio del estimador de Monte Carlo Var(<L_o>_N) = Var(X) / N.
     */
    double average_estimator_variance() const {
        double sum = 0.0;
        for (const auto& p : pixels) {
            sum += p.estimator_variance();
        }
        return sum / static_cast<double>(pixels.size());
    }

    /**
     * @brief Error Cuadrático Medio (MSE) frente a un buffer de referencia de alto muestreo.
     */
    double compute_mse(const RenderBuffer& reference) const {
        double sum_sq_err = 0.0;
        size_t total = pixels.size();
        for (size_t i = 0; i < total; ++i) {
            Color diff = pixels[i].mean - reference.pixels[i].mean;
            sum_sq_err += (diff.x * diff.x + diff.y * diff.y + diff.z * diff.z) / 3.0;
        }
        return sum_sq_err / static_cast<double>(total);
    }
};

/**
 * @brief Motor de Trazado de Caminos (Path Tracer) que resuelve numéricamente
 *        la Ecuación Integral de Renderizado de Kajiya:
 *
 *   L_o(x, \omega_o) = L_e(x, \omega_o) + \int_{\Omega} f_r(x, \omega_i, \omega_o) L_i(x, \omega_i) (\omega_i \cdot n) d\omega_i
 */
class PathTracer {
public:
    /**
     * @brief Evalúa un camino estocástico de transporte de luz para una muestra k:
     *
     * En cada rebote superficial en x con normal n:
     *   1. Se suma la radiancia emitida: L += throughput * L_e(x, \omega_o)
     *   2. Se genera una dirección aleatoria \omega_k en el hemisferio \Omega según p(\omega_k):
     *      - Muestreo Uniforme:
     *          p(\omega_k) = \frac{1}{2\pi}
     *          peso = \frac{f_r \cdot (\omega_k \cdot n)}{p(\omega_k)}
     *               = \frac{(\rho / \pi) \cos\theta_k}{1 / (2\pi)} = 2 \rho \cos\theta_k
     *      - Muestreo Ponderado por el Coseno:
     *          p(\omega_k) = \frac{\cos\theta_k}{\pi}
     *          peso = \frac{f_r \cdot (\omega_k \cdot n)}{p(\omega_k)}
     *               = \frac{(\rho / \pi) \cos\theta_k}{\cos\theta_k / \pi} = \rho
     *        (¡Cancelación algebraica exacta de \cos\theta_k y \pi!)
     */
    static Color trace_ray(
        const Ray& initial_ray,
        const Scene& scene,
        const RenderConfig& config,
        PCG32& rng
    ) {
        Ray current_ray = initial_ray;
        Color radiance(0.0, 0.0, 0.0);
        Color throughput(1.0, 1.0, 1.0);

        for (int bounce = 0; bounce <= config.max_bounces; ++bounce) {
            HitRecord rec;
            if (!scene.intersect(current_ray, 1e-5, 1e30, rec)) {
                // El rayo escapa al entorno (importante para el Furnace Test donde L_env = 1.0)
                radiance += throughput * scene.environment_radiance;
                break;
            }

            const Material& mat = scene.materials[static_cast<size_t>(rec.material_id)];

            // Término de emisión propia L_e(x, \omega_o) en la Ecuación de Kajiya
            if (mat.is_emissive()) {
                radiance += throughput * mat.emission;
                // Las fuentes de área ideales actúan como emisores puros para evitar varianza espuria
                break;
            }

            // Si alcanzamos el número máximo de rebotes de integral hemisférica, terminamos el camino
            if (bounce == config.max_bounces) {
                break;
            }

            // Muestreo estocástico del hemisferio \Omega orientado según la normal rec.normal
            double xi1 = rng.next_double();
            double xi2 = rng.next_double();

            HemisphereSample sample;
            if (config.strategy == SamplingStrategy::UniformHemisphere) {
                sample = sample_hemisphere_uniform(rec.normal, xi1, xi2);
                // Evaluación explícita del estimador de Monte Carlo:
                //   f_r(x, \omega_k, \omega_o) * (\omega_k \cdot n) / p(\omega_k)
                Color brdf = mat.evaluate_brdf(); // \rho / \pi
                double cos_theta = std::max(0.0, sample.cos_theta);
                Color mc_weight = (brdf * cos_theta) / sample.pdf; // = 2 * \rho * \cos\theta
                throughput *= mc_weight;
            } else {
                sample = sample_hemisphere_cosine(rec.normal, xi1, xi2);
                // Con p(\omega_k) = \cos\theta_k / \pi y f_r = \rho / \pi:
                //   \frac{(\rho / \pi) \cos\theta_k}{\cos\theta_k / \pi} = \rho
                // Cancelación algebraica exacta del término coseno (\omega_k \cdot n) = \cos\theta_k
                throughput *= mat.albedo;
            }

            // Ruleta Rusa insesgada para terminación estocástica de caminos profundos
            if (config.enable_russian_roulette && bounce >= 3) {
                double survival_prob = std::clamp(throughput.max_component(), 0.05, 0.95);
                if (rng.next_double() > survival_prob) {
                    break;
                }
                throughput /= survival_prob;
            }

            // Nuevo rayo secundario desde el punto de impacto x + \epsilon * n en dirección \omega_k
            current_ray = Ray(rec.point + rec.normal * 1e-5, sample.direction_world);
        }

        return radiance;
    }

    /**
     * @brief Ejecuta un pase de N muestras por píxel (o acumula 1 muestra adicional si incremental = true)
     *        distribuyendo bloques (tiles) de 16x16 píxeles entre todos los núcleos de CPU con std::thread.
     * @return Tiempo de cómputo en milisegundos.
     */
    static double render_pass(
        const Scene& scene,
        const Camera& camera,
        const RenderConfig& config,
        RenderBuffer& buffer,
        int sample_offset = 0,
        int samples_to_add = -1
    ) {
        int spp_count = (samples_to_add > 0) ? samples_to_add : config.spp;
        int num_threads = config.num_threads;
        if (num_threads <= 0) {
            num_threads = static_cast<int>(std::thread::hardware_concurrency());
            if (num_threads <= 0) num_threads = 4;
        }

        constexpr int TILE_SIZE = 16;
        int tiles_x = (buffer.width + TILE_SIZE - 1) / TILE_SIZE;
        int tiles_y = (buffer.height + TILE_SIZE - 1) / TILE_SIZE;
        int total_tiles = tiles_x * tiles_y;

        std::atomic<int> next_tile(0);
        auto t_start = std::chrono::high_resolution_clock::now();

        auto worker = [&]() {
            while (true) {
                int tile_idx = next_tile.fetch_add(1, std::memory_order_relaxed);
                if (tile_idx >= total_tiles) break;

                int tx = tile_idx % tiles_x;
                int ty = tile_idx / tiles_x;
                int x0 = tx * TILE_SIZE;
                int y0 = ty * TILE_SIZE;
                int x1 = std::min(x0 + TILE_SIZE, buffer.width);
                int y1 = std::min(y0 + TILE_SIZE, buffer.height);

                for (int y = y0; y < y1; ++y) {
                    for (int x = x0; x < x1; ++x) {
                        PixelAccumulator& acc = buffer.at(x, y);
                        uint64_t pixel_id = static_cast<uint64_t>(y * buffer.width + x);

                        for (int s = 0; s < spp_count; ++s) {
                            uint64_t sample_idx = static_cast<uint64_t>(sample_offset + s);
                            // Semilla determinista única por (píxel, muestra, estrategia)
                            uint64_t seed = (pixel_id * 1099511628211ULL) ^
                                            (sample_idx * 14695981039346656037ULL) ^
                                            (config.strategy == SamplingStrategy::CosineWeightedHemisphere ? 0x9E3779B97F4A7C15ULL : 0xBF58476D1CE4E5B9ULL);
                            PCG32 rng(seed, pixel_id + 1ULL);

                            double jx = config.jitter_primary_rays ? rng.next_double() : 0.5;
                            double jy = config.jitter_primary_rays ? rng.next_double() : 0.5;

                            double u = (static_cast<double>(x) + jx) / static_cast<double>(buffer.width);
                            // Invertimos y para que y=0 sea la fila superior de la imagen
                            double v = (static_cast<double>(buffer.height - 1 - y) + jy) / static_cast<double>(buffer.height);

                            Ray primary_ray = camera.get_ray(u, v);
                            Color L = trace_ray(primary_ray, scene, config, rng);
                            acc.add_sample(L);
                        }
                    }
                }
            }
        };

        std::vector<std::thread> threads;
        threads.reserve(static_cast<size_t>(num_threads));
        for (int t = 0; t < num_threads; ++t) {
            threads.emplace_back(worker);
        }
        for (auto& th : threads) {
            th.join();
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = t_end - t_start;
        return elapsed.count();
    }
};

} // namespace kajiya
