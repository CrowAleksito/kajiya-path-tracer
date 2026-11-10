#pragma once

#include "math/vec3.hpp"
#include <cstdint>
#include <cmath>

namespace kajiya {

/**
 * @brief Esquema de muestreo direccional sobre el hemisferio \Omega.
 */
enum class SamplingStrategy {
    UniformHemisphere,        // p(\omega) = 1 / (2\pi)
    CosineWeightedHemisphere  // p(\omega) = \cos\theta / \pi
};

/**
 * @brief Generador de números pseudoaleatorios PCG32 (Permuted Congruential Generator)
 *        implementado desde cero. Determinista, rápido y sin contención entre hilos.
 */
struct PCG32 {
    uint64_t state;
    uint64_t inc;

    explicit constexpr PCG32(uint64_t seed = 0x853c49e6748fea9bULL, uint64_t seq = 0xda3e39cb94b95bdbULL)
        : state(0ULL), inc((seq << 1u) | 1u) {
        next_uint32();
        state += seed;
        next_uint32();
    }

    constexpr uint32_t next_uint32() {
        uint64_t oldstate = state;
        state = oldstate * 6364136223846793005ULL + inc;
        uint32_t xorshifted = static_cast<uint32_t>(((oldstate >> 18u) ^ oldstate) >> 27u);
        uint32_t rot = static_cast<uint32_t>(oldstate >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

    /**
     * @brief Devuelve una variable aleatoria uniforme canónica \xi \in [0, 1).
     */
    constexpr double next_double() {
        return static_cast<double>(next_uint32()) * (1.0 / 4294967296.0);
    }
};

/**
 * @brief Resultado de una muestra direccional en el hemisferio \Omega.
 */
struct HemisphereSample {
    Vec3 direction_world; // Dirección unitaria \omega_k en coordenadas globales
    double cos_theta;     // \omega_k \cdot n = \cos\theta \in [0, 1]
    double pdf;           // Densidad de probabilidad p(\omega_k) respecto al ángulo sólido d\omega
};

/**
 * @brief 1. Muestreo Hemisférico Uniforme:
 *
 * Densidad respecto al ángulo sólido d\omega = \sin\theta d\theta d\phi:
 *   p(\omega) = \frac{1}{2\pi}
 *
 * En coordenadas esféricas:
 *   p(\theta, \phi) = p(\omega) \sin\theta = \frac{\sin\theta}{2\pi}
 *
 * Por el método de la transformada inversa con \xi_1, \xi_2 \sim U(0, 1):
 *   F_\Theta(\theta) = \int_0^\theta \sin\theta' d\theta' = 1 - \cos\theta = \xi_1
 *   => \cos\theta = 1 - \xi_1  (o equivalentemente \xi_1), \sin\theta = \sqrt{1 - \cos^2\theta}
 *   \phi = 2\pi \xi_2
 */
inline HemisphereSample sample_hemisphere_uniform(const Vec3& normal, double xi1, double xi2) {
    double cos_theta = xi1;
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double phi = TWO_PI * xi2;

    Vec3 local_dir(
        sin_theta * std::cos(phi),
        sin_theta * std::sin(phi),
        cos_theta
    );

    OrthonormalBasis onb(normal);
    Vec3 world_dir = normalize(onb.to_world(local_dir));

    return HemisphereSample{
        world_dir,
        cos_theta,
        1.0 / TWO_PI // p(\omega) = 1 / (2\pi)
    };
}

/**
 * @brief 2. Muestreo Hemisférico Ponderado por el Coseno (Cosine-Weighted Importance Sampling):
 *
 * Densidad respecto al ángulo sólido d\omega = \sin\theta d\theta d\phi:
 *   p(\omega) = \frac{\cos\theta}{\pi}
 *
 * Verificación de normalización en el hemisferio \Omega:
 *   \int_0^{2\pi} \int_0^{\pi/2} \frac{\cos\theta}{\pi} \sin\theta d\theta d\phi
 *   = 2 \int_0^{\pi/2} \cos\theta \sin\theta d\theta = 1.
 *
 * Por el método de la transformada inversa con \xi_1, \xi_2 \sim U(0, 1):
 *   p(\theta) = 2 \cos\theta \sin\theta
 *   F_\Theta(\theta) = \int_0^\theta 2 \cos\theta' \sin\theta' d\theta' = \sin^2\theta = \xi_1
 *   => \sin\theta = \sqrt{\xi_1}, \quad \cos\theta = \sqrt{1 - \xi_1}
 *   \phi = 2\pi \xi_2
 *
 * Equivalencia geométrica (Ley de Malus): muestrear uniformemente el disco unitario (r = \sqrt{\xi_1}, \phi = 2\pi\xi_2)
 * y proyectar ortogonalmente hacia la cúpula hemisférica z = \sqrt{1 - r^2}.
 */
inline HemisphereSample sample_hemisphere_cosine(const Vec3& normal, double xi1, double xi2) {
    double cos_theta = std::sqrt(std::max(0.0, 1.0 - xi1));
    double sin_theta = std::sqrt(std::max(0.0, xi1));
    double phi = TWO_PI * xi2;

    Vec3 local_dir(
        sin_theta * std::cos(phi),
        sin_theta * std::sin(phi),
        cos_theta
    );

    OrthonormalBasis onb(normal);
    Vec3 world_dir = normalize(onb.to_world(local_dir));

    return HemisphereSample{
        world_dir,
        cos_theta,
        cos_theta * INV_PI // p(\omega) = \cos\theta / \pi
    };
}

} // namespace kajiya
