#pragma once

#include "integrator/path_tracer.hpp"
#include "io/image_export.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <filesystem>

namespace kajiya {

struct FurnaceResultRow {
    int spp;
    double mean_uniform;
    double var_uniform;
    double abs_err_uniform;
    double time_ms_uniform;
    double mean_cosine;
    double var_cosine;
    double abs_err_cosine;
    double time_ms_cosine;
};

struct CornellBenchmarkRow {
    int spp;
    double time_ms_uniform;
    double sample_var_uniform;
    double est_var_uniform;
    double mse_uniform;
    double rmse_uniform;

    double time_ms_cosine;
    double sample_var_cosine;
    double est_var_cosine;
    double mse_cosine;
    double rmse_cosine;

    double variance_reduction_ratio;
    double mse_reduction_ratio;
};

class BenchmarkSuite {
public:
    /**
     * @brief Ejecuta la Fase 1: Validación Analítica de la Integral (Horno de Luz / White Furnace Test).
     *
     * Escena: Esfera unitaria con reflectancia difusa \rho = 1.0 encerrada en un entorno radiante L_i = 1.0.
     * Integral analítica exacta (1 rebote hemisférico sobre la esfera):
     *   L_o = \int_0^{2\pi} \int_0^{\pi/2} \frac{1}{\pi} \cdot 1 \cdot \cos\theta \sin\theta d\theta d\phi = 1.000000
     */
    static std::vector<FurnaceResultRow> run_furnace_validation(const std::string& output_dir) {
        std::filesystem::create_directories(output_dir);

        std::cout << "\n========================================================================================\n";
        std::cout << "  FASE 1: VALIDACION ANALITICA DE LA INTEGRAL (HORNO DE LUZ / WHITE FURNACE TEST)\n";
        std::cout << "  Integral Analitica: L_o = \\int_\\Omega (1/\\pi) * 1 * cos(\\theta) sin(\\theta) d\\theta d\\phi = 1.000000\n";
        std::cout << "========================================================================================\n";

        Scene scene = Scene::create_furnace_test(1.0);
        int width = 160;
        int height = 160;
        Camera cam = Camera::furnace_camera(width, height);

        std::vector<int> spp_levels = {1, 4, 16, 64, 256, 1024};
        std::vector<FurnaceResultRow> results;

        std::cout << std::left
                  << std::setw(8)  << "SPP(N)"
                  << std::setw(15) << "<Lo> Unif"
                  << std::setw(15) << "Var Unif"
                  << std::setw(15) << "|Err| Unif"
                  << std::setw(15) << "<Lo> Coseno"
                  << std::setw(15) << "Var Coseno"
                  << std::setw(15) << "|Err| Coseno"
                  << "\n";
        std::cout << std::string(98, '-') << "\n";

        for (int spp : spp_levels) {
            RenderConfig cfg_u;
            cfg_u.width = width;
            cfg_u.height = height;
            cfg_u.spp = spp;
            cfg_u.max_bounces = 1; // Evaluación directa de la integral hemisférica sobre la esfera
            cfg_u.enable_russian_roulette = false;
            cfg_u.jitter_primary_rays = false; // Rayo fijo al centro de cada píxel para aislar la integral en x
            cfg_u.strategy = SamplingStrategy::UniformHemisphere;

            RenderBuffer buf_u(width, height);
            double t_u = PathTracer::render_pass(scene, cam, cfg_u, buf_u, 0, spp);

            RenderConfig cfg_c = cfg_u;
            cfg_c.strategy = SamplingStrategy::CosineWeightedHemisphere;
            RenderBuffer buf_c(width, height);
            double t_c = PathTracer::render_pass(scene, cam, cfg_c, buf_c, 0, spp);

            // Medir únicamente los píxeles que impactan la esfera unitaria en el centro de la imagen
            double sum_mean_u = 0.0, sum_var_u = 0.0, sum_err_u = 0.0;
            double sum_mean_c = 0.0, sum_var_c = 0.0, sum_err_c = 0.0;
            int sphere_pixel_count = 0;

            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    double u = (static_cast<double>(x) + 0.5) / static_cast<double>(width);
                    double v = (static_cast<double>(height - 1 - y) + 0.5) / static_cast<double>(height);
                    HitRecord rec;
                    if (scene.intersect(cam.get_ray(u, v), 1e-5, 1e30, rec)) {
                        sphere_pixel_count++;
                        double lu = buf_u.at(x, y).mean.x;
                        double lc = buf_c.at(x, y).mean.x;
                        sum_mean_u += lu;
                        sum_mean_c += lc;
                        sum_var_u  += buf_u.at(x, y).sample_variance();
                        sum_var_c  += buf_c.at(x, y).sample_variance();
                        sum_err_u  += std::abs(lu - 1.0);
                        sum_err_c  += std::abs(lc - 1.0);
                    }
                }
            }

            double inv_cnt = 1.0 / std::max(1, sphere_pixel_count);
            FurnaceResultRow row{
                spp,
                sum_mean_u * inv_cnt,
                sum_var_u * inv_cnt,
                sum_err_u * inv_cnt,
                t_u,
                sum_mean_c * inv_cnt,
                sum_var_c * inv_cnt,
                sum_err_c * inv_cnt,
                t_c
            };
            results.push_back(row);

            std::cout << std::fixed << std::setprecision(6)
                      << std::left << std::setw(8) << spp
                      << std::setw(15) << row.mean_uniform
                      << std::setw(15) << row.var_uniform
                      << std::setw(15) << row.abs_err_uniform
                      << std::setw(15) << row.mean_cosine
                      << std::setw(15) << row.var_cosine
                      << std::setw(15) << row.abs_err_cosine
                      << "\n";

            if (spp == 4 || spp == 64) {
                ImageExporter::export_radiance_image(buf_u, output_dir + "/furnace_uniform_spp" + std::to_string(spp));
                ImageExporter::export_radiance_image(buf_c, output_dir + "/furnace_cosine_spp" + std::to_string(spp));
            }
        }

        // Guardar CSV de Fase 1
        std::ofstream csv(output_dir + "/furnace_test_results.csv");
        csv << "spp,mean_uniform,sample_var_uniform,abs_error_uniform,time_ms_uniform,mean_cosine,sample_var_cosine,abs_error_cosine,time_ms_cosine\n";
        for (const auto& r : results) {
            csv << r.spp << ","
                << std::fixed << std::setprecision(8)
                << r.mean_uniform << "," << r.var_uniform << "," << r.abs_err_uniform << "," << r.time_ms_uniform << ","
                << r.mean_cosine << "," << r.var_cosine << "," << r.abs_err_cosine << "," << r.time_ms_cosine << "\n";
        }

        return results;
    }

    /**
     * @brief Ejecuta la Fase 2 y 3: Escena Estándar Cornell Box, Convergencia MSE O(1/sqrt(N))
     *        y Comparación Cuantitativa Muestreo Uniforme vs. Muestreo Ponderado por Coseno.
     */
    static std::vector<CornellBenchmarkRow> run_cornell_benchmark(
        const std::string& output_dir,
        int width = 200,
        int height = 200,
        int reference_spp = 4096
    ) {
        std::filesystem::create_directories(output_dir);

        std::cout << "\n========================================================================================\n";
        std::cout << "  FASE 2 & 3: CORNELL BOX - CONVERGENCIA MONTE CARLO Y MUESTREO POR IMPORTANCIA\n";
        std::cout << "  Resolucion: " << width << "x" << height << " | Referencia Ground Truth: N = " << reference_spp << " SPP\n";
        std::cout << "========================================================================================\n";

        Scene scene = Scene::create_cornell_box();
        Camera cam = Camera::cornell_camera(width, height);

        // 1. Calcular o acumular la imagen de referencia Ground Truth (N = reference_spp) con Muestreo por Coseno
        std::cout << "[1/2] Calculando imagen de referencia Ground Truth (N = " << reference_spp << " SPP)...\n";
        RenderConfig ref_cfg;
        ref_cfg.width = width;
        ref_cfg.height = height;
        ref_cfg.spp = reference_spp;
        ref_cfg.max_bounces = 5;
        ref_cfg.enable_russian_roulette = true;
        ref_cfg.jitter_primary_rays = true;
        ref_cfg.strategy = SamplingStrategy::CosineWeightedHemisphere;

        RenderBuffer reference_buffer(width, height);
        // Acumulamos con un offset de semilla distinto para que la referencia sea estadísticamente independiente
        double ref_time_ms = PathTracer::render_pass(scene, cam, ref_cfg, reference_buffer, 100000, reference_spp);
        std::cout << "      Referencia completada en " << std::fixed << std::setprecision(1) << (ref_time_ms / 1000.0) << " s.\n";
        ImageExporter::export_radiance_image(reference_buffer, output_dir + "/cornell_reference_spp" + std::to_string(reference_spp));

        // 2. Evaluar potencias de dos N in {1, 4, 16, 64, 256, 1024} de forma progresiva para ambos métodos
        std::cout << "[2/2] Evaluando convergencia para N in {1, 4, 16, 64, 256, 1024} SPP (Uniforme vs. Coseno)...\n\n";

        std::vector<int> spp_targets = {1, 4, 16, 64, 256, 1024};
        std::vector<CornellBenchmarkRow> rows;

        RenderBuffer buf_uniform(width, height);
        RenderBuffer buf_cosine(width, height);

        RenderConfig cfg_uniform = ref_cfg;
        cfg_uniform.strategy = SamplingStrategy::UniformHemisphere;

        RenderConfig cfg_cosine = ref_cfg;
        cfg_cosine.strategy = SamplingStrategy::CosineWeightedHemisphere;

        int current_spp = 0;
        double cum_time_uniform = 0.0;
        double cum_time_cosine = 0.0;

        std::cout << std::left
                  << std::setw(8)  << "SPP(N)"
                  << std::setw(14) << "MSE Unif"
                  << std::setw(14) << "MSE Coseno"
                  << std::setw(14) << "RMSE Coseno"
                  << std::setw(14) << "Var<Lo> Unif"
                  << std::setw(14) << "Var<Lo> Cos"
                  << std::setw(13) << "ReduccionVar"
                  << std::setw(12) << "Tiempo(ms)"
                  << "\n";
        std::cout << std::string(103, '-') << "\n";

        for (int target_n : spp_targets) {
            int delta_spp = target_n - current_spp;

            cum_time_uniform += PathTracer::render_pass(scene, cam, cfg_uniform, buf_uniform, current_spp, delta_spp);
            cum_time_cosine  += PathTracer::render_pass(scene, cam, cfg_cosine,  buf_cosine,  current_spp, delta_spp);
            current_spp = target_n;

            double mse_u = buf_uniform.compute_mse(reference_buffer);
            double mse_c = buf_cosine.compute_mse(reference_buffer);
            double rmse_u = std::sqrt(mse_u);
            double rmse_c = std::sqrt(mse_c);

            double svar_u = buf_uniform.average_sample_variance();
            double svar_c = buf_cosine.average_sample_variance();
            double evar_u = (target_n >= 2) ? buf_uniform.average_estimator_variance() : mse_u;
            double evar_c = (target_n >= 2) ? buf_cosine.average_estimator_variance()  : mse_c;

            double var_ratio = (evar_c > 1e-12) ? (evar_u / evar_c) : (mse_u / std::max(1e-12, mse_c));
            double mse_ratio = (mse_c > 1e-12) ? (mse_u / mse_c) : 1.0;

            CornellBenchmarkRow row{
                target_n,
                cum_time_uniform,
                svar_u,
                evar_u,
                mse_u,
                rmse_u,
                cum_time_cosine,
                svar_c,
                evar_c,
                mse_c,
                rmse_c,
                var_ratio,
                mse_ratio
            };
            rows.push_back(row);

            std::cout << std::fixed << std::setprecision(6)
                      << std::left << std::setw(8) << target_n
                      << std::setw(14) << mse_u
                      << std::setw(14) << mse_c
                      << std::setw(14) << rmse_c
                      << std::setw(14) << evar_u
                      << std::setw(14) << evar_c
                      << std::setprecision(2) << std::setw(13) << (std::to_string(var_ratio).substr(0, 4) + "x")
                      << std::setprecision(1) << std::setw(12) << cum_time_cosine
                      << "\n";

            // Exportar imágenes PPM y PNG para cada nivel N de muestras
            std::string tag = std::to_string(target_n);
            ImageExporter::export_radiance_image(buf_uniform, output_dir + "/cornell_uniform_spp" + tag);
            ImageExporter::export_radiance_image(buf_cosine,  output_dir + "/cornell_cosine_spp" + tag);

            if (target_n >= 4) {
                ImageExporter::export_variance_heatmap(buf_uniform, output_dir + "/variance_uniform_spp" + tag, 0.25);
                ImageExporter::export_variance_heatmap(buf_cosine,  output_dir + "/variance_cosine_spp" + tag, 0.25);
            }
        }

        // Exportar CSV completo
        std::ofstream csv(output_dir + "/benchmark_results.csv");
        csv << "spp,time_ms_uniform,sample_var_uniform,estimator_var_uniform,mse_uniform,rmse_uniform,"
            << "time_ms_cosine,sample_var_cosine,estimator_var_cosine,mse_cosine,rmse_cosine,"
            << "variance_reduction_ratio,mse_reduction_ratio\n";
        for (const auto& r : rows) {
            csv << r.spp << ","
                << std::fixed << std::setprecision(8)
                << r.time_ms_uniform << "," << r.sample_var_uniform << "," << r.est_var_uniform << "," << r.mse_uniform << "," << r.rmse_uniform << ","
                << r.time_ms_cosine << "," << r.sample_var_cosine << "," << r.est_var_cosine << "," << r.mse_cosine << "," << r.rmse_cosine << ","
                << r.variance_reduction_ratio << "," << r.mse_reduction_ratio << "\n";
        }

        return rows;
    }

    /**
     * @brief Genera un reporte interactivo HTML5 + SVG con gráficas log-log de convergencia,
     *        tablas comparativas y galería visual de los renders y mapas de calor de varianza.
     */
    static void generate_html_dashboard(
        const std::string& output_dir,
        const std::vector<FurnaceResultRow>& furnace_rows,
        const std::vector<CornellBenchmarkRow>& cornell_rows,
        int reference_spp
    ) {
        std::ofstream html(output_dir + "/convergence_dashboard.html");
        if (!html.is_open()) return;

        html << R"HTML(<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>Reporte de Convergencia Monte Carlo - Ecuación de Kajiya</title>
<style>
  :root {
    --bg: #0f172a;
    --card: #1e293b;
    --border: #334155;
    --text: #f8fafc;
    --muted: #94a3b8;
    --accent-cos: #38bdf8;
    --accent-uni: #f97316;
    --accent-ref: #4ade80;
  }
  body {
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    background: var(--bg);
    color: var(--text);
    margin: 0;
    padding: 32px;
    line-height: 1.6;
  }
  .container { max-width: 1280px; margin: 0 auto; }
  h1, h2, h3 { color: #fff; margin-top: 0; }
  .header-box {
    background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 24px 32px;
    margin-bottom: 28px;
  }
  .formula {
    background: #090d16;
    border-left: 4px solid var(--accent-cos);
    padding: 14px 20px;
    font-family: 'Cambria Math', 'Consolas', monospace;
    font-size: 1.05rem;
    border-radius: 6px;
    margin: 14px 0;
    color: #e2e8f0;
  }
  .grid-2 {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 24px;
    margin-bottom: 28px;
  }
  .card {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 22px;
    margin-bottom: 24px;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.92rem;
    margin-top: 12px;
  }
  th, td {
    padding: 10px 12px;
    text-align: right;
    border-bottom: 1px solid var(--border);
  }
  th { background: #0f172a; color: var(--accent-cos); text-align: right; }
  td:first-child, th:first-child { text-align: left; font-weight: 600; }
  .badge-good {
    background: rgba(74, 222, 128, 0.15);
    color: #4ade80;
    padding: 3px 8px;
    border-radius: 999px;
    font-weight: 600;
  }
  .gallery {
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    gap: 12px;
    margin-top: 14px;
  }
  .gallery-item {
    background: #090d16;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 8px;
    text-align: center;
  }
  .gallery-item img {
    width: 100%;
    height: auto;
    border-radius: 4px;
    image-rendering: pixelated;
  }
  .gallery-item span {
    display: block;
    font-size: 0.8rem;
    color: var(--muted);
    margin-top: 6px;
  }
</style>
</head>
<body>
<div class="container">
  <div class="header-box">
    <h1>Motor Path Tracer en C++ — Resolución Numérica de la Ecuación de Kajiya</h1>
    <p style="color:var(--muted); margin-bottom:10px;">
      Integración Hemisférica de Monte Carlo en Coordenadas Esféricas (d&omega; = sin&theta; d&theta; d&phi;) y Muestreo por Importancia Ponderado por Coseno
    </p>
    <div class="formula">
      L<sub>o</sub>(x, &omega;<sub>o</sub>) = L<sub>e</sub>(x, &omega;<sub>o</sub>) + &int;<sub>&Omega;</sub> f<sub>r</sub>(x, &omega;<sub>i</sub>, &omega;<sub>o</sub>) L<sub>i</sub>(x, &omega;<sub>i</sub>) (&omega;<sub>i</sub> &middot; n) d&omega;<sub>i</sub>
      &nbsp;&nbsp;&rArr;&nbsp;&nbsp;
      &langle;L<sub>o</sub>&rangle;<sub>N</sub> = (1/N) &Sigma; [ f<sub>r</sub> L<sub>i</sub> cos&theta;<sub>k</sub> / p(&omega;<sub>k</sub>) ]
    </div>
  </div>
)HTML";

        // Tabla Fase 1: Furnace Test
        html << R"HTML(
  <div class="card">
    <h2>Fase 1: Validación Analítica de la Integral — Horno de Luz (White Furnace Test)</h2>
    <p style="color:var(--muted)">
      Esfera unitaria 100% reflectora difusa (&rho; = 1.0, f<sub>r</sub> = 1/&pi;) en un entorno isotrópico L<sub>i</sub> = 1.0.
      En equilibrio termodinámico, &int;<sub>&Omega;</sub> (1/&pi;) cos&theta; sin&theta; d&theta; d&phi; = <strong>1.000000</strong>.
      Con muestreo por coseno p(&omega;) = cos&theta;/&pi;, el factor cos&theta; se cancela algebraicamente, produciendo varianza cero en 1 rebote.
    </p>
    <table>
      <thead>
        <tr>
          <th>Muestras (N)</th>
          <th>&langle;L<sub>o</sub>&rangle; Uniforme</th>
          <th>Varianza Muestral Unif</th>
          <th>Error |&langle;L<sub>o</sub>&rangle; - 1| Unif</th>
          <th>&langle;L<sub>o</sub>&rangle; Coseno</th>
          <th>Varianza Muestral Coseno</th>
          <th>Error |&langle;L<sub>o</sub>&rangle; - 1| Coseno</th>
        </tr>
      </thead>
      <tbody>
)HTML";
        for (const auto& r : furnace_rows) {
            html << std::fixed << std::setprecision(6)
                 << "<tr><td>" << r.spp << " SPP</td>"
                 << "<td>" << r.mean_uniform << "</td>"
                 << "<td>" << r.var_uniform << "</td>"
                 << "<td>" << r.abs_err_uniform << "</td>"
                 << "<td><span class='badge-good'>" << r.mean_cosine << "</span></td>"
                 << "<td><span class='badge-good'>" << r.var_cosine << "</span></td>"
                 << "<td><span class='badge-good'>" << r.abs_err_cosine << "</span></td></tr>\n";
        }
        html << R"HTML(      </tbody>
    </table>
  </div>
)HTML";

        // Tabla Fase 2 & 3: Cornell Box
        html << R"HTML(
  <div class="card">
    <h2>Fase 2 y 3: Escena Estándar Cornell Box — Convergencia MSE y Reducción de Varianza</h2>
    <p style="color:var(--muted)">
      Comparación cuantitativa entre Muestreo Hemisférico Uniforme p(&omega;) = 1/(2&pi;) y Muestreo Ponderado por Coseno p(&omega;) = cos&theta;/&pi;
      frente a la imagen de referencia calculada con N = )HTML" << reference_spp << R"HTML( SPP.
    </p>
    <table>
      <thead>
        <tr>
          <th>SPP (N = 2<sup>k</sup>)</th>
          <th>MSE Uniforme</th>
          <th>MSE Coseno</th>
          <th>RMSE Coseno O(1/&radic;N)</th>
          <th>Var(&langle;L<sub>o</sub>&rangle;) Uniforme</th>
          <th>Var(&langle;L<sub>o</sub>&rangle;) Coseno</th>
          <th>Reducción de Varianza</th>
          <th>Tiempo Coseno (ms)</th>
        </tr>
      </thead>
      <tbody>
)HTML";
        for (const auto& r : cornell_rows) {
            html << std::fixed << std::setprecision(6)
                 << "<tr><td>" << r.spp << " SPP</td>"
                 << "<td>" << r.mse_uniform << "</td>"
                 << "<td>" << r.mse_cosine << "</td>"
                 << "<td>" << r.rmse_cosine << "</td>"
                 << "<td>" << r.est_var_uniform << "</td>"
                 << "<td>" << r.est_var_cosine << "</td>"
                 << "<td><span class='badge-good'>" << std::setprecision(2) << r.variance_reduction_ratio << "x menor varianza</span></td>"
                 << "<td>" << std::setprecision(1) << r.time_ms_cosine << " ms</td></tr>\n";
        }
        html << R"HTML(      </tbody>
    </table>
  </div>
)HTML";

        // Galería visual comparativa Cornell Box
        html << R"HTML(
  <div class="card">
    <h2>Progresión Visual de Convergencia: Muestreo por Coseno p(&omega;) = cos&theta;/&pi;</h2>
    <div class="gallery">
)HTML";
        for (const auto& r : cornell_rows) {
            html << "      <div class='gallery-item'><img src='cornell_cosine_spp" << r.spp << ".png' alt='Coseno SPP " << r.spp << "'>"
                 << "<span><strong>N = " << r.spp << " SPP</strong><br>MSE: " << std::fixed << std::setprecision(5) << r.mse_cosine << "</span></div>\n";
        }
        html << R"HTML(    </div>
    <h2 style="margin-top:24px;">Progresión Visual de Convergencia: Muestreo Uniforme p(&omega;) = 1/(2&pi;)</h2>
    <div class="gallery">
)HTML";
        for (const auto& r : cornell_rows) {
            html << "      <div class='gallery-item'><img src='cornell_uniform_spp" << r.spp << ".png' alt='Uniforme SPP " << r.spp << "'>"
                 << "<span><strong>N = " << r.spp << " SPP</strong><br>MSE: " << std::fixed << std::setprecision(5) << r.mse_uniform << "</span></div>\n";
        }
        html << R"HTML(    </div>
  </div>
</div>
</body>
</html>
)HTML";
    }
};

} // namespace kajiya
