#include "benchmark/benchmark_suite.hpp"
#include "gui/win32_viewer.hpp"
#include <iostream>
#include <string>
#include <thread>

int main(int argc, char* argv[]) {
    bool launch_gui = false;
    bool run_bench = true;
    int width = 200;
    int height = 200;
    int reference_spp = 4096;
    std::string output_dir = "output";
    std::string custom_obj_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--gui") {
            launch_gui = true;
            run_bench = false;
        } else if (arg == "--all") {
            launch_gui = true;
            run_bench = true;
        } else if (arg == "--benchmark") {
            run_bench = true;
            launch_gui = false;
        } else if (arg == "--obj" && i + 1 < argc) {
            custom_obj_path = argv[++i];
            run_bench = false;
        } else if (arg == "--ref-spp" && i + 1 < argc) {
            reference_spp = std::stoi(argv[++i]);
        } else if (arg == "--size" && i + 1 < argc) {
            width = std::stoi(argv[++i]);
            height = width;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Uso: kajiya_pathtracer [opciones]\n"
                      << "  --benchmark       Ejecuta Fase 1 (Furnace Test) y Fase 2/3 (Cornell Box) y exporta imagenes/CSV/HTML (por defecto)\n"
                      << "  --obj <ruta.obj>  Importa un modelo 3D Wavefront (.obj) dentro de la Cornell Box con aceleracion BVH y lo renderiza\n"
                      << "  --gui             Abre la ventana interactiva en tiempo real (Win32 GUI)\n"
                      << "  --all             Ejecuta el benchmark completo y luego abre la interfaz grafica\n"
                      << "  --size <px>       Resolucion cuadrada del benchmark (por defecto: 200)\n"
                      << "  --ref-spp <N>     Muestras por pixel para la referencia Ground Truth (por defecto: 4096)\n";
            return 0;
        }
    }

    std::cout << "========================================================================================\n";
    std::cout << "  MOTOR DE RENDERIZADO PATH TRACER DESDE CERO EN C++ (ECUACION INTEGRAL DE KAJIYA)\n";
    std::cout << "  Hilos paralelos detectados (std::thread): " << std::thread::hardware_concurrency() << " nucleos logicos\n";
    std::cout << "========================================================================================\n";

    if (run_bench) {
        auto furnace_results = kajiya::BenchmarkSuite::run_furnace_validation(output_dir);
        auto cornell_results = kajiya::BenchmarkSuite::run_cornell_benchmark(output_dir, width, height, reference_spp);
        kajiya::BenchmarkSuite::generate_html_dashboard(output_dir, furnace_results, cornell_results, reference_spp);

        std::cout << "\n[OK] Exportacion finalizada en el directorio '" << output_dir << "/':\n"
                  << "  - Imagenes de radiancia (.ppm y .png) para N in {1, 4, 16, 64, 256, 1024} y referencia N = " << reference_spp << "\n"
                  << "  - Mapas de calor de varianza (variance_uniform_spp*.png / variance_cosine_spp*.png)\n"
                  << "  - Datos cuantitativos CSV: " << output_dir << "/furnace_test_results.csv y " << output_dir << "/benchmark_results.csv\n"
                  << "  - Dashboard interactivo HTML: " << output_dir << "/convergence_dashboard.html\n";
    }

    if (!custom_obj_path.empty()) {
        std::filesystem::create_directories(output_dir);
        kajiya::ObjLoadStats stats;
        kajiya::Scene obj_scene = kajiya::Scene::create_cornell_with_obj(custom_obj_path, &stats);
        if (!stats.success) {
            std::cerr << "[ERROR] " << stats.error_message << "\n";
            return 1;
        }
        std::cout << "\n[OBJ IMPORTADO] Archivo: " << custom_obj_path << "\n"
                  << "  - Vertices leidos: " << stats.vertex_count << "\n"
                  << "  - Triangulos del modelo: " << stats.triangle_count << "\n"
                  << "  - Total triangulos en escena (con paredes y pedestal): " << obj_scene.triangles.size() << "\n"
                  << "  - Nodos del arbol BVH construidos: " << obj_scene.bvh.nodes.size() << "\n";

        kajiya::Camera cam = kajiya::Camera::cornell_camera(width, height);
        kajiya::RenderConfig cfg;
        cfg.width = width;
        cfg.height = height;
        cfg.spp = 256;
        cfg.max_bounces = 5;
        cfg.strategy = kajiya::SamplingStrategy::CosineWeightedHemisphere;

        kajiya::RenderBuffer buf_cos(width, height);
        double t_ms = kajiya::PathTracer::render_pass(obj_scene, cam, cfg, buf_cos, 0, cfg.spp);

        kajiya::ImageExporter::export_radiance_image(buf_cos, output_dir + "/cornell_obj_cosine_spp256");
        kajiya::ImageExporter::export_variance_heatmap(buf_cos, output_dir + "/variance_obj_cosine_spp256", 0.15);

        std::cout << "[OK] Modelo OBJ renderizado a " << cfg.spp << " SPP en " << t_ms << " ms.\n"
                  << "     Guardado en: " << output_dir << "/cornell_obj_cosine_spp256.png (.ppm)\n";
    }

#ifdef _WIN32
    if (launch_gui) {
        std::cout << "\nIniciando visor interactivo en tiempo real (Win32 GUI)...\n";
        kajiya::Win32InteractiveViewer viewer(custom_obj_path.empty() ? "assets/models/star_sculpture.obj" : custom_obj_path);
        if (!custom_obj_path.empty()) viewer.active_scene_mode = 2;
        return viewer.run();
    }
#endif

    return 0;
}
