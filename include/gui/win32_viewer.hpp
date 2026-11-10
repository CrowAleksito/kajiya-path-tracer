#pragma once

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#include "integrator/path_tracer.hpp"
#include "io/image_export.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace kajiya {

#ifdef _WIN32

struct HistoryPoint {
    int spp;
    double est_variance;
    double cumulative_ms;
};

class Win32InteractiveViewer {
public:
    int render_w = 240;
    int render_h = 240;
    int win_w = 1120;
    int win_h = 620;

    Scene cornell_scene;
    Scene furnace_scene;
    bool use_furnace = false;
    bool show_variance_heatmap = false;
    bool paused = false;

    RenderConfig cfg_cosine;
    RenderConfig cfg_uniform;
    SamplingStrategy active_strategy = SamplingStrategy::CosineWeightedHemisphere;

    RenderBuffer buf_cosine;
    RenderBuffer buf_uniform;

    int spp_cosine = 0;
    int spp_uniform = 0;
    double time_ms_cosine = 0.0;
    double time_ms_uniform = 0.0;

    std::vector<HistoryPoint> hist_cosine;
    std::vector<HistoryPoint> hist_uniform;

    std::vector<uint32_t> canvas_pixels; // BGRA 32-bit para StretchDIBits

    Win32InteractiveViewer()
        : cornell_scene(Scene::create_cornell_box()),
          furnace_scene(Scene::create_furnace_test(1.0)),
          buf_cosine(render_w, render_h),
          buf_uniform(render_w, render_h),
          canvas_pixels(static_cast<size_t>(win_w * win_h), 0xFF0F172Au)
    {
        cfg_cosine.width = render_w;
        cfg_cosine.height = render_h;
        cfg_cosine.max_bounces = 5;
        cfg_cosine.strategy = SamplingStrategy::CosineWeightedHemisphere;

        cfg_uniform = cfg_cosine;
        cfg_uniform.strategy = SamplingStrategy::UniformHemisphere;
    }

    void reset_buffers() {
        buf_cosine.clear();
        buf_uniform.clear();
        spp_cosine = 0;
        spp_uniform = 0;
        time_ms_cosine = 0.0;
        time_ms_uniform = 0.0;
        hist_cosine.clear();
        hist_uniform.clear();
    }

    void step_simulation() {
        if (paused || (spp_cosine >= 1024 && spp_uniform >= 1024)) return;

        const Scene& active_scene = use_furnace ? furnace_scene : cornell_scene;
        Camera cam = use_furnace
            ? Camera::furnace_camera(render_w, render_h)
            : Camera::cornell_camera(render_w, render_h);

        cfg_cosine.max_bounces = use_furnace ? 1 : 5;
        cfg_uniform.max_bounces = use_furnace ? 1 : 5;

        // Avanzamos 1 muestra en ambos integradores (o hasta 1024 SPP) para graficar en vivo la comparación
        int step = (spp_cosine < 16) ? 1 : ((spp_cosine < 64) ? 2 : 4);

        if (spp_cosine < 1024) {
            time_ms_cosine += PathTracer::render_pass(active_scene, cam, cfg_cosine, buf_cosine, spp_cosine, step);
            spp_cosine += step;
            if (spp_cosine >= 2) {
                hist_cosine.push_back({spp_cosine, buf_cosine.average_estimator_variance(), time_ms_cosine});
            }
        }
        if (spp_uniform < 1024) {
            time_ms_uniform += PathTracer::render_pass(active_scene, cam, cfg_uniform, buf_uniform, spp_uniform, step);
            spp_uniform += step;
            if (spp_uniform >= 2) {
                hist_uniform.push_back({spp_uniform, buf_uniform.average_estimator_variance(), time_ms_uniform});
            }
        }
    }

    inline void put_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b) {
        if (x < 0 || x >= win_w || y < 0 || y >= win_h) return;
        canvas_pixels[static_cast<size_t>(y * win_w + x)] =
            (0xFFu << 24) | (static_cast<uint32_t>(r) << 16) | (static_cast<uint32_t>(g) << 8) | static_cast<uint32_t>(b);
    }

    void fill_rect(int x0, int y0, int w, int h, uint8_t r, uint8_t g, uint8_t b) {
        for (int y = y0; y < y0 + h; ++y) {
            for (int x = x0; x < x0 + w; ++x) {
                put_pixel(x, y, r, g, b);
            }
        }
    }

    void draw_line(int x0, int y0, int x1, int y1, uint8_t r, uint8_t g, uint8_t b, int thickness = 2) {
        int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
        int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
        int err = dx + dy;
        while (true) {
            for (int ty = -thickness / 2; ty <= thickness / 2; ++ty)
                for (int tx = -thickness / 2; tx <= thickness / 2; ++tx)
                    put_pixel(x0 + tx, y0 + ty, r, g, b);
            if (x0 == x1 && y0 == y1) break;
            int e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; }
            if (e2 <= dx) { err += dx; y0 += sy; }
        }
    }

    void render_to_canvas() {
        // Fondo general oscuro
        fill_rect(0, 0, win_w, win_h, 15, 23, 42);

        const RenderBuffer& active_buf =
            (active_strategy == SamplingStrategy::CosineWeightedHemisphere) ? buf_cosine : buf_uniform;

        // 1. Dibujar el viewport 2D escalado a 480x480 en (24, 96)
        int vp_x = 24, vp_y = 96, vp_size = 480;
        fill_rect(vp_x - 2, vp_y - 2, vp_size + 4, vp_size + 4, 56, 189, 248);

        for (int py = 0; py < vp_size; ++py) {
            int sy = (py * render_h) / vp_size;
            for (int px = 0; px < vp_size; ++px) {
                int sx = (px * render_w) / vp_size;
                const PixelAccumulator& acc = active_buf.at(sx, sy);
                uint8_t r, g, b;
                if (!show_variance_heatmap) {
                    r = to_srgb_byte(acc.mean.x);
                    g = to_srgb_byte(acc.mean.y);
                    b = to_srgb_byte(acc.mean.z);
                } else {
                    double v = acc.estimator_variance();
                    double norm = std::sqrt(std::clamp(v / 0.25, 0.0, 1.0));
                    Color hc = heatmap_color(norm);
                    r = static_cast<uint8_t>(hc.x * 255.0);
                    g = static_cast<uint8_t>(hc.y * 255.0);
                    b = static_cast<uint8_t>(hc.z * 255.0);
                }
                put_pixel(vp_x + px, vp_y + py, r, g, b);
            }
        }

        // 2. Gráfico 1: Varianza del Estimador Var(<L_o>) vs. SPP en (535, 96, 555, 225)
        int gx = 535, gy = 96, gw = 555, gh = 225;
        fill_rect(gx, gy, gw, gh, 22, 32, 56);
        // Cuadrícula horizontal
        for (int i = 1; i < 5; ++i) {
            int ly = gy + (gh * i) / 5;
            draw_line(gx + 40, ly, gx + gw - 15, ly, 45, 58, 84, 1);
        }

        if (hist_uniform.size() >= 2 && hist_cosine.size() >= 2) {
            double max_var = std::max(hist_uniform.front().est_variance, hist_cosine.front().est_variance);
            if (max_var < 1e-9) max_var = 1.0;
            int max_spp = std::max(16, hist_cosine.back().spp);

            auto plot_curve = [&](const std::vector<HistoryPoint>& hist, uint8_t cr, uint8_t cg, uint8_t cb) {
                for (size_t i = 1; i < hist.size(); ++i) {
                    double x0_t = static_cast<double>(hist[i - 1].spp) / static_cast<double>(max_spp);
                    double x1_t = static_cast<double>(hist[i].spp) / static_cast<double>(max_spp);
                    // Escala raíz para apreciar toda la caída de varianza
                    double y0_t = std::sqrt(std::clamp(hist[i - 1].est_variance / max_var, 0.0, 1.0));
                    double y1_t = std::sqrt(std::clamp(hist[i].est_variance / max_var, 0.0, 1.0));

                    int px0 = gx + 45 + static_cast<int>(x0_t * (gw - 65));
                    int py0 = gy + gh - 25 - static_cast<int>(y0_t * (gh - 50));
                    int px1 = gx + 45 + static_cast<int>(x1_t * (gw - 65));
                    int py1 = gy + gh - 25 - static_cast<int>(y1_t * (gh - 50));
                    draw_line(px0, py0, px1, py1, cr, cg, cb, 2);
                }
            };

            // Naranja = Uniforme p(\omega) = 1/(2\pi), Celeste = Coseno p(\omega) = \cos\theta/\pi
            plot_curve(hist_uniform, 249, 115, 22);
            plot_curve(hist_cosine, 56, 189, 248);
        }

        // 3. Gráfico 2: Tiempo de cómputo acumulado (ms) vs. SPP en (535, 350, 555, 225)
        int tx = 535, ty = 350, tw = 555, th = 225;
        fill_rect(tx, ty, tw, th, 22, 32, 56);
        for (int i = 1; i < 5; ++i) {
            int ly = ty + (th * i) / 5;
            draw_line(tx + 40, ly, tx + tw - 15, ly, 45, 58, 84, 1);
        }

        if (hist_cosine.size() >= 2) {
            double max_ms = std::max(hist_cosine.back().cumulative_ms, hist_uniform.back().cumulative_ms);
            if (max_ms < 1.0) max_ms = 1.0;
            int max_spp = std::max(16, hist_cosine.back().spp);

            for (size_t i = 1; i < hist_cosine.size(); ++i) {
                double x0_t = static_cast<double>(hist_cosine[i - 1].spp) / static_cast<double>(max_spp);
                double x1_t = static_cast<double>(hist_cosine[i].spp) / static_cast<double>(max_spp);
                double y0_t = hist_cosine[i - 1].cumulative_ms / max_ms;
                double y1_t = hist_cosine[i].cumulative_ms / max_ms;

                int px0 = tx + 45 + static_cast<int>(x0_t * (tw - 65));
                int py0 = ty + th - 25 - static_cast<int>(y0_t * (th - 50));
                int px1 = tx + 45 + static_cast<int>(x1_t * (tw - 65));
                int py1 = ty + th - 25 - static_cast<int>(y1_t * (th - 50));
                draw_line(px0, py0, px1, py1, 74, 222, 128, 2);
            }
        }
    }

    void draw_overlay_text(HDC hdc) {
        SetBkMode(hdc, TRANSPARENT);
        HFONT hFontTitle = CreateFontA(20, 0, 0, 0, FW_BOLD, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
            OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_SWISS, "Segoe UI");
        HFONT hFontSub = CreateFontA(15, 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
            OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_SWISS, "Consolas");

        HFONT oldFont = (HFONT)SelectObject(hdc, hFontTitle);
        SetTextColor(hdc, RGB(248, 250, 252));
        std::string title = "Motor Path Tracer C++ (Ecuacion de Kajiya) - Visualizador de Convergencia Monte Carlo";
        TextOutA(hdc, 24, 16, title.c_str(), static_cast<int>(title.size()));

        SelectObject(hdc, hFontSub);
        SetTextColor(hdc, RGB(56, 189, 248));
        std::string controls = "Controles: [1] Muestreo Coseno  [2] Muestreo Uniforme  [V] Ver Mapa Varianza/RGB  [F] Cambiar Escena (Cornell/Furnace)  [R] Reiniciar  [S] Guardar PNG/PPM";
        TextOutA(hdc, 24, 44, controls.c_str(), static_cast<int>(controls.size()));

        double var_u = buf_uniform.average_estimator_variance();
        double var_c = buf_cosine.average_estimator_variance();
        double ratio = (var_c > 1e-12) ? (var_u / var_c) : 1.0;

        std::ostringstream status;
        status << "Escena: " << (use_furnace ? "Fase 1 Horno de Luz (Furnace)" : "Fase 2 Cornell Box")
               << "  |  Vista: " << (active_strategy == SamplingStrategy::CosineWeightedHemisphere ? "COSENO p(w)=cos(theta)/pi" : "UNIFORME p(w)=1/(2*pi)")
               << (show_variance_heatmap ? " [MAPA DE VARIANZA]" : " [RADIANCIA RGB]")
               << "  |  SPP: " << spp_cosine
               << "  |  Var(Unif): " << std::fixed << std::setprecision(5) << var_u
               << "  |  Var(Cos): " << var_c
               << " (" << std::setprecision(2) << ratio << "x menor varianza)";

        SetTextColor(hdc, RGB(74, 222, 128));
        std::string st = status.str();
        TextOutA(hdc, 24, 68, st.c_str(), static_cast<int>(st.size()));

        // Etiquetas de los gráficos
        SetTextColor(hdc, RGB(248, 250, 252));
        std::string g1 = "Varianza del Estimador Var(<Lo>) vs. SPP  [Naranja: Uniforme | Celeste: Coseno]";
        TextOutA(hdc, 548, 104, g1.c_str(), static_cast<int>(g1.size()));

        std::ostringstream g2;
        g2 << "Tiempo de Computo Acumulado vs. SPP: " << std::fixed << std::setprecision(1)
           << time_ms_cosine << " ms (" << std::thread::hardware_concurrency() << " hilos CPU)";
        std::string g2s = g2.str();
        TextOutA(hdc, 548, 358, g2s.c_str(), static_cast<int>(g2s.size()));

        SelectObject(hdc, oldFont);
        DeleteObject(hFontTitle);
        DeleteObject(hFontSub);
    }

    static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
        Win32InteractiveViewer* viewer = reinterpret_cast<Win32InteractiveViewer*>(GetWindowLongPtrA(hwnd, GWLP_USERDATA));
        switch (msg) {
            case WM_CREATE: {
                CREATESTRUCTA* cs = reinterpret_cast<CREATESTRUCTA*>(lParam);
                SetWindowLongPtrA(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(cs->lpCreateParams));
                return 0;
            }
            case WM_KEYDOWN: {
                if (!viewer) break;
                if (wParam == '1') {
                    viewer->active_strategy = SamplingStrategy::CosineWeightedHemisphere;
                } else if (wParam == '2') {
                    viewer->active_strategy = SamplingStrategy::UniformHemisphere;
                } else if (wParam == 'V' || wParam == 'v') {
                    viewer->show_variance_heatmap = !viewer->show_variance_heatmap;
                } else if (wParam == 'F' || wParam == 'f') {
                    viewer->use_furnace = !viewer->use_furnace;
                    viewer->reset_buffers();
                } else if (wParam == 'R' || wParam == 'r') {
                    viewer->reset_buffers();
                } else if (wParam == 'S' || wParam == 's') {
                    ImageExporter::export_radiance_image(viewer->buf_cosine, "output/snapshot_cosine");
                    ImageExporter::export_radiance_image(viewer->buf_uniform, "output/snapshot_uniform");
                } else if (wParam == VK_ESCAPE) {
                    PostQuitMessage(0);
                }
                InvalidateRect(hwnd, nullptr, FALSE);
                return 0;
            }
            case WM_PAINT: {
                PAINTSTRUCT ps;
                HDC hdc = BeginPaint(hwnd, &ps);
                if (viewer) {
                    BITMAPINFO bmi{};
                    bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
                    bmi.bmiHeader.biWidth = viewer->win_w;
                    bmi.bmiHeader.biHeight = -viewer->win_h; // Top-down DIB
                    bmi.bmiHeader.biPlanes = 1;
                    bmi.bmiHeader.biBitCount = 32;
                    bmi.bmiHeader.biCompression = BI_RGB;

                    StretchDIBits(
                        hdc,
                        0, 0, viewer->win_w, viewer->win_h,
                        0, 0, viewer->win_w, viewer->win_h,
                        viewer->canvas_pixels.data(),
                        &bmi,
                        DIB_RGB_COLORS,
                        SRCCOPY
                    );
                    viewer->draw_overlay_text(hdc);
                }
                EndPaint(hwnd, &ps);
                return 0;
            }
            case WM_DESTROY:
                PostQuitMessage(0);
                return 0;
        }
        return DefWindowProcA(hwnd, msg, wParam, lParam);
    }

    int run() {
        HINSTANCE hInstance = GetModuleHandleA(nullptr);
        WNDCLASSEXA wc{};
        wc.cbSize = sizeof(WNDCLASSEXA);
        wc.style = CS_HREDRAW | CS_VREDRAW;
        wc.lpfnWndProc = WndProc;
        wc.hInstance = hInstance;
        wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
        wc.lpszClassName = "KajiyaPathTracerWindowClass";
        RegisterClassExA(&wc);

        RECT rect = {0, 0, win_w, win_h};
        AdjustWindowRect(&rect, WS_OVERLAPPEDWINDOW, FALSE);

        HWND hwnd = CreateWindowExA(
            0,
            wc.lpszClassName,
            "Motor Path Tracer C++ - Ecuacion Integral de Kajiya (Monte Carlo & Muestreo por Coseno)",
            WS_OVERLAPPEDWINDOW | WS_VISIBLE,
            CW_USEDEFAULT, CW_USEDEFAULT,
            rect.right - rect.left,
            rect.bottom - rect.top,
            nullptr, nullptr, hInstance, this
        );

        MSG msg{};
        bool running = true;
        while (running) {
            while (PeekMessageA(&msg, nullptr, 0, 0, PM_REMOVE)) {
                if (msg.message == WM_QUIT) {
                    running = false;
                    break;
                }
                TranslateMessage(&msg);
                DispatchMessageA(&msg);
            }
            if (!running) break;

            step_simulation();
            render_to_canvas();
            InvalidateRect(hwnd, nullptr, FALSE);
            UpdateWindow(hwnd);
        }
        return 0;
    }
};

#endif // _WIN32

} // namespace kajiya
