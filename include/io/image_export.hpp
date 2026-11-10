#pragma once

#include "integrator/path_tracer.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace kajiya {

/**
 * @brief Mapeo de tonos y corrección gamma (\gamma = 2.2) para visualización sRGB.
 */
inline uint8_t to_srgb_byte(double linear_val) {
    double clamped = std::clamp(linear_val, 0.0, 1.0);
    double gamma_corrected = std::pow(clamped, 1.0 / 2.2);
    return static_cast<uint8_t>(std::clamp( std::floor(gamma_corrected * 255.0 + 0.5), 0.0, 255.0 ));
}

/**
 * @brief Mapa de color perceptual (estilo Inferno) para visualizar la varianza por píxel.
 */
inline Color heatmap_color(double normalized_t) {
    double t = std::clamp(normalized_t, 0.0, 1.0);
    // Interpolación suave entre azul oscuro -> púrpura -> rojo/naranja -> amarillo brillante
    if (t < 0.25) {
        double u = t / 0.25;
        return Color(0.05 + 0.30 * u, 0.02 + 0.05 * u, 0.20 + 0.45 * u);
    } else if (t < 0.50) {
        double u = (t - 0.25) / 0.25;
        return Color(0.35 + 0.45 * u, 0.07 + 0.13 * u, 0.65 - 0.25 * u);
    } else if (t < 0.75) {
        double u = (t - 0.50) / 0.25;
        return Color(0.80 + 0.18 * u, 0.20 + 0.40 * u, 0.40 - 0.30 * u);
    } else {
        double u = (t - 0.75) / 0.25;
        return Color(0.98 + 0.02 * u, 0.60 + 0.38 * u, 0.10 + 0.65 * u);
    }
}

/**
 * @brief Exportador directo de matrices RGB a archivos PPM (sin compresión) y PNG (RFC 2083 desde cero).
 */
class ImageExporter {
public:
    /**
     * @brief Guarda una matriz de píxeles RGB (3 bytes por píxel) en formato PPM binario (P6).
     */
    static bool save_ppm_raw(const std::string& filename, int width, int height, const std::vector<uint8_t>& rgb_bytes) {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) return false;
        out << "P6\n" << width << " " << height << "\n255\n";
        out.write(reinterpret_cast<const char*>(rgb_bytes.data()), static_cast<std::streamsize>(rgb_bytes.size()));
        return true;
    }

    /**
     * @brief Guarda una matriz de píxeles RGB en formato PNG estándar codificado desde cero
     *        (bloques DEFLATE sin compresión + Adler-32 + CRC-32, cero librerías externas).
     */
    static bool save_png_raw(const std::string& filename, int width, int height, const std::vector<uint8_t>& rgb_bytes) {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) return false;

        auto crc32 = [](const uint8_t* data, size_t len) -> uint32_t {
            uint32_t c = 0xFFFFFFFFu;
            for (size_t i = 0; i < len; ++i) {
                c ^= data[i];
                for (int k = 0; k < 8; ++k) {
                    c = (c & 1u) ? (0xEDB88320u ^ (c >> 1u)) : (c >> 1u);
                }
            }
            return c ^ 0xFFFFFFFFu;
        };

        auto write_u32_be = [&](std::vector<uint8_t>& buf, uint32_t v) {
            buf.push_back(static_cast<uint8_t>((v >> 24) & 0xFF));
            buf.push_back(static_cast<uint8_t>((v >> 16) & 0xFF));
            buf.push_back(static_cast<uint8_t>((v >> 8) & 0xFF));
            buf.push_back(static_cast<uint8_t>(v & 0xFF));
        };

        auto write_chunk = [&](const char type[4], const std::vector<uint8_t>& data) {
            std::vector<uint8_t> chunk;
            write_u32_be(chunk, static_cast<uint32_t>(data.size()));
            chunk.push_back(static_cast<uint8_t>(type[0]));
            chunk.push_back(static_cast<uint8_t>(type[1]));
            chunk.push_back(static_cast<uint8_t>(type[2]));
            chunk.push_back(static_cast<uint8_t>(type[3]));
            chunk.insert(chunk.end(), data.begin(), data.end());
            uint32_t crc = crc32(chunk.data() + 4, 4 + data.size());
            write_u32_be(chunk, crc);
            out.write(reinterpret_cast<const char*>(chunk.data()), static_cast<std::streamsize>(chunk.size()));
        };

        // Firma PNG de 8 bytes
        const uint8_t png_sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
        out.write(reinterpret_cast<const char*>(png_sig), 8);

        // Chunk IHDR (13 bytes)
        std::vector<uint8_t> ihdr;
        write_u32_be(ihdr, static_cast<uint32_t>(width));
        write_u32_be(ihdr, static_cast<uint32_t>(height));
        ihdr.push_back(8); // Bit depth: 8 bits por canal
        ihdr.push_back(2); // Color type: 2 (Truecolor RGB)
        ihdr.push_back(0); // Compression method
        ihdr.push_back(0); // Filter method
        ihdr.push_back(0); // Interlace method (no interlace)
        write_chunk("IHDR", ihdr);

        // Construir scanlines con byte de filtro 0 al inicio de cada fila
        std::vector<uint8_t> raw_scanlines;
        raw_scanlines.reserve(static_cast<size_t>(height) * (static_cast<size_t>(width) * 3 + 1));
        for (int y = 0; y < height; ++y) {
            raw_scanlines.push_back(0); // Filtro None
            size_t row_offset = static_cast<size_t>(y) * static_cast<size_t>(width) * 3;
            raw_scanlines.insert(
                raw_scanlines.end(),
                rgb_bytes.begin() + static_cast<std::ptrdiff_t>(row_offset),
                rgb_bytes.begin() + static_cast<std::ptrdiff_t>(row_offset + static_cast<size_t>(width) * 3)
            );
        }

        // Flujo zlib con bloques DEFLATE tipo 0 (sin compresión) de hasta 65535 bytes
        std::vector<uint8_t> zlib_stream;
        zlib_stream.push_back(0x78); // CMF
        zlib_stream.push_back(0x01); // FLG

        size_t pos = 0;
        size_t total = raw_scanlines.size();
        uint32_t s1 = 1, s2 = 0;

        while (pos < total) {
            size_t block_len = std::min<size_t>(65535, total - pos);
            bool is_final = (pos + block_len == total);
            zlib_stream.push_back(is_final ? 0x01 : 0x00);
            uint16_t len16 = static_cast<uint16_t>(block_len);
            uint16_t nlen16 = ~len16;
            zlib_stream.push_back(static_cast<uint8_t>(len16 & 0xFF));
            zlib_stream.push_back(static_cast<uint8_t>((len16 >> 8) & 0xFF));
            zlib_stream.push_back(static_cast<uint8_t>(nlen16 & 0xFF));
            zlib_stream.push_back(static_cast<uint8_t>((nlen16 >> 8) & 0xFF));

            for (size_t i = 0; i < block_len; ++i) {
                uint8_t b = raw_scanlines[pos + i];
                zlib_stream.push_back(b);
                s1 = (s1 + b) % 65521u;
                s2 = (s2 + s1) % 65521u;
            }
            pos += block_len;
        }

        uint32_t adler = (s2 << 16) | s1;
        write_u32_be(zlib_stream, adler);

        write_chunk("IDAT", zlib_stream);
        write_chunk("IEND", {});
        return true;
    }

    /**
     * @brief Exporta la radiancia media acumulada <L_o> de un RenderBuffer a PPM y PNG.
     */
    static void export_radiance_image(
        const RenderBuffer& buffer,
        const std::string& base_path_without_ext
    ) {
        std::vector<uint8_t> rgb(static_cast<size_t>(buffer.width * buffer.height * 3));
        for (size_t i = 0; i < buffer.pixels.size(); ++i) {
            const Color& c = buffer.pixels[i].mean;
            rgb[i * 3 + 0] = to_srgb_byte(c.x);
            rgb[i * 3 + 1] = to_srgb_byte(c.y);
            rgb[i * 3 + 2] = to_srgb_byte(c.z);
        }
        save_ppm_raw(base_path_without_ext + ".ppm", buffer.width, buffer.height, rgb);
        save_png_raw(base_path_without_ext + ".png", buffer.width, buffer.height, rgb);
    }

    /**
     * @brief Exporta el mapa de calor de la varianza del estimador Var(<L_o>) por píxel.
     */
    static void export_variance_heatmap(
        const RenderBuffer& buffer,
        const std::string& base_path_without_ext,
        double max_variance_scale = 0.5
    ) {
        std::vector<uint8_t> rgb(static_cast<size_t>(buffer.width * buffer.height * 3));
        for (size_t i = 0; i < buffer.pixels.size(); ++i) {
            double var = buffer.pixels[i].estimator_variance();
            // Escala logarítmica/raíz para resaltar la estructura espacial del ruido
            double norm = std::sqrt(std::clamp(var / max_variance_scale, 0.0, 1.0));
            Color hc = heatmap_color(norm);
            rgb[i * 3 + 0] = static_cast<uint8_t>(std::clamp(hc.x * 255.0, 0.0, 255.0));
            rgb[i * 3 + 1] = static_cast<uint8_t>(std::clamp(hc.y * 255.0, 0.0, 255.0));
            rgb[i * 3 + 2] = static_cast<uint8_t>(std::clamp(hc.z * 255.0, 0.0, 255.0));
        }
        save_ppm_raw(base_path_without_ext + ".ppm", buffer.width, buffer.height, rgb);
        save_png_raw(base_path_without_ext + ".png", buffer.width, buffer.height, rgb);
    }
};

} // namespace kajiya
