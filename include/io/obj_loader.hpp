#pragma once

#include "geometry/triangle.hpp"
#include "geometry/aabb.hpp"
#include "math/mat4.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

namespace kajiya {

struct ObjLoadStats {
    bool success = false;
    size_t vertex_count = 0;
    size_t triangle_count = 0;
    AABB raw_bounds;
    std::string error_message;
};

/**
 * @brief Importador de mallas tridimensionales Wavefront (.obj) implementado desde cero.
 *
 * Características:
 * - Lee posiciones de vértices `v x y z`.
 * - Lee caras poligonales `f v1 v2 v3 ...` en cualquiera de las variantes del estándar OBJ:
 *     `f v`, `f v/vt`, `f v//vn`, `f v/vt/vn`, incluyendo índices negativos (relativos).
 * - Triangulación automática en abanico (Fan Triangulation) para cuadriláteros y polígonos de N lados.
 * - Normalización opcional a caja unitaria (centra el modelo en X/Z, apoya la base en Y=0
 *   y escala la dimensión máxima a 1.0) antes de aplicar la matriz `Mat4` de transformación.
 */
class ObjLoader {
public:
    static ObjLoadStats load_from_file(
        const std::string& filepath,
        std::vector<Triangle>& out_triangles,
        int material_id,
        const Mat4& transform = Mat4::identity(),
        bool normalize_to_unit_box = true
    ) {
        ObjLoadStats stats;
        std::ifstream file(filepath);
        if (!file.is_open()) {
            stats.error_message = "No se pudo abrir el archivo OBJ: " + filepath;
            return stats;
        }

        std::vector<Point3> raw_vertices;
        std::vector<std::vector<int>> raw_faces;

        std::string line;
        while (std::getline(file, line)) {
            // Eliminar comentarios '#'
            size_t comment_pos = line.find('#');
            if (comment_pos != std::string::npos) {
                line = line.substr(0, comment_pos);
            }

            std::istringstream iss(line);
            std::string prefix;
            if (!(iss >> prefix)) continue;

            if (prefix == "v") {
                double x = 0.0, y = 0.0, z = 0.0;
                if (iss >> x >> y >> z) {
                    Point3 p(x, y, z);
                    raw_vertices.push_back(p);
                    stats.raw_bounds.expand(p);
                }
            } else if (prefix == "f") {
                std::vector<int> face_indices;
                std::string token;
                while (iss >> token) {
                    int v_idx = parse_vertex_index(token, static_cast<int>(raw_vertices.size()));
                    if (v_idx >= 0 && v_idx < static_cast<int>(raw_vertices.size())) {
                        face_indices.push_back(v_idx);
                    }
                }
                if (face_indices.size() >= 3) {
                    raw_faces.push_back(std::move(face_indices));
                }
            }
        }

        if (raw_vertices.empty() || raw_faces.empty()) {
            stats.error_message = "El archivo OBJ no contiene vertices o caras validas.";
            return stats;
        }

        // Construir matriz compuesta: Transformación del usuario * Normalización a caja unitaria
        Mat4 final_transform = transform;
        if (normalize_to_unit_box) {
            Vec3 ext = stats.raw_bounds.max_pt - stats.raw_bounds.min_pt;
            double max_dim = std::max({ext.x, ext.y, ext.z, 1e-9});
            double inv_scale = 1.0 / max_dim;

            double center_x = 0.5 * (stats.raw_bounds.min_pt.x + stats.raw_bounds.max_pt.x);
            double min_y    = stats.raw_bounds.min_pt.y; // Apoyar la base del modelo en y = 0
            double center_z = 0.5 * (stats.raw_bounds.min_pt.z + stats.raw_bounds.max_pt.z);

            Mat4 norm_mat =
                Mat4::scale(Vec3(inv_scale, inv_scale, inv_scale)) *
                Mat4::translate(Vec3(-center_x, -min_y, -center_z));

            final_transform = transform * norm_mat;
        }

        // Transformar todos los vértices una sola vez
        std::vector<Point3> transformed_vertices(raw_vertices.size());
        for (size_t i = 0; i < raw_vertices.size(); ++i) {
            transformed_vertices[i] = final_transform.transform_point(raw_vertices[i]);
        }

        // Triangulación en abanico (Fan Triangulation) para caras de 3, 4 o N vértices
        size_t initial_tri_count = out_triangles.size();
        for (const auto& face : raw_faces) {
            const Point3& v0 = transformed_vertices[static_cast<size_t>(face[0])];
            for (size_t i = 1; i + 1 < face.size(); ++i) {
                const Point3& v1 = transformed_vertices[static_cast<size_t>(face[i])];
                const Point3& v2 = transformed_vertices[static_cast<size_t>(face[i + 1])];
                // Descartar triángulos degenerados de área cero
                if (cross(v1 - v0, v2 - v0).length_squared() > 1e-18) {
                    out_triangles.emplace_back(v0, v1, v2, material_id);
                }
            }
        }

        stats.vertex_count = raw_vertices.size();
        stats.triangle_count = out_triangles.size() - initial_tri_count;
        stats.success = (stats.triangle_count > 0);
        return stats;
    }

private:
    static int parse_vertex_index(const std::string& token, int total_vertices) {
        size_t slash_pos = token.find('/');
        std::string idx_str = (slash_pos == std::string::npos) ? token : token.substr(0, slash_pos);
        if (idx_str.empty()) return -1;
        int idx = std::stoi(idx_str);
        if (idx > 0) {
            return idx - 1; // Índices OBJ son base 1
        } else if (idx < 0) {
            return total_vertices + idx; // Índices negativos son relativos al final
        }
        return -1;
    }
};

} // namespace kajiya
