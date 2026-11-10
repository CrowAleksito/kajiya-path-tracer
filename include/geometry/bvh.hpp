#pragma once

#include "geometry/aabb.hpp"
#include "geometry/triangle.hpp"
#include <vector>
#include <algorithm>

namespace kajiya {

/**
 * @brief Nodo plano de la Jerarquía de Volúmenes Delimitadores (Bounding Volume Hierarchy - BVH).
 *        Permite acelerar la intersección de modelos 3D importados (.obj) de O(M) a O(log_2 M).
 */
struct BVHNode {
    AABB bounds;
    int left_child = -1;   // Índice del hijo izquierdo (-1 si es nodo hoja)
    int right_child = -1;  // Índice del hijo derecho
    int first_tri_idx = 0; // Índice inicial en el arreglo de triángulos ordenados
    int tri_count = 0;     // Número de triángulos en la hoja (> 0 si es hoja)

    constexpr bool is_leaf() const {
        return tri_count > 0;
    }
};

/**
 * @brief Árbol BVH construido desde cero sobre el conjunto de triángulos de la escena.
 */
class BVHAccel {
public:
    std::vector<BVHNode> nodes;
    std::vector<int> tri_indices;

    void build(const std::vector<Triangle>& triangles) {
        nodes.clear();
        tri_indices.resize(triangles.size());
        for (size_t i = 0; i < triangles.size(); ++i) {
            tri_indices[i] = static_cast<int>(i);
        }
        if (triangles.empty()) return;

        nodes.reserve(triangles.size() * 2);
        build_recursive(triangles, 0, static_cast<int>(triangles.size()));
    }

    /**
     * @brief Recorrido iterativo libre de recursión con pila estática para máxima velocidad en CPU.
     */
    bool intersect(
        const std::vector<Triangle>& triangles,
        const Ray& r,
        double t_min,
        double t_max,
        HitRecord& rec
    ) const {
        if (nodes.empty()) return false;

        Vec3 inv_dir(
            1.0 / (std::abs(r.direction.x) > 1e-15 ? r.direction.x : (r.direction.x >= 0 ? 1e-15 : -1e-15)),
            1.0 / (std::abs(r.direction.y) > 1e-15 ? r.direction.y : (r.direction.y >= 0 ? 1e-15 : -1e-15)),
            1.0 / (std::abs(r.direction.z) > 1e-15 ? r.direction.z : (r.direction.z >= 0 ? 1e-15 : -1e-15))
        );

        bool hit_anything = false;
        double closest_so_far = t_max;
        HitRecord temp_rec;

        int stack[64];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0; // Raíz del BVH

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = nodes[static_cast<size_t>(node_idx)];

            if (!node.bounds.intersect(r, inv_dir, t_min, closest_so_far)) {
                continue;
            }

            if (node.is_leaf()) {
                for (int i = 0; i < node.tri_count; ++i) {
                    int tri_id = tri_indices[static_cast<size_t>(node.first_tri_idx + i)];
                    if (triangles[static_cast<size_t>(tri_id)].intersect(r, t_min, closest_so_far, temp_rec)) {
                        hit_anything = true;
                        closest_so_far = temp_rec.t;
                        rec = temp_rec;
                    }
                }
            } else {
                if (node.right_child >= 0 && stack_ptr < 63) stack[stack_ptr++] = node.right_child;
                if (node.left_child >= 0 && stack_ptr < 63)  stack[stack_ptr++] = node.left_child;
            }
        }

        return hit_anything;
    }

private:
    int build_recursive(const std::vector<Triangle>& triangles, int start, int end) {
        int node_idx = static_cast<int>(nodes.size());
        nodes.emplace_back();

        AABB node_bounds;
        AABB centroid_bounds;
        for (int i = start; i < end; ++i) {
            AABB tb = AABB::from_triangle(triangles[static_cast<size_t>(tri_indices[static_cast<size_t>(i)])]);
            node_bounds.expand(tb);
            centroid_bounds.expand(tb.centroid());
        }

        int count = end - start;
        constexpr int MAX_TRIS_PER_LEAF = 4;

        if (count <= MAX_TRIS_PER_LEAF) {
            nodes[static_cast<size_t>(node_idx)].bounds = node_bounds;
            nodes[static_cast<size_t>(node_idx)].first_tri_idx = start;
            nodes[static_cast<size_t>(node_idx)].tri_count = count;
            return node_idx;
        }

        int axis = centroid_bounds.longest_axis();
        int mid = (start + end) / 2;

        // Partición por mediana de centroides en el eje más largo (std::nth_element O(N))
        std::nth_element(
            tri_indices.begin() + start,
            tri_indices.begin() + mid,
            tri_indices.begin() + end,
            [&](int a, int b) {
                Point3 ca = AABB::from_triangle(triangles[static_cast<size_t>(a)]).centroid();
                Point3 cb = AABB::from_triangle(triangles[static_cast<size_t>(b)]).centroid();
                return ca[axis] < cb[axis];
            }
        );

        int left_idx = build_recursive(triangles, start, mid);
        int right_idx = build_recursive(triangles, mid, end);

        nodes[static_cast<size_t>(node_idx)].bounds = node_bounds;
        nodes[static_cast<size_t>(node_idx)].left_child = left_idx;
        nodes[static_cast<size_t>(node_idx)].right_child = right_idx;
        nodes[static_cast<size_t>(node_idx)].tri_count = 0;
        return node_idx;
    }
};

} // namespace kajiya
