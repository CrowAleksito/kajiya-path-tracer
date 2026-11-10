#pragma once

#include "core/ray.hpp"
#include "geometry/triangle.hpp"
#include <algorithm>
#include <cmath>

namespace kajiya {

/**
 * @brief Caja Delimitadora Alineada a los Ejes (Axis-Aligned Bounding Box - AABB)
 *        implementada desde cero mediante el método de losas (Slab Method de Kay & Kajiya).
 *
 * Para cada eje i \in {x, y, z}, el intervalo de intersección del rayo r(t) = o + t*d
 * con los planos x_i = min_i y x_i = max_i es:
 *   t_{0,i} = (min_i - o_i) / d_i,   t_{1,i} = (max_i - o_i) / d_i
 * El rayo intersecta la caja si y solo si max_i(min(t_{0,i}, t_{1,i})) < min_i(max(t_{0,i}, t_{1,i})).
 */
struct AABB {
    Point3 min_pt;
    Point3 max_pt;

    constexpr AABB()
        : min_pt( 1e30,  1e30,  1e30),
          max_pt(-1e30, -1e30, -1e30) {}

    constexpr AABB(const Point3& a, const Point3& b)
        : min_pt(a), max_pt(b) {}

    void expand(const Point3& p) {
        min_pt.x = std::min(min_pt.x, p.x);
        min_pt.y = std::min(min_pt.y, p.y);
        min_pt.z = std::min(min_pt.z, p.z);

        max_pt.x = std::max(max_pt.x, p.x);
        max_pt.y = std::max(max_pt.y, p.y);
        max_pt.z = std::max(max_pt.z, p.z);
    }

    void expand(const AABB& box) {
        expand(box.min_pt);
        expand(box.max_pt);
    }

    Point3 centroid() const {
        return (min_pt + max_pt) * 0.5;
    }

    int longest_axis() const {
        Vec3 ext = max_pt - min_pt;
        if (ext.x >= ext.y && ext.x >= ext.z) return 0;
        if (ext.y >= ext.z) return 1;
        return 2;
    }

    static AABB from_triangle(const Triangle& tri) {
        AABB box;
        box.expand(tri.v0);
        box.expand(tri.v1);
        box.expand(tri.v2);
        // Margen numérico para triángulos coplanares a un eje
        constexpr double PAD = 1e-5;
        box.min_pt -= Vec3(PAD, PAD, PAD);
        box.max_pt += Vec3(PAD, PAD, PAD);
        return box;
    }

    /**
     * @brief Prueba rápida de intersección Rayo-AABB (Slab test).
     */
    inline bool intersect(const Ray& r, const Vec3& inv_dir, double t_min, double t_max) const {
        for (int a = 0; a < 3; ++a) {
            double t0 = (min_pt[a] - r.origin[a]) * inv_dir[a];
            double t1 = (max_pt[a] - r.origin[a]) * inv_dir[a];
            if (inv_dir[a] < 0.0) std::swap(t0, t1);
            t_min = t0 > t_min ? t0 : t_min;
            t_max = t1 < t_max ? t1 : t_max;
            if (t_max < t_min) return false;
        }
        return true;
    }
};

} // namespace kajiya
