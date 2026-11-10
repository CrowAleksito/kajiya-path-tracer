#pragma once

#include "geometry/hittable.hpp"
#include "math/mat4.hpp"
#include <cmath>
#include <vector>

namespace kajiya {

/**
 * @brief Primitiva de triángulo con intersección analítica mediante el
 *        algoritmo de Möller–Trumbore (1997).
 *
 * Ecuación baricéntrica de un punto en el triángulo (v0, v1, v2):
 *        T(u, v) = (1 - u - v) v0 + u v1 + v v2 = v0 + u E1 + v E2
 * donde E1 = v1 - v0 y E2 = v2 - v0, con u >= 0, v >= 0, u + v <= 1.
 *
 * Igualando r(t) = o + t d = v0 + u E1 + v E2:
 *        [-d, E1, E2] [t, u, v]^T = o - v0 = T_vec
 *
 * Resolviendo mediante la Regla de Cramer y el producto triple escalar:
 *        P = d x E2
 *        det = E1 . P
 *        u = (T_vec . P) / det
 *        Q = T_vec x E1
 *        v = (d . Q) / det
 *        t = (E2 . Q) / det
 */
struct Triangle {
    Point3 v0, v1, v2;
    Vec3 edge1, edge2;
    Vec3 geometric_normal;
    int material_id;
    bool single_sided;

    Triangle() : v0(), v1(), v2(), edge1(), edge2(), geometric_normal(0, 1, 0), material_id(0), single_sided(false) {}

    Triangle(const Point3& a, const Point3& b, const Point3& c, int mat_id, bool one_sided = false)
        : v0(a), v1(b), v2(c), material_id(mat_id), single_sided(one_sided) {
        edge1 = v1 - v0;
        edge2 = v2 - v0;
        geometric_normal = normalize(cross(edge1, edge2));
    }

    bool intersect(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
        constexpr double MT_EPS = 1e-9;

        Vec3 pvec = cross(r.direction, edge2);
        double det = dot(edge1, pvec);

        if (single_sided) {
            if (det < MT_EPS) return false;
        } else {
            if (std::abs(det) < MT_EPS) return false;
        }

        double inv_det = 1.0 / det;
        Vec3 tvec = r.origin - v0;

        // Coordenada baricéntrica u
        double u = dot(tvec, pvec) * inv_det;
        if (u < 0.0 || u > 1.0) return false;

        // Coordenada baricéntrica v
        Vec3 qvec = cross(tvec, edge1);
        double v = dot(r.direction, qvec) * inv_det;
        if (v < 0.0 || (u + v) > 1.0) return false;

        // Parámetro t a lo largo del rayo
        double t = dot(edge2, qvec) * inv_det;
        if (t <= t_min || t >= t_max) return false;

        rec.t = t;
        rec.u = u;
        rec.v = v;
        rec.point = r.at(t);
        rec.set_face_normal(r, geometric_normal);
        rec.material_id = material_id;
        return true;
    }
};

/**
 * @brief Añade un cuadrilátero (dividido en 2 triángulos orientados en sentido antihorario)
 *        transformado opcionalmente por una matriz 4x4.
 */
inline void append_quad(
    std::vector<Triangle>& tris,
    const Point3& p0, const Point3& p1, const Point3& p2, const Point3& p3,
    int mat_id,
    const Mat4& transform = Mat4::identity()
) {
    Point3 a = transform.transform_point(p0);
    Point3 b = transform.transform_point(p1);
    Point3 c = transform.transform_point(p2);
    Point3 d = transform.transform_point(p3);
    tris.emplace_back(a, b, c, mat_id);
    tris.emplace_back(a, c, d, mat_id);
}

/**
 * @brief Construye un paralelepípedo (caja 3D de 12 triángulos) entre p_min y p_max
 *        y aplica una matriz de transformación homogénea 4x4 (traslación/rotación/escalado).
 */
inline void append_transformed_box(
    std::vector<Triangle>& tris,
    const Point3& p_min,
    const Point3& p_max,
    int mat_id,
    const Mat4& transform
) {
    Point3 p000(p_min.x, p_min.y, p_min.z);
    Point3 p100(p_max.x, p_min.y, p_min.z);
    Point3 p101(p_max.x, p_min.y, p_max.z);
    Point3 p001(p_min.x, p_min.y, p_max.z);

    Point3 p010(p_min.x, p_max.y, p_min.z);
    Point3 p110(p_max.x, p_max.y, p_min.z);
    Point3 p111(p_max.x, p_max.y, p_max.z);
    Point3 p011(p_min.x, p_max.y, p_max.z);

    // Cara frontal (+Z)
    append_quad(tris, p001, p101, p111, p011, mat_id, transform);
    // Cara trasera (-Z)
    append_quad(tris, p100, p000, p010, p110, mat_id, transform);
    // Cara derecha (+X)
    append_quad(tris, p101, p100, p110, p111, mat_id, transform);
    // Cara izquierda (-X)
    append_quad(tris, p000, p001, p011, p010, mat_id, transform);
    // Cara superior (+Y)
    append_quad(tris, p011, p111, p110, p010, mat_id, transform);
    // Cara inferior (-Y)
    append_quad(tris, p000, p100, p101, p001, mat_id, transform);
}

} // namespace kajiya
