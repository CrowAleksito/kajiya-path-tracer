#pragma once

#include "geometry/hittable.hpp"
#include <cmath>

namespace kajiya {

/**
 * @brief Primitiva geométrica de esfera definida analíticamente por:
 *        ||p - center||^2 = radius^2
 *
 * Sustituyendo la ecuación del rayo r(t) = o + t*d:
 *        (d . d) t^2 + 2 (d . (o - c)) t + (||o - c||^2 - R^2) = 0
 * Sea h = -(d . (c - o)) = d . (o - c). Entonces:
 *        a t^2 + 2 h t + c_term = 0
 * con discriminante reducido \Delta = h^2 - a * c_term.
 */
struct Sphere {
    Point3 center;
    double radius;
    int material_id;

    constexpr Sphere() : center(), radius(1.0), material_id(0) {}
    constexpr Sphere(const Point3& c, double r, int mat_id)
        : center(c), radius(r), material_id(mat_id) {}

    bool intersect(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
        Vec3 oc = r.origin - center;
        double a = r.direction.length_squared();
        double half_b = dot(oc, r.direction);
        double c = oc.length_squared() - radius * radius;

        double discriminant = half_b * half_b - a * c;
        if (discriminant < 0.0) {
            return false;
        }

        double sqrtd = std::sqrt(discriminant);

        // Encontrar la raíz más cercana dentro del intervalo aceptable [t_min, t_max]
        double root = (-half_b - sqrtd) / a;
        if (root <= t_min || root >= t_max) {
            root = (-half_b + sqrtd) / a;
            if (root <= t_min || root >= t_max) {
                return false;
            }
        }

        rec.t = root;
        rec.point = r.at(rec.t);
        Vec3 outward_normal = (rec.point - center) / radius;
        rec.set_face_normal(r, outward_normal);
        rec.material_id = material_id;
        return true;
    }
};

} // namespace kajiya
