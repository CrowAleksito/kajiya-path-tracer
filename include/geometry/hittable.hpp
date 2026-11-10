#pragma once

#include "core/ray.hpp"
#include "core/material.hpp"

namespace kajiya {

/**
 * @brief Información geométrica y radiométrica del punto de impacto x.
 */
struct HitRecord {
    Point3 point;       // Punto de impacto x = r(t)
    Vec3 normal;        // Normal unitaria orientada hacia el hemisferio del rayo incidente
    double t = 0.0;     // Parámetro de distancia a lo largo del rayo
    double u = 0.0;     // Coordenada baricéntrica / superficial u
    double v = 0.0;     // Coordenada baricéntrica / superficial v
    bool front_face = true;
    int material_id = 0;

    inline void set_face_normal(const Ray& r, const Vec3& outward_normal) {
        front_face = dot(r.direction, outward_normal) < 0.0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

} // namespace kajiya
