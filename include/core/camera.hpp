#pragma once

#include "core/ray.hpp"
#include <cmath>

namespace kajiya {

/**
 * @brief Cámara estenopeica (Pinhole Perspective Camera) configurable.
 */
struct Camera {
    Point3 origin;
    Point3 lower_left_corner;
    Vec3 horizontal;
    Vec3 vertical;
    Vec3 u, v, w;

    Camera(
        const Point3& lookfrom,
        const Point3& lookat,
        const Vec3& vup,
        double vfov_degrees,
        double aspect_ratio
    ) {
        double theta = vfov_degrees * PI / 180.0;
        double h = std::tan(theta * 0.5);
        double viewport_height = 2.0 * h;
        double viewport_width = aspect_ratio * viewport_height;

        w = normalize(lookfrom - lookat);
        u = normalize(cross(vup, w));
        v = cross(w, u);

        origin = lookfrom;
        horizontal = viewport_width * u;
        vertical = viewport_height * v;
        lower_left_corner = origin - horizontal * 0.5 - vertical * 0.5 - w;
    }

    /**
     * @brief Genera un rayo primario a través de las coordenadas normalizadas del plano de imagen (s, t) en [0, 1]^2.
     */
    Ray get_ray(double s, double t) const {
        Vec3 dir = lower_left_corner + s * horizontal + t * vertical - origin;
        return Ray(origin, normalize(dir));
    }

    static Camera cornell_camera(int width, int height) {
        double aspect = static_cast<double>(width) / static_cast<double>(height);
        return Camera(
            Point3(0.0, 1.0, 3.65),
            Point3(0.0, 1.0, 0.0),
            Vec3(0.0, 1.0, 0.0),
            39.5,
            aspect
        );
    }

    static Camera furnace_camera(int width, int height) {
        double aspect = static_cast<double>(width) / static_cast<double>(height);
        return Camera(
            Point3(0.0, 0.0, 3.2),
            Point3(0.0, 0.0, 0.0),
            Vec3(0.0, 1.0, 0.0),
            45.0,
            aspect
        );
    }
};

} // namespace kajiya
