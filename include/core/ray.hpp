#pragma once

#include "math/vec3.hpp"

namespace kajiya {

/**
 * @brief Representación paramétrica de un rayo en R^3:
 *        r(t) = origin + t * direction, con t >= 0.
 */
struct Ray {
    Point3 origin;
    Vec3 direction;

    constexpr Ray() : origin(), direction(0.0, 0.0, -1.0) {}
    constexpr Ray(const Point3& o, const Vec3& d) : origin(o), direction(d) {}

    constexpr Point3 at(double t) const {
        return origin + t * direction;
    }
};

} // namespace kajiya
