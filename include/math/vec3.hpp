#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>

namespace kajiya {

constexpr double PI = 3.14159265358979323846;
constexpr double INV_PI = 0.31830988618379067154;
constexpr double TWO_PI = 6.28318530717958647692;
constexpr double EPSILON = 1e-6;

/**
 * @brief Estructura de vector tridimensional escrita desde cero para representar
 *        puntos 3D, direcciones (vectores unitarios) y espectros de radiancia RGB.
 */
struct Vec3 {
    double x, y, z;

    constexpr Vec3() : x(0.0), y(0.0), z(0.0) {}
    constexpr explicit Vec3(double v) : x(v), y(v), z(v) {}
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator-() const { return Vec3(-x, -y, -z); }

    constexpr double operator[](int i) const {
        return (i == 0) ? x : ((i == 1) ? y : z);
    }

    double& operator[](int i) {
        return (i == 0) ? x : ((i == 1) ? y : z);
    }

    Vec3& operator+=(const Vec3& v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& v) {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    Vec3& operator*=(double t) {
        x *= t; y *= t; z *= t;
        return *this;
    }

    Vec3& operator*=(const Vec3& v) {
        x *= v.x; y *= v.y; z *= v.z;
        return *this;
    }

    Vec3& operator/=(double t) {
        double inv = 1.0 / t;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }

    constexpr double length_squared() const {
        return x * x + y * y + z * z;
    }

    double length() const {
        return std::sqrt(length_squared());
    }

    /**
     * @brief Luminancia relativa CIE (Rec. 709) para análisis escalar de varianza y radiancia.
     */
    constexpr double luminance() const {
        return 0.2126 * x + 0.7152 * y + 0.0722 * z;
    }

    constexpr double max_component() const {
        return (x > y) ? ((x > z) ? x : z) : ((y > z) ? y : z);
    }
};

using Point3 = Vec3;
using Color = Vec3;

inline constexpr Vec3 operator+(const Vec3& u, const Vec3& v) {
    return Vec3(u.x + v.x, u.y + v.y, u.z + v.z);
}

inline constexpr Vec3 operator-(const Vec3& u, const Vec3& v) {
    return Vec3(u.x - v.x, u.y - v.y, u.z - v.z);
}

inline constexpr Vec3 operator*(const Vec3& u, const Vec3& v) {
    return Vec3(u.x * v.x, u.y * v.y, u.z * v.z);
}

inline constexpr Vec3 operator*(double t, const Vec3& v) {
    return Vec3(t * v.x, t * v.y, t * v.z);
}

inline constexpr Vec3 operator*(const Vec3& v, double t) {
    return t * v;
}

inline constexpr Vec3 operator/(const Vec3& v, double t) {
    return (1.0 / t) * v;
}

/**
 * @brief Producto punto escalar en R^3: u . v = u_x v_x + u_y v_y + u_z v_z
 */
inline constexpr double dot(const Vec3& u, const Vec3& v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

/**
 * @brief Producto cruz vectorial en R^3: u x v
 */
inline constexpr Vec3 cross(const Vec3& u, const Vec3& v) {
    return Vec3(
        u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x
    );
}

/**
 * @brief Normalización de un vector: v / ||v||
 */
inline Vec3 normalize(const Vec3& v) {
    double len = v.length();
    if (len < 1e-15) return Vec3(0.0, 1.0, 0.0);
    return v / len;
}

/**
 * @brief Base Ortonormal Local (ONB) construida alrededor de la normal de superficie n.
 *        Permite transformar direcciones muestreadas en el hemisferio canónico (z >= 0)
 *        al espacio de coordenadas global orientado según la normal n del punto de impacto x.
 */
struct OrthonormalBasis {
    Vec3 tangent;
    Vec3 bitangent;
    Vec3 normal;

    explicit OrthonormalBasis(const Vec3& n) {
        normal = normalize(n);
        // Construcción robusta de Duff et al. (2017) / Frisvad sin singularidades
        double sign = (normal.z >= 0.0) ? 1.0 : -1.0;
        double a = -1.0 / (sign + normal.z);
        double b = normal.x * normal.y * a;
        tangent = Vec3(1.0 + sign * normal.x * normal.x * a, sign * b, -sign * normal.x);
        bitangent = Vec3(b, sign + normal.y * normal.y * a, -normal.y);
    }

    /**
     * @brief Transforma un vector local (x_l, y_l, z_l) del hemisferio al espacio global.
     */
    constexpr Vec3 to_world(const Vec3& local_dir) const {
        return local_dir.x * tangent + local_dir.y * bitangent + local_dir.z * normal;
    }
};

} // namespace kajiya
