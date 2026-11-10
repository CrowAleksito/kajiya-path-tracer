#pragma once

#include "math/vec3.hpp"
#include <cmath>
#include <cstring>

namespace kajiya {

/**
 * @brief Matriz homogénea 4x4 implementada desde cero para transformaciones afines
 *        (traslación, rotación, escalado) de primitivas geométricas, puntos, vectores y normales.
 */
struct Mat4 {
    double m[4][4];

    Mat4() {
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                m[r][c] = (r == c) ? 1.0 : 0.0;
    }

    static Mat4 identity() {
        return Mat4();
    }

    static Mat4 translate(const Vec3& t) {
        Mat4 res;
        res.m[0][3] = t.x;
        res.m[1][3] = t.y;
        res.m[2][3] = t.z;
        return res;
    }

    static Mat4 scale(const Vec3& s) {
        Mat4 res;
        res.m[0][0] = s.x;
        res.m[1][1] = s.y;
        res.m[2][2] = s.z;
        return res;
    }

    static Mat4 rotate_x(double radians) {
        Mat4 res;
        double c = std::cos(radians);
        double s = std::sin(radians);
        res.m[1][1] = c;  res.m[1][2] = -s;
        res.m[2][1] = s;  res.m[2][2] = c;
        return res;
    }

    static Mat4 rotate_y(double radians) {
        Mat4 res;
        double c = std::cos(radians);
        double s = std::sin(radians);
        res.m[0][0] = c;  res.m[0][2] = s;
        res.m[2][0] = -s; res.m[2][2] = c;
        return res;
    }

    static Mat4 rotate_z(double radians) {
        Mat4 res;
        double c = std::cos(radians);
        double s = std::sin(radians);
        res.m[0][0] = c;  res.m[0][1] = -s;
        res.m[1][0] = s;  res.m[1][1] = c;
        return res;
    }

    Mat4 operator*(const Mat4& rhs) const {
        Mat4 res;
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                res.m[r][c] =
                    m[r][0] * rhs.m[0][c] +
                    m[r][1] * rhs.m[1][c] +
                    m[r][2] * rhs.m[2][c] +
                    m[r][3] * rhs.m[3][c];
            }
        }
        return res;
    }

    /**
     * @brief Transforma un punto 3D en coordenadas homogéneas (w = 1).
     */
    Point3 transform_point(const Point3& p) const {
        double xp = m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3];
        double yp = m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3];
        double zp = m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3];
        double wp = m[3][0] * p.x + m[3][1] * p.y + m[3][2] * p.z + m[3][3];
        if (std::abs(wp - 1.0) > 1e-12 && std::abs(wp) > 1e-15) {
            return Point3(xp / wp, yp / wp, zp / wp);
        }
        return Point3(xp, yp, zp);
    }

    /**
     * @brief Transforma un vector direccional 3D en coordenadas homogéneas (w = 0).
     */
    Vec3 transform_vector(const Vec3& v) const {
        return Vec3(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
        );
    }

    /**
     * @brief Transpuesta de la submatriz 3x3 superior izquierda.
     */
    Mat4 transpose_3x3() const {
        Mat4 res;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                res.m[r][c] = m[c][r];
        return res;
    }

    /**
     * @brief Inversa de la submatriz 3x3 superior izquierda para transformación exacta de normales
     *        mediante (M^{-1})^T.
     */
    Mat4 inverse_transpose_3x3() const {
        double a = m[0][0], b = m[0][1], c = m[0][2];
        double d = m[1][0], e = m[1][1], f = m[1][2];
        double g = m[2][0], h = m[2][1], i = m[2][2];

        double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
        Mat4 res;
        if (std::abs(det) < 1e-15) return res;
        double inv_det = 1.0 / det;

        // Cofactores (que ya equivalen a la transpuesta de la adjunta, es decir (M^{-1})^T)
        res.m[0][0] =  (e * i - f * h) * inv_det;
        res.m[0][1] = -(d * i - f * g) * inv_det;
        res.m[0][2] =  (d * h - e * g) * inv_det;

        res.m[1][0] = -(b * i - c * h) * inv_det;
        res.m[1][1] =  (a * i - c * g) * inv_det;
        res.m[1][2] = -(a * h - b * g) * inv_det;

        res.m[2][0] =  (b * f - c * e) * inv_det;
        res.m[2][1] = -(a * f - c * d) * inv_det;
        res.m[2][2] =  (a * e - b * d) * inv_det;
        return res;
    }

    Vec3 transform_normal(const Vec3& n) const {
        Mat4 inv_t = inverse_transpose_3x3();
        return normalize(inv_t.transform_vector(n));
    }
};

} // namespace kajiya
