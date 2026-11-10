#pragma once

#include "math/vec3.hpp"

namespace kajiya {

/**
 * @brief Material físicamente basado con BRDF Lambertiana ideal y emisión propia.
 *
 * Para una superficie difusa Lambertiana con reflectancia hemisférica (albedo) \rho \in [0, 1]^3:
 *   f_r(x, \omega_i, \omega_o) = \frac{\rho}{\pi}
 *
 * Cumple conservación de energía:
 *   \int_{\Omega} f_r(x, \omega_i, \omega_o) (\omega_i \cdot n) d\omega_i
 *   = \frac{\rho}{\pi} \int_0^{2\pi} \int_0^{\pi/2} \cos\theta \sin\theta d\theta d\phi
 *   = \frac{\rho}{\pi} \cdot \pi = \rho \le 1.
 */
struct Material {
    Color albedo;     // Reflectancia difusa \rho \in [0, 1]^3
    Color emission;   // Radiancia emitida L_e(x, \omega_o)

    constexpr Material() : albedo(0.8, 0.8, 0.8), emission(0.0, 0.0, 0.0) {}
    constexpr Material(const Color& rho, const Color& le = Color(0.0, 0.0, 0.0))
        : albedo(rho), emission(le) {}

    /**
     * @brief Evalúa la Función de Distribución de Reflectancia Bidireccional (BRDF) f_r.
     */
    constexpr Color evaluate_brdf() const {
        return albedo * INV_PI;
    }

    constexpr bool is_emissive() const {
        return emission.x > 0.0 || emission.y > 0.0 || emission.z > 0.0;
    }
};

} // namespace kajiya
