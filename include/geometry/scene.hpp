#pragma once

#include "geometry/hittable.hpp"
#include "geometry/sphere.hpp"
#include "geometry/triangle.hpp"
#include <vector>
#include <string>

namespace kajiya {

/**
 * @brief Contenedor de escena tridimensional que agrupa materiales, esferas,
 *        mallas de triángulos (evaluadas con Möller–Trumbore) y radiancia ambiental L_env.
 */
struct Scene {
    std::string name;
    std::vector<Material> materials;
    std::vector<Sphere> spheres;
    std::vector<Triangle> triangles;
    Color environment_radiance; // L_i(x, \omega_i) cuando el rayo escapa al infinito

    Scene() : name("Unnamed"), environment_radiance(0.0, 0.0, 0.0) {}

    int add_material(const Material& m) {
        materials.push_back(m);
        return static_cast<int>(materials.size() - 1);
    }

    /**
     * @brief Calcula la intersección más cercana r(t) contra todas las esferas y triángulos.
     */
    bool intersect(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
        HitRecord temp_rec;
        bool hit_anything = false;
        double closest_so_far = t_max;

        for (const auto& s : spheres) {
            if (s.intersect(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }

        for (const auto& tri : triangles) {
            if (tri.intersect(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }

        return hit_anything;
    }

    /**
     * @brief Fase 1: Escena de Validación Analítica de la Integral (Horno de Luz / White Furnace Test).
     *
     * Una esfera unitaria con reflectancia difusa \rho = albedo (por defecto 1.0, 100% reflectora)
     * encerrada dentro de un entorno radiante isotrópico homogéneo con L_i = 1.0.
     *
     * En equilibrio termodinámico para \rho = 1:
     *   L_o = \int_{\Omega} \frac{1}{\pi} \cdot 1 \cdot \cos\theta \, d\omega = 1.000000
     */
    static Scene create_furnace_test(double rho = 1.0) {
        Scene scene;
        scene.name = "Fase 1: White Furnace Test (Horno de Luz)";
        scene.environment_radiance = Color(1.0, 1.0, 1.0);

        int mat_white = scene.add_material(Material(Color(rho, rho, rho), Color(0.0, 0.0, 0.0)));
        scene.spheres.emplace_back(Point3(0.0, 0.0, 0.0), 1.0, mat_white);
        return scene;
    }

    /**
     * @brief Fase 2: Escena Estándar Cornell Box.
     *
     * Habitación cerrada en [-1, 1]^3 con:
     * - Pared izquierda roja (\rho = (0.65, 0.05, 0.05))
     * - Pared derecha verde (\rho = (0.12, 0.45, 0.15))
     * - Piso, techo y pared frontal/trasera blancos (\rho = (0.73, 0.73, 0.73))
     * - Fuente de luz de área en el techo con radiancia emitida L_e = (15.0, 15.0, 15.0)
     * - Bloque alto rotado con matriz 4x4 (evaluado con Möller–Trumbore)
     * - Bloque corto rotado con matriz 4x4 + esfera difusa sobre la escena
     */
    static Scene create_cornell_box() {
        Scene scene;
        scene.name = "Fase 2: Cornell Box Estandar";
        scene.environment_radiance = Color(0.0, 0.0, 0.0);

        int mat_white = scene.add_material(Material(Color(0.73, 0.73, 0.73)));
        int mat_red   = scene.add_material(Material(Color(0.65, 0.05, 0.05)));
        int mat_green = scene.add_material(Material(Color(0.12, 0.45, 0.15)));
        int mat_light = scene.add_material(Material(Color(0.78, 0.78, 0.78), Color(7.5, 7.5, 7.5)));
        int mat_cream = scene.add_material(Material(Color(0.85, 0.82, 0.75)));

        // Coordenadas de la habitación: X en [-1, 1], Y en [0, 2], Z en [-1, 1]
        // Piso (y = 0)
        append_quad(scene.triangles,
            Point3(-1.0, 0.0,  1.0),
            Point3( 1.0, 0.0,  1.0),
            Point3( 1.0, 0.0, -1.0),
            Point3(-1.0, 0.0, -1.0),
            mat_white);

        // Techo (y = 2)
        append_quad(scene.triangles,
            Point3(-1.0, 2.0, -1.0),
            Point3( 1.0, 2.0, -1.0),
            Point3( 1.0, 2.0,  1.0),
            Point3(-1.0, 2.0,  1.0),
            mat_white);

        // Pared trasera (z = -1)
        append_quad(scene.triangles,
            Point3(-1.0, 0.0, -1.0),
            Point3( 1.0, 0.0, -1.0),
            Point3( 1.0, 2.0, -1.0),
            Point3(-1.0, 2.0, -1.0),
            mat_white);

        // Pared izquierda ROJA (x = -1)
        append_quad(scene.triangles,
            Point3(-1.0, 0.0,  1.0),
            Point3(-1.0, 0.0, -1.0),
            Point3(-1.0, 2.0, -1.0),
            Point3(-1.0, 2.0,  1.0),
            mat_red);

        // Pared derecha VERDE (x = +1)
        append_quad(scene.triangles,
            Point3( 1.0, 0.0, -1.0),
            Point3( 1.0, 0.0,  1.0),
            Point3( 1.0, 2.0,  1.0),
            Point3( 1.0, 2.0, -1.0),
            mat_green);

        // Fuente de luz de área en el techo (y = 1.998 para evitar z-fighting con el techo)
        append_quad(scene.triangles,
            Point3(-0.48, 1.998, -0.48),
            Point3( 0.48, 1.998, -0.48),
            Point3( 0.48, 1.998,  0.48),
            Point3(-0.48, 1.998,  0.48),
            mat_light);

        // Bloque alto izquierdo transformado con matriz homogénea 4x4 (rotación en Y + traslación)
        Mat4 tall_box_transform =
            Mat4::translate(Vec3(-0.33, 0.0, -0.28)) *
            Mat4::rotate_y(18.0 * PI / 180.0);
        append_transformed_box(
            scene.triangles,
            Point3(-0.28, 0.0, -0.28),
            Point3( 0.28, 1.18,  0.28),
            mat_white,
            tall_box_transform
        );

        // Bloque corto derecho transformado con matriz homogénea 4x4
        Mat4 short_box_transform =
            Mat4::translate(Vec3(0.35, 0.0, 0.25)) *
            Mat4::rotate_y(-16.0 * PI / 180.0);
        append_transformed_box(
            scene.triangles,
            Point3(-0.28, 0.0, -0.28),
            Point3( 0.28, 0.58,  0.28),
            mat_white,
            short_box_transform
        );

        // Esfera difusa sobre el bloque corto para ejercitar simultáneamente Rayo-Esfera y Möller–Trumbore
        scene.spheres.emplace_back(Point3(0.35, 0.86, 0.25), 0.28, mat_cream);

        return scene;
    }
};

} // namespace kajiya
