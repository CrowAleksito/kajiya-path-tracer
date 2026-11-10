#pragma once

#include "geometry/hittable.hpp"
#include "geometry/sphere.hpp"
#include "geometry/triangle.hpp"
#include "geometry/bvh.hpp"
#include "io/obj_loader.hpp"
#include <vector>
#include <string>
#include <iostream>

namespace kajiya {

/**
 * @brief Contenedor de escena tridimensional que agrupa materiales, esferas,
 *        mallas de triángulos (evaluadas con Möller–Trumbore + aceleración BVH)
 *        y radiancia ambiental L_env.
 */
struct Scene {
    std::string name;
    std::vector<Material> materials;
    std::vector<Sphere> spheres;
    std::vector<Triangle> triangles;
    BVHAccel bvh;
    Color environment_radiance; // L_i(x, \omega_i) cuando el rayo escapa al infinito

    Scene() : name("Unnamed"), environment_radiance(0.0, 0.0, 0.0) {}

    int add_material(const Material& m) {
        materials.push_back(m);
        return static_cast<int>(materials.size() - 1);
    }

    /**
     * @brief Construye la Jerarquía de Volúmenes Delimitadores (BVH) sobre todos los triángulos
     *        de la escena para acelerar las intersecciones de O(M) a O(log_2 M).
     */
    void build_bvh() {
        bvh.build(triangles);
    }

    /**
     * @brief Importa un archivo Wavefront .obj dentro de la escena aplicando una matriz Mat4
     *        y reconstruye automáticamente el árbol BVH.
     */
    ObjLoadStats import_obj_model(
        const std::string& filepath,
        int material_id,
        const Mat4& transform = Mat4::identity(),
        bool normalize_to_unit_box = true
    ) {
        ObjLoadStats stats = ObjLoader::load_from_file(
            filepath, triangles, material_id, transform, normalize_to_unit_box
        );
        if (stats.success) {
            build_bvh();
        }
        return stats;
    }

    /**
     * @brief Calcula la intersección más cercana r(t) contra todas las esferas y triángulos (vía BVH).
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

        if (!bvh.nodes.empty()) {
            if (bvh.intersect(triangles, r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        } else {
            for (const auto& tri : triangles) {
                if (tri.intersect(r, t_min, closest_so_far, temp_rec)) {
                    hit_anything = true;
                    closest_so_far = temp_rec.t;
                    rec = temp_rec;
                }
            }
        }

        return hit_anything;
    }

    /**
     * @brief Fase 1: Escena de Validación Analítica de la Integral (Horno de Luz / White Furnace Test).
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
     * @brief Construye únicamente la habitación vacía de la Cornell Box (paredes + luz de área).
     */
    static void add_cornell_room_walls(Scene& scene, int mat_white, int mat_red, int mat_green, int mat_light) {
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

        // Fuente de luz de área en el techo (y = 1.998)
        append_quad(scene.triangles,
            Point3(-0.48, 1.998, -0.48),
            Point3( 0.48, 1.998, -0.48),
            Point3( 0.48, 1.998,  0.48),
            Point3(-0.48, 1.998,  0.48),
            mat_light);
    }

    /**
     * @brief Fase 2: Escena Estándar Cornell Box.
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

        add_cornell_room_walls(scene, mat_white, mat_red, mat_green, mat_light);

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

        // Esfera difusa sobre el bloque corto
        scene.spheres.emplace_back(Point3(0.35, 0.86, 0.25), 0.28, mat_cream);

        scene.build_bvh();
        return scene;
    }

    /**
     * @brief Crea una escena Cornell Box e importa un modelo 3D externo en formato Wavefront (.obj)
     *        colocado sobre un pedestal y acelerado mediante BVH.
     */
    static Scene create_cornell_with_obj(const std::string& obj_filepath, ObjLoadStats* out_stats = nullptr) {
        Scene scene;
        scene.name = "Cornell Box + Modelo OBJ (" + obj_filepath + ")";
        scene.environment_radiance = Color(0.0, 0.0, 0.0);

        int mat_white = scene.add_material(Material(Color(0.73, 0.73, 0.73)));
        int mat_red   = scene.add_material(Material(Color(0.65, 0.05, 0.05)));
        int mat_green = scene.add_material(Material(Color(0.12, 0.45, 0.15)));
        int mat_light = scene.add_material(Material(Color(0.78, 0.78, 0.78), Color(7.5, 7.5, 7.5)));
        int mat_model = scene.add_material(Material(Color(0.88, 0.80, 0.65))); // Tono marfil/dorado difuso

        add_cornell_room_walls(scene, mat_white, mat_red, mat_green, mat_light);

        // Pedestal central bajo (altura 0.32) donde descansa el modelo .obj importado
        Mat4 pedestal_transform =
            Mat4::translate(Vec3(0.0, 0.0, 0.05)) *
            Mat4::rotate_y(12.0 * PI / 180.0);
        append_transformed_box(
            scene.triangles,
            Point3(-0.42, 0.0, -0.42),
            Point3( 0.42, 0.32,  0.42),
            mat_white,
            pedestal_transform
        );

        // Transformación para posicionar el modelo .obj normalizado encima del pedestal (y = 0.32)
        // con tamaño 1.05 unidades y rotación estética de 22 grados
        Mat4 model_transform =
            Mat4::translate(Vec3(0.0, 0.32, 0.05)) *
            Mat4::rotate_y(22.0 * PI / 180.0) *
            Mat4::scale(Vec3(1.05, 1.05, 1.05));

        ObjLoadStats stats = scene.import_obj_model(obj_filepath, mat_model, model_transform, true);
        if (out_stats) {
            *out_stats = stats;
        }
        scene.build_bvh();
        return scene;
    }
};

} // namespace kajiya
