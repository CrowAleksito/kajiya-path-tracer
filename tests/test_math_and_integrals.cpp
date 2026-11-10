#include "math/vec3.hpp"
#include "math/mat4.hpp"
#include "geometry/sphere.hpp"
#include "geometry/triangle.hpp"
#include "integrator/sampling.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace kajiya;

void test_vector_and_onb() {
    Vec3 a(1.0, 2.0, 3.0);
    Vec3 b(4.0, -5.0, 6.0);
    // 1*4 + 2*(-5) + 3*6 = 4 - 10 + 18 = 12
    assert(std::abs(dot(a, b) - 12.0) < 1e-12);

    Vec3 c = cross(Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
    assert(std::abs(c.x) < 1e-12 && std::abs(c.y) < 1e-12 && std::abs(c.z - 1.0) < 1e-12);

    OrthonormalBasis onb(normalize(Vec3(0.3, -0.7, 0.5)));
    assert(std::abs(dot(onb.tangent, onb.normal)) < 1e-12);
    assert(std::abs(dot(onb.bitangent, onb.normal)) < 1e-12);
    assert(std::abs(dot(onb.tangent, onb.bitangent)) < 1e-12);
    std::cout << "[PASS] Algebra vectorial 3D y Base Ortonormal (ONB)\n";
}

void test_ray_sphere_and_moller_trumbore() {
    Sphere s(Point3(0.0, 0.0, -5.0), 1.0, 0);
    Ray r(Point3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, -1.0));
    HitRecord rec;
    bool hit = s.intersect(r, 1e-5, 100.0, rec);
    assert(hit);
    assert(std::abs(rec.t - 4.0) < 1e-9);
    assert(std::abs(rec.normal.z - 1.0) < 1e-9);

    Triangle tri(
        Point3(-1.0, -1.0, -3.0),
        Point3( 1.0, -1.0, -3.0),
        Point3( 0.0,  1.0, -3.0),
        0
    );
    HitRecord trec;
    bool thit = tri.intersect(r, 1e-5, 100.0, trec);
    assert(thit);
    assert(std::abs(trec.t - 3.0) < 1e-9);
    assert(trec.u >= 0.0 && trec.v >= 0.0 && (trec.u + trec.v) <= 1.0);
    std::cout << "[PASS] Interseccion analitica Rayo-Esfera y Rayo-Triangulo (Moller-Trumbore)\n";
}

void test_monte_carlo_hemisphere_integral() {
    // Evaluamos numéricamente I = \int_\Omega (1/\pi) \cos\theta d\omega = 1.0
    // y demostramos que con muestreo uniforme Var[X] = 1/3 \approx 0.333333,
    // mientras que con muestreo ponderado por coseno Var[X] = 0.000000.
    PCG32 rng(42ULL, 1ULL);
    const int N = 100000;
    double sum_u = 0.0, sum_sq_u = 0.0;
    double sum_c = 0.0, sum_sq_c = 0.0;
    Vec3 normal(0.0, 1.0, 0.0);

    for (int i = 0; i < N; ++i) {
        double xi1 = rng.next_double();
        double xi2 = rng.next_double();

        auto su = sample_hemisphere_uniform(normal, xi1, xi2);
        double val_u = (INV_PI * su.cos_theta) / su.pdf; // = 2 * cos_theta
        sum_u += val_u;
        sum_sq_u += val_u * val_u;

        auto sc = sample_hemisphere_cosine(normal, xi1, xi2);
        double val_c = (INV_PI * sc.cos_theta) / sc.pdf; // = 1.0 exacto
        sum_c += val_c;
        sum_sq_c += val_c * val_c;
    }

    double mean_u = sum_u / N;
    double var_u = (sum_sq_u / N) - (mean_u * mean_u);
    double mean_c = sum_c / N;
    double var_c = (sum_sq_c / N) - (mean_c * mean_c);

    assert(std::abs(mean_u - 1.0) < 0.01);
    assert(std::abs(var_u - (1.0 / 3.0)) < 0.01);
    assert(std::abs(mean_c - 1.0) < 1e-12);
    assert(var_c < 1e-24);

    std::cout << "[PASS] Estimador Monte Carlo Hemisferico: <Lo>_unif = " << mean_u
              << " (Var = " << var_u << ", teorica 1/3 = 0.333333) | <Lo>_cos = "
              << mean_c << " (Var = " << var_c << ")\n";
}

int main() {
    std::cout << "Ejecutando pruebas unitarias del nucleo matematico y geometrico...\n";
    test_vector_and_onb();
    test_ray_sphere_and_moller_trumbore();
    test_monte_carlo_hemisphere_integral();
    std::cout << "Todas las pruebas matematicas pasaron exitosamente.\n";
    return 0;
}
