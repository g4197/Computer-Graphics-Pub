#include "utils.hpp"

std::default_random_engine generator;
std::uniform_real_distribution<double> distr(0.0, 1.0);

double max(const Vector3f &v) {
    return (v.x() > v.y() && v.x() > v.z()) ? v.x() : 
           (v.y() > v.z()) ? v.y() : v.z();
}

double erand48() {
    // generate a uniform-distribution random double from 0 to 1
    return distr(generator);
}

void getXYbyZ(Vector3f &x, Vector3f &y, Vector3f &z) {
    // get two axes on the surface with normal z
    Vector3f base = abs(z.x()) > 0.1 ? Vector3f::UP : Vector3f::RIGHT; // don't be linear dependent
    x = Vector3f::cross(base, z).normalized();
    y = Vector3f::cross(x, z).normalized();
}

Vector3f minVec(const Vector3f &v1, const Vector3f &v2) {
    return Vector3f(
        std::min(v1[0], v2[0]),
        std::min(v1[1], v2[1]),
        std::min(v1[2], v2[2])
    );
}

Vector3f maxVec(const Vector3f &v1, const Vector3f &v2) {
    return Vector3f(
        std::max(v1[0], v2[0]),
        std::max(v1[1], v2[1]),
        std::max(v1[2], v2[2])
    );
}