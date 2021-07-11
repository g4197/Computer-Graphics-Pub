#include "triangle.hpp"

Triangle::Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c, Material *m) : Object3D(m) {
    vertices_[0] = a, vertices_[1] = b, vertices_[2] = c;
    normal_ = Vector3f::cross(a - b, a - c);
    normal_.normalize();
}

bool Triangle::intersect(const Ray &ray, Hit &hit, float tmin) {
    Vector3f e_1 = vertices_[0] - vertices_[1];
    Vector3f e_2 = vertices_[0] - vertices_[2];
    Vector3f s = vertices_[0] - ray.getOrigin();
    const Vector3f &r_d = ray.getDirection();
    float det_r_e = Matrix3f(r_d, e_1, e_2).determinant();
    if (abs(det_r_e) < 1e-8) return false;
    const Vector3f solution = Vector3f(
        Matrix3f(s, e_1, e_2).determinant(),
        Matrix3f(r_d, s, e_2).determinant(),
        Matrix3f(r_d, e_1, s).determinant()
    ) / det_r_e;
    float t = solution[0], beta = solution[1], gamma = solution[2];
    if (t <= 0 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1 || beta + gamma > 1) {
        return false;
    }
    if (tmin < t && t < hit.getT()) {
        hit.set(t, material, normal_);
        return true;
    }
    return false;
}