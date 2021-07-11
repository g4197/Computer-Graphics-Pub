#include "triangle.hpp"

Triangle::Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c, Material *m) :
    Triangle(a, b, c, m, Vector2f(), Vector2f(), Vector2f()) {}

Triangle::Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c, Material *m, 
                   const Vector2f &text_a, const Vector2f &text_b, const Vector2f &text_c) :
    Triangle(a, b, c, Vector3f::cross(a - b, a - c).normalized(), m, text_a, text_b, text_c) {}

Triangle::Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c, const Vector3f &n, Material *m,
                   const Vector2f &text_a, const Vector2f &text_b, const Vector2f &text_c) : Object3D(m) {
    vertices_[0] = a, vertices_[1] = b, vertices_[2] = c;
    texture_pos_[0] = text_a, texture_pos_[1] = text_b, texture_pos_[2] = text_c;
    normal_ = n;
    double min_arr[3], max_arr[3];
    for (int i = 0; i < 3; ++i) {
        min_arr[i] = min(min(a[i], b[i]), c[i]);
        max_arr[i] = max(max(a[i], b[i]), c[i]);
    }
    Vector3f l(min_arr[0], min_arr[1], min_arr[2]), r(max_arr[0], max_arr[1], max_arr[2]);
    setAABB(l, r);
}

bool Triangle::intersect(const Ray &ray, Hit &hit, double tmin) {
    Vector3f e_1 = vertices_[0] - vertices_[1];
    Vector3f e_2 = vertices_[0] - vertices_[2];
    Vector3f s = vertices_[0] - ray.origin();
    const Vector3f &r_d = ray.direction();
    const Vector3f solution = Vector3f(
        Matrix3f(s, e_1, e_2).determinant(),
        Matrix3f(r_d, s, e_2).determinant(),
        Matrix3f(r_d, e_1, s).determinant()
    ) / Matrix3f(r_d, e_1, e_2).determinant();
    double t = solution[0], beta = solution[1], gamma = solution[2];
    if (t <= 0 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1 || beta + gamma > 1) {
        return false;
    }
    if (tmin < t && t < hit.t()) {
        hit.set(t, material_, normal_, true, 
               (1 - beta - gamma) * texture_pos_[0] + beta * texture_pos_[1] + gamma * texture_pos_[2]);
        return true;
    }
    return false;
}