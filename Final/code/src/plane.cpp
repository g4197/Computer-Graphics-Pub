#include "plane.hpp"

Plane::Plane() : normal_(Vector3f::ZERO), d_(0.0) {}

Plane::Plane(const Vector3f &normal, double d, Material *m) : 
             Object3D(m), normal_(normal), d_(d) {}

bool Plane::intersect(const Ray &r, Hit &h, double tmin) {
    double normal_direction_dot = Vector3f::dot(normal_, r.direction());
    double t = (d_ - Vector3f::dot(normal_, r.origin())) / (normal_direction_dot);
    if (t <= 0) {
        return false;
    }
    Vector3f pos = r.pointAtParameter(t);
    Vector2f proj = normal_.y() > 0.1 ? pos.xz() : 
                   (normal_.x() > 0.1 ? pos.yz() : pos.xy());
    if (tmin < t && t < h.t()) {
        h.set(t, this->material_, normal_direction_dot > 0 ? -normal_ : normal_, true, proj / kPlaneTextureScale);
        return true;
    }
    return false;
}