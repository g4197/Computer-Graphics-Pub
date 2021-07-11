#include "plane.hpp"

Plane::Plane() : normal_(Vector3f::ZERO), d_(0.0) {}

Plane::Plane(const Vector3f &normal, float d, Material *m) : 
             Object3D(m), normal_(normal), d_(d) {}

bool Plane::intersect(const Ray &r, Hit &h, float tmin) {
    float normal_direction_dot = Vector3f::dot(normal_, r.getDirection());
    float t = (d_ - Vector3f::dot(normal_, r.getOrigin())) / (normal_direction_dot);
    if (t <= 0) {
        return false;
    }
    if (tmin < t && t < h.getT()) {
        h.set(t, this->material, normal_direction_dot > 0 ? -normal_ : normal_);
        return true;
    }
    return false;
}