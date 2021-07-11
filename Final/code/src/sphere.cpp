#include "sphere.hpp"

Sphere::Sphere() : center_(Vector3f::ZERO), radius_(1.0) {}

Sphere::Sphere(const Vector3f &center, double radius, Material *material) : 
    Object3D(material, center - radius * 1.415, center + radius * 1.415), 
    center_(center), radius_(radius) {}

Vector3f Sphere::center() const {
    return center_;
}

double Sphere::radius() const {
    return radius_;
}

bool Sphere::setHit(Hit &h, double tmin, double t, const Vector3f &intersect_pnt, const Ray &r) {
    if (tmin < t && t < h.t()) {
        Vector3f normal = (intersect_pnt - center_).normalized();
        if (Vector3f::dot(normal, r.direction()) > 0) {
            h.set(t, this->material_, -normal, false);
        } else {
            h.set(t, this->material_, normal, true);
        }
        return true;
    }
    return false;
}

bool Sphere::intersect(const Ray &r, Hit &h, double tmin) {
    Vector3f origin_to_center = center_ - r.origin();
    double len_r = r.direction().length();

    double t = tmin - 100;
    double proj_t = origin_to_center.dot(r.direction()) / len_r;

    double det = proj_t * proj_t - origin_to_center.dot(origin_to_center) + radius_ * radius_;

    if (det >= 0) det = sqrt(det);
    else return false;

    t = (proj_t - det) > tmin ? (proj_t - det) : 
       ((proj_t + det) > tmin ? (proj_t + det) : tmin - 100);
    t /= len_r;

    return setHit(h, tmin, t, r.pointAtParameter(t), r);
}