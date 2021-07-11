#include "aabb.hpp"

using namespace std;

AABB::AABB() : l_(-kVecInf), r_(kVecInf) {}

AABB::AABB(const Vector3f &l, const Vector3f &r) : l_(l), r_(r) {}

bool AABB::intersect(const Ray &ray) {
    if (l_ == -kVecInf && r_ == kVecInf) return true;   
    const Vector3f &origin = ray.origin();
    const Vector3f &rev_dir = ray.revDirection();
    Vector3f v1 = (l_ - origin) * rev_dir;
    Vector3f v2 = (r_ - origin) * rev_dir;
    double max_min = max(max(min(v1[0], v2[0]), min(v1[1], v2[1])), min(v1[2], v2[2]));
    double min_max = min(min(max(v1[0], v2[0]), max(v1[1], v2[1])), max(v1[2], v2[2]));
    if (min_max < 0 || max_min > min_max) return false;
    return true;
}

void AABB::print() {
    printf("A bounding box with left");
    l_.print();
    printf("and right");
    r_.print();
}

void AABB::set(const Vector3f &l, const Vector3f &r) {
    this->l_ = l, this->r_ = r;
}

void AABB::fit(const Vector3f &v) {
    set(minVec(l_, v), maxVec(r_, v));
}

Vector3f AABB::l() const {
    return this->l_;
}

Vector3f AABB::r() const {
    return this->r_;
}

void AABB::reset() {
    this->l_ = kVecInf;
    this->r_ = -kVecInf;
}