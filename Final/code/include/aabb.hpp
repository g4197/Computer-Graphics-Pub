#ifndef AABB_H_
#define AABB_H_

#include "vecmath.h"
#include "ray.hpp"

class AABB {
public:
    AABB();
    AABB(const Vector3f &l, const Vector3f &r);
    bool intersect(const Ray &r);
    Vector3f l() const;
    Vector3f r() const;
    void print();
    void set(const Vector3f &l, const Vector3f &r);
    void fit(const Vector3f &v);
    void reset();
private:
    Vector3f l_;
    Vector3f r_;
}; //can't set T for intersect, for ray may in AABB, the T will be wrong.

#endif //AABB_H_