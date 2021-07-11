#ifndef OBJECT3D_H
#define OBJECT3D_H

#include "aabb.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include "material.hpp"
#include "utils.hpp"

// Base class for all 3d entities.
class Object3D {
public:
    Object3D() : material_(nullptr), aabb_(-kVecInf, kVecInf) {}

    explicit Object3D(Material *material) : material_(material), aabb_(-kVecInf, kVecInf) {}

    explicit Object3D(Material *material, const Vector3f &l, const Vector3f &r) : material_(material), aabb_(l, r) {}

    virtual ~Object3D() = default;

    Material *material() const { return material_; }
    const AABB &aabb() const { return aabb_; }

    void setAABB(const Vector3f &l, const Vector3f &r) { aabb_.set(l, r); }
    virtual bool intersect(const Ray &r, Hit &h, double tmin) = 0;
    virtual void print() {}
protected:
    AABB aabb_;
    Material *material_;
};

#endif

