#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement functions and add more fields as necessary

class Sphere : public Object3D {
public:
    enum VectorPlace {kInside, kOnEdge, kOutside};
    Sphere();

    Sphere(const Vector3f &center, float radius, Material *material);

    ~Sphere() override = default;

    VectorPlace getVectorPlace(const Vector3f &v);

    bool intersect(const Ray &r, Hit &h, float tmin) override;

    bool setHit(Hit &h, float tmin, float t, const Vector3f &intersect_pnt, VectorPlace origin_place);

    void print() override {
        printf("A sphere with radius %f and center ", radius_);
        center_.print();
        material->getDiffuseColor().print();
    }

protected:
    Vector3f center_;
    float radius_;
};


#endif
