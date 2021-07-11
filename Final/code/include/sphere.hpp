#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement functions and add more fields as necessary

class Sphere : public Object3D {
public:
    Sphere();

    Sphere(const Vector3f &center, double radius, Material *material);

    ~Sphere() override = default;

    Vector3f center() const;
    
    double radius() const;

    bool intersect(const Ray &r, Hit &h, double tmin) override;

    bool setHit(Hit &h, double tmin, double t, const Vector3f &intersect_pnt, const Ray &r);

    void print() override {
        printf("A sphere with radius %f and center ", radius_);
        center_.print();
        material()->print();
    }

protected:
    Vector3f center_;
    double radius_;
};


#endif
