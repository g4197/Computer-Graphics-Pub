#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane();

    Plane(const Vector3f &normal, float d, Material *m);

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override;

    void print() override {
        printf("A plane with d %f and normal ", d_);
        normal_.print();
    }

protected:
    Vector3f normal_;
    float d_;
};

#endif //PLANE_H
		

