#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>
using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle : public Object3D {

public:
    Triangle() = delete;

    // a b c are three vertex positions of the triangle
    Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c, Material *m);
    Vector3f &normal() {return normal_;}
    bool intersect(const Ray &ray, Hit &hit, float tmin) override;
    void print() override {
        printf("A triangle with vertices\n");
        vertices_[0].print(), vertices_[1].print(), vertices_[2].print();
    }
    Vector3f normal_;
    Vector3f vertices_[3];
protected:
};

#endif //TRIANGLE_H
