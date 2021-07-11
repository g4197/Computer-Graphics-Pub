#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vecmath.h>
#include "object3d.hpp"

// transforms a 3D point using a matrix, returning a 3D point
static Vector3f transformPoint(const Matrix4f &mat, const Vector3f &point) {
    return (mat * Vector4f(point, 1)).xyz();
}

// transform a 3D directino using a matrix, returning a direction
static Vector3f transformDirection(const Matrix4f &mat, const Vector3f &dir) {
    return (mat * Vector4f(dir, 0)).xyz();
}

// TODO: implement this class so that the intersect function first transforms the ray
class Transform : public Object3D {
public:
    Transform() {}

    Transform(const Matrix4f &m, Object3D *obj) : o(obj) {
        transform = m.inverse();
    }

    ~Transform() {
    }

    virtual bool intersect(const Ray &r, Hit &h, double tmin) {
        Vector3f trSource = transformPoint(transform, r.origin());
        Vector3f trDirection = transformDirection(transform, r.direction());
        Ray tr(trSource, trDirection);
        bool inter = o->intersect(tr, h, tmin);
        if (inter) {
            h.set(h.t(), h.material(), transformDirection(transform.transposed(), h.normal()).normalized(), h.into(), h.texturePos());
        }
        return inter;
    }

    void print() override {
        o->print();
    }

protected:
    Object3D *o; //un-transformed object
    Matrix4f transform;
};

#endif //TRANSFORM_H
