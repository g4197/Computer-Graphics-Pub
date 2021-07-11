#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include "utils.hpp"
#include <tuple>

class RevSurface : public Object3D {

    Curve *curve_;

public:
    RevSurface(Curve *pCurve, Material* material);

    ~RevSurface() override {
        delete curve_;
    }

    bool intersect(const Ray &r, Hit &h, double tmin) override;

    void print() override {
        printf("A revsurface with %lu control points\n", curve_->getControls().size());
    }

    bool intersectDyIsZero(const Ray &r, Hit &h, double tmin);

    bool intersectDyNotZero(const Ray &r, Hit &h, double tmin);
};

#endif //REVSURFACE_HPP
