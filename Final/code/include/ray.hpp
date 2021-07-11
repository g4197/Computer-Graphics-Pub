#ifndef RAY_H
#define RAY_H

#include <cassert>
#include <iostream>
#include <Vector3f.h>
#include "utils.hpp"

// Ray class mostly copied from Peter Shirley and Keith Morley
class Ray {
public:

    Ray() = delete;
    Ray(const Vector3f &orig, const Vector3f &dir) {
        origin_ = orig;
        direction_ = dir;
        for (int i = 0; i < 3; ++i) {
            rev_direction_[i] = dir[i] == 0 ? kInf : 1 / dir[i];
        }
    }

    Ray(const Ray &r) {
        origin_ = r.origin_;
        direction_ = r.direction_;
        rev_direction_ = r.rev_direction_;
    }

    const Vector3f &origin() const {
        return origin_;
    }

    const Vector3f &direction() const {
        return direction_;
    }

    const Vector3f &revDirection() const {
        return rev_direction_;
    }

    Vector3f pointAtParameter(double t) const {
        return origin_ + direction_ * t;
    }

    void set(const Vector3f &o, const Vector3f &d) {
        this->origin_ = o, this->direction_ = d;
        for (int i = 0; i < 3; ++i) {
            rev_direction_[i] = d[i] == 0 ? kInf : 1 / d[i];
        }
    }

private:

    Vector3f origin_;
    Vector3f direction_;
    Vector3f rev_direction_;
};

inline std::ostream &operator<<(std::ostream &os, const Ray &r) {
    os << "Ray <" << r.origin() << ", " << r.direction() << ">";
    return os;
}

#endif // RAY_H
