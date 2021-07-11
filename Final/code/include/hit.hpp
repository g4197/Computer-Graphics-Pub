#ifndef HIT_H
#define HIT_H

#include <vecmath.h>
#include "ray.hpp"

class Material;

class Hit {
public:

    // constructors
    Hit() {
        material_ = nullptr;
        t_ = 1e38;
    }

    Hit(double t, Material *m, const Vector3f &n, bool into = true, const Vector2f &texture_pos = Vector2f()) {
        t_ = t;
        material_ = m;
        normal_ = n;
        into_ = into;
        texture_pos_ = texture_pos;
    }

    Hit(const Hit &h) {
        t_ = h.t_;
        material_ = h.material_;
        normal_ = h.normal_;
    }

    // destructor
    ~Hit() = default;

    double t() const {
        return t_;
    }

    Material *material() const {
        return material_;
    }

    const Vector3f &normal() const {
        return normal_;
    }

    const Vector2f &texturePos() const {
        return texture_pos_;
    }

    bool into() const {
        return into_;
    }

    void set(double t, Material *m, const Vector3f &n, bool into = true, const Vector2f &texture_pos = Vector2f()) {
        t_ = t;
        material_ = m;
        normal_ = n;
        into_ = into;
        texture_pos_ = texture_pos;
    }

private:
    double t_;
    Material *material_;
    Vector3f normal_;
    bool into_;
    Vector2f texture_pos_;
};

inline std::ostream &operator<<(std::ostream &os, const Hit &h) {
    os << "Hit <" << h.t() << ", " << h.normal() << ">";
    return os;
}

#endif // HIT_H
