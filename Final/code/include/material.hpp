#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>
#include "texture.hpp"

#include "ray.hpp"
#include "hit.hpp"
#include <iostream>

// TODO: Implement Shade function that computes Phong introduced in class.
class Material {
public:

    enum ReflectType { kDiffuse, kSpecular, kRefractive };

    explicit Material(const char *texture_filename, 
                      const Vector3f &ambient_color, const Vector3f &diffuse_color,
                      const Vector3f &specular_color = Vector3f::ZERO, 
                      ReflectType type = kDiffuse, double refractive_index = 1.0) :
            ambient_color_(ambient_color), diffuse_color_(diffuse_color), specular_color_(specular_color),
            refl_type_(type), refractive_index_(refractive_index) {
        if (texture_filename[0] != '\0') texture_ = new Texture(texture_filename);
        else texture_ = nullptr;
    }

    virtual ~Material() = default;

    virtual Vector3f ambientColor() const {
        return ambient_color_;
    }

    virtual Vector3f diffuseColor(const Hit &h) const {
        if (texture_) return texture_->pixelAt(h.texturePos());
        return diffuse_color_;
    }

    virtual ReflectType reflectType() const {
        return refl_type_;
    }

    virtual double refractiveIndex() const {
        return refractive_index_;
    }

    Texture *texture() {
        return texture_;
    }

    void print() const {
        printf("with material with reflect type %d\n", refl_type_);
        printf("ambient color ");
        ambient_color_.print();
        printf("diffuse color ");
        diffuse_color_.print();
        if (texture_) printf("With texture \n");
        printf("refractive index %f\n", refractive_index_);
    }

protected:
    Vector3f ambient_color_;
    Vector3f diffuse_color_;
    Vector3f specular_color_;
    ReflectType refl_type_;
    Texture *texture_;
    double refractive_index_;
};


#endif // MATERIAL_H
