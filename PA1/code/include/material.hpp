#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include <iostream>

// TODO: Implement Shade function that computes Phong introduced in class.
class Material {
public:

    explicit Material(const Vector3f &d_color, const Vector3f &s_color = Vector3f::ZERO, float s = 0) :
            diffuseColor(d_color), specularColor(s_color), shininess(s) {

    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }


    Vector3f Shade(const Ray &ray, const Hit &hit,
                   const Vector3f &dirToLight, const Vector3f &lightColor) {
        Vector3f shaded = Vector3f::ZERO;
        Vector3f normalized_to_light = dirToLight.normalized();
        float normal_dot_light = Vector3f::dot(normalized_to_light, hit.getNormal());
        if (normal_dot_light >= 0) shaded += normal_dot_light * diffuseColor;
        Vector3f reflect_dir = 2 * (normal_dot_light) * hit.getNormal() - normalized_to_light;
        reflect_dir.normalize();

        //negative because when reflect light is to our eye, we can see most light
        float ray_dot_reflect_dir = -Vector3f::dot(ray.getDirection(), reflect_dir);
        if (ray_dot_reflect_dir >= 0) shaded += pow(ray_dot_reflect_dir, shininess) * specularColor;
        shaded = shaded * lightColor;
        return shaded;
    }

protected:
    Vector3f diffuseColor;
    Vector3f specularColor;
    float shininess;
};


#endif // MATERIAL_H
