#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <glut.h>

// TODO (PA2): Copy from PA1

class Sphere : public Object3D {
public:
    enum VectorPlace {kInside, kOnEdge, kOutside};
    Sphere() : center(Vector3f::ZERO), radius(1.0) {
        // unit ball at the center
    }

    Sphere(const Vector3f &center, float radius, Material *material) : Object3D(material), center(center), radius(radius) {
        //
    }

    ~Sphere() override = default;

    VectorPlace getVectorPlace(const Vector3f &v) {
        float distance = (v - center).length();
        return distance < radius ? kInside : 
            (distance > radius ? kOutside : kOnEdge);
    }

    bool setHit(Hit &h, float tmin, float t, const Vector3f &intersect_pnt, VectorPlace origin_place) {
        if (tmin < t && t < h.getT()) {
            Vector3f normal = (intersect_pnt - center).normalized();
            if (origin_place == kInside) normal = -normal;
            h.set(t, this->material, normal);
            return true;
        }
        return false;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        VectorPlace origin_place = getVectorPlace(r.getOrigin());
        Vector3f origin_to_center = center - r.getOrigin();
        float len_r = r.getDirection().length();

        // proj_t means the dot of origin_to_center and unit r
        // then the square_length is 1 * 1 * proj_t * proj_t
        float proj_t = Vector3f::dot(origin_to_center, r.getDirection()) / len_r;
        float proj_square_len = proj_t * proj_t;
        float origin_to_center_square_len = origin_to_center.squaredLength();
        if (origin_place != kInside && proj_t < 0) {
            return false;
        }

        float square_dis = origin_to_center_square_len - proj_square_len;
        float proj_intersect_square_dis = radius * radius - square_dis;
        if (proj_intersect_square_dis < 0) {
            return false;
        }

        float proj_intersect_dis = sqrt(proj_intersect_square_dis);
        float t = proj_t;
        if (origin_place == kOutside) {
            t -= proj_intersect_dis;
        } else {
            t += proj_intersect_dis;
        }
        t /= len_r;
        return setHit(h, tmin, t, r.pointAtParameter(t), origin_place);
    }

    void drawGL() override {
        Object3D::drawGL();
        glMatrixMode(GL_MODELVIEW); glPushMatrix();
        glTranslatef(center.x(), center.y(), center.z());
        glutSolidSphere(radius, 80, 80);
        glPopMatrix();
    }
protected:
    Vector3f center;
    float radius;

};


#endif
