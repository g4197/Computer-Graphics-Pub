#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>

using namespace std;

// TODO (PA2): Copy from PA1
class Triangle: public Object3D
{

public:
    Triangle() = delete;
        ///@param a b c are three vertex positions of the triangle

	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) : Object3D(m) {
        vertices[0] = a, vertices[1] = b, vertices[2] = c;
        normal = Vector3f::cross(a - b, a - c);
        normal.normalize();
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
        Vector3f e_1 = vertices[0] - vertices[1];
        Vector3f e_2 = vertices[0] - vertices[2];
        Vector3f s = vertices[0] - ray.getOrigin();
        const Vector3f &r_d = ray.getDirection();
        float det_r_e = Matrix3f(r_d, e_1, e_2).determinant();
        if (abs(det_r_e) < 1e-8) return false;
        const Vector3f solution = Vector3f(
            Matrix3f(s, e_1, e_2).determinant(),
            Matrix3f(r_d, s, e_2).determinant(),
            Matrix3f(r_d, e_1, s).determinant()
        ) / det_r_e;
        float t = solution[0], beta = solution[1], gamma = solution[2];
        if (t <= 0 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1 || beta + gamma > 1) {
            return false;
        }
        if (tmin < t && t < hit.getT()) {
            hit.set(t, material, normal);
            return true;
        }
        return false;
	}

    void drawGL() override {
        Object3D::drawGL();
        glBegin(GL_TRIANGLES);
        glNormal3fv(normal);
        glVertex3fv(vertices[0]); glVertex3fv(vertices[1]); glVertex3fv(vertices[2]);
        glEnd();
    }
    Vector3f normal;
	Vector3f vertices[3];

protected:
};

#endif //TRIANGLE_H
