#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <vecmath.h>
#include <float.h>
#include <cmath>


class Camera {
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up);
		this->horizontal.normalize();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// TODO: Implement Perspective camera
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {

public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle) : Camera(center, direction, up, imgW, imgH) {
        // angle is both horizontal and vertical FOV, according to simialrity of triangle, we have the equation below.
        float double_tan_half_angle = 2 * tan(angle / 2);
        fx_ = width / double_tan_half_angle;
        fy_ = height / double_tan_half_angle;
    }

    float centerX() {
        return width / 2;
    }

    float centerY() {
        return height / 2;
    }

    Ray generateRay(const Vector2f &point) override {
        Vector3f cam_axis_ray((point.x() - centerX()) / fx_, (centerY() - point.y()) / fy_, 1);
        cam_axis_ray.normalize();

        //linear combination
        Vector3f normal_axis_ray = cam_axis_ray.x() * horizontal - 
                                   cam_axis_ray.y() * up + 
                                   cam_axis_ray.z() * direction;
        Ray ray(center, normal_axis_ray);
        return ray;
    }
protected:
    float fx_, fy_;
};

#endif //CAMERA_H
