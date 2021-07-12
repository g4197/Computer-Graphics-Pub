#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <vecmath.h>
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
            const Vector3f &up, int imgW, int imgH, double angle_w,
            double focal, double aperture) : Camera(center, direction, up, imgW, imgH) {
        // angle is both horizontal and vertical FOV, according to simialrity of triangle, we have the equation below.
        double double_tan_half_angle = 2 * tan(angle_w / 2);
        fx_ = fy_ = width / double_tan_half_angle;
        focal_ = focal, aperture_ = aperture;
    }

    double centerX() {
        return width / 2;
    }

    double centerY() {
        return height / 2;
    }

    Ray generateRay(const Vector2f &point) override {
        Vector3f cam_axis_ray((point.x() - centerX()) / fx_, (centerY() - point.y()) / fy_, 1);
        cam_axis_ray.normalize();
 
        Vector3f delta = aperture_ * Vector3f(2 * (erand48() - 0.5), 2 * (erand48() - 0.5), 0);
        Vector3f dir = focal_ * cam_axis_ray + delta;

        //linear combination
        Vector3f normal_dir = (dir.x() * horizontal - 
                               dir.y() * up + 
                               dir.z() * direction).normalized(); // generate image on z=aperture_
        


        Ray ray(center - (delta.x() * horizontal - delta.y() * up), normal_dir);
        return ray;
    }
protected:
    double fx_, fy_;
    double focal_, aperture_;
};

#endif //CAMERA_H
