#ifndef TEXTURE_H_
#define TEXTURE_H_

#include <iostream>
#include "utils.hpp"
#include "vecmath.h"
#include <image.hpp>

using namespace std;

class Texture {
public:
    Texture(const char *filename) {
        image_ = Image::LoadBMP(filename);
    }

    Vector3f pixelAt(const Vector2f &x) {
        int w = image_->Width(), h = image_->Height();
        double texture_x = x[0] - int(x[0]), texture_y = x[1] - int(x[1]);
        texture_x = texture_x > 0.0 ? texture_x : texture_x + 1.0;
        texture_y = texture_y > 0.0 ? texture_y : texture_y + 1.0;
        return image_->GetPixel(texture_x * w, texture_y * h);
    }
private:
    Image *image_;
};

#endif // TEXTURE_H_