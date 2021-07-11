#ifndef UTILS_H_
#define UTILS_H_

#include <random>
#include "vecmath.h"

extern std::default_random_engine generator;
extern std::uniform_real_distribution<double> distr;

const double kInf = 1e30;
const double kEps = 1e-4;
const double kNewtonIterEps = 1e-4;
const double kTmin = 1e-4;
const double kEqualEps = 1e-8;
const double kPi = 3.141592653589793;
const double k1Pi = 1 / kPi;
const double kNGlass = 1.5;
const double kNWater = 1.33;
const double kNDiamond = 2.5;
const double kNAir = 1;
const double kAlpha = 0.7;
const double kPlaneTextureScale = 5; // in plane, 5 * 5 will be one texture
const int kPPMTurn = 20; // turns for ppm (for depth and tent filter)
const int kMaxDepth = 8;
const int kMaxIteration = 64;
const int kDeriveShrink = 32;
const int kMaxRandom = 128;
const int kMaxThread = 16;
const int kHPointPoolSize = 1 << 23 | 1;

const Vector3f kVecInf = Vector3f(kInf, kInf, kInf);

inline double DegreesToRadians(double x) {
    return (kPi * x) / 180.0;
}

double max(const Vector3f &v);

double erand48();

void getXYbyZ(Vector3f &x, Vector3f &y, Vector3f &z);

inline Vector3f clamp(const Vector3f &v) {
    Vector3f ret;
    for (int i = 0; i <= 2; ++i) {
        ret[i] = v[i] < 0 ? 0 : (v[i] > 1 ? 1 : v[i]);
    }
    return ret;
}

inline double clamp(double x) {
    return (x > 1.0) ? 1.0 : ((x < -1.0) ? -1.0 : x);
}

Vector3f minVec(const Vector3f &v1, const Vector3f &v2);

Vector3f maxVec(const Vector3f &v1, const Vector3f &v2);

#endif //UTILS_H