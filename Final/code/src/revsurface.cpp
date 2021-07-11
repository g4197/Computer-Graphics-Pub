#include "revsurface.hpp"

RevSurface::RevSurface(Curve *pCurve, Material* material) : curve_(pCurve), Object3D(material) {
    // Check flat.
    Vector3f l = kVecInf, r = -kVecInf;
    for (const auto &cp : pCurve->getControls()) {
        if (cp.z() != 0.0) {
            printf("Profile of revSurface must be flat on xy plane.\n");
            exit(0);
        }
        // l = -max abs(x), miny, -max abs(x); r = max abs(x), max y, max abs(x)
        Vector3f v(abs(cp.x()), cp.y(), abs(cp.z()));
        l[0] = min(l[0], -v[0]), l[1] = min(l[1], v[1]), l[2] = min(l[2], -v[0]);
        r[0] = max(r[0], v[0]), r[1] = max(r[1], v[1]), r[2] = max(r[2], v[0]);
    }
    this->setAABB(l, r);
}

bool RevSurface::intersect(const Ray &r, Hit &h, double tmin) {
    // Force init it
    // t of ray, theta and place on curve
    if (!this->aabb_.intersect(r)) return false;
    if (r.direction().y() != 0) {
        return intersectDyNotZero(r, h, tmin);
    } else {
        return intersectDyIsZero(r, h, tmin);
    }
}

bool RevSurface::intersectDyIsZero(const Ray &r, Hit &h, double tmin) {
    // Intersect when d.y() == 0
    // equation: y - o.y = 0
    // derive: y'
    // (r_t * d_x + o_x) ** 2 + (r_t * d_z + o_z) ** 2 = x_t ** 2
    bool ret = false;
    const Vector3f &o = r.origin(), &d = r.direction(), &rd = r.revDirection();
    double start = curve_->start(), end = curve_->end();
    for (int i = 0; i < kMaxRandom; ++i) {
        double cur_place = start + (end - start) * double(i) / kMaxRandom;
        for (int j = 0; j < kMaxIteration; ++j) {
            auto point = curve_->pointAtParam(cur_place);
            const Vector3f &v = point.V, &t = point.T;
            double solution = v.y() - o.y();
            if (solution > -kNewtonIterEps && solution < kNewtonIterEps) {
                // solve the r_t equation
                double a = 1.0; // d.x() * d.x() + d.z() * d.z(), for d.y = 0
                double half_b = d.x() * o.x() + d.z() * o.z();
                double c = o.x() * o.x() + o.z() * o.z() - v.x() * v.x();
                double r_t = (-half_b - sqrt(half_b * half_b - c)); // a = 1
                if (tmin < r_t && r_t < h.t() - kEps) {
                    double cos_theta = v.x() == 0 ? 0 : (r_t * d.x() + o.x()) / v.x();
                    double sin_theta = v.x() == 0 ? 1 : (r_t * d.z() + o.z()) / v.x();
                    Vector3f n = Vector3f::cross(t, -Vector3f::FORWARD);
                    Vector3f n_rot = Vector3f(
                        cos_theta * n.x() - sin_theta * n.z(),
                        n.y(),
                        sin_theta * n.x() + cos_theta * n.z()
                    );
                    // assume that normal at x < 0 all to x < 0
                    bool into = true;
                    if (Vector3f::dot(n_rot, d) > 0) {
                        n_rot = -n_rot;
                        into = false;
                    }
                    h.set(r_t, material_, n_rot, into, 
                          Vector2f((cur_place - start) / (end - start), acos(clamp(cos_theta)) / (2 * kPi)));
                    ret = true;
                }
                break;
            } // solved
            double derive = t.y();
            cur_place -= solution / (kDeriveShrink * derive);
            if (isnan(cur_place) || cur_place < 0 || cur_place > 1) break;
        }
    }
    return ret;
}

bool RevSurface::intersectDyNotZero(const Ray &r, Hit &h, double tmin) {
    // Intersect when d.y() != 0
    // equation:
    // ((y - o.y) * d.x / d.y + o.x) ** 2 + ((y - o.y) * d.z / d.y + o.z) ** 2 - x ** 2 = 0
    // derive:
    // 2((y - o.y) * d.x / d.y + o.x) * (d.x / d.y) * y' +
    // 2((y - o.y) * d.z / d.y + o.z) * (d.z / d.y) * y' -
    // 2xx'
    bool ret = false;
    const Vector3f &o = r.origin(), &d = r.direction(), &rd = r.revDirection();
    double d_x_y = d.x() * rd.y(), d_z_y = d.z() * rd.y();
    double start = curve_->start(), end = curve_->end();
    for (int i = 0; i < kMaxRandom; ++i) {
        double cur_place = start + (end - start) * double(i) / kMaxRandom;
        for (int j = 0; j < kMaxIteration; ++j) {
            auto point = curve_->pointAtParam(cur_place);
            const Vector3f &v = point.V, &t = point.T;
            double x_cos_theta = (v.y() - o.y()) * d_x_y + o.x();
            double x_sin_theta = (v.y() - o.y()) * d_z_y + o.z();
            double solution = x_cos_theta * x_cos_theta + x_sin_theta * x_sin_theta - v.x() * v.x();
            if (solution > -kNewtonIterEps && solution < kNewtonIterEps) {
                double r_t = max((x_cos_theta - o.x()) * rd.x(), (x_sin_theta - o.z()) * rd.z()); // guarantee no 0 when d.x() == 0 (on Y axis)
                if (tmin < r_t && r_t < h.t() - kEps) {
                    double cos_theta = v.x() == 0 ? 0 : x_cos_theta / v.x();
                    double sin_theta = v.x() == 0 ? 1 : x_sin_theta / v.x();
                    Vector3f n = Vector3f::cross(t, -Vector3f::FORWARD);
                    Vector3f n_rot = Vector3f(
                        cos_theta * n.x() - sin_theta * n.z(),
                        n.y(),
                        sin_theta * n.x() + cos_theta * n.z()
                    );
                    // assume that normal at x < 0 all to x < 0
                    bool into = true;
                    if (Vector3f::dot(n_rot, d) > 0) {
                        n_rot = -n_rot;
                        into = false;
                    }
                    h.set(r_t, material_, n_rot, into, 
                          Vector2f((cur_place - start) / (end - start), acos(clamp(cos_theta)) / (2 * kPi)));
                    ret = true;
                }
                break;
            } // solved
            double derive = 2 * x_cos_theta * d_x_y * t.y() + 
                            2 * x_sin_theta * d_z_y * t.y() - 2 * v.x() * t.x();
            cur_place -= solution / (kDeriveShrink * derive);
            if (isnan(cur_place) || cur_place < 0 || cur_place > 1) break;
        }
    }
    return ret;
}