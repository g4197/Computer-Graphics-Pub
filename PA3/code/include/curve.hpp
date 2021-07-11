#ifndef CURVE_HPP
#define CURVE_HPP

#include "object3d.hpp"
#include <vecmath.h>
#include <vector>
#include <utility>

#include <algorithm>

// TODO (PA3): Implement Bernstein class to compute spline basis function.
//       You may refer to the python-script for implementation.

// The CurvePoint object stores information about a point on a curve
// after it has been tesselated: the vertex (V) and the tangent (T)
// It is the responsiblility of functions that create these objects to fill in all the data.

typedef double real;
const real kEps = 1e-8;

struct CurvePoint {
    Vector3f V; // Vertex
    Vector3f T; // Tangent  (unit)
};

class Curve : public Object3D {
protected:
    std::vector<Vector3f> controls;
public:
    explicit Curve(std::vector<Vector3f> points) : controls(std::move(points)) {}

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        return false;
    }

    std::vector<Vector3f> &getControls() {
        return controls;
    }

    virtual void discretize(int resolution, std::vector<CurvePoint>& data) = 0;

    void drawGL() override {
        Object3D::drawGL();
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glColor3f(1, 1, 0);
        glBegin(GL_LINE_STRIP);
        for (auto & control : controls) { glVertex3fv(control); }
        glEnd();
        glPointSize(4);
        glBegin(GL_POINTS);
        for (auto & control : controls) { glVertex3fv(control); }
        glEnd();
        std::vector<CurvePoint> sampledPoints;
        discretize(30, sampledPoints);
        glColor3f(1, 1, 1);
        glBegin(GL_LINE_STRIP);
        for (auto & cp : sampledPoints) { glVertex3fv(cp.V); }
        glEnd();
        glPopAttrib();
    }
};

class BezierCurve : public Curve {
public:
    explicit BezierCurve(const std::vector<Vector3f> &points) : Curve(points) {
        if (points.size() < 4 || points.size() % 3 != 1) {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
        int deg = points.size() - 1;
        // there are sz - 1 degrees.
        for (int i = 0; i <= deg; ++i) {
            comb_n_.push_back(combination(deg, i));
        }
        int deg_minus_1 = deg - 1;
        for (int i = 0; i <= deg_minus_1; ++i) {
            comb_n_minus_1_.push_back(combination(deg_minus_1, i));
        }
    }

    inline int fact(int a) {
        int ret = 1;
        for (int i = 2; i <= a; ++i) {
            ret *= i;
        }
        return ret;
    }

    real combination(int b, int a) {
        // return C_b^a
        return static_cast<real>(fact(b)) / (fact(a) * fact(b - a));
    }

    CurvePoint pointAtParam(real t) {
        Vector3f tangent, vertex;
        tangent = vertex = Vector3f::ZERO;
        int n = controls.size() - 1;
        for (int i = 0; i <= n; ++i) {
            vertex += comb_n_[i] * pow(1 - t, n - i) * pow(t, i) * controls[i];
        }
        tangent -= n * comb_n_minus_1_[0] * pow(1 - t, n - 1) * controls[0];
        for (int i = 1; i <= n - 1; ++i) {
            tangent += n * (
                pow(1 - t, n - 1 - i) * pow(t, i - 1) * 
                (comb_n_minus_1_[i - 1] * (1 - t) - comb_n_minus_1_[i] * t)
                // comb_n_minus_1_[i - 1] * pow(1 - t, n - i) * pow(t, i - 1) - 
                // comb_n_minus_1_[i] * pow(1 - t, n - 1 - i) * pow(t, i)
            ) * controls[i];
        }
        tangent += n * comb_n_minus_1_[n - 1] * pow(t, n - 1) * controls[n];
        tangent.normalize();
        CurvePoint p;
        p.T = tangent, p.V = vertex;
        return p;
    }

    void discretize(int resolution, std::vector<CurvePoint>& data) override {
        data.clear();
        int sum_points = resolution * (1 + 2 * controls.size());
        real t = 0, unit = 1.0 / sum_points;
        
        for (int i = 0; i < sum_points; ++i, t += unit) {          
            data.push_back(pointAtParam(t));
        }
        // TODO (PA3): fill in data vector
    }

protected:
    std::vector<real> comb_n_;
    std::vector<real> comb_n_minus_1_;
};

class BsplineCurve : public Curve {
public:
    BsplineCurve(const std::vector<Vector3f> &points, int deg = 3) : Curve(points), deg_(deg) {
        if (points.size() < 4) {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
        int knot_num = points.size() + deg; // n + k + 1
        real unit = 1.0 / knot_num;
        real p = 0.0;
        for (int i = 0; i < knot_num + 1; ++i, p += unit) {
            knot_.push_back(p);
        }
        assert(knot_.size() == knot_num + 1); // n + k + 2
    }

    CurvePoint pointAtParam(real t) {
        CurvePoint pnt;
        Vector3f tangent, vertex;
        tangent = vertex = Vector3f::ZERO;
        int base_num = knot_.size() - 1; // n + k + 1
        vector<real> spline_base(base_num + 1);
        for (int i = 0; i < base_num; ++i) spline_base[i] = 0;
        int upper = std::upper_bound(knot_.begin(), knot_.end(), t) - knot_.begin();
        // cerr << knot_[upper] << " " << t << " " << knot_[upper - 1] << endl;
        if (upper != 0) spline_base[upper - 1] = 1;
        else spline_base[0] = 1; // maybe wrong
        int p = 0;
        for (p = 1; p <= deg_ - 1; ++p) {
            for (int i = 0; i < base_num - p; ++i) {
                spline_base[i] = 
                ((t - knot_[i]) * spline_base[i] / (knot_[i + p] - knot_[i])) + 
                ((knot_[i + p + 1] - t) * spline_base[i + 1] / (knot_[i + p + 1] - knot_[i + 1]));
                // cerr << p << " " << i << " " << spline_base[i] << endl;
            }
        } // B(i, k - 1) i from 0 to (n + 1)
        for (int i = 0; i < controls.size(); ++i) {
            tangent += deg_ * (
                spline_base[i] / (knot_[i + deg_] - knot_[i]) - 
                spline_base[i + 1] / (knot_[i + deg_ + 1] - knot_[i + 1])
            ) * controls[i];
        }
        tangent.normalize();
        for (int i = 0, p = deg_; i < base_num - p; ++i) {
            spline_base[i] = 
            ((t - knot_[i]) * spline_base[i] / (knot_[i + p] - knot_[i])) + 
            ((knot_[i + p + 1] - t) * spline_base[i + 1] / (knot_[i + p + 1] - knot_[i + 1]));
        }
        for (int i = 0; i < base_num - deg_; ++i) {
            vertex += spline_base[i] * controls[i];
        }
        pnt.V = vertex, pnt.T = tangent;
        return pnt;
    }


    void discretize(int resolution, std::vector<CurvePoint>& data) override {
        data.clear();
        real start = knot_[deg_];
        // cerr << deg_ << " " << knot_.size() - deg_ - 1 << endl;
        real end = knot_[knot_.size() - deg_ - 1]; // n + 1
        real unit = (knot_[1] - knot_[0]) / resolution;
        // cerr << start << " " << end << endl;
        for (real t = start; t <= end + kEps; t += unit) {
            data.push_back(pointAtParam(t));
        }
        // TODO (PA3): fill in data vector
    }

protected:
    std::vector<real> knot_; // size = n + k + 2
    int deg_; // k
};

#endif // CURVE_HPP
