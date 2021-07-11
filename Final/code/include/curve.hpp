#ifndef CURVE_HPP
#define CURVE_HPP

#include "object3d.hpp"
#include <vecmath.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <immintrin.h>

using namespace std;

struct CurvePoint {
    Vector3f V; // Vertex
    Vector3f T; // Tangent  (unit)
    static CurvePoint illegalPoint() {
        return CurvePoint{ Vector3f(NAN), Vector3f(NAN) };
    }
};

class Curve : public Object3D {
protected:
    std::vector<Vector3f> controls;
public:
    explicit Curve(std::vector<Vector3f> points) : controls(std::move(points)) {}

    bool intersect(const Ray &r, Hit &h, double tmin) override {
        return false;
    }

    std::vector<Vector3f> &getControls() {
        return controls;
    }

    virtual double start() = 0;

    virtual double end() = 0;

    virtual CurvePoint pointAtParam(double t) = 0;
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

    double combination(int b, int a) {
        // return C_b^a
        return static_cast<double>(fact(b)) / (fact(a) * fact(b - a));
    }

    double start() override {
        return 0.0;
    }

    double end() override {
        return 1.0;
    }

    CurvePoint pointAtParam(double t) override {
        // illegal
        if (isnan(t) || t < 0 || t > 1) return CurvePoint::illegalPoint();
        Vector3f tangent, vertex;
        tangent = vertex = Vector3f::ZERO;
        int n = controls.size() - 1;
        static std::vector<double> pow_1_t(n + 1), pow_t(n + 1);
        pow_1_t[0] = 1.0, pow_t[0] = 1.0;
        for (int i = 0; i < n; ++i) {
            pow_1_t[i + 1] = pow_1_t[i] * (1 - t);
            pow_t[i + 1] = pow_t[i] * t;
        }
        for (int i = 0; i <= n; ++i) {
            vertex += comb_n_[i] * pow_1_t[n - i] * pow_t[i] * controls[i];
        }
        tangent -= n * comb_n_minus_1_[0] * pow_1_t[n - 1] * controls[0];
        for (int i = 1; i <= n - 1; ++i) {
            tangent += n * (
                pow_1_t[n - 1 - i] * pow_t[i - 1] * 
                (comb_n_minus_1_[i - 1] * (1 - t) - comb_n_minus_1_[i] * t)
                // comb_n_minus_1_[i - 1] * pow(1 - t, n - i) * pow(t, i - 1) - 
                // comb_n_minus_1_[i] * pow(1 - t, n - 1 - i) * pow(t, i)
            ) * controls[i];
        }
        tangent += n * comb_n_minus_1_[n - 1] * pow_t[n - 1] * controls[n];
        tangent.normalize();
        CurvePoint p;
        p.T = tangent, p.V = vertex;
        return p;
    }

protected:
    std::vector<double> comb_n_;
    std::vector<double> comb_n_minus_1_;
};

class BsplineCurve : public Curve {
public:
    BsplineCurve(const std::vector<Vector3f> &points, int deg = 3) : Curve(points), deg_(deg) {
        if (points.size() < 4) {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
        int knot_num = points.size() + deg; // n + k + 1
        double unit = 1.0 / knot_num;
        double p = 0.0;

        knot_size_ = knot_num + 1;
        knot_ = new double[knot_size_];
        for (int i = 0; i < knot_size_; ++i, p += unit) {
            knot_[i] = p;
        }

        // pre calculate rev of knot minus to avoid division
        for (int p = 0; p <= deg; ++p) {
            rev_knot_minus_.push_back(new double[knot_size_ - p]);
            for (int i = 0; i <= knot_size_ - 1 - p; ++i) {
                rev_knot_minus_[p][i] = 1 / (knot_[i + p] - knot_[i]);
            }
        }
    }

    double start() override {
        return knot_[deg_];
    }

    double end() override {
        return knot_[knot_size_ - deg_ - 1];
        // return knot_[knot_.size() - deg_ - 1];
    }

    CurvePoint pointAtParam(double t) override {
#ifndef __AVX2__
        // Naive version
        double start = knot_[deg_], end = knot_[knot_size_ - deg_ - 1];
        if (isnan(t) || t < start || t > end) return CurvePoint::illegalPoint();
        CurvePoint pnt;
        Vector3f tangent, vertex;
        tangent = vertex = Vector3f::ZERO;
        int base_num = knot_size_ - 1; // n + k + 1
        vector<double> spline_base(base_num + 1);
        for (int i = 0; i < base_num; ++i) spline_base[i] = 0;
        int upper = std::upper_bound(knot_, knot_ + knot_size_, t) - knot_;
        // cerr << knot_[upper] << " " << t << " " << knot_[upper - 1] << endl;
        int upper_1 = upper - 1;
        if (upper != 0) spline_base[upper_1] = 1;
        else spline_base[0] = 1; // maybe wrong
        int p = 0;
        for (p = 1; p <= deg_ - 1; ++p) {
            for (int i = upper_1 - p; i <= upper_1; ++i) {
                spline_base[i] = 
                ((t - knot_[i]) * spline_base[i] * rev_knot_minus_[p][i]) + 
                ((knot_[i + p + 1] - t) * spline_base[i + 1] * rev_knot_minus_[p][i + 1]); // (i + 1 + p) - (i + 1)
            }
        } // B(i, k - 1) i from 0 to (n + 1)

        // pre-calc spline_base[i] / (knot_[i + deg_] - knot_[i])
        vector<double> base_knot(controls.size() + 1);
        for (int i = upper_1 - deg_; i <= upper_1; ++i) {
            base_knot[i] = spline_base[i] * rev_knot_minus_[deg_][i];
        }

        for (int i = upper_1 - p; i <= upper_1; ++i) {
            tangent += deg_ * (base_knot[i] - base_knot[i + 1]) * controls[i];
        }
        tangent.normalize();
        for (int i = upper_1 - deg_, p = deg_; i <= upper_1; ++i) {
            spline_base[i] = ((t - knot_[i]) * base_knot[i]) + ((knot_[i + p + 1] - t) * base_knot[i + 1]);
        }
        for (int i = upper_1 - deg_; i <= upper_1; ++i) {
            vertex += spline_base[i] * controls[i];
        }
#else
        // AVX version (assume deg = 3)
        double start = knot_[deg_], end = knot_[knot_size_ - deg_ - 1];
        if (isnan(t) || t < start || t > end) return CurvePoint::illegalPoint();
        CurvePoint pnt;
        int upper = std::upper_bound(knot_, knot_ + knot_size_, t) - knot_;
        int upper_1 = upper - 1;
        int upper_1_3 = upper_1 - 3; // array start
        __m256d t_s = _mm256_set1_pd(t);
        __m256d degs = _mm256_set1_pd(deg_);
        __m256d t_knot_i = _mm256_sub_pd(t_s, _mm256_loadu_pd(knot_ + upper_1_3)); // t - knot_[i]
        __m256d knot_i_p_1_t[4], rev_knot_minus_i[4], rev_knot_minus_i_1[4];
        for (int i = 1; i <= 3; ++i) {
            knot_i_p_1_t[i] = _mm256_sub_pd(_mm256_loadu_pd(knot_ + upper_1_3 + 1 + i), t_s); // knot_[i + p + 1] - t
            rev_knot_minus_i[i] = _mm256_loadu_pd(rev_knot_minus_[i] + upper_1_3);
            rev_knot_minus_i_1[i] = _mm256_loadu_pd(rev_knot_minus_[i] + upper_1_3 + 1);
        }
        __m256d spline_base_simd = _mm256_set_pd(1, 0, 0, 0);
        __m256d spline_base_1 = _mm256_set_pd(0, 1, 0, 0); // spline_base[i] and spline_base[i + 1]
        __m256d spline_base_permute, spline_base_permute_shift;
        spline_base_simd = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(t_knot_i, spline_base_simd), rev_knot_minus_i[1]),
                                         _mm256_mul_pd(_mm256_mul_pd(knot_i_p_1_t[1], spline_base_1), rev_knot_minus_i_1[1]));
        // x1 x2 x3 x4 => x2 x3 x4 0
        spline_base_permute = _mm256_permute_pd(spline_base_simd, 0b0101); // x2 x1 x4 x3
        spline_base_permute_shift = _mm256_permute2f128_pd(spline_base_permute, spline_base_permute, 0x81); // x4 x3 0 0
        spline_base_1 = _mm256_blend_pd(spline_base_permute, spline_base_permute_shift, 0b1010); // x2 x3 x4 0

        spline_base_simd = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(t_knot_i, spline_base_simd), rev_knot_minus_i[2]),
                                         _mm256_mul_pd(_mm256_mul_pd(knot_i_p_1_t[2], spline_base_1), rev_knot_minus_i_1[2]));

        spline_base_permute = _mm256_permute_pd(spline_base_simd, 0b0101); // x2 x1 x4 x3
        spline_base_permute_shift = _mm256_permute2f128_pd(spline_base_permute, spline_base_permute, 0x81); // x4 x3 0 0
        spline_base_1 = _mm256_blend_pd(spline_base_permute, spline_base_permute_shift, 0b1010); // x2 x3 x4 0

        __m256d base_knot_simd = _mm256_mul_pd(spline_base_simd, rev_knot_minus_i[3]);
        __m256d base_knot_1 = _mm256_mul_pd(spline_base_1, rev_knot_minus_i[3]);
        __m256d deg_base_knot = _mm256_mul_pd(degs, _mm256_sub_pd(base_knot_simd, base_knot_1));
        spline_base_simd = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(t_knot_i, spline_base_simd), rev_knot_minus_i[3]),
                                         _mm256_mul_pd(_mm256_mul_pd(knot_i_p_1_t[3], spline_base_1), rev_knot_minus_i_1[3]));
        __m256d controls_simd[3], tangent_simd[3], vertex_simd[3];
        for (int i = 0; i < 3; ++i) {
            controls_simd[i] = _mm256_set_pd(controls[upper_1_3 + 3][i], controls[upper_1_3 + 2][i], 
                                             controls[upper_1_3 + 1][i], controls[upper_1_3][i]);
            tangent_simd[i] = _mm256_mul_pd(controls_simd[i], deg_base_knot);
            vertex_simd[i] = _mm256_mul_pd(controls_simd[i], spline_base_simd);
        }
  
        // to add all the vectors to one vector
        // {c[0] + c[1], d[0] + d[1], c[2] + c[3], d[2] + d[3]}
        // {a[0] + a[1], b[0] + b[1], a[2] + a[3], b[2] + b[3]}
        __m256d tangent_x_y = _mm256_hadd_pd(tangent_simd[0], tangent_simd[1]);
        __m256d tangent_vertex_z = _mm256_hadd_pd(tangent_simd[2], vertex_simd[2]);
        __m256d vertex_x_y = _mm256_hadd_pd(vertex_simd[0], vertex_simd[1]);

        // add t_xyz and v_x
        // {a[0] + a[1], b[0] + b[1], c[2] + c[3], d[2] + d[3]}
        __m256d t_v_blend = _mm256_blend_pd(tangent_x_y, tangent_vertex_z, 0b1100);
        // {a[2] + a[3], b[2] + b[3], c[0] + c[1], d[0] + d[1]}
        __m256d t_v_perm = _mm256_permute2f128_pd(tangent_x_y, tangent_vertex_z, 0x21);
        __m256d t_v_sum = _mm256_add_pd(t_v_blend, t_v_perm);
        // add v_xy

        __m128d v_xy_high = _mm256_extractf128_pd(vertex_x_y, 1);
        __m128d v_xy_sum = _mm_add_pd(v_xy_high, _mm256_castpd256_pd128(vertex_x_y));

        double txyz_vz[4], vxy[2];
        _mm256_store_pd(txyz_vz, t_v_sum);
        _mm_store_pd(vxy, v_xy_sum);
        Vector3f vertex(vxy[0], vxy[1], txyz_vz[3]);

        // normalize tangent
        __m256d tangent_xyz = _mm256_blend_pd(t_v_sum, _mm256_setzero_pd(), 0b1000);
        __m256d sq_tan = _mm256_mul_pd(tangent_xyz, tangent_xyz);
        __m256d sum_tan = _mm256_hadd_pd(sq_tan, sq_tan); // {x^2 + y^2, x^2 + y^2, z, z}
        __m256d perm_sum_tan = _mm256_permute2f128_pd(sum_tan, sum_tan, 0x01);
        __m256d tan_len = _mm256_sqrt_pd(_mm256_add_pd(sum_tan, perm_sum_tan));
        __m256d norm_tan = _mm256_div_pd(tangent_xyz, tan_len);
        _mm256_store_pd(txyz_vz, norm_tan);
        Vector3f tangent(txyz_vz[0], txyz_vz[1], txyz_vz[2]);
#endif
        pnt.V = vertex, pnt.T = tangent;
        return pnt;
    }

protected:
    double *knot_; // size = n + k + 2
    int knot_size_;
    std::vector<double *> rev_knot_minus_;
    int deg_; // k
};

#endif // CURVE_HPP
