#pragma once

#include <image.hpp>
#include <algorithm>
#include <queue>
#include <stack>
#include <cstdio>

class Element {
public:
    virtual void draw(Image &img) = 0;
    virtual ~Element() = default;
};

class Line : public Element {

public:
    int xA, yA;
    int xB, yB;
    Vector3f color;
    void draw(Image &img) override {
        // Draw from left to right of the long side
        int long_side = abs(xB - xA) < abs(yB - yA);
        int point_matrix[2][2] = {{xA, xB}, {yA, yB}};

        if (point_matrix[long_side][0] > point_matrix[long_side][1]) {
            std::swap(point_matrix[0][0], point_matrix[0][1]);
            std::swap(point_matrix[1][0], point_matrix[1][1]);
        } // Then point_matrix will be from small long_side to big long_side

        // Use pointer to handle long_side and short_side
        int x = point_matrix[0][0], y = point_matrix[1][0];
        int *long_side_ptr = long_side == LONG_SIDE_X?&x:&y;
        int *short_side_ptr = long_side == LONG_SIDE_X?&y:&x;
        if (!long_side_ptr || !short_side_ptr) return;

        // For loop using Bresenham Algorithm
        int d_long_side = point_matrix[long_side][1] - point_matrix[long_side][0];
        int d_short_side = point_matrix[!long_side][1] - point_matrix[!long_side][0];
        int step_short_side = (d_short_side > 0)?1:-1;
        int abs_d_short_side = d_short_side * step_short_side;
        int e = -d_long_side;
        for (int i = 0; i <= d_long_side; ++i) {
            img.SetPixel(x, y, color);
            ++(*long_side_ptr), e += 2 * abs_d_short_side;
            if (e >= 0) {
                (*short_side_ptr) += step_short_side, e -= 2 * d_long_side;
            }
        }
        printf("Draw a line from (%d, %d) to (%d, %d) using color (%f, %f, %f)\n", xA, yA, xB, yB,
                color.x(), color.y(), color.z());
    }
private:
    enum {LONG_SIDE_X, LONG_SIDE_Y};
};

class Circle : public Element {

public:
    int cx, cy;
    int radius;
    Vector3f color;
    void draw(Image &img) override {
        int x = cx, y = cy + radius, d = 5 - 4 * radius;

        // x - cx will be the new x used by mid_point_circle
        // y - cy will be the new y used by mid_point_circle
        int dc = cy - cx; 
        int d_negative = 12 - 8 * cx, d_positive = 8 * (cy - cx) + 20;
        
        //Mid point circle
        for (drawCirclePoints(x, y, img); x + dc <= y; 
             ++x, drawCirclePoints(x, y, img)) {
            if (d < 0) {
                d += 8 * x + d_negative;
            } else {
                d += 8 * (x - y) + d_positive, --y;
            }
        }
        printf("Draw a circle with center (%d, %d) and radius %d using color (%f, %f, %f)\n", cx, cy, radius,
               color.x(), color.y(), color.z());
    }
private:
    inline void drawCirclePoints(int x, int y, Image &img) {
        int dx = x - cx, dy = y - cy;
        static const int kCircSymCoefficient[4][2] = {{1, 1}, {1, -1}, {-1, -1}, {-1, 1}};
        for (int i = 0; i < 4; ++i) {
            img.SetPixel(kCircSymCoefficient[i][0] * dx + cx, 
                         kCircSymCoefficient[i][1] * dy + cy, color);
            img.SetPixel(kCircSymCoefficient[i][0] * dy + cx, 
                         kCircSymCoefficient[i][1] * dx + cy, color);
        }
    }
};

class Fill : public Element {

public:
    int cx, cy;
    Vector3f color;
    using Point = std::pair<int, int>;
    void draw(Image &img) override {
        Vector3f old_color = img.GetPixel(cx, cy);
        if (old_color == color) return;
        int xl, xr;
        bool span_need_fill;
        std::stack<Point> point_stack;
        while (!point_stack.empty()) point_stack.pop();
        point_stack.push(Point(cx, cy));
        int w = img.Width(), h = img.Height();
        while (!point_stack.empty()) {
            Point pt = point_stack.top();
            point_stack.pop();
            int x = 0, y = pt.second;
            for (x = pt.first; x < w && img.GetPixel(x, y) == old_color; ++x) {
                img.SetPixel(x, y, color);
            }
            xr = x - 1;
            for (x = pt.first - 1; x >= 0 && img.GetPixel(x, y) == old_color; --x) {
                img.SetPixel(x, y, color);
            }
            xl = x + 1;
            int y_arr[2] = {y + 1, y - 1};
            for (int i = 0; i < 2; ++i) {
                if (y_arr[i] < 0 || y_arr[i] >= h) continue;
                for (x = xl, y = y_arr[i]; x <= xr;) {
                    span_need_fill = false;
                    while (x <= xr && img.GetPixel(x, y) == old_color) {
                        span_need_fill = true, ++x;
                    }
                    if (span_need_fill) {
                        point_stack.push(Point(x - 1, y));
                        span_need_fill = false;
                    }
                    while (x <= xr && img.GetPixel(x, y) != old_color) ++x;
                }
            }
        }
        printf("Flood fill source point = (%d, %d) using color (%f, %f, %f)\n", cx, cy,
                color.x(), color.y(), color.z());
    }
};