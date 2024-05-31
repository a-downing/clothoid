#include <string>
#include <print>
#include <functional>

#include "clothoid.h"

template<typename T>
T lerp(T a, T b, T t) {
    return a + t * (b - a);
}

// classes that represent shapes to blend conform to the clothoid::Shape concept (for get and len)
class Line {
    double x0, y0, x1, y1;

public:
    Line(double x0, double y0, double x1, double y1): x0(x0), y0(y0), x1(x1), y1(y1) {

    }

    clothoid::shape_point_t<double> get(double t) const {
        auto dx = x1 - x0;
        auto dy = y1 - y0;
        auto x = x0 + t*dx;
        auto y = y0 + t*dy;
        auto phi = std::atan2(dy, dx);

        if(phi < 0) {
            phi += 2*M_PI;
        }

        return { x, y, 0,  phi };
    }

    double len() const {
        auto dx = x1 - x0;
        auto dy = y1 - y0;
        return std::sqrt(dx*dx + dy*dy);
    }
};

class Arc {
    double r;
    double theta_0;
    double theta_e;
    double xc;
    double yc;

public:
    Arc(double r, double theta_0, double theta_e, double xc, double yc): r(r), theta_0(theta_0), theta_e(theta_e), xc(xc), yc(yc) {

    }

    clothoid::shape_point_t<double> get(double t) const {
        auto theta = lerp(theta_0, theta_e, t);
        auto x = xc + r * cos(theta);
        auto y = yc + r * sin(theta);
        auto phi_t = theta + M_PI/2;
        auto k = 1.0 / r;
        return { x, y, k, phi_t };
    }

    double len() const {
        return r * (theta_e - theta_0);
    }
};

// just prints out x and y coordinates
// a blank line separates shapes for plot.py
int main() {
    auto shape1 = Arc(0.5, 0, M_PI/2, 0, 0.5);
    //auto shape2 = Arc(1, M_PI/2, M_PI, 0, 0);
    auto shape2 = Line(0, 1, -0.5, 0);

    auto t = 0.0;
    auto dt = 1e-3;

    // print first arc
    t = 0;
    while(t < 1.0) {
        auto p = shape1.get(t);
        std::printf("%g %g\n", p.x, p.y);
        t += dt;
    }
    std::printf("\n");

    // print second arc
    t = 0;
    while(t < 1.0) {
        auto p = shape2.get(t);
        std::printf("%g %g\n", p.x, p.y);
        t += dt;
    }
    std::printf("\n");

    auto bc = clothoid::fit_biclothoid(shape1, shape2, 0.9, 1e-8);
    //auto bc = clothoid::fit_biclothoid_tol(shape1, shape2, 10, 0.01, 1e-3, 1e-8);

    auto [d, s] = clothoid::corner_deviation(shape1, bc, 10);

    // print the farthest point
    clothoid::vec2<double> p = { 0, 0 };
    p = bc.p(s);
    std::printf("%g %g\n", p.x, p.y);
    std::printf("\n");

    // print first half of bi-clothoid
    t = 0;
    while(t < 1.0) {
        auto s = lerp(0.0, bc.s1, t);
        auto [x, y] = bc.xy(s);
        std::printf("%g %g\n", x, y);
        t += dt;
    }
    std::printf("\n");

    // print second half
    t = 0;
    while(t < 1.0) {
        auto s = lerp(bc.s1, bc.len(), t);
        auto [x, y] = bc.xy(s);
        std::printf("%g %g\n", x, y);
        t += dt;
    }
    std::printf("\n");

    // print the curve midpoint
    p = bc.p(bc.len() / 2);
    //std::printf("%g %g\n", p.x, p.y);
    //std::printf("\n");

    // print the corner
    auto end = shape1.get(1.0);
    auto corner = clothoid::vec2(end.x, end.y);
    //std::printf("%g %g\n", corner.x, corner.y);
    //std::printf("\n");

    return 0;
}
