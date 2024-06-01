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
    clothoid::vec2<double> p0;
    clothoid::vec2<double> p1;

public:
    Line(const clothoid::vec2<double> &p0, const clothoid::vec2<double> &p1): p0(p0), p1(p1) {

    }

    const clothoid::vec2<double> point(double t) const {
        return p0 + t*(p1 - p0);
    }

    double k(double t) const {
        return 0.0;
    }

    double phi(double t) const {
        auto v = p1 - p0;
        auto phi = std::atan2(v.y, v.x);

        if(phi < 0) {
            phi += 2*M_PI;
        }

        return phi;
    }

    double len() const {
        return (p1 - p0).len();
    }
};

class Arc {
    double r;
    double theta_0;
    double rads;
    clothoid::vec2<double> c;

public:
    // an arc defined by radius, center, starting angle, and positive or negative angle to the end (rads)
    Arc(double r, double theta_0, double rads, const clothoid::vec2<double> &c): r(r), theta_0(theta_0), rads(rads), c(c) {

    }

    double theta(double t) const {
        auto theta = std::fmod(theta_0 + lerp(0.0, rads, t), 2*M_PI);

        if(theta < 0) {
            theta += 2*M_PI;
        }

        return theta;
    }

    const clothoid::vec2<double> point(double t) const {
        return { c.x + r * std::cos(theta(t)), c.y + r * std::sin(theta(t)) };
    }

    double k(double t) const {
        return 1.0 / r;
    }

    double phi(double t) const {
        return std::fmod(theta(t) + M_PI/2, 2*M_PI);
    }

    double len() const {
        return r * std::abs(rads);
    }
};

// just prints out x and y coordinates
// a blank line separates shapes for plot.py
int main() {
    auto shape1 = Arc(0.5, -0.5, M_PI/2 + 0.5, clothoid::vec2(0.0, 0.5));
    //auto shape2 = Arc(1, M_PI/2, M_PI, clothoid::vec2(0.0, 0.0));
    auto shape2 = Line(clothoid::vec2(0.0, 1.0), clothoid::vec2(-0.5, 0.0));

    auto t = 0.0;
    auto dt = 1e-3;

    // print first arc
    t = 0;
    while(t < 1.0) {
        auto p = shape1.point(t);
        std::printf("%g %g\n", p.x, p.y);
        t += dt;
    }
    std::printf("\n");

    // print second arc
    t = 0;
    while(t < 1.0) {
        auto p = shape2.point(t);
        std::printf("%g %g\n", p.x, p.y);
        t += dt;
    }
    std::printf("\n");

    auto fr = clothoid::fit_biclothoid(shape1, shape2, 0.5, 1e-8);
    auto [bc, t1, t2] = fr;
    auto sp = bc.p(0);
    auto mp = bc.p(bc.s1);
    auto ep = bc.p(bc.len());
    // not really viable yet
    //auto bc = clothoid::fit_biclothoid_tol(shape1, shape2, 10, 0.01, 1e-3, 1e-8);

    // print the start, middle, and end of the biclothoid
    std::printf("%g %g\n", sp.x, sp.y);
    std::printf("\n");
    std::printf("%g %g\n", mp.x, mp.y);
    std::printf("\n");
    std::printf("%g %g\n", ep.x, ep.y);
    std::printf("\n");

    // print first half of bi-clothoid
    t = 0;
    while(t < 1.0) {
        auto s = lerp(0.0, bc.s1, t);
        auto p = bc.p(s);
        std::printf("%g %g\n", p.x, p.y);
        t += dt;
    }
    std::printf("\n");

    // print second half
    t = 0;
    while(t < 1.0) {
        auto s = lerp(bc.s1, bc.len(), t);
        auto p = bc.p(s);
        std::printf("%g %g\n", p.x, p.y);
        t += dt;
    }
    std::printf("\n");

    return 0;
}
