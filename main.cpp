#include <string>
#include <print>
#include <functional>

#include "clothoid.h"

auto lerp(auto a, auto b, auto t) {
    return a + t * (b - a);
}

// classes that represent shapes to blend conform to the clothoid::Shape concept (for get and len)
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
    auto arc1 = Arc(0.5, 0, M_PI/2, 0, 0.5);
    auto arc2 = Arc(1, M_PI/2, M_PI, 0, 0);

    auto t = 0.0;
    auto dt = 1e-3;

    // print first arc
    t = 0;
    while(t < 1.0) {
        auto p = arc1.get(t);
        std::println("{} {}", p.x, p.y);
        t += dt;
    }
    std::println();

    // print second arc
    t = 0;
    while(t < 1.0) {
        auto p = arc2.get(t);
        std::println("{} {}", p.x, p.y);
        t += dt;
    }
    std::println();

    auto [s1, s2, c1, c2] = clothoid::fit_biclothoid(arc1, arc2, 0.5, 1e-9);

    std::println(stderr, "s1: {}, s2: {}, c1: {}, c2: {}", s1, s2, c1, c2);

    auto [x0, y0, k0, phi_0] = arc1.get(0.5);

    // print first half of bi-clothoid
    t = 0;
    while(t < 1.0) {
        auto s = lerp(0.0, s1, t);
        auto [x, y] = clothoid::xy(x0, y0, phi_0, k0, c1, s);
        std::println("{} {}", x, y);
        t += dt;
    }
    std::println();

    auto [xc, yc] = clothoid::xy(x0, y0, phi_0, k0, c1, s1);
    auto kc = clothoid::k(k0, c1, s1);
    auto phi_c = clothoid::phi(phi_0, k0, c1, s1);

    // print second half
    t = 0;
    while(t < 1.0) {
        auto s = lerp(0.0, s2, t);
        auto [x, y] = clothoid::xy(xc, yc, phi_c, kc, c2, s);
        std::println("{} {}", x, y);
        t += dt;
    }

    return 0;
}
