#include <string>

#include <common.h>
#include <clothoid.h>
#include <shape.h>

class Line : public clothoid::Shape<double> {
    clothoid::vec2<double> p0;
    clothoid::vec2<double> p1;

public:
    Line(const clothoid::vec2<double> &p0, const clothoid::vec2<double> &p1): p0(p0), p1(p1) {

    }

    clothoid::vec2<double> point(double t) const override {
        return p0 + t*(p1 - p0);
    }

    double k(double t) const override {
        return 0.0;
    }

    double phi(double t) const override {
        auto v = p1 - p0;
        auto phi = std::atan2(v.y, v.x);

        if(phi < 0) {
            phi += 2*M_PI;
        }

        return phi;
    }

    double phi_delta(double t) const override {
        return 0.0;
    }

    double len() const override {
        return (p1 - p0).len();
    }
};

class Arc : public clothoid::Shape<double> {
    double r;
    double theta_0;
    double rads;
    clothoid::vec2<double> c;

public:
    // an arc defined by radius, center, starting angle, and positive or negative angle to the end (rads)
    Arc(double r, double theta_0, double rads, const clothoid::vec2<double> &c): r(r), theta_0(theta_0), rads(rads), c(c) {

    }

    double theta(double t) const {
        auto theta = std::fmod(theta_0 + clothoid::lerp(0.0, rads, t), 2*M_PI);

        if(theta < 0) {
            theta += 2*M_PI;
        }

        return theta;
    }

    clothoid::vec2<double> point(double t) const override {
        return { c.x + r * std::cos(theta(t)), c.y + r * std::sin(theta(t)) };
    }

    double k(double t) const override {
        return 1.0 / r;
    }

    double phi(double t) const override {
        return std::fmod(theta(t) + std::copysign(M_PI/2, rads), 2*M_PI);
    }

    double phi_delta(double t) const override {
        return rads * t;
    }

    double len() const override {
        return r * std::abs(rads);
    }
};

template<typename T>
void print_shape(const clothoid::Shape<T> &q, int resolution = 1000) {
    for(int i = 0; i <= resolution; i++) {
        auto t = 1.0 / resolution * i;
        auto p = q.point(t);
        std::printf("%g %g\n", p.x, p.y);
    }
    std::printf("\n");
}

template<typename T>
void print_biclothoid_half(const clothoid::biclothoid_t<T> &bc, bool first_half, int resolution = 1000) {
    for(int i = 0; i <= resolution; i++) {
        auto s = first_half ? bc.s1 / resolution * i : bc.s1 + bc.s2 / resolution * i;
        auto p = bc.p(s);
        std::printf("%g %g\n", p.x, p.y);
    }
    std::printf("\n");
}

// just prints out x and y coordinates
// a blank line separates shapes for plot.py
int main() {
    //auto shape1 = Arc(0.5, -0.5, M_PI/2 + 0.5, clothoid::vec2(0.0, 0.5));
    //auto shape2 = Arc(1, M_PI/2, M_PI, clothoid::vec2(0.0, 0.0));
    //auto shape2 = Line(clothoid::vec2(0.0, 1.0), clothoid::vec2(-0.5, -1.0));

    // cw
    //auto shape1 = Line(clothoid::vec2(-10.0, 0.0), clothoid::vec2(0.0, 10.0));
    //auto shape2 = Line(clothoid::vec2(1.0, 10.0), clothoid::vec2(10.0, 0.0));

    // ccw
    //auto shape1 = Line(clothoid::vec2(10.0, 0.0), clothoid::vec2(1.0, 10.0));
    //auto shape2 = Line(clothoid::vec2(0.0, 10.0), clothoid::vec2(-10.0, 0.0));

    // Fails:
    //auto shape1 = Line(clothoid::vec2(50.0, 0.0), clothoid::vec2(-50.0, 0.0)); // this ends at -50, 0
    //auto shape2 = Arc(50, 0, M_PI, clothoid::vec2(0.0, 0.0)); // this starts at 50, 0
    // corrected paths, works
    //auto shape1 = Line(clothoid::vec2(50.0, 0.0), clothoid::vec2(-50.0, 0.0)); // this ends at -50, 0
    //auto shape2 = Arc(50, M_PI, -M_PI, clothoid::vec2(0.0, 0.0)); // this starts at -50, 0

    // Fails:
    //auto shape1 = Arc(50, M_PI, 2*M_PI, clothoid::vec2(0.0, 5.0)); // this ends at -50, 5
    //auto shape2 = Arc(50, 0, M_PI, clothoid::vec2(0.0, 0.0)); // this starts at 50, 0
    // corrected paths, but still fails unless the t1 param of fit_biclothoid is near 1
    auto shape1 = Arc(50.0, M_PI, 2*M_PI, clothoid::vec2(0.0, 5.0)); // this ends at -50, 5
    auto shape2 = Arc(50.0, M_PI, -M_PI, clothoid::vec2(0.0, 0.0)); // this starts at -50, 0

    // works now
    //auto shape1 = Line(clothoid::vec2(0.0, -100.0), clothoid::vec2(0.0, 0.0));
    //auto shape2 = Arc(50, M_PI, M_PI, clothoid::vec2(50.0, 0.0));

    // works now
    //auto shape1 = Line(clothoid::vec2(0.0, -100.0), clothoid::vec2(0.0, 0.0));
    //auto shape2 = Arc(50, M_PI, -M_PI, clothoid::vec2(50.0, 0.0));

    print_shape(shape1);
    print_shape(shape2);

    auto fr = clothoid::fit_biclothoid(shape1, shape2, 0.9, 1e-8);
    auto [bc, t1, t2] = fr;
    auto sp = bc.p(0);
    auto mp = bc.p(bc.s1);
    auto ep = bc.p(bc.len());

    // test interpolating across the fitted region of two shapes for the calculation on deviation later
    // constexpr int spaces = 16;
    // for(int i = 0; i <= spaces; i++) {
    //     auto t = 1.0 / spaces * i;
    //     auto [shape, _t] = clothoid::lerp(shape1, shape2, fr, t);
    //     auto p = shape.point(_t);
    //     std::printf("%g %g\n", p.x, p.y);
    //     std::printf("\n");
    // }

    //print the start, middle, and end of the biclothoid
    std::printf("%g %g\n", sp.x, sp.y);
    std::printf("\n");
    std::printf("%g %g\n", mp.x, mp.y);
    std::printf("\n");
    std::printf("%g %g\n", ep.x, ep.y);
    std::printf("\n");

    print_biclothoid_half(bc, true);
    print_biclothoid_half(bc, false);

    return 0;
}
