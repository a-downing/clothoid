#ifndef CLOTHOID_H
#define CLOTHOID_H

#include <cmath>

#include <Eigen/Eigen>

#include <common.h>
#include <gauss-legendre.h>
#include <shape.h>
#include <vec2.h>

namespace clothoid
{
    template<typename T>
    T k(T k0, T c, T s) {
        return k0 + c*s;
    }

    template<typename T>
    T phi(T phi_0, T k0, T c, T s) {
        return phi_0 + k0*s + 0.5*c*(s*s);
    }

    template<typename T>
    T x_prime(T phi_0, T k0, T c, T s) {
        return std::cos(phi_0 + k0*s + 0.5*c*s*s);
    }

    template<typename T>
    T y_prime(T phi_0, T k0, T c, T s) {
        return std::sin(phi_0 + k0*s + 0.5*c*s*s);
    }

    template<typename T>
    T x(T x0, T phi_0, T k0, T c, T s) {
        auto [xs, ws] = gauss_legendre::map_xw(gauss_legendre::P9_X, gauss_legendre::P9_W, T(0), s);
        T x = T(0);

        for(int i = 0; i < xs.size(); i++) {
            x += ws[i] * x_prime(phi_0, k0, c, xs[i]);
        }

        return x0 + x;
    }

    template<typename T>
    T y(T y0, T phi_0, T k0, T c, T s) {
        auto [xs, ws] = gauss_legendre::map_xw(gauss_legendre::P9_X, gauss_legendre::P9_W, T(0), s);
        T _y = T(0);

        for(int i = 0; i < xs.size(); i++) {
            y += ws[i] * y_prime(phi_0, k0, c, xs[i]);
        }

        return y0 + y;
    }

    template<typename T>
    std::array<T, 2> xy(T x0, T y0, T phi_0, T k0, T c, T s) {
        auto [xs, ws] = gauss_legendre::map_xw(gauss_legendre::P9_X, gauss_legendre::P9_W, T(0), s);
        T x = T(0);
        T y = T(0);

        for(int i = 0; i < xs.size(); i++) {
            x += ws[i] * x_prime(phi_0, k0, c, xs[i]);
            y += ws[i] * y_prime(phi_0, k0, c, xs[i]);
        }

        return { x0 + x, y0 + y };
    }

    template<typename T>
    struct biclothoid_t {
        T x0;
        T y0;
        T phi_0;
        T k0;
        T s1;
        T s2;
        T c1;

        void print() const {
            std::fprintf(stderr, "x0: %g\n", x0);
            std::fprintf(stderr, "y0: %g\n", y0);
            std::fprintf(stderr, "phi_0: %g\n", phi_0);
            std::fprintf(stderr, "k0: %g\n", k0);
            std::fprintf(stderr, "s1: %g\n", s1);
            std::fprintf(stderr, "s2: %g\n", s2);
            std::fprintf(stderr, "c1: %g\n", c1);
        }

        T c2() const {
            return -c1;
        }

        T len() const {
            return s1 + s2;
        }

        vec2<T> p(T s) const {
            if(s < s1) {
                auto [x, y] = clothoid::xy(x0, y0, phi_0, k0, c1, s);
                return vec2<T>(x, y);
            }

            auto [xc, yc] = clothoid::xy(x0, y0, phi_0, k0, c1, s1);
            auto phi_c = clothoid::phi(phi_0, k0, c1, s1);
            auto kc = clothoid::k(k0, c1, s1);

            auto [x, y] = clothoid::xy(xc, yc, phi_c, kc, c2(), s - s1);
            return vec2<T>(x, y);
        }
    };

    template<typename T>
    struct fit_result_t {
        biclothoid_t<T> bc;
        T t1;
        T t2;
    };

    template<typename T>
    std::array<T, 4> solve_biclothoid(T phi_delta, T k0, T ke, T sbc) {
        T s1, s2, c1, c2;
        bool ccw = phi_delta > 0 ? true : false;

        std::fprintf(stderr, "solve_biclothoid: phi_delta: %gpi, k0: %g, ke: %g, sbc: %g\n", phi_delta/M_PI, k0, ke, sbc);
        phi_delta = std::abs(phi_delta);

        if(ke == k0) {
            s1 = sbc / 2;
            s2 = sbc - s1;

            // I had to figure this out myself (in clothoid.sage). it seems like it should be simpler though
            //c1 = -2*(k0*sbc + phi_0 - phi_e - sqrt(std::pow(k0, 2)*std::pow(sbc, 2) + std::pow(phi_0, 2) - 2*phi_0*phi_e + std::pow(phi_e, 2) + 2*(k0*phi_0 - k0*phi_e)*sbc))/std::pow(sbc, 2);

            // simplifies with phi_0 being zero and phi_e being phi_delta
            c1 = -2*(k0*sbc - phi_delta - sqrt(std::pow(k0, 2)*std::pow(sbc, 2) + std::pow(phi_delta, 2) + 2*(-k0*phi_delta)*sbc))/std::pow(sbc, 2);
        } else {
            // original equations from paper
            auto SQ = std::pow(phi_delta, 2) - phi_delta*sbc*(k0 + ke) + (std::pow(sbc, 2))/2*(std::pow(k0, 2) + std::pow(ke, 2));
            // equation 11 from the paper
            auto s1a = (phi_delta - ke*sbc - sqrt(SQ))/(k0 - ke);
            auto s1b = (phi_delta - ke*sbc + sqrt(SQ))/(k0 - ke);

            std::fprintf(stderr, "solve_biclothoid: s1a: %.17g, s1b: %.17g\n", s1a, s1b);

            s1 = s1a;

            if(s1 < 0) {
                s1 = (phi_delta - ke*sbc + sqrt(SQ))/(k0 - ke);
            }

            if(is_close(s1, sbc)) {
                s1 = sbc;
            }

            s2 = sbc - s1;
            c1 = (ke - k0)/(s1 - s2);
        }

        if(!ccw) {
            c1 = -c1;
        }

        c2 = -c1;

        std::fprintf(stderr, "solve_biclothoid: s1: %g, s2: %g, c1: %g, c2: %g\n", s1, s2, c1, c2);

        return { s1, s2, c1, c2 };
    }

    template<typename T>
    fit_result_t<T> fit_biclothoid(const Shape<T> &q1, const Shape<T> &q2, T t1, T eps, int max_iter = 100) {
        auto p0 = q1.point(t1);
        auto phi_0 = q1.phi(t1);
        auto k0 = q1.k(t1);
        auto sbc = 2.0; //(T(1) - t1) * q1.len() * 2; // this is tricky
        auto t2 = (1 - t1) * q1.len() / q2.len();

        std::fprintf(stderr, "fit_biclothoid: t1: %g, sbc: %g, t2: %g, p0.x: %g, p0.y: %g\n", t1, sbc, t2, p0.x, p0.y);

        T s1, s2, c1, c2;

        auto f1f2 = [&q1, &q2, &s1, &s2, &c1, &c2, &p0, &k0, &t1, &phi_0](T sbc, T t2) -> std::array<T, 2> {
            auto pe = q2.point(t2);
            auto ke = q2.k(t2);
            auto [_s1, _s2, _c1, _c2] = solve_biclothoid(Shape<T>::phi_delta(q1, q2, t1, t2), k0, ke, sbc);
            s1 = _s1;
            s2 = _s2;
            c1 = _c1;
            c2 = _c2;
            auto [xc, yc] = xy(p0.x, p0.y, phi_0, k0, c1, s1);
            auto kc = k(k0, c1, s1);
            auto phi_c = phi(phi_0, k0, c1, s1);
            auto [xe, ye] = xy(xc, yc, phi_c, kc, c2, s2);
            return { xe - pe.x, ye - pe.y };
        };

        for(int i = 0; i < max_iter; i++) {
            auto [f1, f2] = f1f2(sbc, t2);
            auto [f1_sbc_eps, f2_sbc_eps] = f1f2(sbc + eps, t2);
            auto [f1_t2_eps, f2_t2_eps] = f1f2(sbc, t2 + eps);

            std::fprintf(stderr, "fit_biclothoid loop: i: %d, sbc: %g, t2: %g, f1: %g, f2: %g\n", i, sbc, t2, f1, f2);

            auto df1_dsbc = (f1_sbc_eps - f1) / eps;
            auto df2_dsbc = (f2_sbc_eps - f2) / eps;
            auto df1_dt2 = (f1_t2_eps - f1) / eps;
            auto df2_dt2 = (f2_t2_eps - f2) / eps;

            std::fprintf(stderr, "fit_biclothoid loop: df1_dsbc: %g, df2_dsbc: %g, df1_dt2: %g, df2_dt2: %g\n", df1_dsbc, df2_dsbc, df1_dt2, df2_dt2);

            auto J = Eigen::Matrix<T, 2, 2>({
                {df1_dsbc, df1_dt2},
                {df2_dsbc, df2_dt2}
            });

            auto x = J.inverse() * Eigen::Matrix<T, 2, 1>({{f1}, {f2}});

            sbc -= x(0, 0);
            t2 -= x(1, 0);

            if(std::abs(f1) < eps && std::abs(f2) < eps) {
                return fit_result_t<T> { biclothoid_t<T> { p0.x, p0.y, phi_0, k0, s1, s2, c1 }, t1, t2 };
            }
        }

        return fit_result_t<T> { biclothoid_t<T> { p0.x, p0.y, phi_0, k0, s1, s2, c1 }, t1, t2 };
    }

    template<typename T>
    std::tuple<const Shape<T>&, T> lerp(const Shape<T> &q1, const Shape<T> &q2, const fit_result_t<T> &fr, T t) {
        auto len_q1 = q1.len() * (1.0 - fr.t1);
        auto len_q2 = q2.len() * fr.t2;
        auto len_tot = len_q1 + len_q2;
        auto len_t = lerp(0.0, len_tot, t);

        if(len_t > len_q1) {
            return { q2, (len_t - len_q1) / q2.len() };
        }

        return { q1, len_t / q1.len() + fr.t1 };
    }

    // this isn't good enough. the corner isn't always the part with the largest deviation
    template<typename T>
    std::array<T, 2> corner_deviation(const Shape<T> &q1, const biclothoid_t<T> &bc, int iter) {
        auto end = q1.point(1.0);
        auto corner = vec2(end.x, end.y);
        T s = 0.5;
        T ds = 0.25;
        T d;

        for(int i = 0; i < iter; i++) {
            auto dl = (corner - bc.p(s - ds)).len();
            auto dr = (corner - bc.p(s + ds)).len();

            if(dl < dr) {
                s = s - ds;
                d = dl;
            } else {
                s = s + ds;
                d = dr;
            }

            ds *= 0.5;
        }

        return { d, s };
    }

    template<typename T>
    struct deviation_result_t {
        T d;
        T s;
    };

    template<typename T>
    biclothoid_t<T> fit_biclothoid_tol(const Shape<T> &q1, const Shape<T> &q2, int iter, T tol, T tol_eps, T eps) {
        auto t = 1.0;
        auto [bc, _t1, _t2] = fit_biclothoid(q1, q2, T(1) - t, eps);
        auto [d, _s] = corner_deviation(q1, bc, iter);

        if(d < tol) {
            return bc;
        }

        auto alpha = 0.5;

        for(;;) {
            auto relerr = d / tol;
            t = t / (1.0 + alpha * (relerr - 1.0));
            auto [bc, _t1, _t2] = fit_biclothoid(q1, q2, T(1) - t, eps);
            auto [d, _s] = corner_deviation(q1, bc, iter);

            if(d < tol + tol_eps) {
                return bc;
            }
        }
    }
}

#endif //CLOTHOID_H
