#ifndef CLOTHOID_H
#define CLOTHOID_H

#include <cmath>

#include <Eigen/Eigen>

#include "gauss-legendre.h"

namespace clothoid
{
    auto k(auto k0, auto c, auto s) {
        return k0 + c*s;
    }

    auto phi(auto phi_0, auto k0, auto c, auto s) {
        return phi_0 + k0*s + 0.5*c*(s*s);
    }

    auto x_prime(auto phi_0, auto k0, auto c, auto s) {
        return std::cos(phi_0 + k0*s + 0.5*c*s*s);
    }

    auto y_prime(auto phi_0, auto k0, auto c, auto s) {
        return std::sin(phi_0 + k0*s + 0.5*c*s*s);
    }

    template<typename T>
    T x(T x0, T phi_0, T k0, T c, T s) {
        auto [xs, ws] = gauss_legendre::map_xw(gauss_legendre::P5_X, gauss_legendre::P5_W, T(0), s);
        T x = T(0);

        for(int i = 0; i < xs.size(); i++) {
            x += ws[i] * x_prime(phi_0, k0, c, xs[i]);
        }

        return x0 + x;
    }

    template<typename T>
    T y(T y0, T phi_0, T k0, T c, T s) {
        auto [xs, ws] = gauss_legendre::map_xw(gauss_legendre::P5_X, gauss_legendre::P5_W, T(0), s);
        T _y = T(0);

        for(int i = 0; i < xs.size(); i++) {
            y += ws[i] * y_prime(phi_0, k0, c, xs[i]);
        }

        return y0 + y;
    }

    template<typename T>
    std::array<T, 2> xy(T x0, T y0, T phi_0, T k0, T c, T s) {
        auto [xs, ws] = gauss_legendre::map_xw(gauss_legendre::P5_X, gauss_legendre::P5_W, T(0), s);
        T x = T(0);
        T y = T(0);

        for(int i = 0; i < xs.size(); i++) {
            x += ws[i] * x_prime(phi_0, k0, c, xs[i]);
            y += ws[i] * y_prime(phi_0, k0, c, xs[i]);
        }

        return { x0 + x, y0 + y };
    }

    template<typename T>
    std::array<T, 4> solve_biclothoid(T phi_0, T k0, T phi_e, T ke, T sbc) {
        // I need to look more into the math here
        T s1, s2, c1, c2;

        std::println(stderr, "solve_biclothoid: phi_0: {}, k0: {}, phi_e: {}, ke: {}, sbc: {}", phi_0, k0, phi_e, ke, sbc);

        if(ke == k0) {
            s1 = sbc / 2;
            s2 = sbc - s1;
            c1 = -2*(k0*sbc + phi_0 - phi_e - sqrt(std::pow(k0, 2)*std::pow(sbc, 2) + std::pow(phi_0, 2) - 2*phi_0*phi_e + std::pow(phi_e, 2) + 2*(k0*phi_0 - k0*phi_e)*sbc))/std::pow(sbc, 2);
        } else {
            // // original equations from paper
            auto phi_delta = phi_e - phi_0;
            auto SQ = std::pow(phi_delta, 2) - phi_delta*sbc*(k0 + ke) + (std::pow(sbc, 2))/2*(std::pow(k0, 2) + std::pow(ke, 2));
            // equation 11 from the paper
            s1 = (phi_delta - ke*sbc - sqrt(SQ))/(k0 - ke);

            if(s1 < 0) {
                s1 = (phi_delta - ke*sbc + sqrt(SQ))/(k0 - ke);
            }

            s2 = sbc - s1;
            c1 = (ke - k0)/(s1 - s2);

            //std::println(stderr, "s1: {}", s1);

            // auto SQ = sqrt(2*(std::pow(k0, 2) + std::pow(ke, 2))*std::pow(sbc, 2) + 4*std::pow(phi_0, 2) - 8*phi_0*phi_e + 4*std::pow(phi_e, 2) + 4*((k0 + ke)*phi_0 - (k0 + ke)*phi_e)*sbc);
            // auto s1a = -0.5*(2*ke*sbc + 2*phi_0 - 2*phi_e - SQ)/(k0 - ke);
            // auto c1a = -((k0 + ke)*sbc + 2*phi_0 - 2*phi_e + sqrt(2*(std::pow(k0, 2) + std::pow(ke, 2))*std::pow(sbc, 2) + 4*std::pow(phi_0, 2) - 8*phi_0*phi_e + 4*std::pow(phi_e, 2) + 4*((k0 + ke)*phi_0 - (k0 + ke)*phi_e)*sbc))/std::pow(sbc, 2);
            // auto s1b = -0.5*(2*ke*sbc + 2*phi_0 - 2*phi_e + SQ)/(k0 - ke);
            // auto c1b = -((k0 + ke)*sbc + 2*phi_0 - 2*phi_e - sqrt(2*(std::pow(k0, 2) + std::pow(ke, 2))*std::pow(sbc, 2) + 4*std::pow(phi_0, 2) - 8*phi_0*phi_e + 4*std::pow(phi_e, 2) + 4*((k0 + ke)*phi_0 - (k0 + ke)*phi_e)*sbc))/std::pow(sbc, 2);
            //
            // auto s2a = sbc - s1a;
            // auto s2b = sbc - s1b;
            //
            // std::println(stderr, "s1a: {}, s2a: {}", s1a, s2a);
            // std::println(stderr, "s1b: {}, s2b: {}", s1b, s2b);
            //
            // //auto _c1a = (ke - k0)/(s1a - s2a);
            //
            // s1 = s1a;
            // s2 = s2a;
            // c1 = c1a;
        }

        c2 = -c1;
        return { s1, s2, c1, c2 };
    }

    template<typename T>
    struct shape_point_t {
        T x;
        T y;
        T k;
        T phi;
    };

    template<typename T, typename U>
    concept Shape = requires(T a) {
        { a.get(U(0)) } -> std::same_as<shape_point_t<U>>;
        { a.len() } -> std::same_as<U>;
    };

    template<typename T, Shape<T> Q1, Shape<T> Q2>
    std::array<T, 4> fit_biclothoid(Q1 q1, Q2 q2, T t1, T eps) {
        auto [x0, y0, k0, phi_0] = q1.get(t1);
        auto sbc = q2.len();
        auto t2 = (T(1) - t1) * q1.len() / q2.len();

        T s1, s2, c1, c2;

        auto f1f2 = [&q2, &s1, &s2, &c1, &c2, &x0, &y0, &k0, &phi_0](T sbc, T t2) -> std::array<T, 2> {
            auto [xe, ye, ke, phi_e] = q2.get(t2);
            auto [_s1, _s2, _c1, _c2] = solve_biclothoid(phi_0, k0, phi_e, ke, sbc);
            s1 = _s1;
            s2 = _s2;
            c1 = _c1;
            c2 = _c2;
            auto phi_delta = phi_e - phi_0;
            auto [xc, yc] = xy(x0, y0, phi_0, k0, c1, s1);
            auto kc = k(k0, c1, s1);
            auto phi_c = phi(phi_0, k0, c1, s1);
            auto [_xe, _ye] = xy(xc, yc, phi_c, kc, c2, s2);
            return { _xe - xe, _ye - ye };
        };

        for(;;) {
            auto [f1, f2] = f1f2(sbc, t2);
            auto [f1_sbc_eps, f2_sbc_eps] = f1f2(sbc + eps, t2);
            auto [f1_t2_eps, f2_t2_eps] = f1f2(sbc, t2 + eps);

            auto df1_dsbc = (f1_sbc_eps - f1) / eps;
            auto df2_dsbc = (f2_sbc_eps - f2) / eps;
            auto df1_dt2 = (f1_t2_eps - f1) / eps;
            auto df2_dt2 = (f2_t2_eps - f2) / eps;

            auto J = Eigen::Matrix<T, 2, 2>({
                {df1_dsbc, df1_dt2},
                {df2_dsbc, df2_dt2}
            });

            auto x = J.inverse() * Eigen::Matrix<T, 2, 1>({{f1}, {f2}});

            sbc -= x(0, 0);
            t2 -= x(1, 0);

            //return { s1, s2, c1, c2 };
            if(std::abs(f1) < eps && std::abs(f2) < eps) {
                return { s1, s2, c1, c2 };
            }
        }
    }
}

#endif //CLOTHOID_H
