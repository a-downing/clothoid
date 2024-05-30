#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

#include <array>

namespace clothoid::gauss_legendre {
    constexpr std::array P5_X = {
        -0.906179845938664,
        0.906179845938664,
        -0.538469310105683,
        0.538469310105683,
        0.0
    };

    constexpr std::array P5_W = {
        0.23692688505618867,
        0.23692688505618867,
        0.47862867049936625,
        0.47862867049936625,
        0.5688888888888889
    };

    template <typename T, size_t N>
    std::array<T, N> map_x(std::array<T, N> xs, T a, T b) {
        for(auto &x: xs) {
            x = 0.5 * (b - a) * x + 0.5 * (a + b);
        }

        return xs;
    }

    template <typename T, size_t N>
    std::array<T, N> map_w(std::array<T, N> ws, T a, T b) {
        for(auto &w: ws) {
            w = 0.5 * (b - a) * w;
        }

        return ws;
    }

    template <typename T, size_t N>
    std::array<std::array<T, N>, 2> map_xw(std::array<T, N> xs, std::array<T, N> ws, T a, T b) {
        return { map_x(xs, a, b), map_w(ws, a, b) };
    }
}

#endif //GAUSS_LEGENDRE_H
