#ifndef VEC2_H
#define VEC2_H

#include <cmath>

namespace clothoid
{
    template<typename T>
    struct vec2 {
        T x;
        T y;
        vec2(T x, T y) : x(x), y(y) { }

        vec2& operator=(const vec2 &b) {
            x = b.x;
            y = b.y;
            return *this;
        }

        T len() { return std::sqrt(x*x + y*y); }
    };

    template<typename T>
    static vec2<T> operator-(const vec2<T> &a, const vec2<T> &b) {
        return { a.x - b.x, a.y - b.y };
    }

    template<typename T>
    static vec2<T> operator+(const vec2<T> &a, const vec2<T> &b) {
        return { a.x + b.x, a.y + b.y };
    }

    template<typename T>
    static vec2<T> operator*(const vec2<T> &a, T b) {
        return { a.x*b, a.y*b };
    }

    template<typename T>
    static vec2<T> operator*(T a, const vec2<T> &b) {
        return { a*b.x, a*b.y };
    }
}

#endif //VEC2_H
