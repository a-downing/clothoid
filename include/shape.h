#ifndef SHAPE_H
#define SHAPE_H

#include <vec2.h>
#include <common.h>

namespace clothoid
{
    template<typename T>
    class Shape {
    public:
        virtual ~Shape() = default;
        virtual vec2<T> point(T t) const = 0;
        virtual T k(T t) const = 0;
        virtual T phi(T t) const = 0;
        virtual T len() const = 0;
    };

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
}

#endif //SHAPE_H
