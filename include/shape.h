#ifndef SHAPE_H
#define SHAPE_H

#include <vec2.h>

namespace clothoid
{
    template<typename T>
    class Shape {
    public:
        virtual ~Shape() = default;
        virtual vec2<T> point(T t) const = 0;
        virtual T k(T t) const = 0;
        virtual T phi(T t) const = 0;
        virtual T phi_delta(T t) const = 0;
        virtual T len() const = 0;

        static T phi_delta(const Shape &q1, const Shape &q2, T t1, T t2) {
            auto delta = q2.phi(t2) - q1.phi(t1);

            if(delta > M_PI) {
                delta -= 2*M_PI;
            } else if(delta < -M_PI) {
                delta += 2*M_PI;
            }

            return delta;
        }
    };
}

#endif //SHAPE_H
