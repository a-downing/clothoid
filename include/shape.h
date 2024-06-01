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
        virtual T len() const = 0;
    };
}

#endif //SHAPE_H
