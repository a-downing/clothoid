#ifndef COMMON_H
#define COMMON_H

namespace clothoid
{
    template<typename T>
    T lerp(T a, T b, T t) {
        return a + t * (b - a);
    }
}

#endif //COMMON_H
