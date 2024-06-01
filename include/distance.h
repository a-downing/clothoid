#ifndef DISTANCE_H
#define DISTANCE_H

#include <clothoid.h>
#include <shape.h>

namespace clothoid {
    template<typename Q1, typename Q2, typename T>
    T calc_distance(Q1 q1, Q2 q2, fit_result_t<T> fr) {
        auto [bc, t1, t2] = fr;
        const auto &q = lerp(q1, q2, 1.0);
    }
}

#endif //DISTANCE_H
