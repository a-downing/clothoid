#ifndef COMMON_H
#define COMMON_H

namespace clothoid
{
    constexpr double RELTOL = 1e-6;
    constexpr double ABSTOL = 1e-6;

    template<typename T>
    T lerp(T a, T b, T t) {
        return a + t * (b - a);
    }

    inline bool is_close(double x, double target, double reltol = RELTOL, double abstol = ABSTOL) {
        auto err = std::abs(x - target);

        if(target == 0) {
            return err <= abstol;
        }

        return err <= abstol || err / std::abs(target) <= reltol;
    }

    inline bool near_zero(double x) {
        return is_close(x, 0);
    }
}

#endif //COMMON_H
