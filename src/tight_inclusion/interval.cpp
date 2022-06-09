// An interval object.

#include <tight_inclusion/interval.hpp>

#include <iostream>

namespace ticcd {

    static constexpr uint8_t MAX_DENOM_POWER = 8 * sizeof(uint64_t) - 1;

    // calculate a*(2^b)
    uint64_t power(const uint64_t a, const uint8_t b)
    {
        // The fast bit shifting power trick only works if b is not too larger.
        assert(b < MAX_DENOM_POWER);
        // WARNING: Technically this can still fail with `b < MAX_DENOM_POWER` if `a > 1`.
        return a << b;
    }

    // return power t. n=result*2^t
    uint8_t reduction(const uint64_t n, uint64_t &result)
    {
        uint8_t t = 0;
        result = n;
        while (result != 0 && (result & 1) == 0) {
            result >>= 1;
            t++;
        }
        return t;
    }

    NumCCD::NumCCD(Scalar x)
    {
        // Use bisection to find an upper bound of x.
        assert(x >= 0 && x <= 1);
        NumCCD low(0, 0), high(1, 0), mid;

        // Hard code these cases for better accuracy.
        if (x == 0) {
            *this = low;
            return;
        } else if (x == 1) {
            *this = high;
            return;
        }

        do {
            mid = low + high;
            mid.denom_power++;

            if (mid.denom_power >= MAX_DENOM_POWER) {
                break;
            }

            if (x > mid) {
                low = mid;
            } else if (x < mid) {
                high = mid;
            } else {
                break;
            }
        } while (mid.denom_power < MAX_DENOM_POWER);
        *this = high;
        assert(x <= value());
    }

    NumCCD NumCCD::operator+(const NumCCD &other) const
    {
        const uint64_t &k1 = numerator, &k2 = other.numerator;
        const uint8_t &n1 = denom_power, &n2 = other.denom_power;

        NumCCD result;
        if (n1 == n2) {
            result.denom_power = n2 - reduction(k1 + k2, result.numerator);
        } else if (n2 > n1) {
            result.numerator = k1 * pow2(n2 - n1) + k2;
            assert(result.numerator % 2 == 1);
            result.denom_power = n2;
        } else { // n2 < n1
            result.numerator = k1 + k2 * pow2(n1 - n2);
            assert(result.numerator % 2 == 1);
            result.denom_power = n1;
        }
        return result;
    }

    bool NumCCD::operator<(const NumCCD &other) const
    {
        const uint64_t &k1 = numerator, &k2 = other.numerator;
        const uint8_t &n1 = denom_power, &n2 = other.denom_power;

        uint64_t tmp_k1 = k1, tmp_k2 = k2;
        if (n1 < n2) {
            tmp_k1 = pow2(n2 - n1) * k1;
        } else if (n1 > n2) {
            tmp_k2 = pow2(n1 - n2) * k2;
        }
        assert((value() < other.value()) == (tmp_k1 < tmp_k2));
        return tmp_k1 < tmp_k2;
    }

    bool NumCCD::is_sum_leq_1(const NumCCD &num1, const NumCCD &num2)
    {
        if (num1.denom_power == num2.denom_power) {
            // skip the reduction in num1 + num2
            return num1.numerator + num2.numerator <= num1.denominator();
        }
        NumCCD tmp = num1 + num2;
        return tmp.numerator <= tmp.denominator();
    }

    std::pair<Interval, Interval> Interval::bisect() const
    {
        // interval is [k1/pow2(n1), k2/pow2(n2)]
        NumCCD mid = upper + lower;
        mid.denom_power++; // ÷ 2
        assert(mid.value() > lower.value() && mid.value() < upper.value());
        return std::make_pair(Interval(lower, mid), Interval(mid, upper));
    }

    bool Interval::overlaps(const Scalar r1, const Scalar r2) const
    {
        return upper.value() >= r1 && lower.value() <= r2;
    }

    Array3 width(const Interval3 &x)
    {
        return Array3(
            x[0].upper.value() - x[0].lower.value(),
            x[1].upper.value() - x[1].lower.value(),
            x[2].upper.value() - x[2].lower.value());
    }

} // namespace ticcd
