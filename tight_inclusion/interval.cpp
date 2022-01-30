// An interval object.

#include <tight_inclusion/interval.hpp>

#include <iostream>

namespace ticcd {

    // calculate a*(2^b)
    uint64_t power(const uint64_t a, const uint8_t b)
    {
        // The fast bit shifting power trick only works if b is not too larger.
        assert(b < 8 * sizeof(uint64_t) - 1);
        // WARNING: Technically this can still fail with `b < 8 * sizeof(uint64_t) - 1` if `a > 1`.
        return a << b;
    }

    uint64_t pow2(const uint8_t exponent) { return power(1l, exponent); }

    // return power t. n=result*2^t
    uint64_t reduction(const uint64_t n, uint64_t &result)
    {
        int t = 0;
        result = n;
        while (result != 0 && (result & 1) == 0) {
            result >>= 1;
            t++;
        }
        return t;
    }

    NumCCD::NumCCD(float p_x)
    {
        union {
            double val;
            struct {
                uint64_t mantisa : 23;
                uint16_t exponent : 8;
                uint8_t sign : 1;
            } parts;
        } x;
        x.val = p_x;

        denom_power = abs(x.parts.exponent - 127) + 23;
        numerator = x.parts.mantisa | (1l << 23);
        for (int i = 0; i < 23; i++) {
            if ((numerator & 1) == 0) {
                numerator >>= 1;
                denom_power--;
            } else {
                break;
            }
        }

        assert(value() == p_x);
    }

    NumCCD::NumCCD(double p_x)
    {
        union {
            double val;
            struct {
                uint64_t mantisa : 52;
                uint16_t exponent : 11;
                uint8_t sign : 1;
            } parts;
        } x;
        x.val = p_x;

        denom_power = abs(x.parts.exponent - 1023) + 52;
        numerator = x.parts.mantisa | (1l << 52);
        for (int i = 0; i < 52; i++) {
            if ((numerator & 1) == 0) {
                numerator >>= 1;
                denom_power--;
            } else {
                break;
            }
        }

        assert(value() == p_x);
    }

    NumCCD NumCCD::operator+(const NumCCD &other) const
    {
        const uint64_t &k1 = numerator, &k2 = other.numerator;
        const uint8_t &n1 = denom_power, &n2 = other.denom_power;

        uint64_t k;
        uint8_t n;
        if (n1 == n2) {
            n = n2 - reduction(k1 + k2, k);
        } else if (n2 > n1) {
            k = k1 * pow2(n2 - n1) + k2;
            assert(k % 2 == 1);
            n = n2;
        } else { // n2 < n1
            k = k1 + k2 * pow2(n1 - n2);
            assert(k % 2 == 1);
            n = n1;
        }
        assert(value() + other.value() == NumCCD(k, n).value());
        return NumCCD(k, n);
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
