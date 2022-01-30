// An interval object.
#pragma once

#include <Eigen/Core>
#include <tight_inclusion/types.hpp>

namespace ticcd {
    // calculate a*(2^b)
    uint64_t power(const uint64_t a, const uint8_t b);

    // calculate 2^exponent
    uint64_t pow2(const uint8_t exponent);

    // return power t. n=result*2^t
    long reduction(const long n, long &result);

    //<k,n> pair present a number k/pow(2,n)
    struct NumCCD {
        uint64_t numerator = 0;
        uint8_t denom_power = 0;

        NumCCD() {}

        NumCCD(uint64_t p_numerator, uint8_t p_denom_power)
            : numerator(p_numerator), denom_power(p_denom_power)
        {
        }

        NumCCD(float x);
        NumCCD(double x);

        ~NumCCD() {}

        uint64_t denominator() const { return pow2(denom_power); }

        // convert NumCCD to double number
        Scalar value() const { return Scalar(numerator) / denominator(); }

        NumCCD operator+(const NumCCD &other) const;

        bool operator==(const NumCCD &other) const
        {
            return numerator == other.numerator
                   && denom_power == other.denom_power;
        }
        bool operator!=(const NumCCD &other) const { return !(*this == other); }
        bool operator<(const NumCCD &other) const;
        bool operator<=(const NumCCD &other) const
        {
            return (*this == other) || (*this < other);
        }
        bool operator>=(const NumCCD &other) const { return !(*this < other); }
        bool operator>(const NumCCD &other) const { return !(*this <= other); }

        bool operator<(const Scalar other) const { return value() < other; }

        static bool is_sum_leq_1(const NumCCD &num1, const NumCCD &num2);
    };

    // an interval represented by two double numbers
    struct Interval {
        NumCCD lower;
        NumCCD upper;

        Interval() {}

        Interval(const NumCCD &p_lower, const NumCCD &p_upper)
            : lower(p_lower), upper(p_upper)
        {
        }

        ~Interval() {}

        std::pair<Interval, Interval> bisect() const;

        bool overlaps(const Scalar r1, const Scalar r2) const;

        const NumCCD &operator[](int i) const
        {
            switch (i) {
            case 0:
                return lower;
            case 1:
                return upper;
            default:
                throw std::out_of_range("");
            }
        }
    };

    typedef Eigen::Array<Interval, 3, 1> Interval3;
    Array3 width(const Interval3 &x);

} // namespace ticcd
