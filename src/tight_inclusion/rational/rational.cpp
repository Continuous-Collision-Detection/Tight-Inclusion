#include <tight_inclusion/rational/rational.hpp>

#include <gmp.h>

namespace ticcd {

    // clang-format off
    template <> bool Rational::operator< <Rational>(const Rational &r1)
    // clang-format on
    {
        return mpq_cmp(value, r1.value) < 0;
    }

    template <> bool Rational::operator><Rational>(const Rational &r1)
    {
        return mpq_cmp(value, r1.value) > 0;
    }

    template <> bool Rational::operator<=<Rational>(const Rational &r1)
    {
        return mpq_cmp(value, r1.value) <= 0;
    }

    template <> bool Rational::operator>=<Rational>(const Rational &r1)
    {
        return mpq_cmp(value, r1.value) >= 0;
    }

    template <> bool Rational::operator==<Rational>(const Rational &r1)
    {
        return mpq_equal(value, r1.value);
    }

    template <> bool Rational::operator!=<Rational>(const Rational &r1)
    {
        return !mpq_equal(value, r1.value);
    }

} // namespace ticcd
