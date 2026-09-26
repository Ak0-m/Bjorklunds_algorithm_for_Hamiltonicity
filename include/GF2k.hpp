#pragma once

#include <cstdint>
#include <random>

namespace gf2k
{
using elem = uint64_t;

inline elem add(elem a, elem b)
{
    return a ^ b;
}

elem mul(elem a, elem b);
elem mul_f(elem a, elem b);
elem inv(elem a);
elem pow(elem a, uint64_t e);

elem det(std::vector<std::vector<elem>> M, size_t n);

inline elem random(std::mt19937_64 &rng)
{
    return rng();
}

inline elem random_nonzero(std::mt19937_64 &rng)
{
    elem x;
    do
    {
        x = rng();
    } while (x == 0);
    return x;
}

inline elem one()
{
    return 1;
}
inline elem zero()
{
    return 0;
}
} // namespace gf2k
