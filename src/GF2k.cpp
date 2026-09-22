#include "GF2k.hpp"
#include <cstdint>
#include <immintrin.h>
#include <random>

namespace gf2k
{
elem mul(elem a, elem b)
{
    elem res = 0;
    while (b)
    {
        if (b & 1)
            res ^= a;
        uint64_t top = a >> 63;
        a <<= 1;
        if (top)
            a ^= 0x1B;
        b >>= 1;
    }
    return res;
}

elem mul_f(elem a, elem b)
{
    __m128i va = _mm_cvtsi64_si128(a);
    __m128i vb = _mm_cvtsi64_si128(b);
    __m128i prod = _mm_clmulepi64_si128(va, vb, 0x00);

    uint64_t lo = _mm_extract_epi64(prod, 0);
    uint64_t hi = _mm_extract_epi64(prod, 1);

    const __m128i poly = _mm_cvtsi64_si128(0x1B);

    __m128i red1 = _mm_clmulepi64_si128(_mm_cvtsi64_si128(hi), poly, 0x00);
    lo ^= _mm_extract_epi64(red1, 0);
    hi = _mm_extract_epi64(red1, 1);

    __m128i red2 = _mm_clmulepi64_si128(_mm_cvtsi64_si128(hi), poly, 0x00);
    lo ^= _mm_extract_epi64(red2, 0);

    return lo;
}

elem inv(elem a)
{
    if (a == 0)
    {
        return 0;
    }

    return pow(a, 0xFFFFFFFFFFFFFFFEULL);
}
elem pow(elem a, uint64_t e)
{
    elem r = 1;
    while (e)
    {
        if(e & 1)
        {
            r = mul_f(r, a);
        }
        a = mul(a, a);
        e >>= 1;
    }

    return r;
}
} // namespace gf2k
