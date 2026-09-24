#pragma once

#include "GF2k.hpp"
#include <cstddef>
#include <vector>

class Poly
{
    std::vector<gf2k::elem> coeffs_;
    void trim();

  public:
  
    Poly(std::vector<gf2k::elem> a);

    Poly operator+(const Poly &b) const;
    Poly operator*(const Poly &b) const;
    Poly operator*(const gf2k::elem b) const;
    gf2k::elem operator[](std::size_t i) const;

    inline bool is_zero() const
    {
        return coeffs_.empty();
    }


};