#include "Poly.hpp"

void Poly::trim()
{
    if(coeffs_.size() == 0)
    {
        return;
    }
    for (size_t i = coeffs_.size() - 1; true; --i)
    {
        if (coeffs_[i] == 0)
        {
            coeffs_.pop_back();
        }
        else
        {
            return;
        }
        if(i == 0)
        {
            return;
        }
    }
    return;
}

Poly::Poly(std::vector<gf2k::elem> a) : coeffs_(std::move(a))
{
    trim();
}

Poly Poly::operator+(const Poly &b) const
{
    if (coeffs_.size() > b.coeffs_.size())
    {
        std::vector<gf2k::elem> c = coeffs_;
        for (size_t i = 0; i < b.coeffs_.size(); ++i)
        {
            c[i] = gf2k::add(c[i], b[i]);
        }
        Poly new_poly = Poly(c);
        
        return new_poly;
    }
    else
    {
        std::vector<gf2k::elem> c = b.coeffs_;
        for (size_t i = 0; i < coeffs_.size(); ++i)
        {
            c[i] = gf2k::add(c[i], coeffs_[i]);
        }
        Poly new_poly = Poly(c);

        return new_poly;
    }
}

Poly Poly::operator*(const Poly &b) const
{
    std::vector<gf2k::elem> c(coeffs_.size() + b.coeffs_.size());

    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        for (size_t j = 0; j < b.coeffs_.size(); ++j)
        {
            c[i + j] = gf2k::add(c[i + j], gf2k::mul(coeffs_[i], b.coeffs_[j]));
        }
    }

    Poly new_poly = Poly(c);

    return new_poly;
}

Poly Poly::operator*(const gf2k::elem b) const
{
    std::vector<gf2k::elem> c = coeffs_;

    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        c[i] = gf2k::mul(c[i], b);
    }

    Poly new_poly = Poly(c);

    return new_poly;
}

gf2k::elem Poly::operator[](std::size_t i) const
{
    if (i >= coeffs_.size() || i < 0)
    {
        return 0;
    }
    return coeffs_[i];
}