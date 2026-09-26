#include "Poly.hpp"

void Poly::trim()
{
    if (coeffs_.size() == 0)
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
        if (i == 0)
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

gf2k::elem Poly::eval(gf2k::elem x) const
{
    if (coeffs_.size() == 0)
    {
        return 0;
    }
    gf2k::elem res = gf2k::zero();

    for (size_t i = coeffs_.size() - 1; true; --i)
    {
        if (i == 0)
        {
            res = gf2k::add(gf2k::mul_f(res, x), coeffs_[0]);
            return res;
        }

        res = gf2k::add(gf2k::mul_f(res, x), coeffs_[i]);
    }
}

std::vector<std::vector<gf2k::elem>> eval_matrix(std::vector<std::vector<Poly>> &M, size_t n, gf2k::elem x)
{
    std::vector<std::vector<gf2k::elem>> M_elem(n, std::vector<gf2k::elem>(n));

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            M_elem[i][j] = M[i][j].eval(x);
        }
    }

    return M_elem;
}

Poly poly_det(std::vector<std::vector<Poly>> &M, size_t n, size_t d)
{
    std::vector<gf2k::elem> points = gf2k::generate_numbers(d + 1);

    std::vector<gf2k::elem> dets(d + 1);
    std::vector<std::vector<gf2k::elem>> M_elem;

    for (size_t i = 0; i < d + 1; ++i)
    {
        M_elem = eval_matrix(M, n, points[i]);
        dets[i] = gf2k::det(M_elem, n);
    }

    for (size_t i = 1; i < d + 1; ++i)
    {
        for (size_t j = d; j >= i; --j)
        {
            gf2k::elem num = gf2k::add(dets[j], dets[j - 1]);
            gf2k::elem den = gf2k::add(points[j], points[j - i]);
            dets[j] = gf2k::mul_f(num, gf2k::inv(den));
        }
    }

    Poly det({dets[d]});

    if (d > 0)
    {
        for (size_t i = d - 1; true; --i)
        {
            Poly factor({points[i], gf2k::one()});
            Poly factor2({dets[i]});
            det = det * factor;
            det = det + factor2;

            if (i == 0)
            {
                break;
            }
        }
    }

    return det;
}