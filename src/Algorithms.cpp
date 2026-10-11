#include "Algorithms.hpp"
#include "Poly.hpp"
#include <bits/stdc++.h>
#include <vector>

std::vector<gf2k::elem> fhat(const std::vector<std::vector<gf2k::elem>> &g, size_t v1_elem, std::vector<size_t> &V2)
{
    const size_t N2 = 1 << V2.size();
    std::vector<gf2k::elem> dp(N2 * g.size(), gf2k::zero());

    auto at = [g](size_t mask, size_t v) -> size_t { return (size_t)mask * g.size() + v; };

    for (size_t u = 0; u < g.size(); ++u)
    {
        dp[at(0, u)] = g[u][v1_elem];
    }

    for (size_t mask = 1; mask < N2; ++mask)
    {
        for (size_t u = 0; u < g.size(); ++u)
        {
            gf2k::elem sum = gf2k::zero();
            size_t m = mask;

            while (m)
            {
                int bit = __builtin_ctzll(m);
                m &= m - 1;

                gf2k::elem wz = g[u][V2[bit]];

                if (wz == gf2k::zero())
                {
                    continue;
                }

                size_t rest = mask ^ (1 << bit);
                sum = gf2k::add(sum, gf2k::mul_f(wz, dp[at(rest, V2[bit])]));
            }

            dp[at(mask, u)] = sum;
        }
    }
    return dp;
}

std::vector<std::vector<gf2k::elem>> assign_label_weights(std::vector<std::vector<size_t>> &g, std::vector<size_t> &V2,
                                                          std::mt19937_64 &rng)
{
    std::vector<std::vector<gf2k::elem>> w(g.size(), std::vector<gf2k::elem>(g.size(), gf2k::zero()));

    for (auto i : V2)
    {
        for (auto j : g[i])
        {
            if (w[i][j] != 0)
            {
                continue;
            }
            w[i][j] = gf2k::random_nonzero(rng);

            if (j == 0)
            {
                w[j][i] = gf2k::random_nonzero(rng);
                continue;
            }
            w[j][i] = w[i][j];
        }
    }

    return w;
}
std::vector<std::vector<std::vector<gf2k::elem>>> assign_dummy_label_weights(std::vector<std::vector<size_t>> &g,
                                                                             std::vector<size_t> &V1, size_t mmax,
                                                                             std::mt19937_64 &rng)
{
    std::vector<std::vector<std::vector<gf2k::elem>>> w(
        mmax, std::vector<std::vector<gf2k::elem>>(g.size(), std::vector<gf2k::elem>(g.size(), gf2k::zero())));

    std::vector<bool> inV1(g.size(), 0);
    for (auto i : V1)
    {
        inV1[i] = 1;
    }

    for (size_t m = 0; m < mmax; ++m)
    {
        for (auto i : V1)
        {
            for (auto j : g[i])
            {
                if (inV1[j] == 0)
                {
                    continue;
                }
                if (w[m][i][j] != 0)
                {
                    continue;
                }
                w[m][i][j] = gf2k::random_nonzero(rng);

                if (j == 0 || i == 0)
                {
                    w[m][j][i] = gf2k::random_nonzero(rng);
                    continue;
                }
                w[m][j][i] = w[m][i][j];
            }
        }
    }

    return w;
}

std::vector<std::vector<std::vector<gf2k::elem>>> compute_Mf(gf2k::elem r, size_t m,
                                                             std::vector<std::vector<gf2k::elem>> &fh,
                                                             std::vector<std::vector<std::vector<gf2k::elem>>> &dw,
                                                             std::vector<size_t> &V1, std::vector<size_t> &V2)
{
    size_t L = V2.size() + m;
    size_t L_num = 1 << L;
    size_t V2_mask = (1 << V2.size()) - 1;

    std::vector<std::vector<std::vector<gf2k::elem>>> Mf(
        L_num, std::vector<std::vector<gf2k::elem>>(V1.size(), std::vector<gf2k::elem>(V1.size(), gf2k::zero())));

    for (size_t Z = 1; Z < L_num; ++Z)
    {
        size_t V2_part = Z & V2_mask;
        size_t dummy_part = Z >> V2.size();

        if (V2_part != 0 && dummy_part == 0)
        {
            size_t k = __builtin_popcountll(V2_part);
            gf2k::elem rk = gf2k::pow(r, k);

            for (size_t v = 0; v < V1.size(); ++v)
            {
                for (size_t u = 0; u < V1.size(); ++u)
                {
                    if (u == v)
                    {
                        continue;
                    }
                    Mf[Z][u][v] = gf2k::mul_f(fh[v][V2_part * (V1.size() + V2.size()) + V1[u]], rk);
                }
            }
        }
        else if (V2_part == 0 && __builtin_popcountll(dummy_part) == 1)
        {
            size_t k = __builtin_ctzll(dummy_part);

            for (size_t v = 0; v < V1.size(); ++v)
            {
                for (size_t u = 0; u < V1.size(); ++u)
                {
                    if (u == v)
                    {
                        continue;
                    }
                    Mf[Z][u][v] = gf2k::mul_f(dw[k][V1[u]][V1[v]], r);
                }
            }
        }
    }
    return Mf;
}

void yates_transform(std::vector<std::vector<std::vector<gf2k::elem>>> &g, size_t L)
{
    size_t L_num = 1 << L;
    size_t n1 = g[0].size();

    for (size_t i = 0; i < L; ++i)
    {
        size_t bit = 1 << i;
        #pragma omp parallel for schedule(static)
        for (size_t Y = 0; Y < L_num; ++Y)
        {
            if (Y & bit)
            {
                size_t Y_prev = Y ^ bit;
                for (size_t u = 0; u < n1; ++u)
                {
                    for (size_t v = 0; v < n1; ++v)
                    {
                        g[Y][u][v] = gf2k::add(g[Y][u][v], g[Y_prev][u][v]);
                    }
                }
            }
        }
    }
}

gf2k::elem det_sum(std::vector<std::vector<std::vector<gf2k::elem>>> &g, size_t n1)
{
    gf2k::elem sum = gf2k::zero();

    #pragma omp parallel for reduction(^:sum) schedule(dynamic, 64)
    for (int i = 0; i < g.size(); ++i)
    {
        auto T = g[i];
        sum = sum ^ gf2k::det(T, n1);
    }

    return sum;
}

gf2k::elem run_lccs(gf2k::elem r, size_t m, std::vector<std::vector<gf2k::elem>> &fh,
                    std::vector<std::vector<std::vector<gf2k::elem>>> &dw, std::vector<size_t> &V1,
                    std::vector<size_t> &V2)
{
    auto g = compute_Mf(r, m, fh, dw, V1, V2);
    yates_transform(g, V2.size() + m);

    return det_sum(g, V1.size());
}

gf2k::elem extract_coef(size_t m, std::vector<std::vector<gf2k::elem>> &fh,
                        std::vector<std::vector<std::vector<gf2k::elem>>> &dw, std::vector<size_t> &V1,
                        std::vector<size_t> &V2)
{
    size_t L = V2.size() + m;
    size_t deg = L * V1.size();

    std::vector<gf2k::elem> points(deg + 1), targets(deg + 1);

    for (size_t i = 0; i < deg + 1; ++i)
    {
        gf2k::elem i_p = (gf2k::elem)i;
        gf2k::elem i_t = run_lccs(i_p, m, fh, dw, V1, V2);
        points[i] = i_p;
        targets[i] = i_t;
    }

    Poly p = interpolate(points, targets);

    return p[L];
}

void random_partition(size_t n, std::vector<size_t> &V1, std::vector<size_t> &V2, std::mt19937_64 &rng)
{
    V1.clear();
    V2.clear();

    std::vector<size_t> elems(n - 1);

    for (size_t i = 0; i < elems.size(); ++i)
    {
        elems[i] = i + 1;
    }

    std::shuffle(elems.begin(), elems.end(), rng);
    V1.push_back(0);
    for (size_t i = 0; i < elems.size(); ++i)
    {
        if (i % 2 == 0)
        {
            V2.push_back(elems[i]);
        }
        else
        {
            V1.push_back(elems[i]);
        }
    }
}

bool run_once(std::vector<std::vector<size_t>> &g, size_t mmax, std::mt19937_64 &rng)
{
    std::vector<size_t> V1, V2;
    random_partition(g.size(), V1, V2, rng);

    auto w = assign_label_weights(g, V2, rng);
    auto dw = assign_dummy_label_weights(g, V1, mmax, rng);

    std::vector<std::vector<gf2k::elem>> fh(V1.size(), std::vector<gf2k::elem>(V2.size() * g.size()));

    for (size_t i = 0; i < V1.size(); ++i)
    {
        fh[i] = fhat(w, V1[i], V2);
    }

    for (std::size_t m = 0; m <= mmax; ++m)
    {
        gf2k::elem coeff = extract_coef(m, fh, dw, V1, V2);
        if (coeff != gf2k::zero())
            return true;
    }
    return false;
}

bool decide(std::vector<std::vector<size_t>> &g, size_t trials, std::mt19937_64 &rng)
{
    const std::size_t m_max = (std::size_t)(std::ceil(0.205 * static_cast<double>(g.size())));

    for (std::size_t i = 0; i < trials; ++i)
    {
        if (run_once(g, m_max, rng))
        {
            return true;
        }
    }

    return false;
}
