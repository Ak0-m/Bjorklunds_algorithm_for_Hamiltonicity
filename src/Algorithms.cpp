#include "Algorithms.hpp"
#include <bits/stdc++.h>
#include <cstdint>
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
