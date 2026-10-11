#include "Algorithms.hpp"
#include "GF2k.hpp"
#include "Poly.hpp"

#include <iostream>
#include <vector>
#include <random>

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::size_t T;
    if (!(std::cin >> T))
        return 0;

    for (std::size_t t = 0; t < T; ++t)
    {
        std::size_t n, m;
        std::cin >> n >> m;

        std::vector<std::vector<std::size_t>> g(n);
        for (std::size_t i = 0; i < m; ++i)
        {
            std::size_t u, v;
            std::cin >> u >> v;
            g[u].push_back(v);
            g[v].push_back(u);
        }

        std::mt19937_64 rng(67);
        bool hc = decide(g, 30, rng);

        std::cout << "{";
        std::cout << "\"id\":" << t << ",";
        std::cout << "\"n\":" << n << ",";
        std::cout << "\"m\":" << m << ",";
        std::cout << "\"has_cycle\":" << (hc ? "true" : "false");
        std::cout << "}\n";
    }

    return 0;
}
