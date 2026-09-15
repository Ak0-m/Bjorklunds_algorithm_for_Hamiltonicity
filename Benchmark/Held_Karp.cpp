#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>

bool HK(const std::vector<std::vector<int>> &g, int n, bool cycle, std::vector<int> &out, long long &out_us)
{
    out.clear();
    out_us = 0;

    if (n == 0)
    {
        return false;
    }
    if (n < 3 && cycle)
    {
        return false;
    }

    const int N = 1 << n;

    std::vector<bool> dp((size_t)N * n, 0);
    std::vector<signed char> parent((size_t)N * n, -1);

    auto at = [n](int mask, int v) -> size_t { return (size_t)mask * n + v; };

    auto t0 = std::chrono::steady_clock::now();

    if (cycle)
    {
        dp[1] = 1;
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            dp[at(1 << i, i)] = 1;
        }
    }

    for (int mask = 0; mask < N; ++mask)
    {
        for (int v = 0; v < n; ++v)
        {
            if (!dp[at(mask, v)])
            {
                continue;
            }
            for (int u : g[v])
            {
                if (mask & (1 << u))
                {
                    continue;
                }

                int new_mask = mask | (1 << u);

                if (!dp[at(new_mask, u)])
                {
                    dp[at(new_mask, u)] = 1;
                    parent[at(new_mask, u)] = (signed char)v;
                }
            }
        }
    }

    const int full = N - 1;
    int end = -1;

    if (cycle)
    {
        for (int v : g[0])
        {
            if (dp[at(full, v)])
            {
                end = v;
                break;
            }
        }
    }
    else
    {
        for (int v = 0; v < n; ++v)
        {
            if (dp[at(full, v)])
            {
                end = v;
                break;
            }
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    out_us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

    if (end == -1)
    {
        return false;
    }

    std::vector<int> rev(n);
    int i = 0;
    int mask = full, v = end;

    while (v != -1)
    {
        rev[i] = v;
        ++i;

        int p = parent[at(mask, v)];

        mask ^= (1 << v);
        v = p;
    }

    std::reverse(rev.begin(), rev.end());
    out = std::move(rev);
    return true;
}

static void print_vec(const std::vector<int> &v)
{
    std::cout << "[";
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (i)
        {
            std::cout << ",";
        }
        std::cout << v[i];
    }
    std::cout << "]";
}

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int T;
    if (!(std::cin >> T))
        return 0;

    for (int id=0; id < T; ++id)
    {
        int n, m;
        std::cin >> n >> m;

        std::vector<std::vector<int>> g(n);

        for (int i = 0; i < m; ++i)
        {
            int u, v;
            std::cin >> u >> v;
            g[u].push_back(v);
            g[v].push_back(u);
        }

        std::vector<int> cycle, path;
        long long cycle_us = 0, path_us = 0;

        bool has_cycle = HK(g, n, true, cycle, cycle_us);

        bool has_path = false;
        if (has_cycle)
        {
            has_path = true;
            path = cycle;
            path_us = cycle_us;
        }
        else
        {
            has_path = HK(g, n, false, path, path_us);
        }

        std::string s = "none";

        if (has_cycle)
        {
            s = "cycle";
        }
        else if (has_path)
        {
            s = "path";
        }

        std::cout << "{";
        std::cout << "\"id\":" << id << ",";
        std::cout << "\"n\":" << n << ",";
        std::cout << "\"m\":" << m;

        if (has_cycle)
        {
            std::cout << ",\"kind\":\"cycle\"";
            std::cout << ",\"has_cycle\":true";
            std::cout << ",\"cycle_us\":" << cycle_us;
            std::cout << ",\"cycle_order\":";
            std::vector<int> c = cycle;
            c.push_back(c[0]);
            print_vec(c);

            std::cout << ",\"has_path\":true";
            std::cout << ",\"path_us\":" << path_us;
            std::cout << ",\"path_order\":";
            print_vec(path);
        }
        else if (has_path)
        {
            std::cout << ",\"kind\":\"path\"";
            std::cout << ",\"has_cycle\":false";
            std::cout << ",\"has_path\":true";
            std::cout << ",\"path_us\":" << path_us;
            std::cout << ",\"path_order\":";
            print_vec(path);
        }
        else
        {
            std::cout << ",\"kind\":\"none\"";
            std::cout << ",\"has_cycle\":false";
            std::cout << ",\"has_path\":false";
        }

        std::cout << "}\n";
    }
    return 0;
}