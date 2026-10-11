#pragma once

#include "GF2k.hpp"

std::vector<std::vector<gf2k::elem>> assign_label_weights(std::vector<std::vector<size_t>> &g, std::vector<size_t> &V2,
                                                          std::mt19937_64 &rng);

std::vector<std::vector<std::vector<gf2k::elem>>> assign_dummy_label_weights(std::vector<std::vector<size_t>> &g,
                                                                             std::vector<size_t> &V1, size_t mmax,
                                                                             std::mt19937_64 &rng);

std::vector<gf2k::elem> fhat(const std::vector<std::vector<gf2k::elem>> &g, size_t v1_, std::vector<size_t> &V2);

std::vector<std::vector<std::vector<gf2k::elem>>> compute_Mf(gf2k::elem r, size_t m,
                                                             std::vector<std::vector<gf2k::elem>> &fh,
                                                             std::vector<std::vector<std::vector<gf2k::elem>>> &dw,
                                                             std::vector<size_t> &V1, std::vector<size_t> &V2);

void yates_transform(std::vector<std::vector<std::vector<gf2k::elem>>> &g, size_t L);

gf2k::elem det_sum(std::vector<std::vector<std::vector<gf2k::elem>>> &g, size_t n);

gf2k::elem run_lccs(gf2k::elem r, gf2k::elem m, std::vector<std::vector<gf2k::elem>> &fh,
                    std::vector<std::vector<std::vector<gf2k::elem>>> &dw, std::vector<size_t> &V1,
                    std::vector<size_t> &V2);

gf2k::elem extract_coef(gf2k::elem m, std::vector<std::vector<gf2k::elem>> &fh,
                    std::vector<std::vector<std::vector<gf2k::elem>>> &dw, std::vector<size_t> &V1,
                    std::vector<size_t> &V2);

void random_partition(size_t n, std::vector<size_t> &V1, std::vector<size_t> &V2, std::mt19937_64 &rng);

bool run_once(std::vector<std::vector<size_t>> &g, size_t mmax, std::mt19937_64 &rng);

bool decide(std::vector<std::vector<size_t>> &g, size_t trials, std::mt19937_64 &rng);


