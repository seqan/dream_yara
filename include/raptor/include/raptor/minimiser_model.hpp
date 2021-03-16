#include <raptor/shared.hpp>

#pragma once

std::vector<size_t> precompute_threshold(size_t const pattern_size,
                                         size_t const window_size,
                                         uint8_t const kmer_size,
                                         size_t const errors,
                                         double const tau,
                                         double const fpr,
                                         double const correction_threshold,
                                         bool const do_correction);
void do_cerealisation_out(std::vector<size_t> const & vec,
                          std::filesystem::path const & ibf_file,
                          size_t const pattern_size,
                          size_t const window_size,
                          size_t const kmer_size,
                          size_t const errors,
                          double const tau);
bool do_cerealisation_in(std::vector<size_t> & vec,
                         std::filesystem::path const & ibf_file,
                         size_t const pattern_size,
                         size_t const window_size,
                         size_t const kmer_size,
                         size_t const errors,
                         double const tau);
