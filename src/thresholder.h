// ==========================================================================
//                                 d_bloom_filter.h
// ==========================================================================
// Copyright (c) 2017-2022, Enrico Seiler, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Enrico Seiler or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SEILER OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#pragma once

#include <vector>

#include <raptor/minimiser_model.hpp>

class Thresholder
{
public:
    Thresholder() = default;
    Thresholder(Thresholder const &) = default;
    Thresholder & operator=(Thresholder const &) = default;
    Thresholder(Thresholder &&) = default;
    Thresholder & operator=(Thresholder &&) = default;
    ~Thresholder() = default;

    Thresholder(size_t const pattern_size_,
                size_t const window_size_,
                size_t const kmer_size_,
                size_t const errors_,
                std::filesystem::path const & ibf_file)
        : pattern_size(pattern_size_),
          window_size(window_size_),
          kmer_size(kmer_size_),
          errors(errors_)
    {
        kmer_lemma = pattern_size + 1u > (errors + 1u) * kmer_size ? pattern_size + 1u - (errors + 1u) * kmer_size : 0;
        if (window_size != kmer_size)
        {
            if (!do_cerealisation_in(minimiser_thresholds,
                                     ibf_file,
                                     pattern_size,
                                     window_size,
                                     kmer_size,
                                     errors,
                                     tau))
            {
                minimiser_thresholds = precompute_threshold(pattern_size,
                                                            window_size,
                                                            kmer_size,
                                                            errors,
                                                            tau,
                                                            fpr,
                                                            correction_threshold,
                                                            do_correction);

                do_cerealisation_out(minimiser_thresholds,
                                     ibf_file,
                                     pattern_size,
                                     window_size,
                                     kmer_size,
                                     errors,
                                     tau);
            }
        }
    }

    size_t operator()(size_t const num_minimiser = 0) const noexcept
    {
        return window_size == kmer_size ? kmer_lemma : minimiser_thresholds[num_minimiser] + 1;
    }

    size_t check() const noexcept
    {
        if (window_size == kmer_size)
            return kmer_lemma;
        else
            return *std::max_element(minimiser_thresholds.begin(), minimiser_thresholds.end());
    }

private:
    size_t pattern_size{};
    size_t window_size{};
    size_t kmer_size{};
    size_t errors{};
    static constexpr double tau{0.9999f};
    static constexpr double fpr{0.15f};
    static constexpr double correction_threshold{0.05f};
    static constexpr bool do_correction{false};

    size_t kmer_lemma{};
    std::vector<size_t> minimiser_thresholds{};
};
