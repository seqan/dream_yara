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

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include "bit_chunk.hpp"

namespace seqan
{

template <typename string_t = Dna5String>
class SeqAnBloomFilter
{
private:
    using shape_t = Shape<Dna, SimpleShape>;

    uint8_t kmer_size{};
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf{};

public:
    SeqAnBloomFilter() = default;
    SeqAnBloomFilter(SeqAnBloomFilter const &) = default;
    SeqAnBloomFilter & operator=(SeqAnBloomFilter const &) = default;
    SeqAnBloomFilter(SeqAnBloomFilter &&) = default;
    SeqAnBloomFilter & operator=(SeqAnBloomFilter &&) = default;
    ~SeqAnBloomFilter() = default;

    SeqAnBloomFilter(uint32_t bin_count, uint8_t hash_function_count, uint8_t kmer_size, uint64_t full_size):
        kmer_size(kmer_size),
        ibf(seqan3::bin_count{bin_count},
            seqan3::bin_size{full_size / (((bin_count + 63) >> 6) << 6)},
            seqan3::hash_function_count{hash_function_count})
    {}

    SeqAnBloomFilter(const char * file_name)
    {
        load(file_name);
    }

    void addKmers(string_t const & text, uint32_t const bin_index)
    {
        if (length(text) < kmer_size)
            return;

        shape_t kmer_shape;
        resize(kmer_shape, kmer_size);

        auto it = begin(text);
        hashInit(kmer_shape, it);

        for (uint64_t i = 0; i < length(text) - length(kmer_shape) + 1; ++i, ++it)
        {
            ibf.emplace(hashNext(kmer_shape, it), seqan3::bin_index{bin_index});
        }
    }

    void addFastaFile(CharString const & fasta_file, uint32_t const bin_index)
    {
        CharString id;
        IupacString seq;

        SeqFileIn fin;
        if (!open(fin, toCString(fasta_file)))
        {
            CharString msg = "Unable to open contigs File: ";
            append (msg, fasta_file);
            throw toCString(msg);
        }
        while (!atEnd(fin))
        {
            readRecord(id, seq, fin);
            addKmers(seq, bin_index);
        }
        close(fin);
    }

    void clearBins(std::vector<uint32_t> const & bin_range, uint32_t const threads)
    {
        size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(ibf.bin_size() / threads),
                                                     8u,
                                                     64u);

        auto chunked_view = bin_range | extra::views::bit_chunk(chunk_size);

        auto worker = [&] (auto && chunked_range, auto &&)
        {
            ibf.clear(chunked_range | std::views::transform([] (uint32_t const i) { return seqan3::bin_index{i}; }));
        };

        seqan3::detail::execution_handler_parallel executioner{threads};
        executioner.bulk_execute(std::move(worker), std::move(chunked_view), [](){});
    }

    void whichBins(std::vector<bool> & selected, std::vector<uint64_t> & values, string_t const & text, uint16_t const threshold) const
    {
        size_t const possible = length(text) - kmer_size + 1;
        values.clear();

        shape_t kmer_shape;
        resize(kmer_shape, kmer_size);
        hashInit(kmer_shape, begin(text));

        auto it = begin(text);
        for (size_t i = 0; i < possible; ++i, ++it)
            values.push_back(hashNext(kmer_shape, it));

        auto agent = ibf.counting_agent();
        auto & counts = agent.bulk_count(values);
        assert(ibf.bin_count() == counts.size());
        assert(ibf.bin_count() == selected.size());
        assert(selected.size() == counts.size());

        for (size_t bin_index = 0; bin_index < ibf.bin_count(); ++bin_index)
        {
            if (counts[bin_index] >= threshold)
                selected[bin_index] = true;
        }
    }

    void whichBins(std::vector<bool> & selected, string_t const & text, uint16_t const threshold) const
    {
        size_t const possible = length(text) - kmer_size + 1;
        std::vector<uint64_t> values;
        values.reserve(possible);
        whichBins(selected, values, text, threshold);
    }

    std::vector<bool> whichBins(string_t const & text, uint16_t const threshold) const
    {
        std::vector<bool> selected(ibf.bin_count(), false);
        whichBins(selected, text, threshold);
        return selected;
    }

    size_t getNumberOfBins() const noexcept
    {
        return ibf.bin_count();
    }

    uint8_t getKmerSize() const noexcept
    {
        return kmer_size;
    }

    double size_mb() const
    {
        return sdsl::size_in_mega_bytes(ibf.raw_data());
    }

    void save(const char * file_name)
    {
        std::ofstream os{file_name, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        CEREAL_SERIALIZE_FUNCTION_NAME(oarchive);
    }

    void load(const char * file_name)
    {
        std::ifstream os{file_name, std::ios::binary};
        cereal::BinaryInputArchive iarchive{os};
        CEREAL_SERIALIZE_FUNCTION_NAME(iarchive);
    }

    template <typename archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(ibf);
        archive(kmer_size);
    }
};

} // namespace seqan
