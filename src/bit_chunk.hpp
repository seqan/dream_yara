// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides extra::bit_chunk_view.
 */

#pragma once

#include <seqan3/std/bit>

#include <seqan3/range/views/slice.hpp>

namespace extra
{

template <std::ranges::view urng_t>
class bit_chunk_view : public std::ranges::view_interface<bit_chunk_view<urng_t>>
{
private:
    static_assert(std::ranges::forward_range<urng_t>, "The bit_chunk_view only works on forward_ranges");
    static_assert(std::unsigned_integral<std::remove_cvref_t<std::ranges::range_reference_t<urng_t>>>,
                  "The reference type of the underlying range must model std::unsigned_integral.");

    urng_t urange;

    size_t shift_value{};

    template <bool const_range>
    class basic_iterator;

public:
    bit_chunk_view()                                       = default;
    bit_chunk_view(bit_chunk_view const & rhs)             = default;
    bit_chunk_view(bit_chunk_view && rhs)                  = default;
    bit_chunk_view & operator=(bit_chunk_view const & rhs) = default;
    bit_chunk_view & operator=(bit_chunk_view && rhs)      = default;
    ~bit_chunk_view()                                      = default;

    bit_chunk_view(urng_t urange_, size_t const chunk_size) :
        urange{std::move(urange_)},
        shift_value{static_cast<size_t>(std::countr_zero(chunk_size))}
    {}

    template <typename rng_t>
     requires (!std::same_as<std::remove_cvref_t<rng_t>, bit_chunk_view>) &&
              std::ranges::viewable_range<rng_t> &&
              std::constructible_from<urng_t, std::ranges::ref_view<std::remove_reference_t<rng_t>>>
    bit_chunk_view(rng_t && urange_, size_t const chunk_size) :
        urange{std::views::all(std::forward<rng_t>(urange_))},
        shift_value{static_cast<size_t>(std::countr_zero(chunk_size))}
    {}

    auto begin() noexcept
    {
        return basic_iterator<false>{std::ranges::begin(urange), std::ranges::end(urange), shift_value, *this};
    }

    auto begin() const noexcept
        requires seqan3::const_iterable_range<urng_t>
    {
        return basic_iterator<true>{std::ranges::cbegin(urange), std::ranges::cend(urange), shift_value, *this};
    }

    auto end() noexcept
    {
        return std::ranges::end(urange);
    }

    auto end() const noexcept
        requires seqan3::const_iterable_range<urng_t>
    {
        return std::ranges::cend(urange);
    }
};

template <std::ranges::view urng_t>
template <bool const_range>
class bit_chunk_view<urng_t>::basic_iterator
{
private:
    using it_t = seqan3::detail::maybe_const_iterator_t<const_range, urng_t>;

    using sentinel_t = seqan3::detail::maybe_const_sentinel_t<const_range, urng_t>;

    template <bool other_const_range>
    friend class basic_iterator;

public:
    using difference_type = typename std::iter_difference_t<it_t>;
    using value_type = decltype(std::declval<urng_t>() | seqan3::views::slice(0,0));
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::input_iterator_tag;

    constexpr basic_iterator()                                   = default;
    constexpr basic_iterator(basic_iterator const &)             = default;
    constexpr basic_iterator(basic_iterator &&)                  = default;
    constexpr basic_iterator & operator=(basic_iterator const &) = default;
    constexpr basic_iterator & operator=(basic_iterator &&)      = default;
    ~basic_iterator()                                            = default;

    constexpr basic_iterator(basic_iterator<!const_range> const & it) noexcept
        requires const_range
        : shift_value{std::move(it.shift_value)},
          current_chunk{std::move(it.current_chunk)},
          begin_pos{std::move(it.begin_pos)},
          end_pos{std::move(it.end_pos)},
          chunk_it{std::move(it.chunk_it)},
          rng_end{std::move(it.rng_end)},
          parent{std::move(it.parent)}
    {}

    basic_iterator(it_t it_start, sentinel_t it_end, size_t const shift_value_, bit_chunk_view<urng_t> const & parent_) :
        shift_value{shift_value_}, chunk_it{it_start}, rng_end{it_end}, parent{std::addressof(parent_)}
    {
        current_chunk = *it_start >> shift_value;
        get_chunky();
    }

    friend bool operator==(basic_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return lhs.chunk_it == rhs;
    }

    friend bool operator==(sentinel_t const & lhs, basic_iterator const & rhs) noexcept
    {
        return lhs == rhs.chunk_it;
    }

    friend bool operator!=(basic_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    friend bool operator!=(sentinel_t const & lhs, basic_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    basic_iterator & operator++() noexcept
    {
        get_chunky();
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        get_chunky();
        return tmp;
    }

    value_type operator*() const noexcept
    {
        assert(parent != nullptr);
        return parent->urange | seqan3::views::slice(begin_pos, end_pos);
    }

private:
    size_t shift_value{};
    size_t current_chunk{};
    size_t begin_pos{};
    size_t end_pos{};

    it_t chunk_it;
    sentinel_t rng_end;

    bit_chunk_view<urng_t> const * parent{nullptr};

    void get_chunky()
    {
        begin_pos = end_pos;
        for (; chunk_it != rng_end; ++chunk_it)
        {
            if (*chunk_it >> shift_value == current_chunk)
            {
                ++end_pos;
            }
            else
            {
                current_chunk = *chunk_it >> shift_value;
                break;
            }
        }
    }
};

template <std::ranges::viewable_range rng_t>
bit_chunk_view(rng_t &&, size_t const) -> bit_chunk_view<std::views::all_t<rng_t>>;

struct bit_chunk_fn
{

    constexpr auto operator()(size_t const chunk_size) const
    {
        if (!std::has_single_bit(chunk_size))
            throw std::invalid_argument{"The chunk_size must be a power of two."};

        return seqan3::detail::adaptor_from_functor{*this, chunk_size};
    }

    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, size_t const chunk_size) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::bit_chunk cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::bit_chunk must model std::ranges::forward_range.");
        static_assert(std::unsigned_integral<std::remove_cvref_t<std::ranges::range_reference_t<urng_t>>>,
            "The reference type of the underlying range must model std::unsigned_integral.");

        if (!std::has_single_bit(chunk_size))
            throw std::invalid_argument{"The chunk_size must be a power of two."};

        return bit_chunk_view{std::forward<urng_t>(urange), chunk_size};
    }
};

}

namespace extra::views
{

inline constexpr auto bit_chunk = extra::bit_chunk_fn{};

}
