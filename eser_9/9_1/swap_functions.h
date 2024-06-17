#ifndef __SWAP_FUNCTIONS__
#define __SWAP_FUNCTIONS__

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>

namespace PBC_swap
{
    template <std::size_t SIZE>
    std::size_t PBC(const std::size_t idx, const std::size_t offset = 0)
    {
        const std::size_t NEW_SIZE = SIZE - offset;
        const std::size_t new_idx = idx % NEW_SIZE + offset;
        return new_idx;
    }

    template <std::size_t SIZE, class _ForwardIterator>
    void swap_ranges(_ForwardIterator array_begin, std::size_t start_first_range, std::size_t start_second_range, std::size_t range_lenght, std::size_t offset)
    {
        /// @todo Control that there aren't overlaps between ranges

        for (size_t i = 0; i < range_lenght; i++)
        {
            std::iter_swap(array_begin + PBC_swap::PBC<SIZE>(start_first_range, offset), array_begin + PBC_swap::PBC<SIZE>(start_second_range, offset));
            ++start_first_range;
            ++start_second_range;
        }
    }

    template <std::size_t SIZE, class _ForwardIterator>
    void rotate(_ForwardIterator array_begin, std::size_t first, std::size_t middle, std::size_t last, std::size_t offset)
    {
        // Copied from libcpp forward_rotate, modified in its swap parts
        std::size_t i = middle;
        while (true)
        {
            std::iter_swap(array_begin + PBC_swap::PBC<SIZE>(first, offset), array_begin + PBC_swap::PBC<SIZE>(i, offset));
            ++first;
            if (++i == last)
                break;
            assert((i <= last) && "Middle point greater than last elemet");
            if (first == middle)
                middle = i;
        }
        if (first != middle)
        {
            i = middle;
            while (true)
            {
                std::iter_swap(array_begin + PBC_swap::PBC<SIZE>(first, offset), array_begin + PBC_swap::PBC<SIZE>(i, offset));
                ++first;

                if (++i == last)
                {
                    if (first == middle)
                    {
                        break;
                    }
                    i = middle;
                    assert((i <= last) && "Middle point greater than last elemet, during recomposition");
                }
                else if (first == middle)
                    middle = i;
            }
        }
    }

    template <std::size_t SIZE, class _ForwardIterator>
    void reverse(_ForwardIterator array_begin, std::size_t first, std::size_t last, std::size_t offset)
    {
        if (first != last)
            for (; first < --last; ++first)
                std::iter_swap(array_begin + PBC_swap::PBC<SIZE>(first, offset), array_begin + PBC_swap::PBC<SIZE>(last, offset));
    }

} // namespace PBC_swap

#endif // __SWAP_FUNCTIONS__
