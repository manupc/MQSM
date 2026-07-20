/*
 * Copyright (C) 2026 Manuel Pegalajar Cuellar
 * Department of Computer Science and Artificial Intelligence, University of Granada (Spain)
 * E-mail: manupc@ugr.es
 *
 * This file is part of MQSMSim.
 *
 * MQSMSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MQSMSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MQSMSim.  If not, see <https://www.gnu.org/licenses/>.
 *
 * AI-ASSISTED DEVELOPMENT STATEMENT:
 * This file was modified or optimized with the assistance of 
 * Artificial Intelligence. AI was utilized for code completion, debugging, 
 * documentation generation, and fine-tuning optimizations.
 */
#ifndef __BITMASK__H__
#define __BITMASK__H__

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <ostream>
#include <span>
#include <string>
#include <vector>
#include <type_traits>

#include <private/BasicDataTypes.h>


/**
 * @class BitMask
 * @brief An efficient fixed-size bitset representation using unsigned integer arrays.
 *
 * This class provides a compact and high-performance representation for manipulating
 * bits, which is heavily used for tracking qubits and gate targets/controls.
 * 
 * @tparam N The number of bits the mask will track.
 */
template<unsigned int N>
requires(N > 0 && N <= MAX_QUBITS_ALLOWED)
class BitMask {

private: // ---- REPRESENTATION -----

    static constexpr std::size_t NumBlocks = (N + 63) / 64;

    using StorageType = std::conditional_t<N <= 8,  uint8_t,
                        std::conditional_t<N <= 16, uint16_t,
                        std::conditional_t<N <= 32, uint32_t,
                        std::conditional_t<N <= 64, uint64_t,
                        std::array<uint64_t, NumBlocks>>>>>;

    StorageType data;

    // Helper: Zeroes out bits beyond N-1 to prevent padding pollution
    constexpr void clean_padding() noexcept {
        if constexpr (N <= 64) {
            if constexpr (N < sizeof(StorageType) * 8) {
                data &= (static_cast<StorageType>(1ULL << N) - 1);
            }
        } else {
            if constexpr (N % 64 != 0) {
                data.back() &= (1ULL << (N % 64)) - 1;
            }
        }
    }



public: // ---- CONSTRUCTORS ----
    /**
     * @brief Default constructor. Initializes all bits to 0.
     */
    constexpr BitMask() noexcept : data{0} {}
    
    /**
     * @brief Constructs a bitmask from a span of set bit indices.
     * @param bits The span containing the indices of the bits to set.
     * @throws std::invalid_argument If there are duplicate bits.
     */
    BitMask(std::span<const BitIndex> bits) : data{0} {
        for (const auto& bit : bits) {
            if (isSet(bit)) {
                throw std::invalid_argument("BitMask::BitMask: Duplicate bit in initializer span");
            }
            set(bit);
        }
    }

    /**
     * @brief Constructs a bitmask from an initializer list of bit indices.
     */
    BitMask(std::initializer_list<BitIndex> bits) : BitMask(std::span<const BitIndex>(bits)) {}

    /**
     * @brief Constructs a bitmask from a binary string.
     * @param bit_string A string containing exactly N '0' and '1' characters.
     * @throws std::invalid_argument If the string size doesn't match N or contains invalid characters.
     */
    BitMask(std::string_view bit_string) : data{0} {
        if (bit_string.size() != N) {
            throw std::invalid_argument("BitMask::BitMask: Bit string length does not match template width");
        }
        for (std::size_t i = 0; i < N; ++i) {
            if (bit_string[N - 1 - i] == '1') {
                set(static_cast<BitIndex>(i));
            } else if (bit_string[N - 1 - i] != '0') {
                throw std::invalid_argument("BitMask::BitMask: Bit string must only contain '0' or '1'");
            }
        }
    }
    
    ~BitMask() = default;

    // Copy & Move
    BitMask(const BitMask<N> &) = default;
    BitMask(BitMask<N> &&) noexcept = default;
    BitMask<N> & operator=(const BitMask<N> &) = default;
    BitMask<N> & operator=(BitMask<N> &&) noexcept = default;
    
    BitMask<N> & operator=(std::initializer_list<BitIndex> bits) {
        clear();
        for (const auto& bit : bits) set(bit);
        return *this;
    }

    /**
     * @brief Creates a bitmask with a continuous range of set bits.
     * @param first The starting bit index.
     * @param last The ending bit index.
     * @return A new BitMask with bits in [first, last] set to 1.
     */
    static BitMask<N> rangeSet(BitIndex first, BitIndex last) {
        if (first > last) std::swap(first, last);
        if (last >= N) {
            throw std::invalid_argument("BitMask::rangeSet: Last index is out of bounds");
        }
        BitMask<N> result;
        for (BitIndex i = first; i <= last; ++i) {
            result.set(i);
        }
        return result;
    }

public: // ---- BITWISE OPERATORS ----
    
    BitMask<N> & operator&=(const BitMask<N> &other) noexcept {
        if constexpr (N <= 64) {
            data &= other.data;
        } else {
            for (size_t i = 0; i < data.size(); ++i) data[i] &= other.data[i];
        }
        return *this;
    }
    BitMask<N> operator&(const BitMask<N> &other) const noexcept {
        BitMask<N> result = *this;
        result &= other;
        return result;
    }

    BitMask<N> & operator|=(const BitMask<N> &other) noexcept {
        if constexpr (N <= 64) {
            data |= other.data;
        } else {
            for (size_t i = 0; i < data.size(); ++i) data[i] |= other.data[i];
        }
        return *this;
    }
    BitMask<N> operator|(const BitMask<N> &other) const noexcept {
        BitMask<N> result = *this;
        result |= other;
        return result;
    }

    BitMask<N> & operator^=(const BitMask<N> &other) noexcept {
        if constexpr (N <= 64) {
            data ^= other.data;
        } else {
            for (size_t i = 0; i < data.size(); ++i) data[i] ^= other.data[i];
        }
        return *this;
    }
    BitMask<N> operator^(const BitMask<N> &other) const noexcept {
        BitMask<N> result = *this;
        result ^= other;
        return result;
    }

    BitMask<N> operator~() const noexcept {
        BitMask<N> result;
        if constexpr (N <= 64) {
            result.data = static_cast<StorageType>(~data);
        } else {
            for (size_t i = 0; i < data.size(); ++i) result.data[i] = ~data[i];
        }
        result.clean_padding(); // CRITICAL to mask out excess bits
        return result;
    }

public: // ---- INTERFACE ----
    /**
     * @brief Gets the total capacity (number of bits) of the bitmask.
     */
    constexpr BitIndex numBits() const noexcept { return N; }

    /**
     * @brief Casts the bitmask to a different capacity size.
     * @tparam other_N The new size.
     * @return A new BitMask of size other_N.
     * @throws std::invalid_argument If casting to a smaller size and a set bit is truncated.
     */
    template<unsigned int other_N>
    BitMask<other_N> cast() const {
        if constexpr (N == other_N) {
            return *this;
        }
        if (N > other_N) {
            const BitIndex highest_set_bit = max_bit_set();
            if (highest_set_bit != N && highest_set_bit >= other_N) {
                throw std::invalid_argument("BitMask::cast: Cannot cast to a smaller BitMask with set bits outside the target range");
            }
        }

        BitMask<other_N> result;
        for (BitIndex i = 0; i < N; ++i) {
            if (isSet(i)) result.set(i);
        }
        return result;
    }

    /**
     * @brief Permutes the set bits according to a given mapping.
     * @param permutation The mapping array where the i-th bit goes to permutation[i].
     * @return A new permuted BitMask.
     */
    BitMask<N> permute(std::span<const BitIndex> permutation) const {
        BitMask<N> result;
        for (BitIndex i = 0; i < N; ++i) {
            if (i < permutation.size() && permutation[i] < N && isSet(i)) {
                result.set(permutation[i]);
            }
        }

        // Check count for errors
        if (result.count() != count()) {
            throw std::invalid_argument("BitMask::permute: Invalid permutation mapping; some bits were lost or duplicated.");
        }
        return result;
    }

    BitMask<N> permute(const std::initializer_list<BitIndex> &permutation) const {
        return permute(std::span<const BitIndex>(permutation));
    }

    BitMask<N>& permute_inplace(std::span<const BitIndex> permutation) {
        *this = permute(permutation);
        return *this;
    }

    BitMask<N>& permute_inplace(const std::initializer_list<BitIndex> &permutation) {
        return permute_inplace(std::span<const BitIndex>(permutation));
    }
    
    /**
     * @brief Sets the bit at the specified index to 1.
     * @param bit The index of the bit to set.
     */
    void set(BitIndex bit) noexcept {
        if (bit >= N) return;
        if constexpr (N <= 64) {
            data |= (static_cast<StorageType>(1) << bit);
        } else {
            data[bit / 64] |= (1ULL << (bit % 64));
        }
    }

    /**
     * @brief Sets all bits in the mask to 1.
     */
    void set() noexcept {
        if constexpr (N <= 64) {
            data = static_cast<StorageType>(~0ULL);
        } else {
            for (auto& block : data) block = ~0ULL;
        }
        clean_padding(); // CRITICAL
    }
    
    /**
     * @brief Clears (sets to 0) the bit at the specified index.
     * @param bit The index of the bit to clear.
     */
    void clear(BitIndex bit) noexcept {
        if (bit >= N) return;
        if constexpr (N <= 64) {
            data &= ~(static_cast<StorageType>(1) << bit);
        } else {
            data[bit / 64] &= ~(1ULL << (bit % 64));
        }
    }

    /**
     * @brief Clears all bits in the mask (sets to 0).
     */
    void clear() noexcept {
        if constexpr (N <= 64) {
            data = 0;
        } else {
            for (auto& block : data) block = 0;
        }
    }

    /**
     * @brief Checks if a specific bit is set.
     * @param bit The index to check.
     * @return True if the bit is 1, false if it is 0 or out of bounds.
     */
    bool isSet(BitIndex bit) const noexcept {
        if (bit >= N) return false;
        if constexpr (N <= 64) {
            return (data & (static_cast<StorageType>(1) << bit)) != 0;
        } else {
            return (data[bit / 64] & (1ULL << (bit % 64))) != 0;
        }
    }



    /**
     * @brief Checks if any bit is set in the mask.
     */
    bool any() const noexcept {
        if constexpr (N <= 64) {
            return (data != 0);
        } else {
            for (const auto& block : data) 
                if (block != 0) return true;
            return false;
        }
    }

    /**
     * @brief Checks if the mask is completely empty (no bits set).
     */
    bool empty() const noexcept {return !any();}

    /**
     * @brief Counts the total number of set bits (population count).
     */
    BitIndex count() const noexcept {
        if constexpr (N <= 64) {
            return static_cast<BitIndex>(std::popcount(data));
        } else {
            BitIndex total = 0;
            for (const auto& block : data) total += static_cast<BitIndex>(std::popcount(block));
            return total;
        }
    }

    /**
     * @brief Finds the index of the lowest set bit.
     * @return The index, or N if no bits are set.
     */
    BitIndex min_bit_set() const noexcept {
        if constexpr (N <= 64) {
            if (data == 0) return N;
            return static_cast<BitIndex>(std::countr_zero(data));
        } else {
            BitIndex base_bit = 0;
            for (const auto& block : data) {
                if (block != 0) {
                    BitIndex trailing_zeros = static_cast<BitIndex>(std::countr_zero(block));
                    BitIndex final_idx = base_bit + trailing_zeros;
                    return (final_idx < N) ? final_idx : N;
                }
                base_bit += 64;
            }
            return N; 
        }
    }

    /**
     * @brief Finds the index of the highest set bit.
     * @return The index, or N if no bits are set.
     */
    BitIndex max_bit_set() const noexcept {
        if constexpr (N <= 64) {
            if (data == 0) return N;
            return static_cast<BitIndex>(std::bit_width(static_cast<std::uint64_t>(data)) - 1);
        } else {
            for (std::size_t block_idx = data.size(); block_idx-- > 0;) {
                const uint64_t block = data[block_idx];
                if (block != 0) {
                    BitIndex leading_zeros = static_cast<BitIndex>(std::countl_zero(block));
                    BitIndex final_idx = static_cast<BitIndex>(block_idx * 64) + (63 - leading_zeros);
                    return (final_idx < N) ? final_idx : N;
                }
            }
            return N; 
        }
    }

public: // ---- I/O ----
    /**
     * @brief Converts the bitmask to a vector containing the indices of set bits.
     */
    std::vector<BitIndex> toVector() const {
        std::vector<BitIndex> result;
        result.reserve(count());
        for (const auto bit : *this) {
            result.emplace_back(bit);
        }
            
        
        return result;
    }

    /**
     * @brief Returns a binary string representation of the bitmask.
     */
    std::string toString() const {
        if constexpr (N == 0) return "";
        
        std::string result;
        result.reserve(N);

        for (int i = static_cast<int>(N) - 1; i >= 0; --i) {
            result += (isSet(i) ? '1' : '0');
        }
        return result;
    }

    friend std::ostream& operator<<(std::ostream &os, const BitMask<N>& bm) { 
        return (os << bm.toString()); 
    }



public: // --- COMPARISONS ----
    bool operator==(const BitMask<N>& other) const noexcept {
        return data == other.data;
    }

    bool operator!=(const BitMask<N>& other) const noexcept {
        return !(*this == other);
    }


public: // ---- ITERATOR OVER SET BITS ----
    /**
     * @class SetBitIterator
     * @brief An iterator that efficiently skips zero bits and only yields indices of set bits.
     */
    class SetBitIterator {
    private:
        const BitMask<N>& bm;
        BitIndex current_bit;

        void advance_to_next_set_bit() {
            // 1. Out of bounds check (acts as our end() sentinel)
            if (current_bit >= N) {
                current_bit = N;
                return;
            }

            if constexpr (N <= 64) {
                // Create a mask that zeroes out bits we've already passed.
                uint64_t mask = ~0ULL << current_bit; 
                uint64_t remaining = bm.data & mask;

                if (remaining == 0) {
                    current_bit = N;
                } else {
                    current_bit = static_cast<BitIndex>(std::countr_zero(remaining));
                }
            } else {
                size_t block_idx = current_bit / 64;
                BitIndex bit_in_block = current_bit % 64;
                
                // Mask current block to ignore already-processed bits
                uint64_t mask = ~0ULL << bit_in_block;
                uint64_t remaining = bm.data[block_idx] & mask;
                
                if (remaining != 0) {
                    current_bit = static_cast<BitIndex>((block_idx * 64) + std::countr_zero(remaining));
                    return;
                }
                
                // If the current block had no more set bits, scan the remaining blocks
                for (++block_idx; block_idx < bm.data.size(); ++block_idx) {
                    if (bm.data[block_idx] != 0) {
                        current_bit = static_cast<BitIndex>((block_idx * 64) + std::countr_zero(bm.data[block_idx]));
                        return;
                    }
                }
                
                // No set bits found
                current_bit = N;
            }
        }

    public:
        SetBitIterator(const BitMask<N>& bitmask, BitIndex start_bit) 
            : bm(bitmask), current_bit(start_bit) {
            advance_to_next_set_bit();
        }

        BitIndex operator*() const { return current_bit; }

        SetBitIterator& operator++() {
            ++current_bit; // Move past the current bit
            advance_to_next_set_bit(); // Fast-forward to the next set bit
            return *this;
        }

        bool operator!=(const SetBitIterator& other) const { 
            return current_bit != other.current_bit; 
        }
    };

    SetBitIterator begin() const { return SetBitIterator(*this, 0); }
    SetBitIterator end() const { return SetBitIterator(*this, N); } 
};

#endif