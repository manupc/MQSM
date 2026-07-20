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
#ifndef __SAMPLE_COUNT__H__
#define __SAMPLE_COUNT__H__

#include <algorithm>
#include <cstdint>
#include <format>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <private/BasicDataTypes.h>

using std::map;
using std::string;
using std::ostream;
using std::string_view;
using std::ostringstream;


// Advanced declarations
template<unsigned int max_q>
requires(max_q > 0 && max_q <= MAX_QUBITS_ALLOWED)
class MQSM;



/**
 * @class SampleCount
 * @brief Represents the result of taking multiple measurement shots of a quantum state.
 *
 * This class stores a histogram of measured bitstrings (as strings of '0' and '1')
 * and their corresponding occurrence counts. It provides utilities for filtering
 * and finding the most frequent outcomes.
 */
class SampleCount {

private: // ---------- REPRESENTATION ----------
    map<string, size_t> counts;

    SampleCount(const map<string, size_t>& counts) : counts(counts) {}
    SampleCount(map<string, size_t>&& counts) : counts(std::move(counts)) {}


public: // ---------- CONSTRUCTORS ----------
    /**
     * @brief Default constructor creating an empty sample count map.
     */
    SampleCount() = default;

    /**
     * @brief Default destructor.
     */
    ~SampleCount() = default;

    SampleCount(const SampleCount& other) = default;
    SampleCount(SampleCount&& other) noexcept = default;
    SampleCount& operator=(const SampleCount& other) = default;
    SampleCount& operator=(SampleCount&& other) noexcept = default;


public: // ---------- INTERFACE ----------

    /**
     * @brief Gets the maximum count among all sampled outcomes.
     * @return The highest occurrence count.
     */
    double getMaxCount() const {
        double max_count = 0.0;
        for (const auto& [label, count] : counts) {
            if (count > max_count) {
                max_count = count;
            }
        }
        return max_count;
    }

    /**
     * @brief Gets the bitstring label that occurred most frequently.
     * @return The most frequent bitstring label.
     */
    string getMaxCountKet() const {
        string max_ket;
        size_t max_count = 0;
        for (const auto& [label, count] : counts) {
            if (count > max_count) {
                max_count = count;
                max_ket = label;
            }
        }
        return max_ket;
    }

    /**
     * @brief Filters the sample counts, retaining only those above a minimum threshold.
     * @param min_count The minimum count required to retain an outcome.
     * @return A new SampleCount object with the filtered results.
     */
    SampleCount filterMinimumCount(size_t min_count) const {
        map<string, size_t> filtered_counts;
        for (const auto& [label, count] : counts) {
            if (count >= min_count) {
                filtered_counts[label] = count;
            }
        }
        return SampleCount(std::move(filtered_counts));
    }

    /**
     * @brief Returns a new SampleCount containing only the top N most frequent outcomes.
     * @param top_n The number of top outcomes to retain.
     * @return A new SampleCount object with the most frequent outcomes.
     */
    SampleCount mostFrequent(size_t top_n) const {
        if (top_n == 0) return SampleCount();

        // Create a vector of pairs and sort it by count in descending order
        std::vector<std::pair<string, size_t>> sorted_counts(counts.begin(), counts.end());
        std::sort(sorted_counts.begin(), sorted_counts.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        // Take the top N entries
        map<string, size_t> top_counts;
        for (size_t i = 0; i < std::min(top_n, static_cast<size_t>(sorted_counts.size())); ++i) {
            top_counts[sorted_counts[i].first] = sorted_counts[i].second;
        }

        return SampleCount(std::move(top_counts));
    }

    /**
     * @brief Gets the total number of shots (sum of all counts).
     */
    size_t totalCounts() const {
        size_t total = 0;
        for (const auto& [label, count] : counts) {
            total += count;
        }
        return total;
    }
    /**
     * @brief Gets the count for a specific bitstring label.
     * @param label The bitstring label (e.g., "010").
     * @return The number of times the label was measured.
     */
    size_t getCount(const string& label) const {
        auto it = counts.find(label);
        return (it != counts.end()) ? it->second : 0;
    }
    /**
     * @brief Gets the number of unique outcomes sampled.
     */
    size_t size() const { return counts.size(); }

    /**
     * @brief Checks if the sample count map is empty.
     */
    bool empty() const { return counts.empty(); }


public: // ---------- I/O --------------

    string toString(string_view padding= "") const {
        ostringstream result;
        for (const auto& [label, count] : counts) {
            result << std::format("{}|{}>: {}\n", padding, label, count);
        }
        return result.str();
    }


    friend ostream & operator<<(ostream &os, const SampleCount& sc) { return (os << sc.toString()); }


    

    template<unsigned int max_q>
    requires(max_q > 0 && max_q <= MAX_QUBITS_ALLOWED)
    friend class MQSM;

};


#endif