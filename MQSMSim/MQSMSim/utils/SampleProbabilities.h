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
#ifndef __SAMPLE_PROBABILITIES__H__
#define __SAMPLE_PROBABILITIES__H__

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
using std::ostringstream;
using std::string;
using std::ostream;
using std::string_view;


// Advanced declarations
template<unsigned int max_q>
requires(max_q > 0 && max_q <= MAX_QUBITS_ALLOWED)
class MQSM;


/**
 * @class SampleProbabilities
 * @brief Represents the exact theoretical measurement probabilities of a quantum state.
 *
 * This class stores a map of bitstrings to their corresponding exact probabilities 
 * as computed from the MQSM state vector.
 */
class SampleProbabilities {

private: // ---------- REPRESENTATION ----------
    map<string, double> probs;

    SampleProbabilities(const map<string, double>& probs) : probs(probs) {}
    SampleProbabilities(map<string, double>&& probs) : probs(std::move(probs)) {}


public: // ---------- CONSTRUCTORS ----------
    /**
     * @brief Default constructor creating an empty probability map.
     */
    SampleProbabilities() = default;

    /**
     * @brief Default destructor.
     */
    ~SampleProbabilities() = default;

    SampleProbabilities(const SampleProbabilities& other) = default;
    SampleProbabilities(SampleProbabilities&& other) noexcept = default;
    SampleProbabilities& operator=(const SampleProbabilities& other) = default;
    SampleProbabilities& operator=(SampleProbabilities&& other) noexcept = default;


public: // ---------- INTERFACE ----------

    /**
     * @brief Gets the maximum probability among all outcomes.
     * @return The highest probability value.
     */
    double getMaxProbability() const {
        double max_prob = 0.0;
        for (const auto& [label, prob] : probs) {
            if (prob > max_prob) {
                max_prob = prob;
            }
        }
        return max_prob;
    }

    /**
     * @brief Gets the bitstring label with the highest probability.
     * @return The most probable bitstring label.
     */
    string getMaxProbabilityKet() const {
        string max_ket;
        double max_prob = 0.0;
        for (const auto& [label, prob] : probs) {
            if (prob > max_prob) {
                max_prob = prob;
                max_ket = label;
            }
        }
        return max_ket;
    }

    /**
     * @brief Filters the probabilities, retaining only those above a minimum threshold.
     * @param min_prob The minimum probability required to retain an outcome.
     * @return A new SampleProbabilities object with the filtered results.
     */
    SampleProbabilities filterMinimumProbability(double min_prob) const {
        map<string, double> filtered_probs;
        for (const auto& [label, prob] : probs) {
            if (prob >= min_prob) {
                filtered_probs[label] = prob;
            }
        }
        return SampleProbabilities(std::move(filtered_probs));
    }

    /**
     * @brief Returns a new SampleProbabilities containing only the top N most probable outcomes.
     * @param top_n The number of top outcomes to retain.
     * @return A new SampleProbabilities object with the most probable outcomes.
     */
    SampleProbabilities mostFrequent(size_t top_n) const {
        if (top_n == 0) return SampleProbabilities();

        // Create a vector of pairs and sort it by count in descending order
        std::vector<std::pair<string, double>> sorted_probs(probs.begin(), probs.end());
        std::sort(sorted_probs.begin(), sorted_probs.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        // Take the top N entries
        map<string, double> top_probs;
        for (size_t i = 0; i < std::min(top_n, sorted_probs.size()); ++i) {
            top_probs[sorted_probs[i].first] = sorted_probs[i].second;
        }

        return SampleProbabilities(std::move(top_probs));
    }

    /**
     * @brief Gets the number of unique outcomes.
     */
    size_t size() const { return probs.size(); }
    
    /**
     * @brief Checks if the probability map is empty.
     */
    bool empty() const { return probs.empty(); }
    
    /**
     * @brief Gets the total sum of all probabilities (should typically be close to 1.0).
     */
    double totalProbability() const {
        double total = 0.0;
        for (const auto& [label, prob] : probs) {
            total += prob;
        }
        return total;
    }
    /**
     * @brief Gets the probability for a specific bitstring label.
     * @param label The bitstring label (e.g., "010").
     * @return The computed probability for the label.
     */
    double getProbability(const string& label) const {
        auto it = probs.find(label);
        return (it != probs.end()) ? it->second : 0.0;
    }
    

public: // ---------- I/O --------------

    string toString(string_view padding= "") const {
        ostringstream result;
        result.precision(6);
        for (const auto& [label, prob] : probs) {
            result << padding << "|"<<label<<">: " << prob << "\n";
        }
        return result.str();
    }


    friend ostream & operator<<(ostream &os, const SampleProbabilities& sp) { return (os << sp.toString()); }


    template<unsigned int max_q>
    requires(max_q > 0 && max_q <= MAX_QUBITS_ALLOWED)
    friend class MQSM;

};


#endif