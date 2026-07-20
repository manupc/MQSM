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
#ifndef __QCIRCUIT_LEVEL__H__
#define __QCIRCUIT_LEVEL__H__

#include <private/Numbers.h>
#include <private/BitMask.h>
#include <private/QInstruction.h>
#include <private/BasicDataTypes.h>

#include <algorithm>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>



// Advanced declarations
class MQGTSim;



/**
 * @class QCircuitLevel
 * @brief Represents a single execution level within a quantum circuit.
 *
 * A circuit level acts as a time-slice in a quantum circuit where multiple independent
 * quantum instructions (gates) can be executed simultaneously. It ensures that no two
 * instructions within the same level act on the same qubit.
 *
 * @tparam max_qubits Maximum number of qubits this circuit level can handle.
 */
template<unsigned int max_qubits>
class QCircuitLevel {

private: // ---------- REPRESENTATION ----------
    std::vector<QInstruction<max_qubits>> instructions;
    BitMask<max_qubits> used_qubits; 
    QubitIndex nqb; // Number of qubits in the level (<= max_qubits)

public: // ---------- CONSTRUCTORS ----------
    /**
     * @brief Default constructor creating an empty circuit level.
     */
    QCircuitLevel() = default;

    /**
     * @brief Default destructor.
     */
    ~QCircuitLevel() = default;

    /**
     * @brief Constructs a circuit level for a specific number of qubits.
     * @param num_qubits The number of qubits available in this level.
     * @throws std::invalid_argument if num_qubits exceeds max_qubits.
     */
    QCircuitLevel(QubitIndex num_qubits)  {
        if (num_qubits > max_qubits) {
            throw std::invalid_argument("QCircuitLevel::QCircuitLevel: Number of qubits exceeds max_qubits");
        }
        nqb = num_qubits;
    }

    // Copy & Move
    QCircuitLevel(const QCircuitLevel& other) = default;
    QCircuitLevel(QCircuitLevel&& other) noexcept = default;
    QCircuitLevel& operator=(const QCircuitLevel& other) = default;
    QCircuitLevel& operator=(QCircuitLevel&& other) noexcept = default;

public: // ---------- INTERFACE ----------
    /**
     * @brief Checks if an instruction can be inserted into this level without overlapping.
     * @param instr The quantum instruction to check.
     * @return True if it can be inserted, false otherwise.
     */
    bool can_insert(const QInstruction<max_qubits> &instr) const noexcept {
        BitMask<max_qubits> instr_qubits = BitMask<max_qubits>::rangeSet(instr.getMinimumQubitIndex(), instr.getMaximumQubitIndex());
        if (instr.getMaximumQubitIndex() >= nqb) {
            return false; // Instruction qubits are out of bounds for this level
        }
        return !(used_qubits & instr_qubits).any();
    }

    /**
     * @brief Adds an instruction to the level, ensuring no qubit conflicts.
     * @param instr The instruction to insert.
     * @throws std::invalid_argument If the instruction overlaps or is out of bounds.
     */
    void insert(const QInstruction<max_qubits> &instr) {

        // Error checking
        BitMask<max_qubits> instr_qubits = BitMask<max_qubits>::rangeSet(instr.getMinimumQubitIndex(), instr.getMaximumQubitIndex());
        if (instr.getMaximumQubitIndex() >= nqb) {
            throw std::invalid_argument("QCircuitLevel::insert: Instruction qubit index exceeds level's qubit count.");
        }
        if ((used_qubits & instr_qubits).any()) {
            throw std::invalid_argument("QCircuitLevel::insert: Instruction overlaps with existing instructions in the level.");
        }

        // Update level qubit mask
        used_qubits |= instr_qubits;

        // Insert instruction in sorted order based on minimum qubit index
        auto pos = std::lower_bound(instructions.begin(), instructions.end(), instr.getMinimumQubitIndex(),
            [](const QInstruction<max_qubits> &candidate, QubitIndex qubit) {
                return candidate.getMinimumQubitIndex() < qubit;
            }
        );
        instructions.insert(pos, instr);
    }


    /**
     * @brief Adds an instruction to the level (move semantics), ensuring no qubit conflicts.
     * @param instr The instruction to move-insert.
     * @throws std::invalid_argument If the instruction overlaps or is out of bounds.
     */
    void insert(QInstruction<max_qubits> &&instr) {

        // Error checking
        BitMask<max_qubits> instr_qubits = BitMask<max_qubits>::rangeSet(instr.getMinimumQubitIndex(), instr.getMaximumQubitIndex());
        if (instr.getMaximumQubitIndex() >= nqb) {
            throw std::invalid_argument("QCircuitLevel::insert: Instruction qubit index exceeds level's qubit count.");
        }
        if ((used_qubits & instr_qubits).any()) {
            throw std::invalid_argument("QCircuitLevel::insert: Instruction overlaps with existing instructions in the level.");
        }

        // Update level qubit mask
        used_qubits |= instr_qubits;

        // Insert instruction in sorted order based on minimum qubit index
        auto pos = std::lower_bound(instructions.begin(), instructions.end(), instr.getMinimumQubitIndex(),
            [](const QInstruction<max_qubits> &candidate, QubitIndex qubit) {
                return candidate.getMinimumQubitIndex() < qubit;
            }
        );
        instructions.insert(pos, std::move(instr));
    }

    /**
     * @brief Gets the bitmask of all qubits currently used by instructions in this level.
     */
    const BitMask<max_qubits>& getUsedQubits() const noexcept { 
        return used_qubits; 
    }
    
    /**
     * @brief Gets the list of instructions scheduled in this level.
     */
    const std::vector<QInstruction<max_qubits>>& getInstructions() const noexcept {
        return instructions; 
    }

    /**
     * @brief Gets the number of instructions in this level.
     */
    size_t getInstructionCount() const noexcept { 
        return static_cast<size_t>(instructions.size()); 
    }

    template<unsigned int other_max_qubits>
    QCircuitLevel<other_max_qubits> cast() const {
        QCircuitLevel<other_max_qubits> new_level(nqb);
        new_level.used_qubits = used_qubits.template cast<other_max_qubits>();
        new_level.instructions.reserve(instructions.size());
        for (const auto& instr : instructions) {
            new_level.instructions.push_back(instr.template cast<other_max_qubits>());
        }
        return new_level;
    }


    /**
     * @brief Computes the inverse of all instructions in this level.
     * @return A new QCircuitLevel containing the inverted instructions.
     */
    QCircuitLevel<max_qubits> inverse() const {
        QCircuitLevel<max_qubits> inv_level(nqb);
        inv_level.used_qubits = used_qubits;
        inv_level.instructions.reserve(instructions.size());
        for (const auto& instr : instructions) {
            inv_level.instructions.push_back(instr.inverse());
        }
        return inv_level;
    }

    QCircuitLevel<max_qubits> & inverse_inplace() noexcept {
        for (auto &instr : instructions) {
            instr.inverse_inplace();
        }
        return *this;
    }


    /**
     * @brief In-place permutation of the qubit indices for all instructions in this level.
     * @param permutation The mapping array indicating the new index for each original index.
     * @return A new permuted QCircuitLevel.
     */
    QCircuitLevel<max_qubits> permute(std::span<const QubitIndex> permutation) const {

        QCircuitLevel<max_qubits> new_level(nqb);
        new_level.used_qubits = used_qubits.permute(permutation);
        new_level.instructions.reserve(instructions.size());
        for (const auto& instr : instructions) {
            new_level.instructions.push_back(instr.permute(permutation));
        }
        auto less_by_min_qubit = [](const QInstruction<max_qubits> &a, const QInstruction<max_qubits> &b) {
            return a.getMinimumQubitIndex() < b.getMinimumQubitIndex();
        };
        std::sort(new_level.instructions.begin(), new_level.instructions.end(), less_by_min_qubit);
        return new_level;
    }


    template<unsigned int other_max_qubits>
    QCircuitLevel<other_max_qubits> permute(std::span<const QubitIndex> permutation) const {
        QCircuitLevel<other_max_qubits> new_level(nqb);
        new_level.used_qubits = used_qubits.permute(permutation).template cast<other_max_qubits>();
        new_level.instructions.reserve(instructions.size());
        for (const auto& instr : instructions) {
            new_level.instructions.push_back(instr.permute(permutation).template cast<other_max_qubits>());
        }
        auto less_by_min_qubit = [](const QInstruction<other_max_qubits> &a, const QInstruction<other_max_qubits> &b) {
            return a.getMinimumQubitIndex() < b.getMinimumQubitIndex();
        };
        std::sort(new_level.instructions.begin(), new_level.instructions.end(), less_by_min_qubit);
        return new_level;
    }

    QCircuitLevel<max_qubits> & permute_inplace(std::span<const QubitIndex> permutation) {
        used_qubits.permute_inplace(permutation);
        for (auto &instr : instructions) {
            instr.permute_inplace(permutation);
        }

        auto less_by_min_qubit = [](const QInstruction<max_qubits> &a, const QInstruction<max_qubits> &b) {
            return a.getMinimumQubitIndex() < b.getMinimumQubitIndex();
        };
        std::sort(instructions.begin(), instructions.end(), less_by_min_qubit);
        return *this;
    }



public: // ---------- I/O ----------
    std::string toString(std::string_view prefix = "") const {
        std::string result;
        for (const auto& instr : instructions) {
            result += prefix;
            result += instr.toString();
            result += "\n";
        }
        return result;
    }

    std::ostream & operator<<(std::ostream &os) const { return (os << toString()); }


private: // ------- FRIENDS -------

    friend class MQGTSim;

    template<unsigned int>
    friend class QCircuitLevel;


};

#endif