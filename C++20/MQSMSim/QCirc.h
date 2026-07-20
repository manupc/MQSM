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
#ifndef __QCIRCUIT_H__
#define __QCIRCUIT_H__

#include <private/BasicDataTypes.h>
#include <cstddef>
#include <array>
#include <initializer_list>
#include <iterator>
#include <algorithm>
#include <vector>
#include <span>
#include <stdexcept>
#include <string>
#include <format>
#include <ranges>
#include <type_traits> 

#include <private/QCircuitLevel.h>
#include <private/QInstruction.h>
#include <private/BasicDataTypes.h>

using std::array;
using std::initializer_list;
using std::vector;
using std::span;


// Advanced declarations
class MQGTSim;



/**
 * @class QCircuit
 * @brief Representation of a quantum circuit composed of discrete levels and instructions.
 *
 * QCircuit manages quantum instructions organized into levels (moments).
 * It supports standard gates, appending circuits, controlled operations,
 * and text/LaTeX visualizations.
 *
 * @tparam max_qubits Maximum number of qubits this circuit can support.
 */
template<unsigned int max_qubits>
class QCircuit {

public: // ---------- INTERFACE TYPES & PROXIES ----------

    /**
     * @class ControlQubits
     * @brief Proxy class to seamlessly unify single values, initializer lists, and spans for control qubits.
     *
     * This avoids having to write multiple overloads for every controlled gate.
     */
    class ControlQubits {
    private:
        span<const QubitIndex> m_span;
        QubitIndex m_single; // Local backup storage for single value lifetime tracking

    public:
        // Option 1: Single Value
        ControlQubits(QubitIndex q) noexcept : m_single(q) {
            m_span = span<const QubitIndex>(&m_single, 1);
        }

        // Option 2: Initializer List
        ControlQubits(initializer_list<QubitIndex> il) noexcept : m_span(il.begin(), il.end()), m_single(0) {}

        // Option 3: Contiguous Ranges (span, vector, array, etc.)
        template <std::ranges::contiguous_range R>
        requires std::same_as<std::remove_cvref_t<std::ranges::range_value_t<R>>, QubitIndex>
        ControlQubits(R&& range) noexcept : m_span(std::forward<R>(range)), m_single(0) {}

        span<const QubitIndex> to_span() const noexcept { return m_span; }
    };

private: // ---------- REPRESENTATION ----------
    vector<QCircuitLevel<max_qubits>> levels; // List of levels in the circuit, each containing a set of non-overlapping instructions
    array<size_t, max_qubits> last_used_level_per_qubit; // Cache to track the last level index where each qubit was used, for efficient insertion of new instructions
    QubitIndex nqb; // Number of qubits in the circuit (<= max_qubits)

    vector<QCircuitLevel<max_qubits>>& getLevels() noexcept { return levels; }
    const vector<QCircuitLevel<max_qubits>>& getLevels() const noexcept { return levels; }


private: // ---------- PRIVATE METHODS: QUBIT CACHE AND INSTRUCTIONS ----------
    void initializeLevelIndexCache() noexcept {
        last_used_level_per_qubit.fill(-1);
    }

    void rebuildLevelIndexCache() noexcept {
        initializeLevelIndexCache();
        for (size_t level_idx = 0; level_idx < levels.size(); ++level_idx) {
            for (QubitIndex q : levels[level_idx].getUsedQubits()) {
                last_used_level_per_qubit[q] = level_idx;
            }
        }
    }

    template <unsigned int other_max_qubits>
    void updateLevelIndexCache(const BitMask<other_max_qubits>& touched_qubits, size_t level_idx) noexcept {
        for (const QubitIndex q : touched_qubits) {
            last_used_level_per_qubit[q] = level_idx;
        }
    }

    template <typename T>
    requires std::same_as<std::decay_t<T>, QInstruction<max_qubits>>
    QCircuit<max_qubits> & addInstruction(T&& instr) noexcept {
        const BitMask<max_qubits> instr_qubits = BitMask<max_qubits>::rangeSet(instr.getMinimumQubitIndex(), instr.getMaximumQubitIndex());

        size_t insertion_level_idx = 0;
        for (QubitIndex q : instr_qubits) {
            const size_t last_level = last_used_level_per_qubit[q];
            insertion_level_idx = std::max(insertion_level_idx, last_level + 1);
        }

        if (insertion_level_idx == levels.size()) {
            levels.emplace_back(nqb);
        }

        levels[insertion_level_idx].insert(std::forward<T>(instr));
        updateLevelIndexCache(instr_qubits, insertion_level_idx);
        
        return *this;
    }


    QInstruction<max_qubits> make_instruction(QInstructionType type,QubitIndex target1, QubitIndex target2= 0, const Real &parameter= 0, const ControlQubits &ctrl={}, std::string_view ctrl_values="") {

        BitMask<max_qubits> ctrl_mask;
        BitMask<max_qubits> ctrl_value_mask;
        QubitIndex idx= 0;

        for (const QubitIndex qb : ctrl.to_span()) {
            if (qb >= nqb) {
                throw std::invalid_argument("QCircuit::make_instruction: Control qubit index out of bounds");
            }
            ctrl_mask.set(qb);
            if (ctrl_values.empty() || ctrl_values[idx] == '1') {
                ctrl_value_mask.set(qb);
            }
            ++idx;
        }

        if (target1 >= nqb || (isTwoQubitGate(type) && target2 >= nqb)) {
            throw std::invalid_argument("QCircuit::make_instruction: Target qubit index out of bounds");
        }
        if (isControlledGate(type) && !ctrl_mask.any()) {
            throw std::invalid_argument("QCircuit::make_instruction: Controlled gate must have at least one control qubit");
        }
        if (!isControlledGate(type) && ctrl_mask.any()) {
            throw std::invalid_argument("QCircuit::make_instruction: Non-controlled gate cannot have control qubits");
        }

        switch(type) {
            case QInstructionType::X: return QInstruction<max_qubits>::X(target1);
            case QInstructionType::Y: return QInstruction<max_qubits>::Y(target1);
            case QInstructionType::Z: return QInstruction<max_qubits>::Z(target1);
            case QInstructionType::H: return QInstruction<max_qubits>::H(target1);
            case QInstructionType::S: return QInstruction<max_qubits>::S(target1);
            case QInstructionType::Sdg: return QInstruction<max_qubits>::Sdg(target1);
            case QInstructionType::T: return QInstruction<max_qubits>::T(target1);
            case QInstructionType::Tdg: return QInstruction<max_qubits>::Tdg(target1);
            case QInstructionType::SX: return QInstruction<max_qubits>::SX(target1);
            case QInstructionType::SXdg: return QInstruction<max_qubits>::SXdg(target1);

            case QInstructionType::P: return QInstruction<max_qubits>::P(target1, parameter);
            case QInstructionType::Rz: return QInstruction<max_qubits>::Rz(target1, parameter);
            case QInstructionType::Rx: return QInstruction<max_qubits>::Rx(target1, parameter);
            case QInstructionType::Ry: return QInstruction<max_qubits>::Ry(target1, parameter);

            case QInstructionType::SWAP: return QInstruction<max_qubits>::SWAP(target1, target2);
            case QInstructionType::Rzz: return QInstruction<max_qubits>::Rzz(target1, target2, parameter);
            case QInstructionType::Rxx: return QInstruction<max_qubits>::Rxx(target1, target2, parameter);
            case QInstructionType::Ryy: return QInstruction<max_qubits>::Ryy(target1, target2, parameter);

            case QInstructionType::CX: return QInstruction<max_qubits>::CX(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CY: return QInstruction<max_qubits>::CY(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CZ: return QInstruction<max_qubits>::CZ(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CH: return QInstruction<max_qubits>::CH(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CS: return QInstruction<max_qubits>::CS(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CSdg: return QInstruction<max_qubits>::CSdg(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CT: return QInstruction<max_qubits>::CT(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CTdg: return QInstruction<max_qubits>::CTdg(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CSX: return QInstruction<max_qubits>::CSX(target1, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CSXdg: return QInstruction<max_qubits>::CSXdg(target1, ctrl_mask, ctrl_value_mask);

            case QInstructionType::CP: return QInstruction<max_qubits>::CP(target1, parameter, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CRz: return QInstruction<max_qubits>::CRz(target1, parameter, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CRx: return QInstruction<max_qubits>::CRx(target1, parameter, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CRy: return QInstruction<max_qubits>::CRy(target1, parameter, ctrl_mask, ctrl_value_mask);

            case QInstructionType::CSWAP: return QInstruction<max_qubits>::CSWAP(target1, target2, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CRzz: return QInstruction<max_qubits>::CRzz(target1, target2, parameter, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CRxx: return QInstruction<max_qubits>::CRxx(target1, target2, parameter, ctrl_mask, ctrl_value_mask);
            case QInstructionType::CRyy: return QInstruction<max_qubits>::CRyy(target1, target2, parameter, ctrl_mask, ctrl_value_mask);

            default:
                throw std::invalid_argument("QCircuit::make_instruction: Unsupported instruction type");
        }

    }

    template <unsigned int other_max_qubits>
    void appendLevelsAndRefreshIndex(const vector<QCircuitLevel<other_max_qubits>> &other_levels) {
        if (other_levels.empty()) {
            return;
        }

        const size_t offset = levels.size();
        levels.reserve(levels.size() + other_levels.size());
        if constexpr (other_max_qubits == max_qubits) {
            levels.insert(levels.end(), other_levels.begin(), other_levels.end());
            for (size_t rel_idx = 0; rel_idx < other_levels.size(); ++rel_idx) {
                updateLevelIndexCache(other_levels[rel_idx].getUsedQubits(), offset + rel_idx);
            }
        } else {
            for (const auto& lev : other_levels) {
                levels.push_back(lev.template cast<max_qubits>());
            }
            for (size_t rel_idx = 0; rel_idx < other_levels.size(); ++rel_idx) {
                updateLevelIndexCache(levels[offset + rel_idx].getUsedQubits(), offset + rel_idx);
            }
        }
    }

    template <unsigned int other_max_qubits>
    void appendPermutedLevelsAndRefreshIndex(
        const vector<QCircuitLevel<other_max_qubits>> &other_levels,
        span<const QubitIndex> qubit_mapping
    ) {
        if (other_levels.empty()) {
            return;
        }

        const size_t offset = levels.size();
        levels.reserve(levels.size() + other_levels.size());

        for (const auto& other_level : other_levels) {
            levels.push_back(other_level.template permute<max_qubits>(qubit_mapping));
        }

        for (size_t rel_idx = 0; rel_idx < other_levels.size(); ++rel_idx) {
            updateLevelIndexCache(levels[offset + rel_idx].getUsedQubits(), offset + rel_idx);
        }
    }

public: // ---------- CONSTRUCTORS ----------
    /**
     * @brief Constructs a new quantum circuit.
     *
     * @param num_qubits The number of active qubits in the circuit.
     *
     * @pre num_qubits <= max_qubits
     * @post An empty quantum circuit with the specified number of qubits is created.
     *
     * @par Example:
     * @code
     * QCircuit<5> qc(3); // A circuit with 3 active qubits, max capacity 5.
     * @endcode
     */
    QCircuit(QubitIndex num_qubits) : nqb(num_qubits) { 
        initializeLevelIndexCache(); 
        if (nqb > max_qubits) {
            throw std::invalid_argument("QCircuit::QCircuit: Number of qubits exceeds max_qubits");
        }
    }

    ~QCircuit() = default;

    QCircuit(const QCircuit& other) = default;
    QCircuit(QCircuit&& other) noexcept = default;
    QCircuit& operator=(const QCircuit& other) = default;
    QCircuit& operator=(QCircuit&& other) noexcept = default;


public: // ---------- INTERFACE ----------
    /**
     * @brief Gets the maximum number of qubits allowed.
     */
    constexpr QubitIndex getMaxNumQubits() const noexcept { return max_qubits; }
    
    /**
     * @brief Gets the current number of active qubits.
     */
    QubitIndex getNumQubits() const noexcept { return nqb; }
    
    /**
     * @brief Checks if the circuit is empty (has no instructions or qubits).
     */
    bool empty() const noexcept { return nqb==0 || levels.empty(); }
    size_t getNumLevels() const noexcept { return levels.size(); }
    size_t getTotalInstructionCount() const noexcept {
        size_t total = 0;
        for (const auto& level : levels) {
            total += level.getInstructionCount();
        }
        return total;
    }

    /**
     * @brief Creates a new quantum circuit that is the inverse (adjoint) of this circuit.
     *
     * Reverses the order of levels and inverts all instructions.
     *
     * @return A new QCircuit representing the inverse operation.
     *
     * @pre Circuit must be valid.
     * @post The returned circuit applies the adjoint operation.
     *
     * @par Example:
     * @code
     * QCircuit<3> inv_qc = qc.inverse();
     * @endcode
     */
    QCircuit<max_qubits> inverse() const {
        QCircuit<max_qubits> inv_circuit(nqb);
        inv_circuit.levels.reserve(levels.size());

        for (const auto& level : std::views::reverse(levels)) {
            inv_circuit.levels.push_back(level.inverse());
        }
        inv_circuit.rebuildLevelIndexCache();
        return inv_circuit;
    }

    /**
     * @brief Inverts the current circuit in place.
     *
     * Reverses the order of levels and inverts all instructions within this circuit.
     *
     * @return A reference to this modified circuit.
     *
     * @pre Circuit must be valid.
     * @post This circuit now represents its own adjoint operation.
     */
    QCircuit<max_qubits> & inverse_inplace() noexcept {
        for (auto& level : levels) {
            level.inverse_inplace();
        }

        // Reverse the order of levels in place
        std::reverse(levels.begin(), levels.end());
        rebuildLevelIndexCache();
        return *this;
    }


    /**
     * @brief Creates a controlled version of this entire circuit.
     *
     * @param ctrl The control qubits.
     * @param ctrl_values Optional string of '0' or '1' specifying the control state (default all '1's).
     * @return A new QCircuit where all operations are controlled.
     *
     * @pre Control qubits must be valid and must not overlap with active qubits in the circuit.
     * @post A new controlled circuit is returned.
     *
     * @par Example:
     * @code
     * auto cqc = qc.controlled({0, 1}, "10");
     * @endcode
     */
    QCircuit<max_qubits> controlled(ControlQubits ctrl, std::string_view ctrl_values="") const {

        BitMask<max_qubits> ctrl_mask;
        BitMask<max_qubits> ctrl_value_mask;
        QubitIndex idx= 0;

        for (const QubitIndex qb : ctrl.to_span()) {
            if (qb >= nqb) {
                throw std::invalid_argument("QCircuit::controlled: Control qubit index out of bounds");
            }
            if (last_used_level_per_qubit[qb] >= 0) {
                throw std::invalid_argument("QCircuit::controlled: Control qubit overlaps with existing instructions in the circuit");
            }
            ctrl_mask.set(qb);
            if (ctrl_values.empty() || ctrl_values[idx] == '1') {
                ctrl_value_mask.set(qb);
            }
            ++idx;
        }

        QCircuit<max_qubits> new_circuit;
        new_circuit.levels.reserve(levels.size());
        for (const auto& level : levels) {
            for (const auto& instr : level.getInstructions()) {
                new_circuit.addInstruction(instr.controlled(ctrl_mask, ctrl_value_mask));
            }
        }
        
        return new_circuit;
    }



    /**
     * @brief Appends another quantum circuit to the end of this one.
     *
     * @tparam other_max_qubits Max qubits of the other circuit.
     * @param other The circuit to append.
     * @param qubit_mapping Optional mapping from the other circuit's qubits to this circuit's qubits.
     * @return A reference to this modified circuit.
     *
     * @pre other.getNumQubits() <= this->getNumQubits()
     * @post The levels of the other circuit are appended.
     *
     * @par Example:
     * @code
     * qc1.append(qc2, {2, 3}); // Appends qc2 mapping its q0->q2 and q1->q3
     * @endcode
     */
    template<unsigned int other_max_qubits>
    QCircuit<max_qubits> & append(const QCircuit<other_max_qubits> &other, span<const QubitIndex> qubit_mapping = {}) {
        
        // Error checking
        if (!qubit_mapping.empty() && qubit_mapping.size() != other.nqb) {
            throw std::invalid_argument("QCircuit::append: Qubit mapping size does not match the number of qubits in the other circuit");
        }

        if (nqb < other.nqb) {
            throw std::invalid_argument("QCircuit::append: Current circuit does not have enough qubits (" + std::to_string(nqb) + ") to accommodate the other circuit (" + std::to_string(other.nqb) + ")");
        }

        if (!qubit_mapping.empty() && qubit_mapping.size() != other.nqb) {
            throw std::invalid_argument("QCircuit::append: Qubit mapping size (" + std::to_string(qubit_mapping.size()) + ") does not match the number of qubits in the other circuit (" + std::to_string(other.nqb) + ")");
        }

        if (!qubit_mapping.empty()) {
            BitMask<max_qubits> seen;
            for (QubitIndex mapped : qubit_mapping) {
                if (mapped >= max_qubits) {
                    throw std::invalid_argument("QCircuit::append: Mapped qubit index out of bounds");
                }
                if (seen.isSet(mapped)) {
                    throw std::invalid_argument("QCircuit::append: Qubit mapping contains duplicates");
                }
                seen.set(mapped);
            }

            if (seen.max_bit_set() >= nqb) {
                throw std::invalid_argument("QCircuit::append: Mapped qubit index " + std::to_string(seen.max_bit_set()) + " exceeds the number of qubits in the current circuit (" + std::to_string(nqb) + ")");
            }
        }

        if (qubit_mapping.empty()) {
            appendLevelsAndRefreshIndex(other.levels);
        } else {
            appendPermutedLevelsAndRefreshIndex(other.levels, qubit_mapping);
        }


        return *this;
    }

    /**
     * @brief Appends another quantum circuit using an initializer list for the mapping.
     */
    template<unsigned int other_max_qubits>
    QCircuit<max_qubits> & append(const QCircuit<other_max_qubits> &other, initializer_list<QubitIndex> qubit_mapping) {
        return append<other_max_qubits>(other, span<const QubitIndex>(qubit_mapping));
    }

    /**
     * @name Standard Gate Factories
     * @brief Methods to append standard single and multi-qubit gates to the circuit.
     *
     * Each method adds the corresponding quantum instruction to the circuit.
     * @param target The target qubit index.
     * @return A reference to this circuit to allow chaining.
     *
     * @pre Target qubits must be < getNumQubits().
     * @post The corresponding instruction is appended to the circuit.
     *
     * @par Example:
     * @code
     * qc.h(0).cx({0}, 1).x(2);
     * @endcode
     */
    ///@{
    
    QCircuit<max_qubits> & x(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::X, target)); }
    QCircuit<max_qubits> & y(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::Y, target)); }
    QCircuit<max_qubits> & z(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::Z, target)); }
    QCircuit<max_qubits> & h(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::H, target)); }
    QCircuit<max_qubits> & s(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::S, target)); }
    QCircuit<max_qubits> & sdg(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::Sdg, target)); }
    QCircuit<max_qubits> & t(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::T, target)); }
    QCircuit<max_qubits> & tdg(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::Tdg, target)); }
    QCircuit<max_qubits> & sx(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::SX, target)); }
    QCircuit<max_qubits> & sxdg(QubitIndex target) { return addInstruction(make_instruction(QInstructionType::SXdg, target)); }

    QCircuit<max_qubits> & p(QubitIndex target, Real angle) { return addInstruction(make_instruction(QInstructionType::P, target, 0, angle)); }
    QCircuit<max_qubits> & rz(QubitIndex target, Real angle) { return addInstruction(make_instruction(QInstructionType::Rz, target, 0, angle)); }
    QCircuit<max_qubits> & rx(QubitIndex target, Real angle) { return addInstruction(make_instruction(QInstructionType::Rx, target, 0, angle)); }
    QCircuit<max_qubits> & ry(QubitIndex target, Real angle) { return addInstruction(make_instruction(QInstructionType::Ry, target, 0, angle)); }

    QCircuit<max_qubits> & swap(QubitIndex target1, QubitIndex target2) { return addInstruction(make_instruction(QInstructionType::SWAP, target1, target2)); }
    QCircuit<max_qubits> & rzz(QubitIndex target1, QubitIndex target2, Real angle) { return addInstruction(make_instruction(QInstructionType::Rzz, target1, target2, angle)); }
    QCircuit<max_qubits> & rxx(QubitIndex target1, QubitIndex target2, Real angle) { return addInstruction(make_instruction(QInstructionType::Rxx, target1, target2, angle)); }
    QCircuit<max_qubits> & ryy(QubitIndex target1, QubitIndex target2, Real angle) { return addInstruction(make_instruction(QInstructionType::Ryy, target1, target2, angle)); }
    ///@}

    /**
     * @name Controlled Gate Factories
     * @brief Methods to append controlled standard gates.
     *
     * @param ctrl Control qubits.
     * @param target Target qubit.
     * @param ctrl_values Optional string of control values ('0' or '1').
     * @return A reference to this circuit.
     *
     * @pre Target and Control qubits must be valid and disjoint.
     */
    ///@{
    // Controlled factory methods updated to use ControlQubits proxy
    QCircuit<max_qubits> & cx(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CX, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & cy(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CY, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & cz(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CZ, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & ch(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CH, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & cs(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CS, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & csdg(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CSdg, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & ct(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CT, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & ctdg(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CTdg, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & csx(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CSX, target, 0, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & csxdg(ControlQubits ctrl, QubitIndex target, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CSXdg, target, 0, 0, ctrl, ctrl_values)); }

    QCircuit<max_qubits> & cp(ControlQubits ctrl, QubitIndex target, Real angle, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CP, target, 0, angle, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & crz(ControlQubits ctrl, QubitIndex target, Real angle, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CRz, target, 0, angle, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & crx(ControlQubits ctrl, QubitIndex target, Real angle, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CRx, target, 0, angle, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & cry(ControlQubits ctrl, QubitIndex target, Real angle, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CRy, target, 0, angle, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & cswap(ControlQubits ctrl, QubitIndex target1, QubitIndex target2, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CSWAP, target1, target2, 0, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & crzz(ControlQubits ctrl, QubitIndex target1, QubitIndex target2, Real angle, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CRzz, target1, target2, angle, ctrl, ctrl_values)); }  
    QCircuit<max_qubits> & crxx(ControlQubits ctrl, QubitIndex target1, QubitIndex target2, Real angle, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CRxx, target1, target2, angle, ctrl, ctrl_values)); }
    QCircuit<max_qubits> & cryy(ControlQubits ctrl, QubitIndex target1, QubitIndex target2, Real angle, std::string_view ctrl_values = "") { return addInstruction(make_instruction(QInstructionType::CRyy, target1, target2, angle, ctrl, ctrl_values)); }
    ///@}


public: // ------------ I/O -------------
    /**
     * @brief Generates an ASCII representation of the quantum circuit.
     *
     * @return A string containing the ASCII art of the circuit.
     *
     * @pre None
     * @post Circuit remains unchanged.
     */
    std::string draw() const {
        const size_t qubit_label_width = std::to_string(nqb == 0 ? 0 : nqb - 1).size();
        auto make_row_prefix = [qubit_label_width](QubitIndex q) -> std::string {
            const std::string q_index = std::to_string(q);
            std::string prefix = "q_" + q_index + ":";
            prefix += std::string(qubit_label_width - q_index.size(), ' ');
            prefix += " --";
            return prefix;
        };

        if (levels.empty()) {
            std::string empty_circuit;
            for (QubitIndex q = 0; q < nqb; ++q) {
                empty_circuit += make_row_prefix(q);
                empty_circuit += "-\n";
            }
            return empty_circuit;
        }

        const size_t num_lines = 2 * nqb - 1;
        std::vector<std::string> lines(num_lines);

        for (QubitIndex q = 0; q < nqb; ++q) {
            lines[2 * q] = make_row_prefix(q);
            if (q < nqb - 1) {
                lines[2 * q + 1] = std::string(lines[2 * q].length(), ' ');
            }
        }

        auto format_centered_slot = [](const std::string& content, size_t width, char fill_char = '-') -> std::string {
            if (width <= content.length()) {
                return content;
            }

            const size_t padding = width - content.length();
            const size_t left_padding = padding / 2;
            const size_t right_padding = padding - left_padding;

            return std::string(left_padding, fill_char) + content + std::string(right_padding, fill_char);
        };

        auto format_gate_slot = [&](const std::string& name, size_t width) -> std::string {
            return format_centered_slot("[" + name + "]", width, '-');
        };

        auto format_control_slot = [&](char symbol, size_t width) -> std::string {
            return format_centered_slot(std::string(1, symbol), width, '-');
        };

        auto format_vertical_slot = [&](char symbol, size_t width) -> std::string {
            return format_centered_slot(std::string(1, symbol), width, ' ');
        };

        auto instructionLabel = [](const QInstruction<max_qubits>& instr) -> std::string {
            std::string label= instruction_label(instr.getType());
            if (instr.isParametrized()) {
                return std::format("{}({})", label, instr.getParameter().toString());
            }
            return label;
        };

        for (const auto& level : levels) {
            size_t slot_width = 7;
            for (const auto& instr : level.getInstructions()) {
                const std::string name = instructionLabel(instr);
                slot_width = std::max(slot_width, name.length() + 4);
            }
            if ((slot_width % 2) == 0) {
                ++slot_width;
            }

            std::vector<std::string> col_qubits(nqb, std::string(slot_width, '-'));
            std::vector<std::string> col_spaces(nqb - 1, std::string(slot_width, ' '));

            for (const auto& instr : level.getInstructions()) {
                const std::string name = instructionLabel(instr);
                QubitIndex min_q = nqb;
                QubitIndex max_q = 0;
                const auto& control_values = instr.getControlValueMask();

                if (instr.isTwoQubitGate()) {
                    const auto [target1, target2] = instr.getTargetPair();
                    min_q = std::min(min_q, target1);
                    max_q = std::max(max_q, target1);
                    col_qubits[target1] = format_gate_slot(name, slot_width);

                    min_q = std::min(min_q, target2);
                    max_q = std::max(max_q, target2);
                    col_qubits[target2] = format_gate_slot(name, slot_width);
                } else {
                    const QubitIndex target = instr.getTarget();
                    min_q = std::min(min_q, target);
                    max_q = std::max(max_q, target);
                    col_qubits[target] = format_gate_slot(name, slot_width);
                }

                for (const QubitIndex c : instr.getControlMask()) {
                    min_q = std::min(min_q, c);
                    max_q = std::max(max_q, c);
                    col_qubits[c] = format_control_slot(control_values.isSet(c) ? 'x' : 'o', slot_width);
                }

                if (max_q > min_q) {
                    for (QubitIndex s = min_q; s < max_q; ++s) {
                        col_spaces[s] = format_vertical_slot('|', slot_width);
                    }
                    for (QubitIndex r = min_q + 1; r < max_q; ++r) {
                        if (col_qubits[r] == std::string(slot_width, '-')) {
                            col_qubits[r] = format_control_slot('|', slot_width);
                        }
                    }
                }
            }

            for (QubitIndex q = 0; q < nqb; ++q) {
                lines[2 * q] += col_qubits[q];
                if (q < nqb - 1) {
                    lines[2 * q + 1] += col_spaces[q];
                }
            }
        }

        for (QubitIndex q = 0; q < nqb; ++q) {
            lines[2 * q] += "-";
        }

        std::string final_rendering;
        final_rendering.reserve(lines.size() * lines[0].size() * 2); 
        for (const auto& line : lines) {
            final_rendering += line + "\n";
        }

        return final_rendering;
    }

    /**
     * @brief Generates a LaTeX representation of the circuit using the Quantikz package.
     *
     * @return A string containing the LaTeX code.
     */
    std::string toQuantikz() const {
        if (nqb == 0) {
            return "\\begin{quantikz}\n\\end{quantikz}\n";
        }

        // Initialize a 2D matrix: [nqb][levels + 1]
        // Column 0 is reserved for the qubit start labels (\lstick)
        // All other cells default to standard quantum wires (\qw)
        std::vector<std::vector<std::string>> grid(nqb, std::vector<std::string>(levels.size() + 1, "\\qw"));

        // Setup row prefixes
        for (QubitIndex q = 0; q < nqb; ++q) {
            grid[q][0] = std::format("\\lstick{{q_{}}}", q);
        }

        // Reused lambda to extract formatted labels consistently
        auto instructionLabel = [](const QInstruction<max_qubits>& instr) -> std::string {
            std::string label = instruction_label(instr.getType());
            if (instr.isParametrized()) {
                return std::format("{}({})", label, instr.getParameter().toString());
            }
            return label;
        };

        // Populate the matrix level by level
        for (size_t level_idx = 0; level_idx < levels.size(); ++level_idx) {
            const size_t col_idx = level_idx + 1;
            const auto& level = levels[level_idx];

            for (const auto& instr : level.getInstructions()) {
                const std::string name = instructionLabel(instr);
                QubitIndex main_target = 0;

                // 1. Process Targets
                if (instr.isTwoQubitGate()) {
                    const auto [t1, t2] = instr.getTargetPair();
                    const QubitIndex top = std::min(t1, t2);
                    const QubitIndex bottom = std::max(t1, t2);
                    main_target = top;
                    
                    // Draw the gate on both targets and link them with a vertical wire (\vqw)
                    grid[top][col_idx] = std::format("\\gate{{{}}} \\vqw{{{}}}", name, bottom - top);
                    grid[bottom][col_idx] = std::format("\\gate{{{}}}", name);
                } else {
                    main_target = instr.getTarget();
                    grid[main_target][col_idx] = std::format("\\gate{{{}}}", name);
                }

                // 2. Process Controls (Closed and Open)
                const auto& control_values = instr.getControlValues();
                for (const QubitIndex c : instr.getControls()) {
                    // Offset handles both upward (negative) and downward (positive) pointing controls
                    const int offset = static_cast<int>(main_target) - static_cast<int>(c);
                    const std::string ctrl_macro = control_values.isSet(c) ? "\\ctrl" : "\\octrl";
                    
                    grid[c][col_idx] = std::format("{}{{{}}}", ctrl_macro, offset);
                }
            }
        }

        // Compile the matrix into a final LaTeX formatted string
        std::string result = "\\begin{quantikz}\n";
        for (QubitIndex q = 0; q < nqb; ++q) {
            for (size_t col = 0; col <= levels.size(); ++col) {
                result += grid[q][col];
                if (col < levels.size()) {
                    result += " & ";
                }
            }
            // End of row syntax for TikZ/Quantikz
            if (q < nqb - 1) {
                result += " \\\\\n";
            } else {
                result += "\n";
            }
        }
        result += "\\end{quantikz}\n";

        return result;
    }

    /**
     * @brief Serializes the circuit level by level into a string.
     */
    std::string toString() const {
        std::string result;
        result.reserve(levels.size() * 48); 
        
        QubitIndex level_index = 0;
        for (const auto& level : levels) {
            std::format_to(std::back_inserter(result), "Level {}:\n{}\n", level_index++, level.toString("  "));
        }
        return result;
    }

    friend std::ostream & operator<<(std::ostream &os, const QCircuit<max_qubits> &circuit) { 
        return (os << circuit.toString()); 
    }

    template <unsigned int>
    friend class QCircuit;

    friend class MQGTSim;
};



typedef QCircuit<MAX_QUBITS_ALLOWED> QCirc;
typedef QCircuit<64> QCirc64;
typedef QCircuit<32> QCirc32;
typedef QCircuit<16> QCirc16;
typedef QCircuit<8> QCirc8;



#endif