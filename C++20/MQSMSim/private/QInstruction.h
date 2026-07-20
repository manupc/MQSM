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
#ifndef __QINSTRUCTION__H__
#define __QINSTRUCTION__H__

#include <stdexcept>
#include <span>
#include <format>
#include <iterator>
#include <array>
#include <utility>
#include <vector>
#include <string>
#include <algorithm>
#include <private/BitMask.h>
#include <private/BasicDataTypes.h>

using std::span;
using std::string;
using std::pair;
using std::vector;
using std::array;


/**
 * @enum QInstructionType
 * @brief Defines the supported quantum gate types.
 *
 * Each gate type is encoded as a bitmask where specific bits indicate properties
 * like being a controlled gate, a two-qubit gate, or a parameterized gate.
 */
enum class QInstructionType : uint8_t {
    I    = 0b00010000,
    Z    = 0b00010001,
    S    = 0b00010010,
    Sdg  = 0b00010011,
    T    = 0b00010100,
    Tdg  = 0b00010101,
    X    = 0b00100000,
    Y    = 0b00100001,
    P    = 0b00011000,
    Rz   = 0b00011001,
    Rx   = 0b00001000,
    Ry   = 0b00001001,
    H    = 0b00000000,
    SX   = 0b00000001,
    SXdg = 0b00000010,
    SWAP = 0b01000000,
    Rzz  = 0b01011000,
    Rxx  = 0b01001000,
    Ryy  = 0b01001001,

    CZ   = 0b10010001,
    CS   = 0b10010010,
    CSdg = 0b10010011,
    CT   = 0b10010100,
    CTdg = 0b10010101,
    CX   = 0b10100000,
    CY   = 0b10100001,
    CP   = 0b10011000,
    CRz  = 0b10011001,
    CRx  = 0b10001000,
    CRy  = 0b10001001,
    CH   = 0b10000000,
    CSX  = 0b10000001,
    CSXdg= 0b10000010,
    CSWAP= 0b11000000,
    CRzz = 0b11011000,
    CRxx = 0b11001000,
    CRyy = 0b11001001
};

using QInstructionTypeNumber = std::underlying_type<QInstructionType>::type; // Get the underlying type of the enum class

constexpr QInstructionTypeNumber ControlledGateMask    = 0b10000000; // Mask to identify controlled gates
constexpr QInstructionTypeNumber TwoQubitGateMask      = 0b01000000; // Mask to identify two-qubit gates
constexpr QInstructionTypeNumber ParameterizedGateMask = 0b00001000; // Mask to identify parameterized gates

inline bool isAntidiagonalGate(QInstructionType gate) {
    return (static_cast<QubitIndex>(gate) & 0b00100000)!= 0;
}

inline bool isDiagonalGate(QInstructionType gate) {
    return (static_cast<QubitIndex>(gate) & 0b00010000)!= 0;
}


inline bool isControlledGate(QInstructionType gate) {
    return (static_cast<QInstructionTypeNumber>(gate) & ControlledGateMask) != 0;
}

inline bool isTwoQubitGate(QInstructionType gate) {
    return (static_cast<QInstructionTypeNumber>(gate) & TwoQubitGateMask) != 0;
}

inline QInstructionType baseGateType(QInstructionType gate) {
    return static_cast<QInstructionType>(static_cast<QInstructionTypeNumber>(gate) & ~ControlledGateMask);
}

inline QInstructionType controlledGateType(QInstructionType gate) {
    return static_cast<QInstructionType>(static_cast<QInstructionTypeNumber>(gate) | ControlledGateMask);
}

inline bool isParameterizedGate(QInstructionType gate) {
    return (static_cast<QInstructionTypeNumber>(gate) & ParameterizedGateMask) != 0;
}

inline string instruction_label(QInstructionType type, bool include_control_info = false) {
    switch (include_control_info? type : baseGateType(type)) {
        case QInstructionType::I:    return "I";
        case QInstructionType::X:    return "X";
        case QInstructionType::Y:    return "Y";
        case QInstructionType::Z:    return "Z";
        case QInstructionType::H:    return "H";
        case QInstructionType::S:    return "S";
        case QInstructionType::Sdg:  return "Sdg";
        case QInstructionType::T:    return "T";
        case QInstructionType::Tdg:  return "Tdg";
        case QInstructionType::SX:   return "SX";
        case QInstructionType::SXdg: return "SXdg";
        case QInstructionType::P:    return "P";
        case QInstructionType::Rz:   return "Rz";
        case QInstructionType::Rx:   return "Rx";
        case QInstructionType::Ry:   return "Ry";
        case QInstructionType::SWAP: return "SWAP";
        case QInstructionType::Rzz:  return "Rzz";
        case QInstructionType::Rxx:  return "Rxx";
        case QInstructionType::Ryy:  return "Ryy";

        case QInstructionType::CX:    return "CX";
        case QInstructionType::CY:    return "CY";
        case QInstructionType::CZ:    return "CZ";
        case QInstructionType::CH:    return "CH";
        case QInstructionType::CS:    return "CS";
        case QInstructionType::CSdg:  return "CSdg";
        case QInstructionType::CT:    return "CT";
        case QInstructionType::CTdg:  return "CTdg";
        case QInstructionType::CSX:   return "CSX";
        case QInstructionType::CSXdg: return "CSXdg";
        case QInstructionType::CP:    return "CP";
        case QInstructionType::CRz:   return "CRz";
        case QInstructionType::CRx:   return "CRx";
        case QInstructionType::CRy:   return "CRy";
        case QInstructionType::CSWAP: return "CSWAP";
        case QInstructionType::CRzz:  return "CRzz";
        case QInstructionType::CRxx:  return "CRxx";
        case QInstructionType::CRyy:  return "CRyy";
    }
    throw std::invalid_argument("instruction_label: Unknown QInstructionType");
};




// Advanced declarations
class MQGTSim;



/**
 * @class QInstruction
 * @brief Represents a single quantum instruction (gate) in a quantum circuit.
 *
 * This class stores the gate type, target qubits, control qubits, and parameters.
 * It provides methods for inversion, permutation, and controlled gate generation.
 *
 * @tparam max_qubits Maximum number of qubits this instruction can target or control.
 */
template<unsigned int max_qubits>
class QInstruction {

public:
    class ControlsView;

private: // ---------- REPRESENTATION ----------
    BitMask<max_qubits> control_mask; // Control mask in case of controlled gates
    BitMask<max_qubits> control_values;  // Control value mask in case of controlled gates
    Real parameter; // Parameter in [0, 2pi] in case of parameterized gates
    QubitIndex t1, t2; // Target qubit indices (t1 for single-qubit gates, t1 and t2 for two-qubit gates)
    QInstructionType type; // Instruction type
    
    

private: // ---------- PRIVATE METHODS -------------
    QInstruction(QInstructionType t, const BitMask<max_qubits>& c, const BitMask<max_qubits>& cv, 
                    QubitIndex target1, QubitIndex target2, Real p)
        : control_mask(c), control_values(cv), parameter(p), t1(target1), t2(target2), type(t) {

        if (isTwoQubitGate()) {

            if (t1 == t2) {
                throw std::invalid_argument("QInstruction::QInstruction: Two-qubit gate cannot have the same target qubit.");
            }
            if (t1 >= max_qubits || t2 >= max_qubits) {
                throw std::invalid_argument("QInstruction::QInstruction: Target qubit index is out of bounds.");
            }
            if (t1 > t2) {
                std::swap(t1, t2); // Ensure t1 < t2 for consistency
            }

            // Check controls & targets do not overlap
            if (control_mask.isSet(t1) || control_mask.isSet(t2)) {
                throw std::invalid_argument("QInstruction::QInstruction: Control qubits cannot overlap with target qubits.");
            }
        } else {

            if (t1 >= max_qubits) {
                throw std::invalid_argument("QInstruction::QInstruction: Target qubit index is out of bounds.");
            }

            // Check controls & targets do not overlap
            if (control_mask.isSet(t1)) {
                throw std::invalid_argument("QInstruction::QInstruction: Control qubits cannot overlap with target qubits.");
            }
        }

        // Check controls
        if (!control_mask.any() && isControlledGate(type)) {
            throw std::invalid_argument("QInstruction::QInstruction: Controlled gate must have at least one control qubit.");
        } else if (control_mask.any() && !isControlledGate(type)) {
            throw std::invalid_argument("QInstruction::QInstruction: Non-controlled gate cannot have control qubits.");
        }
    }

    // Creates a controlled gate instruction
    static QInstruction<max_qubits> make_controlled(QInstructionType t, const Real &p, 
                                                    QubitIndex target1, QubitIndex target2, 
                                                    BitMask<max_qubits> ctrl_qubits, 
                                                    BitMask<max_qubits> ctrl_vals) {

        // Create the instruction
        return QInstruction<max_qubits>(controlledGateType(t), ctrl_qubits, ctrl_vals, target1, target2, p);
    }


public: // ---------- CONSTRUCTORS ----------
    /**
     * @brief Default constructor creating an Identity (I) instruction.
     */
    QInstruction() : control_mask(), control_values(), parameter(0.0), t1(0), t2(0), type(QInstructionType::I) {}

    /**
     * @brief Default destructor.
     */
    ~QInstruction() = default;

    // Copy & Move
    QInstruction(const QInstruction<max_qubits>& other) = default;
    QInstruction(QInstruction<max_qubits>&& other) noexcept = default;
    QInstruction& operator=(const QInstruction<max_qubits>& other) = default;
    QInstruction& operator=(QInstruction<max_qubits>&& other) noexcept = default;

    // ----- FACTORY METHODS -----
    static QInstruction<max_qubits> X(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::X, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> Y(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::Y, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> Z(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::Z, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> H(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::H, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> S(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::S, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> Sdg(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::Sdg, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> T(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::T, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> Tdg(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::Tdg, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> SX(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::SX, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }
    static QInstruction<max_qubits> SXdg(QubitIndex target) { return QInstruction<max_qubits>(QInstructionType::SXdg, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, 0); }

    static QInstruction<max_qubits> P(QubitIndex target, Real angle) { return QInstruction<max_qubits>(QInstructionType::P, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, angle); }
    static QInstruction<max_qubits> Rz(QubitIndex target, Real angle) { return QInstruction<max_qubits>(QInstructionType::Rz, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, angle); }
    static QInstruction<max_qubits> Rx(QubitIndex target, Real angle) { return QInstruction<max_qubits>(QInstructionType::Rx, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, angle); }
    static QInstruction<max_qubits> Ry(QubitIndex target, Real angle) { return QInstruction<max_qubits>(QInstructionType::Ry, BitMask<max_qubits>(), BitMask<max_qubits>(), target, 0, angle); }

    static QInstruction<max_qubits> SWAP(QubitIndex target1, QubitIndex target2) { return QInstruction<max_qubits>(QInstructionType::SWAP, BitMask<max_qubits>(), BitMask<max_qubits>(), target1, target2, 0); }
    static QInstruction<max_qubits> Rzz(QubitIndex target1, QubitIndex target2, Real angle) { return QInstruction<max_qubits>(QInstructionType::Rzz, BitMask<max_qubits>(), BitMask<max_qubits>(), target1, target2, angle); }
    static QInstruction<max_qubits> Rxx(QubitIndex target1, QubitIndex target2, Real angle) { return QInstruction<max_qubits>(QInstructionType::Rxx, BitMask<max_qubits>(), BitMask<max_qubits>(), target1, target2, angle); }
    static QInstruction<max_qubits> Ryy(QubitIndex target1, QubitIndex target2, Real angle) { return QInstruction<max_qubits>(QInstructionType::Ryy, BitMask<max_qubits>(), BitMask<max_qubits>(), target1, target2, angle); }

    // Controlled factory methods
    static QInstruction<max_qubits> CX(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CX, 0, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CY(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CY, 0, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CZ(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CZ, 0, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CH(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CH, 0, target, 0, ctrl, ctrl_values);}
    static QInstruction<max_qubits> CS(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CS, 0, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CSdg(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CSdg, 0, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CT(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CT, 0, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CTdg(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CTdg, 0, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CSX(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CSX, 0, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CSXdg(QubitIndex target, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CSXdg, 0, target, 0, ctrl, ctrl_values); }
    
    static QInstruction<max_qubits> CP(QubitIndex target, Real angle, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CP, angle, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CRz(QubitIndex target, Real angle, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CRz, angle, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CRx(QubitIndex target, Real angle, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CRx, angle, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CRy(QubitIndex target, Real angle, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CRy, angle, target, 0, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CSWAP(QubitIndex target1, QubitIndex target2, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CSWAP, Real(), target1, target2, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CRzz(QubitIndex target1, QubitIndex target2, Real angle, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CRzz, angle, target1, target2, ctrl, ctrl_values); }  
    static QInstruction<max_qubits> CRxx(QubitIndex target1, QubitIndex target2, Real angle, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CRxx, angle, target1, target2, ctrl, ctrl_values); }
    static QInstruction<max_qubits> CRyy(QubitIndex target1, QubitIndex target2, Real angle, const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) { return make_controlled(QInstructionType::CRyy, angle, target1, target2, ctrl, ctrl_values); }


public: // --------- INTERFACE ----------

    /**
     * @brief Creates a new instruction with permuted qubit indices.
     * @param permutation The mapping array for the permutation.
     * @return A new permuted QInstruction.
     */
    QInstruction<max_qubits> permute(span<const QubitIndex> permutation) const {
        auto map_qubit = [&](QubitIndex q) -> QubitIndex {
            if (q < permutation.size() && permutation[q] < max_qubits) {
                return permutation[q];
            }
            throw std::invalid_argument("QInstruction::permute: Invalid permutation mapping for qubit index " + std::to_string((int)q));
        };

        const QubitIndex mapped_t1 = static_cast<QubitIndex>(map_qubit(t1));
        const QubitIndex mapped_t2 = isTwoQubitGate() ? static_cast<QubitIndex>(map_qubit(t2)) : (QubitIndex)0;

        
        BitMask<max_qubits> new_control_mask = control_mask.permute(permutation);
        BitMask<max_qubits> new_control_values = control_values.permute(permutation);
        
        // Error checking
        if (new_control_mask.count() != control_mask.count()) {
            throw std::logic_error("QInstruction::permute: Overlapping qubits detected in control mask after permutation.");
        }

        return QInstruction<max_qubits>(type, new_control_mask, new_control_values,
                                        mapped_t1, mapped_t2, parameter);
    }

    /**
     * @brief Permutes the qubit indices of this instruction in-place.
     * @param permutation The mapping array for the permutation.
     * @return Reference to the permuted instruction.
     */
    QInstruction<max_qubits> & permute_inplace(span<const QubitIndex> permutation) {
        *this = permute(permutation);
        return *this;
    }

    /**
     * @brief Gets a bitmask of all qubits (targets and controls) used by this instruction.
     */
    BitMask<max_qubits> getQubitMask() const {
        BitMask<max_qubits> mask= control_mask;
        mask.set(t1);
        if (isTwoQubitGate()) {
            mask.set(t2);
        }

        if (mask.count() != control_mask.count() + (isTwoQubitGate() ? 2 : 1)) {
            throw std::logic_error("QInstruction::getQubitMask: Overlapping qubits detected in control and target masks.");
        }

        return mask;
    }

    QInstructionType getType() const noexcept { return type; }
    bool isControlled() const noexcept { return isControlledGate(type); }
    const BitMask<max_qubits>& getControlMask() const noexcept { return control_mask; }
    const BitMask<max_qubits>& getControlValueMask() const noexcept { return control_values; }
    pair<vector<QubitIndex>, string> getControlAndValues() const {

        const QubitIndex num_controls= control_mask.count();
        vector<QubitIndex> indices;
        indices.reserve(num_controls);
        string values(num_controls, '0');
        QubitIndex i = 0;
        for (const auto bit : control_mask) {
            indices.push_back(bit);
            values[i++] = control_values.isSet(bit) ? '1' : '0';
        }
        return {indices, values};
    }

    ControlsView controls() const noexcept { return ControlsView(control_mask, control_values); }

    bool isParametrized() const noexcept {
        return isParameterizedGate(type);
    }

    Real getParameter() const noexcept { return parameter; }
    

    bool isTwoQubitGate() const noexcept {
        return ::isTwoQubitGate(type);
    }

    QubitIndex getTarget() const { 
        if (isTwoQubitGate()) {
            throw std::invalid_argument("QInstruction::getTarget: Instruction does not have a single target"); 
        }
        return t1; 
    }

    array<QubitIndex,2> getTargets() const { 
        return {t1, t2};
    }

    
    pair<QubitIndex, QubitIndex> getTargetPair() const { 
        if (!isTwoQubitGate()) {
            throw std::invalid_argument("QInstruction::getTargetPair: Instruction does not have two targets"); 
        }
        return {t1, t2}; 
    }


    

    /**
     * @brief Computes the inverse of this instruction.
     * @return A new QInstruction representing the inverse operation.
     */
    QInstruction<max_qubits> inverse() const {
        switch (type) {
            case QInstructionType::S:    return QInstruction<max_qubits>(QInstructionType::Sdg, control_mask, control_values, t1, t2, parameter);
            case QInstructionType::Sdg:  return QInstruction<max_qubits>(QInstructionType::S,   control_mask, control_values, t1, t2, parameter);
            case QInstructionType::SX:   return QInstruction<max_qubits>(QInstructionType::SXdg,control_mask, control_values, t1, t2, parameter);
            case QInstructionType::SXdg: return QInstruction<max_qubits>(QInstructionType::SX,  control_mask, control_values, t1, t2, parameter);
            case QInstructionType::T:    return QInstruction<max_qubits>(QInstructionType::Tdg, control_mask, control_values, t1, t2, parameter);
            case QInstructionType::Tdg:  return QInstruction<max_qubits>(QInstructionType::T,   control_mask, control_values, t1, t2, parameter);
            case QInstructionType::P:    return QInstruction<max_qubits>(QInstructionType::P,   control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::Rz:   return QInstruction<max_qubits>(QInstructionType::Rz,  control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::Rx:   return QInstruction<max_qubits>(QInstructionType::Rx,  control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::Ry:   return QInstruction<max_qubits>(QInstructionType::Ry,  control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::Rzz:  return QInstruction<max_qubits>(QInstructionType::Rzz, control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::Rxx:  return QInstruction<max_qubits>(QInstructionType::Rxx, control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::Ryy:  return QInstruction<max_qubits>(QInstructionType::Ryy, control_mask, control_values, t1, t2, -parameter);

            case QInstructionType::CS:    return QInstruction<max_qubits>(QInstructionType::CSdg, control_mask, control_values, t1, t2, parameter);
            case QInstructionType::CSdg:  return QInstruction<max_qubits>(QInstructionType::CS,   control_mask, control_values, t1, t2, parameter);
            case QInstructionType::CSX:   return QInstruction<max_qubits>(QInstructionType::CSXdg,control_mask, control_values, t1, t2, parameter);
            case QInstructionType::CSXdg: return QInstruction<max_qubits>(QInstructionType::CSX,  control_mask, control_values, t1, t2, parameter);
            case QInstructionType::CT:    return QInstruction<max_qubits>(QInstructionType::CTdg, control_mask, control_values, t1, t2, parameter);
            case QInstructionType::CTdg:  return QInstruction<max_qubits>(QInstructionType::CT,   control_mask, control_values, t1, t2, parameter);
            case QInstructionType::CP:    return QInstruction<max_qubits>(QInstructionType::CP,   control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::CRz:   return QInstruction<max_qubits>(QInstructionType::CRz,  control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::CRx:   return QInstruction<max_qubits>(QInstructionType::CRx,  control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::CRy:   return QInstruction<max_qubits>(QInstructionType::CRy,  control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::CRzz:  return QInstruction<max_qubits>(QInstructionType::CRzz, control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::CRxx:  return QInstruction<max_qubits>(QInstructionType::CRxx, control_mask, control_values, t1, t2, -parameter);
            case QInstructionType::CRyy:  return QInstruction<max_qubits>(QInstructionType::CRyy, control_mask, control_values, t1, t2, -parameter);

            default:                     return *this; // For gates that are their own inverse (X, Y, Z, H, SWAP)
        }
    }

    
    /**
     * @brief Computes the inverse of this instruction in-place.
     * @return Reference to the inverted instruction.
     */
    QInstruction<max_qubits> & inverse_inplace() {
        switch (type) {
            case QInstructionType::S:    type= QInstructionType::Sdg; break;
            case QInstructionType::Sdg:  type= QInstructionType::S; break;
            case QInstructionType::SX:   type= QInstructionType::SXdg; break;
            case QInstructionType::SXdg: type= QInstructionType::SX; break;
            case QInstructionType::T:    type= QInstructionType::Tdg; break;
            case QInstructionType::Tdg:  type= QInstructionType::T; break;
            case QInstructionType::P:    parameter = -parameter; break;
            case QInstructionType::Rz:   parameter = -parameter; break;
            case QInstructionType::Rx:   parameter = -parameter; break;
            case QInstructionType::Ry:   parameter = -parameter; break;
            case QInstructionType::Rzz:  parameter = -parameter; break;
            case QInstructionType::Rxx:  parameter = -parameter; break;
            case QInstructionType::Ryy:  parameter = -parameter; break;

            case QInstructionType::CS:    type= QInstructionType::CSdg; break;
            case QInstructionType::CSdg:  type= QInstructionType::CS;   break;
            case QInstructionType::CSX:   type= QInstructionType::CSXdg; break;
            case QInstructionType::CSXdg: type= QInstructionType::CSX;  break;
            case QInstructionType::CT:    type= QInstructionType::CTdg; break;
            case QInstructionType::CTdg:  type= QInstructionType::CT;   break;
            case QInstructionType::CP:    parameter= -parameter;   break;
            case QInstructionType::CRz:   parameter = -parameter;  break;
            case QInstructionType::CRx:   parameter = -parameter;  break;
            case QInstructionType::CRy:   parameter = -parameter;  break;
            case QInstructionType::CRzz:  parameter = -parameter; break;
            case QInstructionType::CRxx:  parameter = -parameter; break;
            case QInstructionType::CRyy:  parameter = -parameter; break;
            default:                    break; // For gates that are their own inverse (X, Y, Z, H, SWAP)
        }
        return *this;
    }


    /**
     * @brief Gets the minimum qubit index used by this instruction (targets and controls).
     */
    QubitIndex getMinimumQubitIndex() const {
        QubitIndex min_index = t1;
        if (isControlled()) {
            min_index = std::min(min_index, control_mask.min_bit_set()); 
        }
        return min_index;
    }

    /**
     * @brief Gets the maximum qubit index used by this instruction (targets and controls).
     */
    QubitIndex getMaximumQubitIndex() const {
        QubitIndex max_index = isTwoQubitGate() ? t2 : t1;
        if (isControlled()) {
            max_index = std::max(max_index, control_mask.max_bit_set());
        }
        return max_index;
    }

    // Defaulted equality operator natively implements field-by-field == and !=
    bool operator==(const QInstruction<max_qubits>&) const = default;

    // Custom strict weak ordering based on minimum qubit index
    bool operator<(const QInstruction<max_qubits>& other) const {
        return getMinimumQubitIndex() < other.getMinimumQubitIndex();
    }

    template<unsigned int other_max_qubits>
    QInstruction<other_max_qubits> cast() const {
        if (getMaximumQubitIndex() >= other_max_qubits) {
            throw std::invalid_argument("QInstruction::cast: Cannot cast QInstruction: target or control qubit index exceeds new size.");
        }

        return QInstruction<other_max_qubits>(type, control_mask.template cast<other_max_qubits>(),
                                                control_values.template cast<other_max_qubits>(),
                                                t1, t2, parameter);
    }

    QInstruction<max_qubits> & controlled_inplace(const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) {
        for (QubitIndex qb : ctrl) {
            if (qb >= max_qubits) {
                throw std::invalid_argument("QInstruction::controlled_inplace: Control qubit index is out of bounds.");
            }
            if (qb == t1 || (isTwoQubitGate() && qb == t2)) {
                throw std::invalid_argument("QInstruction::controlled_inplace: Control qubit index cannot be the same as target qubit index.");
            }
            if (control_mask.isSet(qb)) {
                throw std::invalid_argument("QInstruction::controlled_inplace: Duplicate control qubit index.");
            }

            control_mask.set(qb);
            if (ctrl_values.isSet(qb)) {
                control_values.set(qb);
            }
        }
        type= controlledGateType(type);
        return *this;
    }

    QInstruction<max_qubits> controlled(const BitMask<max_qubits> &ctrl, const BitMask<max_qubits> &ctrl_values) const {
        QInstruction<max_qubits> new_instr = *this;
        new_instr.controlled_inplace(ctrl, ctrl_values);
        return new_instr;
    }


public: // ------ I/O ------
    string toString() const {
        string result;
        string name= instruction_label(type, true);

        if (isParametrized()) {
            result = std::format("{}({})", name, parameter.toString());
        } else {
            result = string(name);
        }

        if (isTwoQubitGate()) {
            std::format_to(std::back_inserter(result), "[{},{}]", t1, t2);
        } else {
            std::format_to(std::back_inserter(result), "[{}]", t1);
        }
        
        if (isControlled()) {
            result += " {Controls: ";
            bool first = true;
            for (const auto [bit, value] : controls()) {
                if (!first) result += ",";
                std::format_to(std::back_inserter(result), "{}/{}", bit, value ? 1 : 0);
                first = false;
            }
            result += "}"; // Fixed bug: Was missing a closing brace
        }

        return result;
    }

    friend std::ostream & operator<<(std::ostream &os, const QInstruction<max_qubits> &instr) { 
        return (os << instr.toString()); 
    }


public: // ---------- ITERATOR OVER CONTROL QUBITS ----------
    class ControlIterator {
    private:
        typename BitMask<max_qubits>::SetBitIterator bit_iterator;
        const BitMask<max_qubits>* control_values;

    public:
        using value_type = pair<QubitIndex, bool>;

        ControlIterator(typename BitMask<max_qubits>::SetBitIterator iterator, const BitMask<max_qubits>& values) noexcept
            : bit_iterator(iterator), control_values(&values) {}

        value_type operator*() const {
            const QubitIndex qubit = *bit_iterator;
            return {qubit, control_values->isSet(qubit)};
        }

        ControlIterator& operator++() {
            ++bit_iterator;
            return *this;
        }

        bool operator!=(const ControlIterator& other) const {
            return bit_iterator != other.bit_iterator;
        }
    };

    class ControlsView {
    private:
        const BitMask<max_qubits>* control_mask;
        const BitMask<max_qubits>* control_values;

    public:
        ControlsView(const BitMask<max_qubits>& controls, const BitMask<max_qubits>& values) noexcept
            : control_mask(&controls), control_values(&values) {}

        ControlIterator begin() const noexcept {
            return ControlIterator(control_mask->begin(), *control_values);
        }

        ControlIterator end() const noexcept {
            return ControlIterator(control_mask->end(), *control_values);
        }
    };


private: // ------- FRIENDS ------- 

    friend class MQGTSim;

    template<unsigned int>
    friend class QInstruction;
};



#endif  // QINSTR_H