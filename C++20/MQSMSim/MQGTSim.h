/*
 * Copyright (C) 2026 Manuel Pegalajar Cuellar
 * Department of Computer Science and Artificial Intelligence, University of
 * Granada (Spain) E-mail: manupc@ugr.es
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
#ifndef __MQGTSIM_H__
#define __MQGTSIM_H__

#include <MQSMQuantumState.h>
#include <QCirc.h>
#include <UnitaryOp.h>
#include <private/MQGT.h>
#include <private/MQSM.h>
#include <private/Numbers.h>
#include <private/QCircuitLevel.h>

#include <array>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using std::array;
using std::pair;
using std::string;
using std::vector;

/**
 * @class MQGTSim
 * @brief Simulator for Quantum Circuits using Multi-Qubit Gate Transducers
 * (MQGT).
 *
 * This class provides the main interface for evolving quantum states or
 * computing unitary matrices from quantum circuits. It relies on the MQGT
 * representation for efficient simulation.
 */
class MQGTSim {

private: // ---------- REPRESENTATION ----------
  /**
   * @brief Aligns the target and control qubits for a controlled gate.
   *
   * Reorders the MQGT instruction levels so that the local sorted structure
   * matches the physical qubit indices provided in controls and targets.
   *
   * @tparam max_qubits Maximum number of qubits supported.
   * @param mqgt_instr The base MQGT instruction to align.
   * @param controls List of control qubit indices.
   * @param targets Array of target qubit indices.
   * @param is_two_qubit_gate True if the base gate is a two-qubit gate, false
   * otherwise.
   * @return An MQGT instruction correctly aligned to the target and control
   * qubits.
   *
   * @pre max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED
   * @post The returned MQGT has its internal levels permuted to match the
   * requested indices.
   *
   * @par Example:
   * @code
   * // For internal use only:
   * auto aligned_gate = alignControlledGate(base_cx, {0}, {1, 0}, false);
   * @endcode
   */
  template <unsigned int max_qubits>
    requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
  static MQGT<max_qubits> alignControlledGate(
      MQGT<max_qubits> mqgt_instr, const std::vector<QubitIndex> &controls,
      const std::array<QubitIndex, 2> &targets, bool is_two_qubit_gate) {
    std::vector<QubitIndex> qubit_order = controls;
    qubit_order.push_back(targets[0]);
    if (is_two_qubit_gate) {
      qubit_order.push_back(targets[1]);
    }

    std::vector<QubitIndex> sorted_qubits = qubit_order;
    std::sort(sorted_qubits.begin(), sorted_qubits.end());

    std::vector<QubitIndex> local_to_sorted(qubit_order.size(), 0);
    for (size_t i = 0; i < qubit_order.size(); ++i) {
      auto it = std::lower_bound(sorted_qubits.begin(), sorted_qubits.end(),
                                 qubit_order[i]);
      if (it == sorted_qubits.end() || *it != qubit_order[i]) {
        throw std::runtime_error(
            "MQGTSim: Failed to align controlled gate qubits.");
      }
      local_to_sorted[i] =
          static_cast<QubitIndex>(std::distance(sorted_qubits.begin(), it));
    }

    // Apply the full permutation with swaps until each local index reaches its
    // sorted slot.
    for (QubitIndex i = 0; i < local_to_sorted.size(); ++i) {
      while (local_to_sorted[i] != i) {
        const QubitIndex j = local_to_sorted[i];
        mqgt_instr = mqgt_instr.swapLevels(i, j);
        std::swap(local_to_sorted[i], local_to_sorted[j]);
      }
    }

    return mqgt_instr;
  }

  /**
   * @brief Converts a standard QInstruction into an MQGT representation.
   *
   * Parses the instruction type, parameters, and controls to construct the
   * equivalent MQGT state machine that applies the operation.
   *
   * @tparam max_qubits Maximum number of qubits supported.
   * @param instr The quantum instruction to convert.
   * @return The corresponding MQGT instruction.
   *
   * @pre max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED
   * @post A valid MQGT representing the instruction is returned.
   *
   * @par Example:
   * @code
   * // For internal use only:
   * QInstruction<5> hadamard(QInstructionType::H, {0});
   * MQGT<5> h_mqgt = instructionToMQGT(hadamard);
   * @endcode
   */
  template <unsigned int max_qubits>
    requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
  MQGT<max_qubits> instructionToMQGT(const QInstruction<max_qubits> &instr) {

    MQGT<max_qubits> mqgt_instr;
    vector<QubitIndex> control_qubits;
    string control_values;
    QubitIndex num_controls = 0;

    if (instr.isControlled()) {
      std::tie(control_qubits, control_values) = instr.getControlAndValues();
      num_controls = control_qubits.size();
    }

    // Create an MQGT representation of the instruction
    switch (instr.getType()) {
    case QInstructionType::X:
      mqgt_instr = MQGT<max_qubits>::X();
      break;
    case QInstructionType::Y:
      mqgt_instr = MQGT<max_qubits>::Y();
      break;
    case QInstructionType::Z:
      mqgt_instr = MQGT<max_qubits>::Z();
      break;
    case QInstructionType::H:
      mqgt_instr = MQGT<max_qubits>::H();
      break;
    case QInstructionType::S:
      mqgt_instr = MQGT<max_qubits>::S();
      break;
    case QInstructionType::Sdg:
      mqgt_instr = MQGT<max_qubits>::Sdg();
      break;
    case QInstructionType::T:
      mqgt_instr = MQGT<max_qubits>::T();
      break;
    case QInstructionType::Tdg:
      mqgt_instr = MQGT<max_qubits>::Tdg();
      break;
    case QInstructionType::SX:
      mqgt_instr = MQGT<max_qubits>::SX();
      break;
    case QInstructionType::SXdg:
      mqgt_instr = MQGT<max_qubits>::SXdg();
      break;
    case QInstructionType::P:
      mqgt_instr = MQGT<max_qubits>::P(instr.getParameter());
      break;
    case QInstructionType::Rz:
      mqgt_instr = MQGT<max_qubits>::RZ(instr.getParameter());
      break;
    case QInstructionType::Rx:
      mqgt_instr = MQGT<max_qubits>::RX(instr.getParameter());
      break;
    case QInstructionType::Ry:
      mqgt_instr = MQGT<max_qubits>::RY(instr.getParameter());
      break;
    case QInstructionType::SWAP:
      mqgt_instr = MQGT<max_qubits>::SWAP();
      break;
    case QInstructionType::Rzz:
      mqgt_instr = MQGT<max_qubits>::RZZ(instr.getParameter());
      break;
    case QInstructionType::Rxx:
      mqgt_instr = MQGT<max_qubits>::RXX(instr.getParameter());
      break;
    case QInstructionType::Ryy:
      mqgt_instr = MQGT<max_qubits>::RYY(instr.getParameter());
      break;

    case QInstructionType::CX:
      mqgt_instr = MQGT<max_qubits>::CX(num_controls, control_values);
      break;
    case QInstructionType::CY:
      mqgt_instr = MQGT<max_qubits>::CY(num_controls, control_values);
      break;
    case QInstructionType::CZ:
      mqgt_instr = MQGT<max_qubits>::CZ(num_controls, control_values);
      break;
    case QInstructionType::CH:
      mqgt_instr = MQGT<max_qubits>::CH(num_controls, control_values);
      break;
    case QInstructionType::CS:
      mqgt_instr = MQGT<max_qubits>::CS(num_controls, control_values);
      break;
    case QInstructionType::CSdg:
      mqgt_instr = MQGT<max_qubits>::CSdg(num_controls, control_values);
      break;
    case QInstructionType::CT:
      mqgt_instr = MQGT<max_qubits>::CT(num_controls, control_values);
      break;
    case QInstructionType::CTdg:
      mqgt_instr = MQGT<max_qubits>::CTdg(num_controls, control_values);
      break;
    case QInstructionType::CSX:
      mqgt_instr = MQGT<max_qubits>::CSX(num_controls, control_values);
      break;
    case QInstructionType::CSXdg:
      mqgt_instr = MQGT<max_qubits>::CSXdg(num_controls, control_values);
      break;
    case QInstructionType::CP:
      mqgt_instr = MQGT<max_qubits>::CP(instr.getParameter(), num_controls,
                                        control_values);
      break;
    case QInstructionType::CRz:
      mqgt_instr = MQGT<max_qubits>::CRZ(instr.getParameter(), num_controls,
                                         control_values);
      break;
    case QInstructionType::CRx:
      mqgt_instr = MQGT<max_qubits>::CRX(instr.getParameter(), num_controls,
                                         control_values);
      break;
    case QInstructionType::CRy:
      mqgt_instr = MQGT<max_qubits>::CRY(instr.getParameter(), num_controls,
                                         control_values);
      break;
    case QInstructionType::CSWAP:
      mqgt_instr = MQGT<max_qubits>::CSWAP(num_controls, control_values);
      break;
    case QInstructionType::CRzz:
      mqgt_instr = MQGT<max_qubits>::CRZZ(instr.getParameter(), num_controls,
                                          control_values);
      break;
    case QInstructionType::CRxx:
      mqgt_instr = MQGT<max_qubits>::CRXX(instr.getParameter(), num_controls,
                                          control_values);
      break;
    case QInstructionType::CRyy:
      mqgt_instr = MQGT<max_qubits>::CRYY(instr.getParameter(), num_controls,
                                          control_values);
      break;
    default:
      throw std::invalid_argument(
          "MQGTSim::instructionToMQGT: Unsupported instruction type" +
          instruction_label(instr.getType()) + " in MQGTSim.");
    }

    // Create controlled gate if needed
    array<QubitIndex, 2> targets = instr.getTargets();
    if (instr.isControlled()) {
      mqgt_instr = alignControlledGate(std::move(mqgt_instr), control_qubits,
                                       targets, instr.isTwoQubitGate());
    }
    return mqgt_instr;
  }

public: // ---------- CONSTRUCTORS ----------
  /**
   * @brief Default constructor.
   */
  MQGTSim() = default;

  /**
   * @brief Default destructor.
   */
  ~MQGTSim() = default;

  MQGTSim(const MQGTSim &other) = default;
  MQGTSim(MQGTSim &&other) noexcept = default;
  MQGTSim &operator=(const MQGTSim &other) = default;
  MQGTSim &operator=(MQGTSim &&other) noexcept = default;

public: // ---------- INTERFACE ----------
  /**
   * @brief Computes the full unitary matrix of a given quantum circuit.
   *
   * Iterates through all levels and instructions of the quantum circuit and
   * composes them into a single UnitaryOp (MQGT).
   *
   * @tparam max_qubits Maximum number of qubits supported.
   * @param qc The quantum circuit to convert into a unitary matrix.
   * @return A UnitaryOp representing the full evolution of the circuit.
   *
   * @pre max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED
   * @post The returned UnitaryOp accurately represents the combined unitary of
   * the circuit.
   *
   * @par Example:
   * @code
   * MQGTSim sim;
   * QCircuit<3> qc(3);
   * qc.h(0);
   * qc.cx(0, 1);
   * UnitaryOp<3> unitary = sim.getUnitary(qc);
   * @endcode
   */
  template <unsigned int max_qubits>
    requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
  UnitaryOp<max_qubits> getUnitary(const QCircuit<max_qubits> &qc) {

    MQGT<max_qubits> mqgt = MQGT<max_qubits>::I(qc.getNumQubits());
    size_t level_idx = 0;
    for (const auto &level : qc.getLevels()) { // Iterate over circuit levels

      ++level_idx;
      for (const auto &instr :
           level.getInstructions()) { // Iterate over level instructions

        // Apply the instruction to the MQGT
        mqgt.evolve_inplace(instructionToMQGT(instr), instr.getQubitMask());
      }
    }

    return UnitaryOp<max_qubits>(std::move(mqgt));
  }

  /**
   * @brief Evolves an initial quantum state according to the given quantum
   * circuit.
   *
   * Applies all instructions in the circuit level by level to the initial MQSM
   * state. If the initial state is empty, it assumes the |0...0> basis state.
   *
   * @tparam max_qubits Maximum number of qubits supported.
   * @param qc The quantum circuit specifying the evolution.
   * @param initial_state The initial quantum state (defaults to |0...0>).
   * @return The final MQSMQuantumState after applying all instructions.
   *
   * @pre max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED
   * @pre The circuit's active qubits must not exceed max_qubits.
   * @post The returned quantum state reflects the application of the circuit.
   *
   * @par Example:
   * @code
   * MQGTSim sim;
   * QCircuit<5> qc(5);
   * qc.h(0);
   * qc.cx(0, 1);
   * MQSMQuantumState<5> state = sim.evolveMQSM(qc);
   * @endcode
   */
  template <unsigned int max_qubits>
    requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
  MQSMQuantumState<max_qubits>
  evolveMQSM(const QCircuit<max_qubits> &qc,
             const MQSMQuantumState<max_qubits> &initial_state =
                 MQSMQuantumState<max_qubits>()) {

    MQSM<max_qubits> mqsm =
        initial_state.empty()
            ? MQSM<max_qubits>::fromBasis(0, qc.getNumQubits())
            : initial_state.getMQSM();

    for (const auto &level : qc.getLevels()) { // Iterate over circuit levels

      // Get the mqgt for the level
      MQGT<max_qubits> mqgt;
      BitMask<max_qubits> active_qubits_mask;

      for (const auto &instr :
           level.getInstructions()) { // Iterate over level instructions
        if (mqgt.empty()) {
          mqgt = instructionToMQGT(instr);
          active_qubits_mask = instr.getQubitMask();
        } else {
          mqgt ^= instructionToMQGT(instr);
          active_qubits_mask |= instr.getQubitMask();
        }
      }

      // Apply the level MQGT to the MQSM
      mqsm.evolve_inplace(mqgt, active_qubits_mask);
    }

    return MQSMQuantumState(std::move(mqsm));
  }
};

#endif