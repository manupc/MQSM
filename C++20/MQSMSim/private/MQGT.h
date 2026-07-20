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
#ifndef __MQGT_H__
#define __MQGT_H__

#include <bit>
#include <private/FSMDataTypes.h>
#include <private/QInstruction.h>

/**
 * @class MQGT
 * @brief Represents a Mealy Quantum Gate Transducer (MQGT).
 *
 * An MQGT is a specialized Finite State Machine (FSM) that models the unitary
 * evolution of a quantum gate or circuit across multiple qubits. It serves as
 * the core algebraic structure for multiplying and tensor-producting quantum
 * gates.
 *
 * @tparam max_qubits Maximum number of qubits this MQGT can operate on.
 */
template <unsigned int max_qubits>
  requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
class MQGT {

private:                  // ---------- REPRESENTATION ----------
  FSM<max_qubits, 4> fsm; // FSM representation of the MQGT
  MQGT(FSM<max_qubits, 4> &&fsm_) noexcept : fsm(std::move(fsm_)) {}

private: // ---------- HELPERS ----------
  static void add_1qb(QInstructionType type, FSM<max_qubits, 4> &fsm,
                      const Real &theta = 0.0) {

    if (fsm.getNumLevels() < 1) {
      throw std::invalid_argument(
          "MQGT:add_1qb: FSM cannot accommodate a single qubit gate qubits.");
    }
    const QubitIndex start_level = fsm.getNumLevels() - 1;

    // Reserve state and transitions for the new gate
    auto &states = fsm.level_states[start_level];
    states.reserve(states.size() + 1);
    MachineIndex state_index = states.size();
    states.emplace_back();
    auto &state = states[state_index];
    auto &transitions = fsm.level_transitions[start_level];
    transitions.reserve(
        transitions.size() +
        (isAntidiagonalGate(type) || isDiagonalGate(type) ? 2 : 4));

    // Add sub-machine for the new qubit based on the gate type
    if (isAntidiagonalGate(type)) {

      state.T[1] = transitions.size();
      state.T[2] = transitions.size() + 1;
      if (type == QInstructionType::X) {
        transitions.emplace_back(1, 0);
        transitions.emplace_back(1, 0);
      } else { // Y gate
        transitions.emplace_back(Complex::minus_i(), 0);
        transitions.emplace_back(Complex::i(), 0);
      }
    } else if (isDiagonalGate(type)) {

      state.T[0] = transitions.size();
      state.T[3] = transitions.size() + 1;
      Complex up_weight = (type == QInstructionType::Rz)
                              ? Complex::fromPolar(1, -theta / 2.0)
                              : 1;
      Complex down_weight;
      switch (type) {
      case QInstructionType::I:
        down_weight = 1;
        break;
      case QInstructionType::Z:
        down_weight = -1;
        break;
      case QInstructionType::S:
        down_weight = Complex::i();
        break;
      case QInstructionType::Sdg:
        down_weight = Complex::minus_i();
        break;
      case QInstructionType::T:
        down_weight = Complex::fromPolar(1, Real::pi() / 4.0);
        break;
      case QInstructionType::Tdg:
        down_weight = Complex::fromPolar(1, -Real::pi() / 4.0);
        break;
      case QInstructionType::P:
        down_weight = Complex::fromPolar(1, theta);
        break;
      case QInstructionType::Rz:
        down_weight = Complex::fromPolar(1, theta / 2.0);
        break;
      default:
        throw std::invalid_argument(
            "MQGT:add_1qb: Unsupported diagonal gate type.");
      }

      transitions.emplace_back(up_weight, 0);
      transitions.emplace_back(down_weight, 0);

    } else {
      for (int i = 0; i < 4; i++)
        state.T[i] = transitions.size() + i;

      Complex w00, w01, w10, w11;

      switch (type) {
      case QInstructionType::H: {
        Complex inv_sqrt2 = Complex::inv_sqrt2();
        w00 = inv_sqrt2;
        w01 = inv_sqrt2;
        w10 = inv_sqrt2;
        w11 = -inv_sqrt2;
      } break;
      case QInstructionType::SX: {
        Complex one_plus_i_div_2 = Complex(0.5, 0.5);
        Complex one_minus_i_div_2 = Complex(0.5, -0.5);
        w00 = one_plus_i_div_2;
        w01 = one_minus_i_div_2;
        w10 = one_minus_i_div_2;
        w11 = one_plus_i_div_2;
      } break;

      case QInstructionType::SXdg: {
        Complex one_plus_i_div_2 = Complex(0.5, 0.5);
        Complex one_minus_i_div_2 = Complex(0.5, -0.5);
        w00 = one_minus_i_div_2;
        w01 = one_plus_i_div_2;
        w10 = one_plus_i_div_2;
        w11 = one_minus_i_div_2;
      } break;

      case QInstructionType::Rx: {
        Complex cos_theta_div_2 = Complex(Real::cos(theta / 2.0), 0);
        Complex minus_i_sin_theta_div_2 = Complex(0, -Real::sin(theta / 2.0));
        w00 = cos_theta_div_2;
        w01 = minus_i_sin_theta_div_2;
        w10 = minus_i_sin_theta_div_2;
        w11 = cos_theta_div_2;
      } break;
      case QInstructionType::Ry: {
        Complex cos_theta_div_2 = Complex(Real::cos(theta / 2.0), 0);
        Complex sin_theta_div_2 = Complex(Real::sin(theta / 2.0), 0);
        w00 = cos_theta_div_2;
        w01 = -sin_theta_div_2;
        w10 = sin_theta_div_2;
        w11 = cos_theta_div_2;
        // CHECK: w00= cos_theta_div_2; w01= sin_theta_div_2; w10=
        // -sin_theta_div_2; w11= cos_theta_div_2;
      } break;
      default:
        throw std::invalid_argument("MQGT:add_1qb: Unsupported gate type.");
      }

      transitions.emplace_back(w00, 0);
      transitions.emplace_back(w01, 0);
      transitions.emplace_back(w10, 0);
      transitions.emplace_back(w11, 0);
    }
  }

  static MQGT make_1qb(QInstructionType type, const Real &theta = 0.0) {
    FSM<max_qubits, 4> fsm(1);
    add_1qb(type, fsm, theta);
    return MQGT(std::move(fsm));
  }

  static MQGT make_I(QubitIndex tensor_n) {
    FSM<max_qubits, 4> fsm(tensor_n);

    for (QubitIndex level = 0; level < tensor_n; ++level) {
      auto &states = fsm.level_states[level];
      auto &transitions = fsm.level_transitions[level];
      states.reserve(1);
      transitions.reserve(2);
      auto &state = states.emplace_back();
      state.T[0] = transitions.size();
      state.T[3] = transitions.size() + 1;
      transitions.emplace_back(1, 0);
      transitions.emplace_back(1, 0);
    }
    return MQGT(std::move(fsm));
  }

  static void add_2qb(QInstructionType type, FSM<max_qubits, 4> &fsm,
                      const Real &theta = 0.0) {
    if (fsm.getNumLevels() < 2) {
      throw std::invalid_argument(
          "MQGT:add_2qb: FSM cannot accommodate a two-qubit gate.");
    }

    const QubitIndex start_level = fsm.getNumLevels() - 2;

    // Reserve state and transitions for the new gate
    auto &states0 = fsm.level_states[start_level];
    auto &states1 = fsm.level_states[start_level + 1];
    auto &transitions0 = fsm.level_transitions[start_level];
    auto &transitions1 = fsm.level_transitions[start_level + 1];

    const MachineIndex states_0_start = states0.size();

    const MachineIndex num_transitions_0 =
        (type == QInstructionType::Rzz) ? 2 : 4;
    const MachineIndex transitions_0_start = transitions0.size();

    const MachineIndex num_states_1 = (type == QInstructionType::SWAP) ? 4 : 2;
    const MachineIndex states_1_start = states1.size();

    const MachineIndex num_transitions_1 = 4;
    const MachineIndex transitions_1_start = transitions1.size();

    states0.reserve(states_0_start + 1);
    states1.reserve(states_1_start + num_states_1);
    transitions0.reserve(transitions_0_start + num_transitions_0);
    transitions1.reserve(transitions_1_start + num_transitions_1);

    if (type == QInstructionType::SWAP) {

      // Level 0
      auto &state0 = states0.emplace_back();
      for (MachineIndex idx = 0; idx < 4; ++idx) {
        state0.T[idx] = transitions_0_start + idx;
        transitions0.emplace_back(1, states_1_start + idx);
      }

      // Level 1
      for (MachineIndex idx = 0; idx < 4; ++idx) {
        auto &state1 = states1.emplace_back();
        if (idx == 1)
          state1.T[2] = transitions_1_start + idx;
        else if (idx == 2)
          state1.T[1] = transitions_1_start + idx;
        else
          state1.T[idx] = transitions_1_start + idx;
        transitions1.emplace_back(1, 0);
      }

    } else if (type == QInstructionType::Rzz) {
      Complex m_phase = Complex::fromPolar(1, -theta / 2.0);
      Complex p_phase = Complex::fromPolar(1, theta / 2.0);

      // Level 0
      auto &state0 = states0.emplace_back();
      state0.T[0] = transitions_0_start;
      state0.T[3] = transitions_0_start + 1;
      for (MachineIndex idx = 0; idx < 2; ++idx) {
        transitions0.emplace_back(1, states_1_start + idx);
      }

      // Level 1
      auto &state10 = states1.emplace_back();
      state10.T[0] = transitions_1_start;
      state10.T[3] = transitions_1_start + 1;
      auto &state11 = states1.emplace_back();
      state11.T[0] = transitions_1_start + 2;
      state11.T[3] = transitions_1_start + 3;

      transitions1.emplace_back(m_phase, 0);
      transitions1.emplace_back(p_phase, 0);
      transitions1.emplace_back(p_phase, 0);
      transitions1.emplace_back(m_phase, 0);

    } else {

      Complex cos_v = Complex(Real::cos(theta / 2.0), 0);
      Complex minus_i_sin_v = Complex(0, -Real::sin(theta / 2.0));
      Complex i_sin_v = Complex(0, Real::sin(theta / 2.0));

      // Level 0
      auto &state0 = states0.emplace_back();
      state0.T[0] = transitions_0_start;
      state0.T[1] = transitions_0_start + 1;
      state0.T[2] = transitions_0_start + 2;
      state0.T[3] = transitions_0_start + 3;
      transitions0.emplace_back(1, states_1_start);
      transitions0.emplace_back(1, states_1_start + 1);
      transitions0.emplace_back((type == QInstructionType::Rxx ? 1 : -1),
                                states_1_start + 1);
      transitions0.emplace_back(1, states_1_start);

      // Level 1
      auto &state10 = states1.emplace_back();
      state10.T[0] = transitions_1_start;
      state10.T[3] = transitions_1_start + 1;
      auto &state11 = states1.emplace_back();
      state11.T[1] = transitions_1_start + 2;
      state11.T[2] = transitions_1_start + 3;

      transitions1.emplace_back(cos_v, 0);
      transitions1.emplace_back(cos_v, 0);
      transitions1.emplace_back(
          (type == QInstructionType::Rxx ? minus_i_sin_v : i_sin_v), 0);
      transitions1.emplace_back(minus_i_sin_v, 0);
    }
  }

  static MQGT make_2qb(QInstructionType type, const Real &theta = 0.0) {
    FSM<max_qubits, 4> fsm(2);
    add_2qb(type, fsm, theta);
    return MQGT(std::move(fsm));
  }

  static void add_control_skeleton(FSM<max_qubits, 4> &fsm,
                                   QubitIndex num_controls,
                                   string_view ctrl_states) {

    if (num_controls == 0) {
      throw std::invalid_argument("MQGT:add_control_skeleton: Number of "
                                  "controls must be greater than zero.");
    }
    if (ctrl_states.size() != 0 && ctrl_states.size() != num_controls) {
      throw std::invalid_argument(
          "MQGT:add_control_skeleton: Control states string length must match "
          "number of controls.");
    }

    // First control level
    auto &state0 = fsm.level_states[0].emplace_back();
    QubitIndex control_qubit =
        ctrl_states.empty() || ctrl_states[0] == '1' ? 3 : 0;
    state0.T[3 - control_qubit] = 1;
    state0.T[control_qubit] = 0;
    auto &transition0 = fsm.level_transitions[0];
    transition0.reserve(2);
    transition0.emplace_back(1, 0);
    transition0.emplace_back(1, 1);

    // Fill in the remaining control qubits' levels
    for (QubitIndex i = 1; i < num_controls; ++i) {
      auto &states = fsm.level_states[i];
      auto &transitions = fsm.level_transitions[i];

      // Add states
      states.reserve(2);
      auto &stateok = states.emplace_back();
      auto &statefail = states.emplace_back();
      control_qubit = ctrl_states.empty() || ctrl_states[i] == '1' ? 3 : 0;
      stateok.T[3 - control_qubit] = 1;
      stateok.T[control_qubit] = 0;
      statefail.T[0] = 1;
      statefail.T[3] = 1;

      // Add transitions
      transitions.reserve(4);
      const MachineIndex next_state_index =
          (i == num_controls - 1) ? fsm.level_states[i + 1].size() : 1;
      transitions.emplace_back(1, 0);                // Transition for stateok
      transitions.emplace_back(1, next_state_index); // Transition for stateok
      transitions.emplace_back(1, next_state_index); // Transition for statefail
      transitions.emplace_back(1, next_state_index); // Transition for statefail
    }

    // Fill in the target levels
    for (QubitIndex i = num_controls; i < fsm.getNumLevels(); ++i) {
      auto &states = fsm.level_states[i];
      auto &transitions = fsm.level_transitions[i];

      // Add I state
      const MachineIndex transitions_start = transitions.size();
      states.reserve(states.size() + 1);
      auto &state = states.emplace_back();
      state.T[0] = transitions_start;
      state.T[3] = transitions_start + 1;

      transitions.reserve(transitions.size() + 2);
      const MachineIndex next_state_index =
          (i == fsm.getNumLevels() - 1) ? 0 : fsm.level_states[i + 1].size();
      transitions.emplace_back(1, next_state_index); // Transition for I
      transitions.emplace_back(1, next_state_index); // Transition for I
    }
  }

  static MQGT make_controlled_gate(QInstructionType type,
                                   QubitIndex num_controls,
                                   const Real &theta = 0.0,
                                   string_view ctrl_states = "") {

    type = baseGateType(type); // Ensure we are working with the base gate type
    const QubitIndex num_targets = isTwoQubitGate(type) ? 2 : 1;

    if (num_controls + num_targets > max_qubits) {
      throw std::invalid_argument("MQGT:make_controlled_gate: Total number of "
                                  "qubits exceeds maximum allowed.");
    }

    FSM<max_qubits, 4> fsm(num_controls + num_targets);

    // Append the gate to the FSM
    if (num_targets == 1)
      add_1qb(type, fsm, theta);
    else
      add_2qb(type, fsm, theta);

    // Append the control skeleton
    add_control_skeleton(fsm, num_controls, ctrl_states);

    return MQGT(std::move(fsm));
  }

  static MQGT make_controlled_MQGT(const MQGT &mqgt, QubitIndex num_controls,
                                   string_view ctrl_states = "") {

    const QubitIndex num_targets = mqgt.getNumQubits();

    FSM<max_qubits, 4> fsm(num_controls + num_targets);

    // Append the MQGT to the FSM
    const QubitIndex start_level = fsm.getNumLevels() - num_targets;
    std::copy_n(mqgt.fsm.level_states.begin(), num_targets,
                fsm.level_states.begin() + start_level);
    std::copy_n(mqgt.fsm.level_transitions.begin(), num_targets,
                fsm.level_transitions.begin() + start_level);

    // Append the control skeleton
    add_control_skeleton(fsm, num_controls, ctrl_states);

    return MQGT(std::move(fsm)).reduce();
  }

  template <typename T>
    requires(std::is_same_v<T, std::complex<double>> ||
             std::is_same_v<T, Complex>)
  vector<T> toMatrixImpl() const {

    QubitIndex num_qubits = getNumQubits();
    size_t rows = 1ULL << num_qubits;
    vector<T> U(rows * rows, 0.0);

    for (size_t r = 0; r < rows; ++r) {
      for (size_t c = 0; c < rows; ++c) {
        size_t i = r * rows + c;
        MachineIndex current_node = 0;
        T amplitude = 1;
        for (QubitIndex current_qubit = 0; current_qubit < num_qubits;
             ++current_qubit) {

          QubitIndex shift = num_qubits - 1 - current_qubit;
          QubitIndex y = (c >> shift) & 1;
          QubitIndex x = (r >> shift) & 1;
          MachineInputSymbol qubit_value = (x << 1) | y;

          const auto &state = fsm.level_states[current_qubit][current_node];
          if (!state.hasTransition(qubit_value)) {
            amplitude = 0.0;
            break;
          }
          const MachineIndex transition_index =
              state.getTransitionIndex(qubit_value);
          const auto &transition =
              fsm.level_transitions[current_qubit][transition_index];
          amplitude *= transition.weight.toStdComplex();
          current_node = transition.next_state;
        }
        U[i] = amplitude;
      }
    }
    return U;
  }

public: // ---------- CONSTRUCTORS ----------
  /**
   * @brief Default constructor creating an empty MQGT.
   */
  MQGT() = default;

  /**
   * @brief Default destructor.
   */
  ~MQGT() = default;

  MQGT(const MQGT &other) = default;
  MQGT(MQGT &&other) noexcept = default;
  MQGT &operator=(const MQGT &other) = default;
  MQGT &operator=(MQGT &&other) noexcept = default;

public: // ---------- FACTORY BUILDERS -----------
  /**
   * @name Standard Gate Builders
   * @brief Factory methods for creating standard 1-qubit and 2-qubit MQGT
   * gates.
   */
  ///@{
  static MQGT I(QubitIndex tensor_n = 1) { return make_I(tensor_n); }
  static MQGT X() { return make_1qb(QInstructionType::X); }
  static MQGT Y() { return make_1qb(QInstructionType::Y); }
  static MQGT Z() { return make_1qb(QInstructionType::Z); }
  static MQGT H() { return make_1qb(QInstructionType::H); }
  static MQGT S() { return make_1qb(QInstructionType::S); }
  static MQGT Sdg() { return make_1qb(QInstructionType::Sdg); }
  static MQGT T() { return make_1qb(QInstructionType::T); }
  static MQGT Tdg() { return make_1qb(QInstructionType::Tdg); }
  static MQGT SX() { return make_1qb(QInstructionType::SX); }
  static MQGT SXdg() { return make_1qb(QInstructionType::SXdg); }
  static MQGT P(const Real &theta) {
    return make_1qb(QInstructionType::P, theta);
  }
  static MQGT RZ(const Real &theta) {
    return make_1qb(QInstructionType::Rz, theta);
  }
  static MQGT RY(const Real &theta) {
    return make_1qb(QInstructionType::Ry, theta);
  }
  static MQGT RX(const Real &theta) {
    return make_1qb(QInstructionType::Rx, theta);
  }

  static MQGT SWAP() { return make_2qb(QInstructionType::SWAP); }
  static MQGT RZZ(const Real &theta) {
    return make_2qb(QInstructionType::Rzz, theta);
  }
  static MQGT RXX(const Real &theta) {
    return make_2qb(QInstructionType::Rxx, theta);
  }
  static MQGT RYY(const Real &theta) {
    return make_2qb(QInstructionType::Ryy, theta);
  }
  ///@}

  /**
   * @name Controlled Gate Builders
   * @brief Factory methods for creating controlled MQGT gates.
   */
  ///@{
  static MQGT CX(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CX, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CY(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CY, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CZ(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CZ, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CH(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CH, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CS(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CS, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CSdg(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CSdg, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CT(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CT, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CTdg(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CTdg, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CSX(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CSX, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CSXdg(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CSXdg, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CP(const Real &theta, QubitIndex num_controls = 1,
                 string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CP, num_controls, theta,
                                ctrl_states);
  }
  static MQGT CRX(const Real &theta, QubitIndex num_controls = 1,
                  string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CRx, num_controls, theta,
                                ctrl_states);
  }
  static MQGT CRY(const Real &theta, QubitIndex num_controls = 1,
                  string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CRy, num_controls, theta,
                                ctrl_states);
  }
  static MQGT CRZ(const Real &theta, QubitIndex num_controls = 1,
                  string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CRz, num_controls, theta,
                                ctrl_states);
  }

  static MQGT CSWAP(QubitIndex num_controls = 1, string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CSWAP, num_controls, 0.0,
                                ctrl_states);
  }
  static MQGT CRZZ(const Real &theta, QubitIndex num_controls = 1,
                   string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CRzz, num_controls, theta,
                                ctrl_states);
  }
  static MQGT CRXX(const Real &theta, QubitIndex num_controls = 1,
                   string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CRxx, num_controls, theta,
                                ctrl_states);
  }
  static MQGT CRYY(const Real &theta, QubitIndex num_controls = 1,
                   string_view ctrl_states = "") {
    return make_controlled_gate(QInstructionType::CRyy, num_controls, theta,
                                ctrl_states);
  }
  ///@}

public: // ---------- INTERFACE ----------
  /**
   * @brief Gets the number of qubits the MQGT acts on.
   */
  QubitIndex getNumQubits() const { return fsm.getNumLevels(); }

  /**
   * @brief Gets the total number of states in the underlying FSM.
   */
  size_t getNumStates() const { return fsm.getNumStates(); }

  /**
   * @brief Gets the total number of transitions in the underlying FSM.
   */
  size_t getNumTransitions() const { return fsm.getNumTransitions(); }

  /**
   * @brief Checks if the MQGT is empty (contains no levels or states).
   */
  bool empty() const {
    return fsm.getNumLevels() == 0 || fsm.level_states[0].empty();
  }

  /**
   * @brief In-place tensor product with another MQGT.
   */
  MQGT &operator^=(const MQGT &other) {
    fsm ^= other.fsm;
    return *this;
  }

  /**
   * @brief Tensor product with another MQGT.
   */
  MQGT operator^(const MQGT &other) const { return MQGT(fsm ^ other.fsm); }

  /**
   * @brief Reduces the underlying FSM to eliminate redundant states and
   * transitions.
   * @return Reference to the reduced MQGT.
   */
  MQGT &reduce() {
    fsm.reduce();
    return *this;
  }

  /**
   * @brief Computes the MQGT raised to an integer power.
   * @param exponent The power to raise the MQGT to.
   * @return A new MQGT representing the exponentiated gate.
   */
  MQGT power(QubitIndex exponent) const {
    if (exponent == 0)
      return MQGT::I(getNumQubits());
    if (exponent == 1)
      return *this;
    MQGT result = *this;
    for (QubitIndex i = 1; i < exponent; ++i) {
      result.evolve_inplace(*this);
    }
    return result;
  }

  /**
   * @brief Creates a controlled version of this MQGT.
   * @param num_controls The number of control qubits.
   * @param ctrl_states Optional string of '0's and '1's indicating active
   * control states.
   * @return A new controlled MQGT.
   */
  MQGT controlled(QubitIndex num_controls = 1,
                  string_view ctrl_states = "") const {
    return make_controlled_MQGT(*this, num_controls, ctrl_states);
  }

  MQGT &swapLevels(QubitIndex q1, QubitIndex q2) {
    if (q1 >= getNumQubits() || q2 >= getNumQubits()) {
      throw std::invalid_argument("MQGT::swapLevels: Level indices must be "
                                  "within the range of the MQGT.");
    }
    if (q1 == q2)
      return *this;

    QubitIndex low, high;
    if (q1 < q2) {
      low = q1;
      high = q2;
    } else {
      low = q2;
      high = q1;
    }

    // SWAP qubits low and high by evolving with a SWAP gate on those qubits
    evolve_inplace(MQGT::SWAP(), {low, high});
    return *this;
  }

public: // ----------- EVOLUTION INTERFACE -----------
  /**
   * @brief In-place evolution (matrix multiplication) with another MQGT.
   * @param V The other MQGT to evolve with.
   * @param qubits Optional bitmask specifying target qubits.
   * @return Reference to the evolved MQGT.
   */
  MQGT &evolve_inplace(const MQGT &V, BitMask<max_qubits> qubits = {}) {
    if (this == &V) {
      throw std::invalid_argument(
          "MQGT::evolve_inplace: Cannot evolve in place with itself.");
    }

    fsm.compose_inplace(V.fsm, qubits);
    return *this;
  }

  /**
   * @brief Evolution (matrix multiplication) with another MQGT.
   * @param V The other MQGT to evolve with.
   * @param qubits Optional bitmask specifying target qubits.
   * @return A new MQGT representing the evolved result.
   */
  MQGT evolve(const MQGT &V, BitMask<max_qubits> qubits = {}) const {
    MQGT result = *this;
    return result.evolve_inplace(V, qubits);
  }

public: // ---------- DEBUG -----------
  FSM<max_qubits, 4> getFSM() const { return fsm; }

public: // ---------- I/O ----------
  string toString() const { return fsm.toString("u"); }

  vector<Complex> toMatrix() const { return toMatrixImpl<Complex>(); }

  vector<std::complex<double>> toStdMatrix() const {
    return toMatrixImpl<std::complex<double>>();
  }

  string toMatrixString() const {
    vector<Complex> U = toMatrix();
    ostringstream os;
    QubitIndex num_qubits = getNumQubits();
    size_t rows = 1ULL << num_qubits;

    for (size_t r = 0; r < rows; ++r) {
      for (size_t c = 0; c < rows; ++c) {
        size_t i = r * rows + c;
        if (U[i].hasImag() && U[i].hasReal()) {
          os << "(" << U[i] << ") ";
        } else {
          os << U[i] << " ";
        }
      }
      os << "\n";
    }
    return os.str();
  }

  friend std::ostream &operator<<(std::ostream &os, const MQGT &mqgt) {
    os << mqgt.toString();
    return os;
  }
};

#endif