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
#ifndef __MQSM_H__
#define __MQSM_H__

#include <algorithm>
#include <bit>
#include <private/FSMDataTypes.h>
#include <private/MQGT.h>
#include <random>
#include <string>
#include <utils/SampleCount.h>
#include <utils/SampleProbabilities.h>

using std::string;

template <unsigned int max_qubits>
  requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
class UnitaryOp;

// Advanced declarations
template <unsigned int max_qubits>
  requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
class MQSMQuantumState;

/**
 * @class MQSM
 * @brief Represents a Mealy Quantum State Model (MQSM).
 *
 * An MQSM models a quantum state using a finite state machine structure, where
 * each level corresponds to a qubit. This allows for highly efficient
 * representation and manipulation of quantum states, particularly those with
 * significant structure or entanglement.
 *
 * @tparam max_qubits Maximum number of qubits this MQSM can track.
 */
template <unsigned int max_qubits>
  requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
class MQSM {

private:                  // ---------- REPRESENTATION ----------
  FSM<max_qubits, 2> fsm; // FSM representation of the MQSM

  MQSM(FSM<max_qubits, 2> &&fsm_) noexcept : fsm(std::move(fsm_)) {}

private: // ---------- HELPERS ----------
  template <typename T>
    requires(std::is_same_v<T, std::complex<double>> ||
             std::is_same_v<T, Complex>)
  vector<T> toStatevectorImpl() const {

    QubitIndex num_qubits = getNumQubits();
    size_t sv_size = 1ULL << num_qubits;
    vector<T> sv(sv_size, 0);

    for (size_t i = 0; i < sv_size; ++i) {
      size_t current_node = 0;
      T amplitude = 1;
      for (QubitIndex current_qubit = 0; current_qubit < num_qubits;
           ++current_qubit) {

        QubitIndex shift = num_qubits - 1 - current_qubit;
        QubitIndex qubit_value = (i >> shift) & 1;

        const auto &state = fsm.level_states[current_qubit][current_node];

        if (!state.hasTransition(qubit_value)) {
          amplitude = 0;
          break;
        }
        const auto &transition =
            fsm.level_transitions[current_qubit]
                                 [state.getTransitionIndex(qubit_value)];
        amplitude *= transition.weight.toStdComplex();

        current_node = transition.next_state;
      }
      sv[i] = amplitude;
    }
    return sv;
  }

  template <typename T>
    requires(std::is_same_v<T, std::complex<double>> ||
             std::is_same_v<T, Complex>)
  static MQSM fromStatevectorImpl(const vector<T> &sv) {

    if (sv.empty()) {
      throw std::invalid_argument(
          "MQSM::fromStatevectorImpl: Statevector cannot be empty");
    } else if (!std::has_single_bit(sv.size())) {
      throw std::invalid_argument(
          "MQSM::fromStatevectorImpl: Statevector size must be a power of 2");
    }

    QubitIndex num_qubits = std::countr_zero(sv.size());
    if (num_qubits > max_qubits) {
      throw std::invalid_argument("MQSM::fromStatevectorImpl: Number of qubits "
                                  "exceeds maximum allowed.");
    }
    FSM<max_qubits, 2> fsm(num_qubits);
    fsm.addState(0); // Initial State
    Complex one = 1;

    for (size_t i = 0; i < sv.size(); ++i) {
      if (Complex(sv[i]).NearZero())
        continue;

      MachineIndex current_node = 0;
      for (QubitIndex current_qubit = 0; current_qubit < num_qubits;
           ++current_qubit) {

        auto &states = fsm.level_states[current_qubit];
        auto &transitions = fsm.level_transitions[current_qubit];

        QubitIndex shift = num_qubits - 1 - current_qubit;
        QubitIndex qubit_value = (i >> shift) & 1;

        if (current_qubit == num_qubits - 1) {
          states[current_node].T[qubit_value] = transitions.size();
          transitions.emplace_back(Complex(sv[i]), 0);

        } else {

          if (!fsm.stateHasTransition(current_qubit, current_node,
                                      qubit_value)) {
            MachineIndex new_node_idx = fsm.addState(current_qubit + 1);
            states[current_node].T[qubit_value] = transitions.size();
            transitions.emplace_back(one, new_node_idx);
          }

          current_node =
              transitions[states[current_node].getTransitionIndex(qubit_value)]
                  .next_state;
        }
      }
    }

    // Normalize the FSM and return
    return MQSM(std::move(fsm.reduce()));
  }

public: // ---------- CONSTRUCTORS ----------
  /**
   * @brief Default constructor creating an empty state.
   */
  MQSM() = default;

  /**
   * @brief Default destructor.
   */
  ~MQSM() = default;

  // Copy & move
  MQSM(const MQSM &other) = default;
  MQSM(MQSM &&other) noexcept = default;
  MQSM &operator=(const MQSM &other) = default;
  MQSM &operator=(MQSM &&other) noexcept = default;

  /**
   * @brief Constructs an MQSM from a string label.
   *
   * The string label can contain characters: '0', '1', '+', '-', 'r', 'l'
   * representing standard basis and superposition states for each qubit.
   *
   * @param label The state string label.
   * @return A newly constructed MQSM.
   */
  static MQSM fromLabel(string_view label) {
    QubitIndex num_levels = label.size();

    if (num_levels > max_qubits) {
      throw std::invalid_argument(
          "MQSM::fromLabel: Number of levels exceeds maximum allowed.");
    }
    if (num_levels == 0) {
      throw std::invalid_argument("MQSM::fromLabel: Label must not be empty.");
    }

    Complex inv_sqrt2 = Complex::inv_sqrt2();
    Complex i_inv_sqrt2 = Complex::i_inv_sqrt2();
    FSM<max_qubits, 2> fsm(num_levels);

    for (QubitIndex level = 0; level < num_levels; ++level) {
      auto &states = fsm.level_states[level];
      auto &transitions = fsm.level_transitions[level];
      char c = label[level];
      states.resize(1);
      auto &state = states[0];
      if (c == '0' || c == '1') {
        transitions.resize(1);
        state.T[(QubitIndex)(c - '0')] = 0;
        transitions[0] = {1, 0};

      } else {
        transitions.resize(2);
        state.T[0] = 0;
        state.T[1] = 1;
        transitions[0] = {inv_sqrt2, 0};
        if (c == '+') {
          transitions[1] = {inv_sqrt2, 0};
        } else if (c == '-') {
          transitions[1] = {-inv_sqrt2, 0};
        } else if (c == 'r') {
          transitions[1] = {i_inv_sqrt2, 0};
        } else if (c == 'l') {
          transitions[1] = {-i_inv_sqrt2, 0};
        } else {
          throw std::invalid_argument(
              "MQSM::fromLabel: Label must only contain '0', '1', '+', '-', "
              "'r', or 'l'. Character '" +
              std::string(1, c) + "' is invalid.");
        }
      }
    }
    return MQSM(std::move(fsm));
  }

  /**
   * @brief Constructs an MQSM representing a specific computational basis state
   * (ket).
   * @param ket_index The integer index of the basis state.
   * @param num_qubits The number of qubits in the state.
   * @return A newly constructed MQSM.
   */
  static MQSM fromBasis(BasisKetIndex ket_index, QubitIndex num_qubits) {
    if (num_qubits > max_qubits) {
      throw std::invalid_argument(
          "MQSM::fromBasis: Number of qubits exceeds maximum allowed.");
    }
    if (ket_index >= (1ULL << num_qubits)) {
      throw std::invalid_argument(
          "MQSM::fromBasis: ket_index must be less than 2^num_qubits.");
    }

    FSM<max_qubits, 2> fsm(num_qubits);
    for (QubitIndex level = 0; level < num_qubits; ++level) {
      auto &states = fsm.level_states[level];
      auto &transitions = fsm.level_transitions[level];
      states.resize(1);
      transitions.resize(1);
      QubitIndex bit_value = (ket_index >> (num_qubits - 1 - level)) & 1;
      transitions[0] = {1,
                        0}; // Transition to state 0 in next level with weight 1
      states[0].T[bit_value] = 0; // Transition to next state in the next level
    }
    return MQSM(std::move(fsm));
  }

  /**
   * @brief Constructs an MQSM from a dense statevector (complex amplitudes).
   * @param sv The statevector array.
   * @return A newly constructed MQSM.
   */
  static MQSM fromStatevector(const vector<Complex> &sv) {
    return fromStatevectorImpl<Complex>(sv);
  }

  static MQSM fromStatevector(const vector<std::complex<double>> &sv) {
    return fromStatevectorImpl<std::complex<double>>(sv);
  }

public: // ---------- INTERFACE ----------
  /**
   * @brief Gets the number of qubits the MQSM represents.
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
   * @brief Checks if the MQSM is empty.
   */
  bool empty() const {
    return fsm.getNumLevels() == 0 || fsm.level_states[0].empty();
  }

  /**
   * @brief In-place tensor product with another MQSM state.
   */
  MQSM &operator^=(const MQSM &other) {
    fsm ^= other.fsm;
    return *this;
  }

  /**
   * @brief Tensor product with another MQSM state.
   */
  MQSM operator^(const MQSM &other) const {
    MQSM result(fsm ^ other.fsm);
    return result;
  }

  /**
   * @brief Reduces the underlying FSM to eliminate redundant states and
   * transitions.
   * @return Reference to the reduced MQSM.
   */
  MQSM &reduce() {
    fsm.reduce();
    return *this;
  }

public: // ---------- EVOLUTION INTERFACE ----------
  /**
   * @brief Evolves the state in-place using a unitary MQGT operation.
   * @param U The unitary operation (MQGT) to apply.
   * @param qubits Optional bitmask specifying which qubits the gate applies to.
   * @return Reference to the evolved MQSM.
   */
  MQSM &evolve_inplace(const MQGT<max_qubits> &U,
                       BitMask<max_qubits> qubits = {}) {

    fsm.compose_inplace(U.getFSM(), qubits);
    return *this;
  }

  /**
   * @brief Evolves the state using a unitary MQGT operation.
   * @param U The unitary operation (MQGT) to apply.
   * @param qubits Optional bitmask specifying which qubits the gate applies to.
   * @return A new MQSM representing the evolved state.
   */
  MQSM evolve(const MQGT<max_qubits> &U,
              BitMask<max_qubits> qubits = {}) const {
    MQSM result = *this;
    return result.evolve_inplace(U, qubits);
  }

public: // ------- MEASUREMENT INTERFACE -------
  /**
   * @brief Computes the probabilities of measuring specific basis states for a
   * subset of qubits.
   *
   * @param qubits Span of qubit indices to measure. If empty, computes full
   * state probabilities.
   * @return SampleProbabilities A map from bitstring outcomes to probabilities.
   */
  SampleProbabilities
  sampleProbabilities(span<const QubitIndex> qubits = {}) const {

    map<string, double> probabilities;
    if (empty())
      return probabilities;

    bool full_state = qubits.empty();
    vector<QubitIndex> target_qubits;
    QubitIndex max_target_level = 0;

    if (!full_state) {
      target_qubits.assign(qubits.begin(), qubits.end());
      std::sort(target_qubits.begin(), target_qubits.end());
      if (target_qubits.empty())
        return probabilities;
      max_target_level = target_qubits.back();

      if (max_target_level >= getNumQubits()) {
        throw std::invalid_argument(
            "MQSM::measureProbabilities: Target qubit index exceeds number of "
            "qubits in the MQSM.");
      }

    } else {
      max_target_level = getNumQubits() - 1;
    }

    unordered_map<uint64_t, double> fast_probs;

    // Depth-First Search using real probability accumulation and early stopping
    auto dfs = [&](auto &self, QubitIndex level, MachineIndex node_idx,
                   double current_prob, uint64_t ket_int) -> void {
      if (current_prob <= 0.0)
        return;

      // Algorithmic Optimization: Early stopping for marginalization
      if (!full_state && level > max_target_level) {
        fast_probs[ket_int] += current_prob;
        return;
      }

      // Base case for full state
      if (level == getNumQubits()) {
        fast_probs[ket_int] += current_prob;
        return;
      }

      bool is_target =
          full_state ||
          std::binary_search(target_qubits.begin(), target_qubits.end(), level);
      const auto &state = fsm.level_states[level][node_idx];
      const auto &transitions = fsm.level_transitions[level];

      for (MachineInputSymbol symbol = 0; symbol < 2; ++symbol) {
        if (state.hasTransition(symbol)) {

          const MachineIndex transition_index =
              state.getTransitionIndex(symbol);
          const auto &transition = transitions[transition_index];
          double trans_prob = transition.weight.norm().toDouble();

          if (trans_prob > 0.0) {
            uint64_t new_ket_int = ket_int;
            if (is_target)
              new_ket_int = (ket_int << 1ULL) | symbol;

            const MachineIndex next_node = transition.next_state;
            self(self, level + 1, next_node, current_prob * trans_prob,
                 new_ket_int);
          }
        }
      }
    };

    // Root node starts with a probability mass of 1.0
    dfs(dfs, 0, 0, 1.0, 0ULL);

    // Convert to map of strings
    QubitIndex num_bits = full_state ? getNumQubits() : target_qubits.size();

    // Buffer for string conversion
    std::string buffer(num_bits, '0');
    for (const auto &[ket_int, prob] : fast_probs) {
      uint64_t x = ket_int;
      for (int b = num_bits - 1; b >= 0; --b) {
        buffer[b] = (x & 1) ? '1' : '0';
        x >>= 1;
      }
      probabilities[buffer] = prob;
    }

    return SampleProbabilities(std::move(probabilities));
  }

  SampleProbabilities
  sampleProbabilities(initializer_list<QubitIndex> qubits) const {
    return sampleProbabilities(
        span<const QubitIndex>(qubits.begin(), qubits.end()));
  }

  /**
   * @brief Simulates taking multiple measurement shots of the state.
   *
   * @param shots Number of times to measure the state.
   * @param qubits Span of qubit indices to measure. If empty, measures all
   * qubits.
   * @return SampleCount A map from bitstring outcomes to their occurrence
   * counts.
   */
  SampleCount sample(size_t shots, span<const QubitIndex> qubits = {}) const {
    map<string, size_t> samples;
    if (shots == 0 || empty())
      return SampleCount();

    bool full_state = qubits.empty();
    vector<QubitIndex> target_qubits;

    if (!full_state) {
      target_qubits.assign(qubits.begin(), qubits.end());
      std::sort(target_qubits.begin(), target_qubits.end());

      if (target_qubits.empty())
        return SampleCount();
      if (target_qubits.back() >= getNumQubits()) {
        throw std::invalid_argument("MQSM::sample: Target qubit index exceeds "
                                    "number of qubits in the MQSM.");
      }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::unordered_map<uint64_t, size_t> fast_samples;

    // Probabilistic DAG walk (O(N) per shot)
    for (size_t s = 0; s < shots; ++s) {
      uint64_t outcome_int = 0;

      MachineIndex current_node = 0;

      for (QubitIndex level = 0; level < getNumQubits(); ++level) {
        const auto &state = fsm.level_states[level][current_node];
        const auto &transitions = fsm.level_transitions[level];
        double p0 = 0.0, p1 = 0.0;

        if (state.hasTransition(0)) {
          MachineIndex t0_index = state.getTransitionIndex(0);
          p0 = transitions[t0_index].weight.norm().toDouble();
        }
        if (state.hasTransition(1)) {
          MachineIndex t1_index = state.getTransitionIndex(1);
          p1 = transitions[t1_index].weight.norm().toDouble();
        }

        double total_prob = p0 + p1;
        if (total_prob <= 0.0)
          break; // Failsafe for unnormalized or dead-end states

        // Roll the dice for this branch
        double r = dist(gen) * total_prob;
        MachineInputSymbol choice = (r < p0) ? 0 : 1;
        if (!state.hasTransition(choice)) { // Failsafe for dead-end transitions
          choice = 1 - choice; // Fallback to the other transition if available
        }

        // Only track the qubit if it's targeted
        bool is_target =
            full_state || std::binary_search(target_qubits.begin(),
                                             target_qubits.end(), level);
        if (is_target) {
          outcome_int = (outcome_int << 1ULL) | choice;
        }

        current_node = transitions[state.getTransitionIndex(choice)].next_state;
      }
      fast_samples[outcome_int]++;
    }

    // Convert to map of strings
    QubitIndex num_bits = full_state ? getNumQubits() : target_qubits.size();

    // Buffer for string conversion
    std::string buffer(num_bits, '0');
    for (const auto &[outcome_int, count] : fast_samples) {
      uint64_t x = outcome_int;
      for (int b = num_bits - 1; b >= 0; --b) {
        buffer[b] = (x & 1) ? '1' : '0';
        x >>= 1;
      }
      samples[buffer] = count;
    }

    return SampleCount(std::move(samples));
  }

  SampleCount sample(size_t shots, initializer_list<QubitIndex> qubits) const {
    return sample(shots, span<const QubitIndex>(qubits.begin(), qubits.end()));
  }

public: // ---------- DEBUG -----------
  FSM<max_qubits, 2> getFSM() const { return fsm; }

public: // ---------- I/O ----------
  string toString() const { return fsm.toString(); }

  vector<Complex> toStatevector() const { return toStatevectorImpl<Complex>(); }

  vector<std::complex<double>> toStdStatevector() const {
    return toStatevectorImpl<std::complex<double>>();
  }

  string toStatevectorString() const {
    vector<Complex> sv = toStatevector();
    ostringstream os;
    QubitIndex num_qubits = getNumQubits();
    bool first = true;

    // Memory Optimization: Hoist the string allocation OUTSIDE the loop!
    std::string buffer(num_qubits, '0');

    for (size_t i = 0; i < sv.size(); ++i) {
      if (sv[i].NearZero())
        continue;

      if (!first)
        os << " + ";
      first = false;

      if (sv[i].hasImag() && sv[i].hasReal()) {
        os << "(" << sv[i] << ")";
      } else {
        if (!sv[i].Near(1))
          os << sv[i];
      }

      // Fill buffer
      size_t x = i;
      for (int b = num_qubits - 1; b >= 0; --b) {
        buffer[b] = (x & 1) ? '1' : '0';
        x >>= 1;
      }
      os << "|" << buffer << ">";
    }
    return os.str();
  }

  friend std::ostream &operator<<(std::ostream &os, const MQSM &mqsm) {
    os << mqsm.toString();
    return os;
  }

public: // ---------- FRIENDS ----------
  friend class MQSMQuantumState<max_qubits>;

  friend class UnitaryOp<max_qubits>;
};

#endif