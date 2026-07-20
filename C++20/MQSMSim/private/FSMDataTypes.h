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
#ifndef FSMDATATYPES_H
#define FSMDATATYPES_H

#include <private/BasicDataTypes.h>
#include <private/BitMask.h>

#include <algorithm>
#include <array>
#include <initializer_list>
#include <ostream>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using std::array;
using std::initializer_list;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::set;
using std::span;
using std::string;
using std::string_view;
using std::tuple;
using std::unordered_map;
using std::vector;

// ------- MQSM / MQGT INDEXING TYPES -------

/**
 * @typedef MachineIndex
 * @brief Type for indexing states and transitions within a finite state
 * machine.
 */
typedef uint32_t MachineIndex;

/**
 * @def NullMachineIndex
 * @brief Sentinel value indicating an invalid or null machine index.
 */
const constexpr MachineIndex NullMachineIndex =
    std::numeric_limits<MachineIndex>::max();

/**
 * @typedef BasisKetIndex
 * @brief Type for indexing computational basis states (kets) in a statevector.
 */
typedef uint64_t BasisKetIndex;

/**
 * @typedef MachineInputSymbol
 * @brief Type for indexing machine input symbols (0, 1 for MQSM; 0..3 for
 * MQGT).
 */
typedef uint8_t MachineInputSymbol;

/**
 * @typedef MachineOutputWeight
 * @brief Type representing the complex weight/amplitude output of a transition.
 */
typedef Complex MachineOutputWeight;

// Advanced declarations
template <unsigned int max_q>
  requires(max_q > 0 && max_q <= MAX_QUBITS_ALLOWED)
class MQGT;

template <unsigned int max_q>
  requires(max_q > 0 && max_q <= MAX_QUBITS_ALLOWED)
class MQSM;

/**
 * @class FSM
 * @brief Represents a Finite State Machine structured by levels (e.g., qubits).
 *
 * An FSM consists of a sequence of levels. Each level contains states and
 * transitions. This structure forms the core backend for both MQSM (Multi-Qubit
 * State Machines) and MQGT (Multi-Qubit Gate Transducers).
 *
 * @tparam max_qubits Maximum number of levels (qubits) allowed.
 * @tparam num_input_symbols Number of input symbols per state (2 for MQSM, 4
 * for MQGT).
 */
template <unsigned int max_qubits, unsigned int num_input_symbols>
  requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED &&
           (num_input_symbols == 2 || num_input_symbols == 4))
class FSM {

private: // ---------- REPRESENTATION ----------
  // Transition output
  struct TransitionOutput {
    MachineOutputWeight weight;
    MachineIndex next_state;

    TransitionOutput() : weight(0), next_state(NullMachineIndex) {}
    TransitionOutput(const MachineOutputWeight &w, MachineIndex ns)
        : weight(w), next_state(ns) {}
    void clear() {
      weight = 0;
      next_state = NullMachineIndex;
    }
    bool empty() const {
      return next_state == NullMachineIndex || weight.NearZero();
    }
    bool operator==(const TransitionOutput &other) const {
      return next_state == other.next_state &&
             (weight - other.weight).NearZero();
    }
    bool operator!=(const TransitionOutput &other) const {
      return !(*this == other);
    }
  };

  // State
  struct FSMState {
    array<MachineIndex, num_input_symbols> T; // State transition indices
    FSMState() { T.fill(NullMachineIndex); }
    bool hasTransition(MachineInputSymbol symbol) const {
      return T[symbol] != NullMachineIndex;
    }
    MachineIndex getTransitionIndex(MachineInputSymbol symbol) const {
      return T[symbol];
    }
    bool empty() const {
      return std::all_of(T.begin(), T.end(), [](MachineIndex idx) {
        return idx == NullMachineIndex;
      });
    }
  };

  // States by level
  array<vector<FSMState>, max_qubits>
      level_states; // States for each level of the FSM
  array<vector<TransitionOutput>, max_qubits>
      level_transitions; // Transitions for each level of the FSM
  QubitIndex num_levels; // Number of levels in the FSM

private: // ---------- HELPERS ----------
  bool states_equivalent(QubitIndex level, MachineIndex s1_index,
                         MachineIndex s2_index) const {
    const FSMState &s1 = level_states[level][s1_index];
    const FSMState &s2 = level_states[level][s2_index];
    for (MachineInputSymbol symbol = 0; symbol < num_input_symbols; ++symbol) {
      if (s1.T[symbol] != s2.T[symbol]) {
        MachineIndex t1_index = s1.getTransitionIndex(symbol);
        MachineIndex t2_index = s2.getTransitionIndex(symbol);

        if (t1_index == NullMachineIndex || t2_index == NullMachineIndex) {
          return false; // One state has a transition while the other does not
        }
        const TransitionOutput &t1 = level_transitions[level][t1_index];
        const TransitionOutput &t2 = level_transitions[level][t2_index];
        if (t1 != t2)
          return false; // Transition weights or next states differ
      }
    }
    return true;
  }

  // State normalization for binary input symbols (MQSM)
  template <unsigned int T = num_input_symbols>
    requires(T == 2)
  MachineOutputWeight normalize_state(QubitIndex level,
                                      MachineIndex state_index) {

    MachineOutputWeight global_phase = 1.0;
    MachineOutputWeight w0 = 0, w1 = 0;
    FSMState &state = level_states[level][state_index];
    vector<TransitionOutput> &transitions = level_transitions[level];
    bool has_t0 = state.hasTransition(0), has_t1 = state.hasTransition(1);

    if (has_t0) {
      MachineIndex t_index = state.getTransitionIndex(0);

      if (t_index >= transitions.size())
        throw std::runtime_error(
            "FSM::normalize_state: Transition 0 points to invalid index.");
      TransitionOutput &t = transitions[t_index];
      if (t.empty())
        throw std::runtime_error(
            "FSM::normalize_state: Transition 0 is empty.");

      global_phase = t.weight;
      t.weight = w0 = 1;
    } else if (has_t1) {
      MachineIndex t_index = state.getTransitionIndex(1);

      if (t_index >= transitions.size())
        throw std::runtime_error(
            "FSM::normalize_state: Transition 1 points to invalid index.");
      TransitionOutput &t = transitions[t_index];
      if (t.empty())
        throw std::runtime_error(
            "FSM::normalize_state: Transition 1 is empty.");

      global_phase = t.weight;
      t.weight = w1 = 1;
    }

    if (has_t0 && has_t1) {
      MachineIndex t_index = state.getTransitionIndex(1);
      if (t_index >= transitions.size())
        throw std::runtime_error(
            "FSM::normalize_state: Transition 1 points to invalid index.");

      TransitionOutput &t = transitions[t_index];
      if (t.empty())
        throw std::runtime_error(
            "FSM::normalize_state: Transition 1 is empty.");
      t.weight /= global_phase;
      w1 = t.weight;
    }

    Real norm = (w1.norm() + w0.norm()).sqrt();
    if (norm.NearZero()) {
      return 0.0; // State has zero norm, return zero weight
    }

    if (has_t0)
      transitions[state.getTransitionIndex(0)].weight /= norm;
    if (has_t1)
      transitions[state.getTransitionIndex(1)].weight /= norm;

    MachineOutputWeight state_weight = global_phase * norm;
    return state_weight;
  }

  // State normalization for binary input symbols (MQGT)
  template <unsigned int T = num_input_symbols>
    requires(T == 4)
  MachineOutputWeight normalize_state(QubitIndex level,
                                      MachineIndex state_index,
                                      MachineInputSymbol output_symbol) {

    MachineOutputWeight w0 = 0, w1 = 0;
    FSMState &state = level_states[level][state_index];
    vector<TransitionOutput> &transitions = level_transitions[level];

    output_symbol <<=
        1u; // Shift left to get the correct output symbol for MQGT

    bool has_t0 = state.hasTransition(output_symbol | 0u),
         has_t1 = state.hasTransition(output_symbol | 1u);

    if (has_t0) {
      MachineIndex t_index = state.getTransitionIndex(output_symbol | 0u);

      if (t_index >= transitions.size())
        throw std::runtime_error("FSM::normalize_state: Transition (" +
                                 std::to_string((int)output_symbol) +
                                 ",0) points to invalid index.");
      const TransitionOutput &t = transitions[t_index];
      if (t.empty())
        throw std::runtime_error("FSM::normalize_state: Transition (" +
                                 std::to_string((int)output_symbol) +
                                 ",0) is empty.");
      w0 = t.weight;
    }

    if (has_t1) {
      MachineIndex t_index = state.getTransitionIndex(output_symbol | 1u);
      if (t_index >= transitions.size())
        throw std::runtime_error("FSM::normalize_state: Transition (" +
                                 std::to_string((int)output_symbol) +
                                 ",1) points to invalid index.");

      const TransitionOutput &t = transitions[t_index];
      if (t.empty())
        throw std::runtime_error("FSM::normalize_state: Transition (" +
                                 std::to_string((int)output_symbol) +
                                 ",1) is empty.");
      w1 = t.weight;
    }

    Real norm = (w1.norm() + w0.norm()).sqrt();
    if (norm.NearZero())
      throw std::runtime_error("FSM::normalize_state: State has zero norm.");

    if (has_t0)
      transitions[state.getTransitionIndex(output_symbol | 0u)].weight /= norm;
    if (has_t1)
      transitions[state.getTransitionIndex(output_symbol | 1u)].weight /= norm;

    MachineOutputWeight state_weight = norm;
    return state_weight;
  }

public: // ------- CONSTRUCTORS -------
  /**
   * @brief Default constructor creating an empty FSM with 0 levels.
   */
  FSM() : num_levels(0) {}

  /**
   * @brief Default destructor.
   */
  ~FSM() = default;

  /**
   * @brief Constructs an FSM with a specified number of levels.
   * @param levels The number of levels to initialize.
   * @pre levels <= max_qubits
   */
  FSM(QubitIndex levels) : num_levels(levels) {
    if (levels > max_qubits)
      throw std::invalid_argument(
          "FSM: Number of levels exceeds maximum allowed.");
  }

  // Copy & move
  FSM(const FSM &other) = default;
  FSM(FSM &&other) noexcept = default;
  FSM &operator=(const FSM &other) = default;
  FSM &operator=(FSM &&other) noexcept = default;

  // reserve: List of pairs (num_states, num_transitions) for each level
  // transitions: List of tuples (level, state, symbol, weight, next_state)
  static FSM
  createFSM(span<const pair<MachineIndex, MachineIndex>> reserve,
            span<const tuple<QubitIndex, MachineIndex, MachineInputSymbol,
                             MachineOutputWeight, MachineIndex>>
                transitions,
            unsigned int times = 1) {

    FSM fsm(reserve.size() * times);

    // Reserve states and transitions for each level
    QubitIndex level = 0;
    for (auto [num_states, num_transitions] : reserve) {
      auto &states = fsm.level_states[level];
      auto &transitions = fsm.level_transitions[level];
      states.resize(num_states);
      transitions.reserve(num_transitions);
      if (times > 1) {
        for (unsigned int t = 1; t < times; ++t) {
          auto &states_t = fsm.level_states[level + t * reserve.size()];
          auto &transitions_t =
              fsm.level_transitions[level + t * reserve.size()];
          states_t.resize(num_states);
          transitions_t.reserve(num_transitions);
        }
      }
      ++level;
    }

    // Add transitions
    for (auto [level, state, symbol, weight, next_state] : transitions) {
      fsm.addTransition(level, state, symbol, weight, next_state);
    }

    // Copy the reserved states and transitions for the specified number of
    // times
    for (unsigned int t = 1; t < times; ++t) {
      for (QubitIndex level = 0; level < reserve.size(); ++level) {
        auto &states = fsm.level_states[level + t * reserve.size()];
        auto &transitions = fsm.level_transitions[level + t * reserve.size()];
        states = fsm.level_states[level];
        transitions = fsm.level_transitions[level];
      }
    }
    return fsm;
  }

  static FSM
  createFSM(initializer_list<pair<MachineIndex, MachineIndex>> reserve,
            initializer_list<tuple<QubitIndex, MachineIndex, MachineInputSymbol,
                                   MachineOutputWeight, MachineIndex>>
                transitions,
            unsigned int times = 1) {

    return createFSM(
        span<const pair<MachineIndex, MachineIndex>>(reserve),
        span<const tuple<QubitIndex, MachineIndex, MachineInputSymbol,
                         MachineOutputWeight, MachineIndex>>(transitions),
        times);
  }

public: // ---------- INTERFACE ----------
  /**
   * @brief Gets the number of active levels (qubits) in the FSM.
   */
  QubitIndex getNumLevels() const noexcept { return num_levels; }

  /**
   * @brief Gets the total number of states across all levels.
   */
  size_t getNumStates() const noexcept {
    size_t total_states = 0;
    for (QubitIndex level = 0; level < num_levels; ++level) {
      total_states += level_states[level].size();
    }
    return total_states;
  }

  /**
   * @brief Gets the total number of transitions across all levels.
   */
  size_t getNumTransitions() const noexcept {
    size_t total_transitions = 0;
    for (QubitIndex level = 0; level < num_levels; ++level) {
      total_transitions += level_transitions[level].size();
    }
    return total_transitions;
  }

  /**
   * @brief Adds a new state to a specified level.
   * @param level The index of the level.
   * @return The local index of the newly created state.
   */
  MachineIndex addState(QubitIndex level) {
    if (level >= num_levels)
      throw std::out_of_range("FSM::addState: Level index out of range.");
    level_states[level].emplace_back();
    return static_cast<MachineIndex>(level_states[level].size() - 1);
  }

  /**
   * @brief Adds a transition to an existing state.
   * @param level The level of the source state.
   * @param state The index of the source state.
   * @param symbol The input symbol for the transition.
   * @param weight The complex weight associated with the transition.
   * @param next_state The index of the target state in the next level.
   * @return The local index of the newly created transition.
   */
  MachineIndex addTransition(QubitIndex level, MachineIndex state,
                             MachineInputSymbol symbol,
                             const MachineOutputWeight &weight,
                             MachineIndex next_state) {
    if (level >= num_levels)
      throw std::out_of_range("FSM::addTransition: Level index out of range.");
    if (state >= level_states[level].size())
      throw std::out_of_range("FSM::addTransition: State index out of range.");
    if (symbol >= num_input_symbols)
      throw std::out_of_range(
          "FSM::addTransition: Input symbol index out of range.");
    if (level == num_levels - 1 && next_state != 0)
      throw std::invalid_argument(
          "FSM::addTransition: Next state must be 0 for the last level.");
    if (level < num_levels - 1 && next_state >= level_states[level + 1].size())
      throw std::out_of_range("FSM::addTransition: Next state index out of "
                              "range for the next level.");

    level_transitions[level].emplace_back(weight, next_state);
    MachineIndex transition_index =
        static_cast<MachineIndex>(level_transitions[level].size() - 1);
    level_states[level][state].T[symbol] = transition_index;
    return transition_index;
  }

  /**
   * @brief Checks if a specific state has a valid transition for a given
   * symbol.
   */
  bool stateHasTransition(QubitIndex level, MachineIndex state,
                          MachineInputSymbol symbol) const {
    if (level >= num_levels)
      throw std::out_of_range(
          "FSM::stateHasTransition: Level index out of range.");
    if (state >= level_states[level].size())
      throw std::out_of_range(
          "FSM::stateHasTransition: State index out of range.");
    if (symbol >= num_input_symbols)
      throw std::out_of_range(
          "FSM::stateHasTransition: Input symbol index out of range.");

    auto &state_repr = level_states[level][state];
    if (!state_repr.hasTransition(symbol))
      return false;
    auto &transition =
        level_transitions[level][state_repr.getTransitionIndex(symbol)];
    return !transition.empty();
  }

  FSM &operator^=(const FSM &other) {
    if (getNumLevels() == 0) {
      *this = other;
      return *this;
    } else if (other.getNumLevels() == 0) {
      return *this;
    }
    if (getNumLevels() + other.getNumLevels() > max_qubits) {
      throw std::invalid_argument("FSM::operator^=: Resulting number of levels "
                                  "exceeds maximum allowed.");
    }

    std::copy_n(other.level_states.begin(), other.getNumLevels(),
                level_states.begin() + getNumLevels());
    std::copy_n(other.level_transitions.begin(), other.getNumLevels(),
                level_transitions.begin() + getNumLevels());
    num_levels += other.getNumLevels();
    return *this;
  }

  FSM operator^(const FSM &other) const {
    FSM result(num_levels + other.num_levels);
    std::copy_n(level_states.begin(), num_levels, result.level_states.begin());
    std::copy_n(level_transitions.begin(), num_levels,
                result.level_transitions.begin());
    std::copy_n(other.level_states.begin(), other.num_levels,
                result.level_states.begin() + num_levels);
    std::copy_n(other.level_transitions.begin(), other.num_levels,
                result.level_transitions.begin() + num_levels);
    return result;
  }

  FSM &reduce(QubitIndex bottom_level = max_qubits, QubitIndex top_level = 0) {

    if (bottom_level < top_level)
      std::swap(bottom_level, top_level);

    QubitIndex num_qubits = num_levels;
    if (bottom_level == max_qubits)
      bottom_level = num_qubits - 1;

    struct StateInfo {
      array<MachineIndex, num_input_symbols> next_states;
      array<MachineOutputWeight, num_input_symbols> weights;
      bool operator==(const StateInfo &other) const {
        return next_states == other.next_states && weights == other.weights;
      }

      StateInfo() {
        next_states.fill(NullMachineIndex);
        weights.fill(0);
      }
    };

    struct hashStateInfo {
      std::size_t operator()(const StateInfo &key) const noexcept {
        std::size_t h = 0;
        // Only hash active transitions to prevent heavy collisions on sparse
        // nodes
        for (QubitIndex i = 0; i < num_input_symbols; ++i) {
          if (key.next_states[i] != NullMachineIndex) {
            h ^= std::hash<size_t>{}(key.next_states[i]) + 0x9e3779b9 +
                 (h << 6) + (h >> 2);
            h ^= std::hash<size_t>{}(i) + 0x9e3779b9 + (h << 6) + (h >> 2);
            Complex rounded_w = key.weights[i].toleranceRound();
            std::size_t hw_r = std::hash<double>{}(rounded_w.real().toDouble());
            std::size_t hw_i = std::hash<double>{}(rounded_w.imag().toDouble());
            h ^= hw_r + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= hw_i + 0x9e3779b9 + (h << 6) + (h >> 2);
          }
        }
        return h;
      }
    };

    // Bottom-up strategy. Last level points to a single state 0 that retursn
    // weight=1
    vector<MachineIndex> current_isomorphisms, next_isomorphisms;
    vector<MachineOutputWeight> current_weights, next_weights;

    // Set state to the last level's state 0 with weight 1
    if (bottom_level == num_qubits - 1) {
      next_isomorphisms.emplace_back(0);
      next_weights.emplace_back(1.0);

    } else { // Set state to the next level's unique states with weight 1
      next_isomorphisms.resize(level_states[bottom_level + 1].size());
      std::iota(next_isomorphisms.begin(), next_isomorphisms.end(), 0);
      next_weights.resize(next_isomorphisms.size(), 1.0);
    }

    MachineIndex current_state_count, current_transition_count;

    for (QubitIndex level = bottom_level + 1; level > 0; --level) {
      QubitIndex level_idx = level - 1;

      auto &current_states = level_states[level_idx];
      auto &current_transitions = level_transitions[level_idx];

      current_state_count = current_states.size();
      current_transition_count = current_transitions.size();

      // For each state in the current level, we have to:
      //   0. Detach dead transitions from states
      //   1. Update each transition's weight with the propagated weight from
      //   the next level
      //   2. Set the transitions' next state indices to point to unique states
      //   in the next level
      //   3. Normalize the state and store its propagated weight
      //   4. Identify unique states in the current level
      //   5. Detach transitions from duplicate states
      //   6. Update transitions and states according to the new configuration

      current_isomorphisms.clear();
      current_isomorphisms.resize(current_state_count);
      std::iota(current_isomorphisms.begin(), current_isomorphisms.end(), 0);
      current_weights.clear();
      current_weights.resize(current_state_count, 0.0);

      unordered_map<StateInfo, MachineIndex, hashStateInfo> unique_nodes;
      unique_nodes.reserve(current_state_count);

      for (MachineIndex state_idx = 0; state_idx < current_states.size();
           ++state_idx) {

        auto &state = current_states[state_idx];
        StateInfo s_info;
        bool dead_state = true;

        // Check transitions
        for (MachineInputSymbol symbol = 0; symbol < num_input_symbols;
             ++symbol) {
          if (state.hasTransition(symbol)) {

            auto &transition =
                current_transitions[state.getTransitionIndex(symbol)];

            if (!transition.empty()) {
              // Update transition weight
              transition.weight *= next_weights[transition.next_state];

              if (transition.weight.NearZero()) {
                transition.clear();
                state.T[symbol] = NullMachineIndex; // Detach dead transition
                current_transition_count--;

              } else { // Point to unique state in the next level
                transition.next_state =
                    next_isomorphisms[transition.next_state];
                dead_state = false;
              }
            } else {
              state.T[symbol] = NullMachineIndex; // Detach dead transition
              current_transition_count--;
            }
          }
        }

        // Normalize the state and store its weight for propagation to the
        // previous level ONLY FOR MQSM-type FSM
        if constexpr (num_input_symbols == 2) {
          if (!dead_state)
            current_weights[state_idx] = normalize_state(level_idx, state_idx);
          else
            current_weights[state_idx] = 0.0;
        } else {
          current_weights[state_idx] = (dead_state ? 0.0 : 1.0);
        }

        if (!dead_state) {

          // Update s_info with normalized transition weights and next states
          for (MachineInputSymbol symbol = 0; symbol < num_input_symbols;
               ++symbol) {
            if (state.hasTransition(symbol)) {
              auto &transition =
                  current_transitions[state.getTransitionIndex(symbol)];
              s_info.next_states[symbol] = transition.next_state;
              s_info.weights[symbol] = transition.weight;
            }
          }

          // Check uniqueness of the state in the current level
          auto [it, inserted] =
              unique_nodes.try_emplace(s_info, unique_nodes.size());
          if (!inserted) {
            current_state_count--;

            // Detach transitions from duplicate state
            for (MachineInputSymbol symbol = 0; symbol < num_input_symbols;
                 ++symbol) {
              if (state.hasTransition(symbol)) {
                current_transitions[state.getTransitionIndex(symbol)].clear();
                state.T[symbol] = NullMachineIndex;
                current_transition_count--;
              }
            }
          }
          current_isomorphisms[state_idx] =
              it->second; // Map to the unique state
        } else {
          current_state_count--;
        }
      }

      // Pack the current level's states and transitions to remove duplicates
      // and dead transitions
      vector<FSMState> packed_states;
      vector<TransitionOutput> packed_transitions;

      packed_states.resize(current_state_count);
      packed_transitions.reserve(current_transition_count);

      // Iterate unique nodes
      for (const auto &[s_info, unique_idx] : unique_nodes) {
        auto &state = packed_states[unique_idx];
        for (MachineInputSymbol symbol = 0; symbol < num_input_symbols;
             ++symbol) {
          if (s_info.next_states[symbol] == NullMachineIndex)
            continue; // Skip dead transitions

          state.T[symbol] =
              packed_transitions.size(); // Point to the next transition index
          packed_transitions.emplace_back(s_info.weights[symbol],
                                          s_info.next_states[symbol]);
        }
      }

      if (packed_states.empty()) {
        throw std::runtime_error("FSM::reduce: All states in level " +
                                 std::to_string((int)level_idx) +
                                 " are dead after reduction.");
      }
      if (packed_transitions.empty()) {
        throw std::runtime_error("FSM::reduce: All transitions in level " +
                                 std::to_string((int)level_idx) +
                                 " are dead after reduction.");
      }

      // Update the FSM with packed states and transitions
      level_states[level_idx] = std::move(packed_states);
      level_transitions[level_idx] = std::move(packed_transitions);

      // Prepare next level iteration
      std::swap(current_isomorphisms, next_isomorphisms);
      std::swap(current_weights, next_weights);

      // Check if current level was the last level to reduce
      if (level_idx <= top_level) {
        bool stop_reduction = true;

        // Check if no weights must be propagated
        auto fw = next_weights[0];
        for (const auto &w : next_weights) {
          if (!w.Near(fw)) {
            stop_reduction = false;
            break;
          }
        }
        if (stop_reduction)
          break;
      }
    }

    return *this;
  }

  FSM &compose_inplace(const FSM<max_qubits, 4> &U,
                       BitMask<max_qubits> levels = {}) {

    if (getNumLevels() == 0)
      throw std::runtime_error(
          "FSM::compose_inplace: Cannot compose an empty FSM.");
    if (U.getNumLevels() == 0)
      return *this; // Composing with an empty FSM does nothing

    if (!levels.empty() && levels.count() != U.getNumLevels()) {
      throw std::invalid_argument(
          "FSM::compose_inplace: The number of levels (" +
          std::to_string(levels.count()) +
          ") must match the number of levels in the MQGT (" +
          std::to_string(U.getNumLevels()) + ").");
    }

    // Error checks and target qubit preparation
    bool full_machine = levels.empty();
    vector<QubitIndex> target_levels;

    QubitIndex min_target_level = 0;

    if (!full_machine) {
      target_levels = levels.toVector();
      std::sort(target_levels.begin(), target_levels.end());

    } else {
      target_levels.resize(U.getNumLevels());
      std::iota(target_levels.begin(), target_levels.end(), 0);
    }
    min_target_level = target_levels[0];

    if (target_levels.size() != U.getNumLevels()) {
      throw std::invalid_argument(
          "FSM::compose_inplace: The number of levels in the MQGT must match "
          "the size of the provided levels list.");
    }

    QubitIndex final_level = getNumLevels() - 1;

    // Representation of composite states
    struct StatePair {
      MachineIndex u;
      MachineIndex s;
      MachineOutputWeight weight;

      bool operator==(const StatePair &other) const {
        return u == other.u && s == other.s && weight.Near(other.weight);
      }
      bool operator!=(const StatePair &other) const {
        return !(*this == other);
      }
      bool operator<(const StatePair &other) const {
        if (u != other.u)
          return u < other.u;
        if (s != other.s)
          return s < other.s;
        if (!weight.real().Near(other.weight.real()))
          return weight.real() < other.weight.real();
        return weight.imag() < other.weight.imag() &&
               !weight.imag().Near(other.weight.imag());
      }

      StatePair(MachineIndex u_, MachineIndex s_, MachineOutputWeight w_ = 1.0)
          : u(u_), s(s_), weight(w_) {}
    };

    struct CompositeState {
      vector<StatePair> pairs;

      bool empty() const { return pairs.empty(); }

      void add(MachineIndex u, MachineIndex s, MachineOutputWeight w) {
        if (w.NearZero())
          return;
        for (auto &p : pairs) {
          if (p.u == u && p.s == s) {
            p.weight += w;
            return;
          }
        }
        pairs.emplace_back(u, s, w);
      }

      void normalize(MachineOutputWeight &extracted_weight) {
        vector<StatePair> valid_pairs;
        for (auto &p : pairs) {
          if (!p.weight.NearZero())
            valid_pairs.push_back(p);
        }
        pairs = std::move(valid_pairs);
        if (pairs.empty()) {
          extracted_weight = 0;
          return;
        }
        std::sort(pairs.begin(), pairs.end());
        extracted_weight = pairs[0].weight;
        for (auto &p : pairs) {
          p.weight /= extracted_weight;
        }
      }

      bool operator==(const CompositeState &other) const {
        if (pairs.size() != other.pairs.size())
          return false;
        for (size_t i = 0; i < pairs.size(); ++i) {
          if (pairs[i] != other.pairs[i])
            return false;
        }
        return true;
      }
    };

    // HASH FOR COMPOSITE STATE
    struct hashCompositeState {
      std::size_t operator()(const CompositeState &state) const noexcept {
        std::size_t h = 0;
        for (const auto &p : state.pairs) {
          const std::size_t hu = std::hash<MachineIndex>{}(p.u);
          const std::size_t hs = std::hash<MachineIndex>{}(p.s);
          Complex rounded_w = p.weight.toleranceRound();
          const std::size_t hw_r =
              std::hash<double>{}(rounded_w.real().toDouble());
          const std::size_t hw_i =
              std::hash<double>{}(rounded_w.imag().toDouble());
          const std::size_t pair_hash =
              hu ^ (hs + 0x9e3779b9 + (hu << 6) + (hu >> 2)) ^
              (hw_r + 0x9e3779b9 + (hw_r << 6) + (hw_r >> 2)) ^
              (hw_i + 0x9e3779b9 + (hw_i << 6) + (hw_i >> 2));
          h ^= pair_hash + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
      }
    };

    // Get the FSM source submachine
    vector<vector<FSMState>> level_states_src(final_level - min_target_level +
                                              1);
    vector<vector<TransitionOutput>> level_transitions_src(
        final_level - min_target_level + 1);
    std::move(level_states.begin() + min_target_level,
              level_states.begin() + final_level + 1, level_states_src.begin());
    std::move(level_transitions.begin() + min_target_level,
              level_transitions.begin() + final_level + 1,
              level_transitions_src.begin());

    // Data structures to store composite states and transitions of current
    // level
    unordered_map<CompositeState, MachineIndex, hashCompositeState>
        composite_state_map; // Map to track unique composite states
    unordered_map<CompositeState, MachineIndex, hashCompositeState>
        next_composite_state_map; // Map to track unique composite states

    // Create the composite states for the initial level
    {
      auto &states = level_states_src[0];
      composite_state_map.reserve(states.size());
      for (MachineIndex s_idx = 0; s_idx < states.size(); ++s_idx) {
        CompositeState init_state;
        init_state.pairs.push_back({0, s_idx, 1.0});
        composite_state_map[init_state] = s_idx;
      }
    }

    QubitIndex current_level = min_target_level;
    QubitIndex next_processing_level_idx = 0;

    // Perform breath-first composition
    while (current_level <= final_level) {

      bool is_target_level =
          (next_processing_level_idx < target_levels.size() &&
           current_level == target_levels[next_processing_level_idx]);

      // Get current level's states and transitions from U and the FSM
      const auto &U_states = U.level_states[next_processing_level_idx];
      const auto &U_transitions =
          U.level_transitions[next_processing_level_idx];
      const auto &s_states = level_states_src[current_level - min_target_level];
      const auto &s_transitions =
          level_transitions_src[current_level - min_target_level];

      // Get the target level's states and transitions in the FSM
      auto &states_target = level_states[current_level];
      states_target.resize(composite_state_map.size());

      auto &transitions_target = level_transitions[current_level];
      transitions_target.clear();
      transitions_target.reserve(
          s_transitions
              .size()); // Initial reserve based on the source FSM's transitions

      // Initialize next map of composite states
      next_composite_state_map.clear();
      next_composite_state_map.reserve(composite_state_map.size() *
                                       U_states.size());

      // For each composite state in the level, we need to apply the transitions
      // from U and the FSM
      for (const auto &[composite_state, state_idx] : composite_state_map) {

        // Get the actual state in the target FSM that corresponds to this
        // composite state
        auto &state_target = states_target[state_idx];

        // For the current composite state, we need to get transitions for each
        // input symbol
        for (MachineInputSymbol symbol = 0; symbol < num_input_symbols;
             ++symbol) {

          TransitionOutput
              new_transition; // New transition for the composite state
          CompositeState new_composite_state; // Target composite state

          // For each pair in the composite state, apply the transitions from U
          // and the FSM
          for (const auto &pair : composite_state.pairs) {
            const MachineIndex u_state = pair.u;
            const MachineIndex s_state = pair.s;
            const MachineOutputWeight pair_weight = pair.weight;

            // Composite states may carry dead or out-of-range indices while
            // composing sparse levels.
            if (s_state >= s_states.size())
              continue;

            // Get possible ways that U can output the symbol
            for (MachineInputSymbol x = 0; x < 2; ++x) {
              MachineInputSymbol s_symbol, u_symbol;
              if constexpr (num_input_symbols == 2) {
                s_symbol = x;
                u_symbol =
                    (symbol << 1u) | x; // U's input symbol is the composite of
                                        // FSM's input and output
              } else {
                MachineInputSymbol w = symbol & 1;
                MachineInputSymbol y = symbol >> 1;
                s_symbol = (x << 1u) | w;
                u_symbol = (y << 1u) | x;
              }

              if (!s_states[s_state].hasTransition(s_symbol))
                continue;
              const auto &st =
                  s_transitions[s_states[s_state].getTransitionIndex(s_symbol)];

              // Apply transitions from U or identity (when level is not
              // targeted).
              if (is_target_level) {
                if (u_state >= U_states.size())
                  continue;
                const auto &u_state_repr = U_states[u_state];
                if (!u_state_repr.hasTransition(u_symbol))
                  continue;

                const auto &ut =
                    U_transitions[u_state_repr.getTransitionIndex(u_symbol)];

                new_composite_state.add(ut.next_state, st.next_state,
                                        pair_weight * ut.weight * st.weight);

              } else {
                // Identity MQGT contributes only when output == input.
                MachineInputSymbol y =
                    (num_input_symbols == 2) ? symbol : (symbol >> 1);
                if (x != y)
                  continue;

                new_composite_state.add(u_state, st.next_state,
                                        pair_weight * st.weight);
              }
            }
          }

          MachineOutputWeight extracted_weight = 0;
          new_composite_state.normalize(extracted_weight);
          new_transition.weight = extracted_weight;

          // Check if the new transition is not empty
          if (!new_transition.weight.NearZero() &&
              !new_composite_state.empty()) {

            // On the final level of composition, transitions must point to the
            // target s_state or to 0 if it's the very last level of the entire
            // FSM.
            if (current_level == final_level) {
              if (new_composite_state.empty()) {
                new_transition.next_state = NullMachineIndex;
              } else {
                new_transition.next_state = new_composite_state.pairs[0].s;
              }

            } else {
              // Get composite state index for the next level, or create if it
              // doesn't exist
              new_transition.next_state = next_composite_state_map.size();
              auto [it, inserted] = next_composite_state_map.try_emplace(
                  new_composite_state, new_transition.next_state);
              if (!inserted) {
                new_transition.next_state =
                    it->second; // Already exists, get the existing index
              }
            }

            state_target.T[symbol] =
                transitions_target.size(); // Point to the new transition
            transitions_target.push_back(
                new_transition); // Add the new transition to the target FSM

          } else {
            state_target.T[symbol] = NullMachineIndex; // No valid transition
          }
        }
      }

      if (is_target_level)
        ++next_processing_level_idx; // Move to the next target level

      // Move transitions
      std::swap(composite_state_map, next_composite_state_map);

      ++current_level; // Move to the next level in the FSM
    }

    this->reduce(final_level, 0);
    return *this;
  }

  FSM compose(const FSM<max_qubits, 4> &U,
              BitMask<max_qubits> levels = {}) const {
    FSM result = *this;
    result.compose_inplace(U, levels);
    return result;
  }

public: // --------- DEBUGGING ---------
  bool isValid() const {

    if (num_levels == 0)
      return true; // An empty FSM is valid

    // First level contains a single state
    if (level_states[0].size() != 1)
      return false;

    // Check that all transitions point to valid next states
    // and that all states have valid transitions
    // Also that there are no dead states (states with no outgoing transitions)
    // nor transitions to dead states
    for (QubitIndex level = 0; level < num_levels; ++level) {

      const auto &states = level_states[level];
      const auto &transitions = level_transitions[level];
      set<MachineIndex> reachable_states;
      set<MachineIndex> used_transitions;
      for (MachineIndex state_index = 0; state_index < states.size();
           ++state_index) {
        const FSMState &state = states[state_index];
        if (state.empty())
          return false; // Dead state

        for (MachineInputSymbol symbol = 0; symbol < num_input_symbols;
             ++symbol) {
          if (state.hasTransition(symbol)) {
            MachineIndex t_index = state.getTransitionIndex(symbol);
            if (t_index >= transitions.size())
              return false; // Invalid transition index
            const TransitionOutput &t = transitions[t_index];
            if (t.empty())
              return false; // Dead transition (no next state or zero weight)
            if (level < num_levels - 1 &&
                t.next_state >= level_states[level + 1].size())
              return false; // Invalid next state index
            if (level == num_levels - 1 && t.next_state != 0)
              return false; // Last level must transition to state 0
            reachable_states.insert(t.next_state);
            used_transitions.insert(t_index);
          }
        }
      }

      if (used_transitions.size() != transitions.size())
        return false; // There are unused transitions
      if (level < num_levels - 1) {
        // Check that all states in the next level are reachable
        if (reachable_states.size() != level_states[level + 1].size())
          return false;
      }
    }
    return true;
  }

  vector<string> getValidationLog() const {
    vector<string> log;

    if (num_levels == 0) {
      log.push_back("FSM is empty.");
      return log;
    }

    // First level contains a single state
    if (level_states[0].size() != 1) {
      log.push_back("Level 0 must contain a single state, but contains " +
                    std::to_string(level_states[0].size()) + " states.");
    }

    // Check that all transitions point to valid next states
    // and that all states have valid transitions
    // Also that there are no dead states (states with no outgoing transitions)
    // nor dead transitions
    for (QubitIndex level = 0; level < num_levels; ++level) {

      const auto &states = level_states[level];
      const auto &transitions = level_transitions[level];
      set<MachineIndex> reachable_states;
      set<MachineIndex> used_transitions;
      for (MachineIndex state_index = 0; state_index < states.size();
           ++state_index) {
        const FSMState &state = states[state_index];
        if (state.empty())
          log.push_back("Level " + std::to_string((int)level) + ", State " +
                        std::to_string(state_index) +
                        ": Dead state (no outgoing transitions).");
        else {

          for (MachineInputSymbol symbol = 0; symbol < num_input_symbols;
               ++symbol) {
            if (state.hasTransition(symbol)) {
              MachineIndex t_index = state.getTransitionIndex(symbol);
              if (t_index >= transitions.size()) {
                log.push_back("Level " + std::to_string((int)level) +
                              ", State " + std::to_string(state_index) +
                              ", Symbol " + std::to_string((int)symbol) +
                              ": Invalid transition index " +
                              std::to_string(t_index) + ".");
              } else {
                const TransitionOutput &t = transitions[t_index];
                if (t.empty()) {
                  log.push_back(
                      "Level " + std::to_string((int)level) + ", State " +
                      std::to_string(state_index) + ", Symbol " +
                      std::to_string((int)symbol) + ": Dead transition " +
                      std::to_string(t_index) +
                      " (next state= " + std::to_string(t.next_state) +
                      ", weight= " + t.weight.toString() + ").");

                } else {
                  if (level < num_levels - 1 &&
                      t.next_state >= level_states[level + 1].size()) {
                    log.push_back("Level " + std::to_string((int)level) +
                                  ", State " + std::to_string(state_index) +
                                  ", Symbol " + std::to_string((int)symbol) +
                                  ": Invalid next state index " +
                                  std::to_string(t.next_state) +
                                  " for next level.");
                  }
                  if (level == num_levels - 1 && t.next_state != 0) {
                    log.push_back("Level " + std::to_string((int)level) +
                                  ", State " + std::to_string(state_index) +
                                  ", Symbol " + std::to_string((int)symbol) +
                                  ": Last level must transition to state 0, "
                                  "but next state is " +
                                  std::to_string(t.next_state) + ".");
                  }
                  reachable_states.insert(t.next_state);
                  used_transitions.insert(t_index);
                }
              }
            }
          }
        }
      }

      if (used_transitions.size() != transitions.size()) {
        log.push_back("Level " + std::to_string((int)level) +
                      ": There are unused transitions.");
      }
      if (level < num_levels - 1) {
        // Check that all states in the next level are reachable
        if (reachable_states.size() != level_states[level + 1].size()) {
          log.push_back("Level " + std::to_string((int)level) +
                        ": Not all states in the next level are reachable.");
        }
      }
    }

    return log;
  }

public: // ---------- I/O ----------
  string toString(const string &prefix = "s") const {
    ostringstream os;
    for (QubitIndex level = 0; level < num_levels; ++level) {
      os << "Level " << (int)level << ":\n";
      for (MachineIndex state_index = 0;
           state_index < level_states[level].size(); ++state_index) {
        string state_str = prefix + "(" + std::to_string((int)level) + "," +
                           std::to_string(state_index) + ")";
        os << "  " << state_str << ":\n";
        for (MachineInputSymbol symbol = 0; symbol < num_input_symbols;
             ++symbol) {
          if (level_states[level][state_index].hasTransition(symbol)) {
            MachineIndex t_index =
                level_states[level][state_index].getTransitionIndex(symbol);
            const TransitionOutput &t = level_transitions[level][t_index];

            string next_state_str = prefix + "(" +
                                    std::to_string((int)(level + 1)) + "," +
                                    std::to_string(t.next_state) + ")";
            string symbol_str;
            if constexpr (num_input_symbols == 2) {
              symbol_str = std::to_string((int)symbol);
            } else if constexpr (num_input_symbols == 4) {
              symbol_str = "<" +
                           std::to_string((unsigned int)(symbol >> 1) & 1) +
                           "," + std::to_string((unsigned int)symbol & 1) + ">";
            } else {
              symbol_str = std::to_string(symbol);
            }

            os << "   [" << state_str << " , " << symbol_str << "] -> ["
               << next_state_str << " , " << t.weight << "]";
            os << "\n";
          }
        }
      }
    }
    return os.str();
  }

  friend std::ostream &operator<<(std::ostream &os, const FSM &fsm) {
    os << fsm.toString();
    return os;
  }

public: // ---------- FRIENDS -----------
  friend class MQGT<max_qubits>;

  friend class MQSM<max_qubits>;

  template <unsigned int N, unsigned int M>
    requires(N > 0 && N <= MAX_QUBITS_ALLOWED && (M == 2 || M == 4))
  friend class FSM;
};

#endif