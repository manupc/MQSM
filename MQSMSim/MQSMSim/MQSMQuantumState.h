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
#ifndef __MQSMQUANTUMSTATE__H__
#define __MQSMQUANTUMSTATE__H__

#include <bit>
#include <complex>
#include <cstdint>
#include <initializer_list>
#include <private/BasicDataTypes.h>
#include <private/MQSM.h>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using std::complex;
using std::initializer_list;
using std::invalid_argument;
using std::span;
using std::string;
using std::string_view;
using std::vector;

// Advanced declarations
template <unsigned int max_qubits>
  requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
class UnitaryOp;

class MQGTSim;

/**
 * @class MQSMQuantumState
 * @brief High-level representation of a quantum state using Multi-Qubit State
 * Machines (MQSM).
 *
 * Provides a user-friendly interface to create, manipulate, and measure quantum
 * states relying on the efficient underlying MQSM data structure.
 *
 * @tparam max_qubits The maximum number of qubits this state can hold.
 */
template <unsigned int max_qubits>
  requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
class MQSMQuantumState {

private: // ---------- REPRESENTATION ----------
  MQSM<max_qubits> mqsm;

  MQSMQuantumState(const MQSM<max_qubits> &mqsm) : mqsm(mqsm) {}
  MQSMQuantumState(MQSM<max_qubits> &&mqsm) : mqsm(std::move(mqsm)) {}

public: // ---------- CONSTRUCTORS ----------
  /**
   * @brief Default constructor creating an empty quantum state.
   */
  MQSMQuantumState() = default;

  /**
   * @brief Default destructor.
   */
  ~MQSMQuantumState() = default;

  MQSMQuantumState(const MQSMQuantumState &other) = default;
  MQSMQuantumState(MQSMQuantumState &&other) noexcept = default;
  MQSMQuantumState &operator=(const MQSMQuantumState &other) = default;
  MQSMQuantumState &operator=(MQSMQuantumState &&other) noexcept = default;

  /**
   * @brief Constructs a quantum state initialized to the |0...0> basis state.
   *
   * @param num_qubits The number of active qubits for this state.
   *
   * @pre num_qubits > 0
   * @pre num_qubits <= max_qubits
   * @post The state is initialized to |0>^{\otimes num_qubits}.
   *
   * @par Example:
   * @code
   * QuantumState state(5); // Creates a 5-qubit state |00000>
   * @endcode
   */
  MQSMQuantumState(QubitIndex num_qubits) {
    if (num_qubits <= 0) {
      throw std::invalid_argument(
          "MQSMQuantumState: Number of qubits must be positive.");
    }

    if (num_qubits > max_qubits) {
      throw std::invalid_argument("MQSMQuantumState: Number of qubits exceeds "
                                  "the maximum allowed limit of " +
                                  std::to_string(max_qubits) + ".");
    }

    mqsm = MQSM<max_qubits>::fromBasis(0, num_qubits);
  }

  /**
   * @brief Constructs a quantum state from a string label (e.g., "01+").
   *
   * @param label A string view representing the state of each qubit.
   * @return The constructed MQSMQuantumState.
   *
   * @pre The length of the label must not exceed max_qubits.
   * @pre The label characters must be valid state descriptors ('0', '1', '+',
   * '-', etc.).
   * @post The returned state matches the tensor product of the specified
   * single-qubit states.
   *
   * @par Example:
   * @code
   * auto state = QuantumState::fromLabel("01+");
   * @endcode
   */
  static MQSMQuantumState fromLabel(std::string_view label) {

    return MQSMQuantumState(MQSM<max_qubits>::fromLabel(label));
  }

  /**
   * @brief Constructs a quantum state from a full statevector (std::complex).
   *
   * @param sv A vector of std::complex<double> amplitudes representing the
   * statevector.
   * @return The constructed MQSMQuantumState.
   *
   * @pre sv.size() must be a power of 2, and log2(sv.size()) <= max_qubits.
   * @post The constructed MQSM matches the provided statevector amplitudes.
   *
   * @par Example:
   * @code
   * std::vector<std::complex<double>> sv = {{0.707, 0}, {0}, {0}, {0.707, 0}};
   * auto state = QuantumState::fromStatevector(sv);
   * @endcode
   */
  static MQSMQuantumState
  fromStatevector(const vector<std::complex<double>> &sv) {

    return MQSMQuantumState(MQSM<max_qubits>::fromStatevector(sv));
  }

  /**
   * @brief Constructs a quantum state from a full statevector using the custom
   * Complex type.
   *
   * @param sv A vector of Complex amplitudes representing the statevector.
   * @return The constructed MQSMQuantumState.
   *
   * @pre sv.size() must be a power of 2, and log2(sv.size()) <= max_qubits.
   * @post The constructed MQSM matches the provided statevector amplitudes.
   *
   * @par Example:
   * @code
   * std::vector<Complex> sv = {Complex(0.707, 0), Complex(0), Complex(0),
   * Complex(0.707, 0)}; auto state = QuantumState::fromStatevector(sv);
   * @endcode
   */
  static MQSMQuantumState fromStatevector(const vector<Complex> &sv) {

    return MQSMQuantumState(MQSM<max_qubits>::fromStatevector(sv));
  }

  /**
   * @brief Constructs a quantum state initialized to a specific basis state.
   *
   * @param basis_index The integer index representing the computational basis
   * state (e.g., 5 for |101>).
   * @param num_qubits The total number of active qubits.
   * @return The constructed MQSMQuantumState.
   *
   * @pre basis_index < (1 << num_qubits)
   * @pre num_qubits <= max_qubits
   * @post The state is initialized to the requested basis state.
   *
   * @par Example:
   * @code
   * auto state = QuantumState::fromBasisState(3, 2); // Creates state |11>
   * @endcode
   */
  static MQSMQuantumState fromBasisState(BasisKetIndex basis_index,
                                         QubitIndex num_qubits) {
    if (basis_index >= (1ull << num_qubits)) {
      throw std::invalid_argument(
          "MQSMQuantumState::fromBasisState: Basis index exceeds the maximum "
          "allowed limit for " +
          std::to_string(num_qubits) + " qubits.");
    }
    return MQSMQuantumState(
        MQSM<max_qubits>::fromBasis(basis_index, num_qubits));
  }

public: // ------ INTERFACE ----------
  /**
   * @brief Retrieves the number of active qubits in the state.
   * @return The number of qubits.
   */
  QubitIndex getNumQubits() const noexcept { return mqsm.getNumQubits(); }

  /**
   * @brief Checks if the quantum state is empty (uninitialized).
   * @return True if empty, false otherwise.
   */
  bool empty() const noexcept { return mqsm.empty(); }

  /**
   * @brief Computes the tensor product of this state with another state
   * in-place.
   *
   * @param other The other quantum state to tensor product with.
   * @return A reference to this modified state.
   *
   * @pre this->getNumQubits() + other.getNumQubits() <= max_qubits
   * @post This state becomes the tensor product of the two original states.
   *
   * @par Example:
   * @code
   * state1 ^= state2;
   * @endcode
   */
  MQSMQuantumState &operator^=(const MQSMQuantumState &other) {
    mqsm ^= other.mqsm;
    return *this;
  }

  /**
   * @brief Computes the tensor product of this state with another state.
   *
   * @param other The other quantum state to tensor product with.
   * @return A new MQSMQuantumState representing the combined state.
   *
   * @pre this->getNumQubits() + other.getNumQubits() <= max_qubits
   * @post A new state representing the tensor product is returned.
   *
   * @par Example:
   * @code
   * auto combined_state = state1 ^ state2;
   * @endcode
   */
  MQSMQuantumState operator^(const MQSMQuantumState &other) const {
    return MQSMQuantumState(mqsm ^ other.mqsm);
  }

  /**
   * @brief Calculates the exact probabilities of measuring specific qubits.
   *
   * @param qubits A span of qubit indices to measure. If empty, all qubits are
   * measured.
   * @return A SampleProbabilities object mapping measurement outcomes to their
   * probabilities.
   *
   * @pre All specified qubits must be valid (0 <= q < getNumQubits()).
   * @post The original state remains unchanged.
   *
   * @par Example:
   * @code
   * std::vector<QubitIndex> to_measure = {0, 2};
   * auto probs = state.sampleProbabilities(to_measure);
   * @endcode
   */
  SampleProbabilities
  sampleProbabilities(span<const QubitIndex> qubits = {}) const {
    return mqsm.sampleProbabilities(qubits);
  }

  /**
   * @brief Calculates the exact probabilities of measuring specific qubits
   * (initializer list overload).
   *
   * @param qubits An initializer list of qubit indices to measure.
   * @return A SampleProbabilities object mapping measurement outcomes to their
   * probabilities.
   */
  SampleProbabilities
  measureProbabilities(initializer_list<QubitIndex> qubits) const {
    return mqsm.measureProbabilities(qubits);
  }

  /**
   * @brief Simulates measurement by sampling the state a given number of times.
   *
   * @param shots The number of times to sample the state.
   * @param qubits A span of qubit indices to sample. If empty, samples all
   * qubits.
   * @return A SampleCount object containing a histogram of the measurement
   * outcomes.
   *
   * @pre shots > 0
   * @pre All specified qubits must be valid (0 <= q < getNumQubits()).
   * @post The original state remains unchanged.
   *
   * @par Example:
   * @code
   * auto counts = state.sampleCounts(1024);
   * @endcode
   */
  SampleCount sampleCounts(size_t shots,
                           span<const QubitIndex> qubits = {}) const {
    return mqsm.sample(shots, qubits);
  }

  /**
   * @brief Simulates measurement by sampling specific qubits (initializer list
   * overload).
   *
   * @param shots The number of times to sample the state.
   * @param qubits An initializer list of qubit indices to sample.
   * @return A SampleCount object containing a histogram of the measurement
   * outcomes.
   */
  SampleCount sampleCounts(size_t shots,
                           initializer_list<QubitIndex> qubits) const {
    return mqsm.sample(shots, qubits);
  }

private: // ---------- DEBUG ----------
  /**
   * @brief Exposes the underlying MQSM structure.
   * @return A copy of the internal MQSM.
   */
  MQSM<max_qubits> getMQSM() const { return mqsm; }

public: // ---------- I/O ----------
  /**
   * @brief Returns a string representation of the state (as a statevector
   * string).
   * @return A string showing the amplitudes of all basis states.
   */
  string toString() const { return mqsm.toStatevectorString(); }



  /**
   * @brief Returns a string representation of the underlying Finite State Machine representing the quantum state.
   * @return A string describing the structure of the FSM.
   */
  string toFSMString() const { return mqsm.toString(); }


  /**
   * @brief Converts the MQSM state into a full standard complex statevector.
   *
   * @return A std::vector of std::complex<double> representing the statevector
   * amplitudes.
   *
   * @pre getNumQubits() should be reasonably small to avoid memory overflow
   * (O(2^n)).
   * @post A full 2^N length vector is returned.
   */
  vector<complex<double>> toStatevector() const {
    return mqsm.toStdStatevector();
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const MQSMQuantumState &state) {
    return (os << state.toString());
  }

  friend class MQGTSim;

  friend class UnitaryOp<max_qubits>;
};

typedef MQSMQuantumState<MAX_QUBITS_ALLOWED> QuantumState;
typedef MQSMQuantumState<64> QuantumState64;
typedef MQSMQuantumState<32> QuantumState32;
typedef MQSMQuantumState<16> QuantumState16;
typedef MQSMQuantumState<8> QuantumState8;

#endif