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
#ifndef __UNITARY_OP_H__
#define __UNITARY_OP_H__

#include <MQSMQuantumState.h>
#include <complex>
#include <private/MQGT.h>

using std::complex;
using std::string;

// Advanced declarations
class MQGTSim;

/**
 * @class UnitaryOp
 * @brief Represents a unitary quantum operation using Multi-Qubit Gate
 * Transducers (MQGT).
 *
 * This class abstracts the MQGT structure to provide a mathematical and
 * user-friendly interface for creating, composing, and applying unitary
 * operations to quantum states.
 *
 * @tparam max_qubits The maximum number of qubits this operation can support.
 */
template <unsigned int max_qubits>
  requires(max_qubits > 0 && max_qubits <= MAX_QUBITS_ALLOWED)
class UnitaryOp {

private: // ---------- REPRESENTATION ----------
  MQGT<max_qubits> mqgt;

  UnitaryOp(const MQGT<max_qubits> &mqgt) : mqgt(mqgt) {}
  UnitaryOp(MQGT<max_qubits> &&mqgt) : mqgt(std::move(mqgt)) {}

public: // ---------- CONSTRUCTORS ----------
  /**
   * @brief Default constructor creating an empty or identity unitary operation.
   */
  UnitaryOp() = default;

  /**
   * @brief Default destructor.
   */
  ~UnitaryOp() = default;

  UnitaryOp(const UnitaryOp &other) = default;
  UnitaryOp(UnitaryOp &&other) noexcept = default;
  UnitaryOp &operator=(const UnitaryOp &other) = default;
  UnitaryOp &operator=(UnitaryOp &&other) noexcept = default;

public: // ---------- INTERFACE ----------
  /**
   * @brief Gets the number of qubits the unitary operation acts on.
   * @return Number of active qubits.
   */
  size_t getNumQubits() const noexcept { return mqgt.getNumQubits(); }

  /**
   * @brief In-place tensor product with another unitary operation.
   *
   * @param other The other UnitaryOp to tensor with.
   * @return Reference to this modified operation.
   *
   * @pre this->getNumQubits() + other.getNumQubits() <= max_qubits
   * @post This operation becomes the tensor product of both.
   *
   * @par Example:
   * @code
   * op1 ^= op2;
   * @endcode
   */
  UnitaryOp &operator^=(const UnitaryOp &other) {
    mqgt ^= other.mqgt;
    return *this;
  }

  /**
   * @brief Tensor product with another unitary operation.
   *
   * @param other The other UnitaryOp to tensor with.
   * @return A new UnitaryOp representing the tensor product.
   *
   * @pre this->getNumQubits() + other.getNumQubits() <= max_qubits
   */
  UnitaryOp operator^(const UnitaryOp &other) const {
    UnitaryOp result = *this;
    result.tensoreq(other);
    return result;
  }

  /**
   * @brief Applies the unitary operation to a quantum state.
   *
   * @param state The MQSMQuantumState to evolve.
   * @return A new MQSMQuantumState representing the evolved state.
   *
   * @pre State must not be empty.
   * @pre Number of qubits in state and unitary must match.
   * @post The returned state is the result of U * |psi>.
   */
  MQSMQuantumState<max_qubits>
  operator*(const MQSMQuantumState<max_qubits> &state) const {
    if (state.empty()) {
      throw std::invalid_argument("UnitaryOp: Cannot apply unitary operation "
                                  "to an empty MQSMQuantumState.");
    }
    if (mqgt.getNumQubits() != state.getNumQubits()) {
      throw std::invalid_argument(
          "UnitaryOp: Number of qubits in the unitary operation does not match "
          "the number of qubits in the MQSMQuantumState.");
    }
    return MQSMQuantumState<max_qubits>(state.getMQSM().evolve(mqgt));
  }

  /**
   * @brief In-place matrix multiplication (composition) with another unitary
   * operation.
   *
   * @param other The other UnitaryOp to multiply by.
   * @return Reference to this modified operation.
   */
  UnitaryOp &operator*=(const UnitaryOp &other) {
    FSM<max_qubits, 4> new_fsm = other.mqgt.getFSM();
    new_fsm.compose(mqgt.getFSM());
    mqgt = MQGT<max_qubits>(std::move(new_fsm));
    return *this;
  }

  /**
   * @brief Matrix multiplication (composition) of two unitary operations.
   *
   * @param other The other UnitaryOp.
   * @return A new UnitaryOp representing the composed operations.
   */
  UnitaryOp operator*(const UnitaryOp &other) const {
    UnitaryOp result = *this;
    result *= other;
    return result;
  }

  /**
   * @brief Creates a controlled version of this unitary operation.
   *
   * @param num_controls The number of control qubits.
   * @param ctrl_states Optional string of '0's and '1's indicating active
   * states (default all '1's).
   * @return A new controlled UnitaryOp.
   */
  UnitaryOp controlled(QubitIndex num_controls = 1,
                       string_view ctrl_states = "") const {
    return UnitaryOp(mqgt.controlled(num_controls, ctrl_states));
  }

public: // ---------- FACTORY -------------
  /**
   * @name Standard Gate Unitaries
   * @brief Static factory methods returning standard unitary operations.
   */
  ///@{
  static UnitaryOp I(QubitIndex tensor_n = 1) {
    return UnitaryOp(MQGT<max_qubits>::I(tensor_n));
  }
  static UnitaryOp X() { return UnitaryOp(MQGT<max_qubits>::X()); }
  static UnitaryOp Y() { return UnitaryOp(MQGT<max_qubits>::Y()); }
  static UnitaryOp Z() { return UnitaryOp(MQGT<max_qubits>::Z()); }
  static UnitaryOp H() { return UnitaryOp(MQGT<max_qubits>::H()); }
  static UnitaryOp S() { return UnitaryOp(MQGT<max_qubits>::S()); }
  static UnitaryOp Sdg() { return UnitaryOp(MQGT<max_qubits>::Sdg()); }
  static UnitaryOp T() { return UnitaryOp(MQGT<max_qubits>::T()); }
  static UnitaryOp Tdg() { return UnitaryOp(MQGT<max_qubits>::Tdg()); }
  static UnitaryOp SX() { return UnitaryOp(MQGT<max_qubits>::SX()); }
  static UnitaryOp SXdg() { return UnitaryOp(MQGT<max_qubits>::SXdg()); }
  static UnitaryOp P(const Real &theta) {
    return UnitaryOp(MQGT<max_qubits>::P(theta));
  }
  static UnitaryOp RZ(const Real &theta) {
    return UnitaryOp(MQGT<max_qubits>::Rz(theta));
  }
  static UnitaryOp RY(const Real &theta) {
    return UnitaryOp(MQGT<max_qubits>::Ry(theta));
  }
  static UnitaryOp RX(const Real &theta) {
    return UnitaryOp(MQGT<max_qubits>::Rx(theta));
  }

  static UnitaryOp SWAP() { return UnitaryOp(MQGT<max_qubits>::SWAP()); }
  static UnitaryOp RZZ(const Real &theta) {
    return UnitaryOp(MQGT<max_qubits>::RZZ(theta));
  }
  static UnitaryOp RXX(const Real &theta) {
    return UnitaryOp(MQGT<max_qubits>::RXX(theta));
  }
  static UnitaryOp RYY(const Real &theta) {
    return UnitaryOp(MQGT<max_qubits>::RYY(theta));
  }
  ///@}

  /**
   * @name Controlled Gate Unitaries
   * @brief Static factory methods returning standard controlled unitary
   * operations.
   */
  ///@{
  static UnitaryOp CX(QubitIndex num_controls = 1,
                      string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CX(num_controls, ctrl_states));
  }
  static UnitaryOp CY(QubitIndex num_controls = 1,
                      string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CY(num_controls, ctrl_states));
  }
  static UnitaryOp CZ(QubitIndex num_controls = 1,
                      string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CZ(num_controls, ctrl_states));
  }
  static UnitaryOp CH(QubitIndex num_controls = 1,
                      string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CH(num_controls, ctrl_states));
  }
  static UnitaryOp CS(QubitIndex num_controls = 1,
                      string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CS(num_controls, ctrl_states));
  }
  static UnitaryOp CSdg(QubitIndex num_controls = 1,
                        string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CSdg(num_controls, ctrl_states));
  }
  static UnitaryOp CT(QubitIndex num_controls = 1,
                      string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CT(num_controls, ctrl_states));
  }
  static UnitaryOp CTdg(QubitIndex num_controls = 1,
                        string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CTdg(num_controls, ctrl_states));
  }
  static UnitaryOp CSX(QubitIndex num_controls = 1,
                       string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CSX(num_controls, ctrl_states));
  }
  static UnitaryOp CSXdg(QubitIndex num_controls = 1,
                         string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CSXdg(num_controls, ctrl_states));
  }
  static UnitaryOp CP(const Real &theta, QubitIndex num_controls = 1,
                      string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CP(num_controls, theta, ctrl_states));
  }
  static UnitaryOp CRX(const Real &theta, QubitIndex num_controls = 1,
                       string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CRx(num_controls, theta, ctrl_states));
  }
  static UnitaryOp CRY(const Real &theta, QubitIndex num_controls = 1,
                       string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CRy(num_controls, theta, ctrl_states));
  }
  static UnitaryOp CRZ(const Real &theta, QubitIndex num_controls = 1,
                       string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CRz(num_controls, theta, ctrl_states));
  }

  static UnitaryOp CSWAP(QubitIndex num_controls = 1,
                         string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CSWAP(num_controls, ctrl_states));
  }
  static UnitaryOp CRZZ(const Real &theta, QubitIndex num_controls = 1,
                        string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CRzz(num_controls, theta, ctrl_states));
  }
  static UnitaryOp CRXX(const Real &theta, QubitIndex num_controls = 1,
                        string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CRxx(num_controls, theta, ctrl_states));
  }
  static UnitaryOp CRYY(const Real &theta, QubitIndex num_controls = 1,
                        string_view ctrl_states = "") {
    return UnitaryOp(MQGT<max_qubits>::CRyy(num_controls, theta, ctrl_states));
  }
  ///@}

private: // ---------- DEBUG -----------
  /**
   * @brief Exposes the underlying MQGT structure.
   * @return A copy of the internal MQGT.
   */
  MQGT<max_qubits> getMQGT() const { return mqgt; }

public: // ---------- I/O ----------
  /**
   * @brief Returns an ASCII string representing the sparse matrix form of the
   * unitary.
   */
  string toString() const { return mqgt.toMatrixString(); }


  /**
   * @brief Returns a string representation of the underlying Finite State Machine representing the operation.
   * @return A string describing the structure of the FSM.
   */
  string toFSMString() const { return mqgt.toString(); }



  /**
   * @brief Converts the unitary operation to a full standard dense matrix.
   *
   * @return A flattened std::vector representing the matrix in row-major order.
   *
   * @pre Matrix size (2^N x 2^N) should fit in memory.
   */
  vector<std::complex<double>> toMatrix() const { return mqgt.toStdMatrix(); }

  friend std::ostream &operator<<(std::ostream &os, const UnitaryOp &op) {
    return (os << op.toString());
  }

public: // ---------- FRIENDS ----------
  friend class MQGTSim;
};

typedef UnitaryOp<MAX_QUBITS_ALLOWED> Unitary;
typedef UnitaryOp<64> Unitary64;
typedef UnitaryOp<32> Unitary32;
typedef UnitaryOp<16> Unitary16;
typedef UnitaryOp<8> Unitary8;

#endif // __UNITARY_OP_H__