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
#ifndef QFT_ALG_H
#define QFT_ALG_H

#include <QCirc.h>
#include <cmath>
#include <cstdint>
#include <numbers>
using namespace std;

/**
 * @brief Constructs the standard Quantum Fourier Transform (QFT) circuit.
 * @tparam max_qubits Maximum capacity of the generated circuit.
 * @param num_qubits Number of qubits the QFT acts on.
 * @param do_swaps Whether to append SWAP gates at the end to reverse the qubit order (default true).
 * @return A QCircuit implementing the QFT.
 */
template<unsigned int max_qubits>
QCircuit<max_qubits> templateQFT(QubitIndex num_qubits, bool do_swaps= true) {
 
    double pi= std::numbers::pi_v<double>;

    QCircuit<max_qubits> qc(num_qubits);    
    for (int i = num_qubits-1; i >= 0; --i) {
        qc.h(i);
        for (int j = i - 1; j >= 0; --j) {
            double angle = pi / (std::size_t(1) << (i - j));
            qc.cp(j, i, angle);
        }
    }
    // SWAP final
    if (do_swaps) {
        unsigned int half = num_qubits >> 1;
        for (unsigned int i = 0; i < half; ++i) {
            qc.swap(i, num_qubits - 1 - i);
        }
    }
    return qc;
}

/**
 * @brief Constructs the Inverse Quantum Fourier Transform (QFT) circuit.
 * @tparam max_qubits Maximum capacity of the generated circuit.
 * @param num_qubits Number of qubits the inverse QFT acts on.
 * @param do_swaps Whether to prepend SWAP gates to reverse the qubit order (default true).
 * @return A QCircuit implementing the Inverse QFT.
 */
template<unsigned int max_qubits>
QCircuit<max_qubits> templateQFTinv(QubitIndex num_qubits, bool do_swaps= true) {

    double pi= std::numbers::pi_v<double>;
    QCircuit<max_qubits> qc(num_qubits);
    

    // SWAP final
    if (do_swaps) {
        unsigned int half = num_qubits >> 1;
        for (int i = half-1; i >= 0; --i) {
            qc.swap(i, num_qubits-1-i);
        }
    }

    for (int i = 0; i < (int)num_qubits; ++i) {
        for (int j = 0; j < i; ++j) {
            double angle = pi / (size_t(1) << (i - j));
            qc.cp(j, i, -angle); // Control en j, target en i
        }
        qc.h(i);
    }

    return qc;

}




inline QCirc QFT(QubitIndex num_qubits, bool do_swaps= true) { return templateQFT<64>(num_qubits, do_swaps); }
inline QCirc64 QFT64(QubitIndex num_qubits, bool do_swaps= true) { return templateQFT<64>(num_qubits, do_swaps); }
inline QCirc32 QFT32(QubitIndex num_qubits, bool do_swaps= true) { return templateQFT<32>(num_qubits, do_swaps); }
inline QCirc16 QFT16(QubitIndex num_qubits, bool do_swaps= true) { return templateQFT<16>(num_qubits, do_swaps); }
inline QCirc8 QFT8(QubitIndex num_qubits, bool do_swaps= true) { return templateQFT<8>(num_qubits, do_swaps); }

inline QCirc QFTinv(QubitIndex num_qubits, bool do_swaps= true) { return templateQFTinv<64>(num_qubits, do_swaps); }
inline QCirc64 QFTinv64(QubitIndex num_qubits, bool do_swaps= true) { return templateQFTinv<64>(num_qubits, do_swaps); }
inline QCirc32 QFTinv32(QubitIndex num_qubits, bool do_swaps= true) { return templateQFTinv<32>(num_qubits, do_swaps); }
inline QCirc16 QFTinv16(QubitIndex num_qubits, bool do_swaps= true) { return templateQFTinv<16>(num_qubits, do_swaps); }
inline QCirc8 QFTinv8(QubitIndex num_qubits, bool do_swaps= true) { return templateQFTinv<8>(num_qubits, do_swaps); }





#endif