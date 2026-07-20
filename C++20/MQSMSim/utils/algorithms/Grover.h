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
#ifndef GROVER_ALG_H
#define GROVER_ALG_H



#include <QCirc.h>
#include <numbers>
#include <cmath>
#include <cstdint>
#include <private/BasicDataTypes.h>

using namespace std;


/**
 * @brief Calculates the optimal number of iterations for Grover's algorithm.
 * @param num_qubits The number of qubits in the search space.
 * @param num_solutions The expected number of valid solutions (default 1).
 * @return The optimal integer number of Grover iterations to perform.
 */
inline size_t GroverOptimalIterations(QubitIndex num_qubits, size_t num_solutions= 1) {
    if (num_solutions == 0 || num_qubits <= 1) return 0;
    long double ratio = static_cast<long double>(1ull << num_qubits) / num_solutions;
    long double pi= std::numbers::pi_v<long double>;
    return static_cast<size_t>(std::floor((pi / 4) * std::sqrt(ratio)));
}

/**
 * @brief Constructs a standard Grover diffusion operator circuit.
 * @tparam max_qubits Maximum capacity of the generated circuit.
 * @param num_qubits Number of qubits the diffuser acts on.
 * @return A QCircuit implementing the Grover diffuser.
 */
template<unsigned int max_qubits>
QCircuit<max_qubits> templateGroverDiffuser(QubitIndex num_qubits) {

    // Máscaras para aplicar H y X a todos los qubits excepto el último
    vector<QubitIndex> all_but_last_qubits(num_qubits - 1);
    std::iota(all_but_last_qubits.begin(), all_but_last_qubits.end(), 0);

    QCircuit<max_qubits> qc(num_qubits);

    for (QubitIndex q= 0; q< num_qubits; ++q) {
        qc.h(q).x(q);
    }
    qc.cz(all_but_last_qubits, num_qubits - 1);
    for (QubitIndex q= 0; q< num_qubits; ++q) {
        qc.x(q).h(q);
    }

    return qc;
}

/**
 * @brief Constructs a single Grover iteration (Oracle + Diffuser).
 * @tparam max_qubits Maximum capacity of the generated circuit.
 * @param oracle The oracle circuit that marks the solution states.
 * @param diffuser_qubits Optional explicit set of qubits to apply the diffuser on. If empty, uses oracle.getNumQubits().
 * @return A QCircuit representing one full Grover iteration.
 */
template<unsigned int max_qubits>
QCircuit<max_qubits> templateGroverIteration(const QCircuit<max_qubits>& oracle,
                                         span<const QubitIndex> diffuser_qubits= {}) {
    

    if (max_qubits < diffuser_qubits.size()) {
        throw std::invalid_argument("GroverIteration: Oracle circuit has fewer qubits than the diffuser circuit");
    }
    if (oracle.getNumQubits() > max_qubits) {
        throw std::invalid_argument("GroverIteration: Oracle circuit exceeds max_qubits");
    }
    if (oracle.getNumQubits() < diffuser_qubits.size()) {
        throw std::invalid_argument("GroverIteration: Oracle circuit has fewer qubits than the diffuser qubits");
    }

    QCircuit<max_qubits> qc= oracle;

    qc.append(templateGroverDiffuser<max_qubits>(diffuser_qubits.size()), diffuser_qubits);
    return qc;
}

/**
 * @brief Constructs the full Grover's search algorithm circuit.
 *
 * This function automatically determines the optimal number of iterations and applies
 * the oracle and diffuser repeatedly.
 *
 * @tparam max_qubits Maximum capacity of the generated circuit.
 * @param oracle The oracle circuit marking the solution states.
 * @param expected_solutions Number of expected solutions in the search space (default 1).
 * @param diffuser_qubits Optional explicit set of qubits for the diffuser.
 * @return A QCircuit implementing the full Grover search.
 */
template<unsigned int max_qubits>
QCircuit<max_qubits> templateGroverAlgorithm(const QCircuit<max_qubits>& oracle, 
                                            unsigned int expected_solutions= 1,
                                            span<const QubitIndex> diffuser_qubits= {}) {


    if (oracle.getNumQubits() < diffuser_qubits.size()) {
        throw std::invalid_argument("GroverAlgorithm: Oracle circuit has fewer qubits than the diffuser circuit");
    }
    if (expected_solutions == 0) {
        throw std::invalid_argument("GroverAlgorithm: Expected solutions must be greater than zero");
    }
    if (oracle.getNumQubits() > max_qubits) {
        throw std::invalid_argument("GroverAlgorithm: Oracle circuit exceeds max_qubits");
    }
    if (oracle.getNumQubits() < diffuser_qubits.size()) {
        throw std::invalid_argument("GroverAlgorithm: Oracle circuit has fewer qubits than the diffuser qubits");
    }


    if (!diffuser_qubits.empty() && diffuser_qubits.size() > oracle.getNumQubits()) {
        throw std::invalid_argument("GroverAlgorithm: Diffuser qubits size (" + std::to_string(diffuser_qubits.size()) + ") is greater than the number of qubits in the oracle circuit (" + std::to_string(oracle.getNumQubits()) + ")");
    }

    int max_n_qubits= diffuser_qubits.empty() ? oracle.getNumQubits() : diffuser_qubits.size();
    size_t iterations= GroverOptimalIterations(max_n_qubits, expected_solutions);

    QCircuit<max_qubits> grover_iteration= templateGroverIteration<max_qubits>(oracle, diffuser_qubits);

    // Repeat Grover iteration the optimal number of times
    QCircuit<max_qubits> qc(grover_iteration.getNumQubits());
    for (uint32_t i = 0; i < iterations; ++i) {
        qc.append(grover_iteration);
    }
    return qc;
}





inline QCirc64 GroverDiffuser(QubitIndex num_qubits) { return templateGroverDiffuser<64>(num_qubits); }
inline QCirc64 GroverDiffuser64(QubitIndex num_qubits) { return templateGroverDiffuser<64>(num_qubits); }
inline QCirc32 GroverDiffuser32(QubitIndex num_qubits) { return templateGroverDiffuser<32>(num_qubits); }
inline QCirc16 GroverDiffuser16(QubitIndex num_qubits) { return templateGroverDiffuser<16>(num_qubits); }
inline QCirc8 GroverDiffuser8(QubitIndex num_qubits) { return templateGroverDiffuser<8>(num_qubits); }

inline QCirc64 GroverIteration(const QCirc64& oracle, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverIteration<64>(oracle, diffuser_qubits); }
inline QCirc64 GroverIteration64(const QCirc64& oracle, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverIteration<64>(oracle, diffuser_qubits); }
inline QCirc32 GroverIteration32(const QCirc32& oracle, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverIteration<32>(oracle, diffuser_qubits); }
inline QCirc16 GroverIteration16(const QCirc16& oracle,  span<const QubitIndex> diffuser_qubits= {}) { return templateGroverIteration<16>(oracle, diffuser_qubits); }
inline QCirc8 GroverIteration8(const QCirc8& oracle, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverIteration<8>(oracle, diffuser_qubits); }

inline QCirc64 GroverAlgorithm(const QCirc64& oracle, unsigned int expected_solutions= 1, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverAlgorithm<64>(oracle, expected_solutions, diffuser_qubits); }
inline QCirc64 GroverAlgorithm64(const QCirc64& oracle, unsigned int expected_solutions= 1, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverAlgorithm<64>(oracle, expected_solutions, diffuser_qubits); }
inline QCirc32 GroverAlgorithm32(const QCirc32& oracle, unsigned int expected_solutions= 1, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverAlgorithm<32>(oracle, expected_solutions, diffuser_qubits); }
inline QCirc16 GroverAlgorithm16(const QCirc16& oracle, unsigned int expected_solutions= 1, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverAlgorithm<16>(oracle, expected_solutions, diffuser_qubits); }
inline QCirc8 GroverAlgorithm8(const QCirc8& oracle, unsigned int expected_solutions= 1, span<const QubitIndex> diffuser_qubits= {}) { return templateGroverAlgorithm<8>(oracle, expected_solutions, diffuser_qubits); }








#endif