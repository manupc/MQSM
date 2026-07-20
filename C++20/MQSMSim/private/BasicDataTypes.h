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
#ifndef __BASICDATATYPES__H__
#define __BASICDATATYPES__H__

#include <cstdint>
#include <limits>

#include <private/Numbers.h>

// ------- BIT / QUBIT INDEXING TYPES -------

/**
 * @def MAX_QUBITS_ALLOWED
 * @brief The maximum number of qubits the simulator can handle natively using bitmasks.
 */
const int MAX_QUBITS_ALLOWED = 64; 

/**
 * @typedef BitIndex
 * @brief Type for indexing bits within a BitMask.
 */
typedef uint8_t BitIndex; 

/**
 * @typedef QubitIndex
 * @brief Type for indexing qubits in a quantum circuit or state.
 */
typedef uint8_t QubitIndex; 


/**
 * @typedef Complex
 * @brief Custom complex literal wrapper from the Numbers namespace.
 */
using Complex = Numbers::ComplexLiteral;

/**
 * @typedef Real
 * @brief Custom real literal wrapper from the Numbers namespace.
 */
using Real = Numbers::RealLiteral;





#endif