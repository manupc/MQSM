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
#ifndef __MYNUMBERS__H__
#define __MYNUMBERS__H__

#include <iostream>
#include <cmath>
#include <complex>
#include <numbers>
#include <numeric>
#include <sstream>
#include <string>


namespace Numbers {

using BaseReal= double;
const int decimals= 3;
const double ErrorTolerance= 1e-12;


// Advanced declarations
class RealLiteral;
class ComplexLiteral;


/**
 * @class RealLiteral
 * @brief A wrapper class for real numbers with built-in tolerance for floating-point comparisons.
 *
 * This class ensures that arithmetic operations and comparisons robustly handle 
 * small precision errors common in quantum simulations.
 */
class RealLiteral {
private: // -------------------------- REPRESENTATION -------------------------
    BaseReal value;

public: // -------------------------- CONSTRUCTORS -------------------------
    constexpr RealLiteral() : value(0.0) {}
    ~RealLiteral() = default;
    constexpr RealLiteral(BaseReal real) : value(real) {}

    // Copy & Move
    RealLiteral(const RealLiteral& other)= default;
    RealLiteral(RealLiteral&& other) noexcept= default;
    RealLiteral& operator=(const RealLiteral& other)= default;
    RealLiteral& operator=(RealLiteral&& other) noexcept= default;

    static RealLiteral cos(RealLiteral angle) { return RealLiteral(std::cos(angle.value)); }
    static RealLiteral sin(RealLiteral angle) { return RealLiteral(std::sin(angle.value)); }

public: // -------------------------- ARITHMETICS --------------------------
    RealLiteral operator+(const RealLiteral& other) const { return RealLiteral(value + other.value); }
    RealLiteral operator-(const RealLiteral& other) const { return RealLiteral(value - other.value); }
    RealLiteral operator*(const RealLiteral& other) const { return RealLiteral(value * other.value); }
    RealLiteral operator/(const RealLiteral& other) const { return RealLiteral(value / other.value); }

    // Overloaded operators for scalar BaseReal
    RealLiteral operator+(BaseReal r) const { return RealLiteral(value + r); }
    RealLiteral operator-(BaseReal r) const { return RealLiteral(value - r); }
    RealLiteral operator*(BaseReal r) const { return RealLiteral(value * r); }
    RealLiteral operator/(BaseReal r) const { return RealLiteral(value / r); }

    // Compound assignment operators
    RealLiteral& operator+=(const RealLiteral& other) { value += other.value; return *this; }
    RealLiteral& operator-=(const RealLiteral& other) { value -= other.value; return *this; }
    RealLiteral& operator*=(const RealLiteral& other) { value *= other.value; return *this; }
    RealLiteral& operator/=(const RealLiteral& other) { value /= other.value; return *this; }

    // Compound assignment operators for scalar BaseReal
    RealLiteral& operator+=(BaseReal r) { value += r; return *this; }
    RealLiteral& operator-=(BaseReal r) { value -= r; return *this; }
    RealLiteral& operator*=(BaseReal r) { value *= r; return *this; }
    RealLiteral& operator/=(BaseReal r) { value /= r; return *this; }

    // Unary operators
    RealLiteral operator-() const { return RealLiteral(-value); }
    RealLiteral operator+() const { return *this; }

public: // -------------------------- MATH CONSTANTS --------------------------
    static constexpr RealLiteral one() { return RealLiteral(1.0); }
    static constexpr RealLiteral zero() { return RealLiteral(0.0); }
    static constexpr RealLiteral pi() { return RealLiteral(std::numbers::pi_v<BaseReal>); }
    static constexpr RealLiteral e() { return RealLiteral(std::numbers::e_v<BaseReal>); }
    static constexpr RealLiteral sqrt2() { return RealLiteral(std::numbers::sqrt2_v<BaseReal>); }
    static constexpr RealLiteral inv_sqrt2() { return RealLiteral(((BaseReal)1.0) / std::numbers::sqrt2_v<BaseReal>); }

public: // -------------------------- MATH FUNCTIONS --------------------------
    RealLiteral abs() const { return RealLiteral(std::abs(value)); }
    bool NearZero() const { return std::abs(value) < ErrorTolerance; }
    bool Near(const RealLiteral& other) const { return std::abs(value - other.value) < ErrorTolerance; }
    RealLiteral sqrt() const { return RealLiteral(std::sqrt(value)); }
    RealLiteral sqr() const { return RealLiteral(value * value); }
    RealLiteral exp() const { return RealLiteral(std::exp(value)); }
    RealLiteral logn() const { return RealLiteral(std::log(value)); }
    RealLiteral norm() const { return sqr(); } // abs()^2
    RealLiteral sin() const { return RealLiteral(std::sin(value)); }
    RealLiteral cos() const { return RealLiteral(std::cos(value)); }
    RealLiteral toleranceRound() const { 
        BaseReal rounded(std::round(value / ErrorTolerance) * ErrorTolerance);
        return RealLiteral(rounded); 
    }
public: // -------------------------- COMPARISON OPERATORS --------------------------
    bool operator==(const RealLiteral& other) const { 
        if (value == other.value) return true;
        return std::abs(value - other.value) < ErrorTolerance;
    }
    bool operator!=(const RealLiteral& other) const { return !(*this == other); }
    bool operator<(const RealLiteral& other) const { return value < other.value && !(*this == other); }
    bool operator<=(const RealLiteral& other) const { return value < other.value || *this == other; }
    bool operator>(const RealLiteral& other) const { return value > other.value && !(*this == other); }
    bool operator>=(const RealLiteral& other) const { return value > other.value || *this == other; }


public: // -------------------------- I/O --------------------------

    double toDouble() const { return static_cast<double>(value); }
    float toFloat() const { return static_cast<float>(value); }
    long double toLongDouble() const { return static_cast<long double>(value); }
    
    std::string toString() const {
        std::ostringstream os;
        os.precision(decimals);
        os << value;
        return os.str();
    }

    friend std::ostream & operator<<(std::ostream & os, const RealLiteral & r) { return (os << r.toString()); }

    friend class ComplexLiteral; // Allow ComplexLiteral to access private members for efficient construction
};


/**
 * @class ComplexLiteral
 * @brief A wrapper class for complex numbers with built-in tolerance for floating-point comparisons.
 *
 * It extends std::complex to support robust equality checks, polar/rectangular factory methods,
 * and standard quantum amplitude operations.
 */
class ComplexLiteral {
private: // -------------------------- REPRESENTATION -------------------------
    using BaseComplex= std::complex<BaseReal>;
    BaseComplex value;

public: // -------------------------- CONSTRUCTORS -------------------------
    constexpr ComplexLiteral() : value(BaseReal(0.0), BaseReal(0.0)) {}
    ~ComplexLiteral() = default;
    constexpr ComplexLiteral(RealLiteral real) : value(real.value, BaseReal(0.0)) {}
    constexpr ComplexLiteral(RealLiteral real, RealLiteral imag) : value(real.value, imag.value) {}
    constexpr ComplexLiteral(BaseReal real, BaseReal imag= BaseReal(0.0)) : value(real, imag) {}
    constexpr ComplexLiteral(const BaseComplex& c) : value(c) {}

    // Copy & Move
    ComplexLiteral(const ComplexLiteral& other)= default;
    ComplexLiteral(ComplexLiteral&& other) noexcept= default;
    ComplexLiteral& operator=(const ComplexLiteral& other)= default;
    ComplexLiteral& operator=(ComplexLiteral&& other) noexcept= default;

    // From polar coordinates
    static ComplexLiteral fromPolar(RealLiteral magnitude, RealLiteral angle) {
        return ComplexLiteral(std::polar(magnitude.value, angle.value));
    }

    static ComplexLiteral fromPolar(BaseReal magnitude, BaseReal angle) {
        return ComplexLiteral(std::polar(magnitude, angle));
    }

    static ComplexLiteral from_rectangular(RealLiteral real, RealLiteral imag) {
        return ComplexLiteral(real.value, imag.value);
    }

    static ComplexLiteral from_rectangular(BaseReal real, BaseReal imag) {
        return ComplexLiteral(real, imag);
    }

public: // -------------------------- MATH CONSTANTS --------------------------
    static constexpr ComplexLiteral one() { return ComplexLiteral(BaseReal(1.0)); }
    static constexpr ComplexLiteral zero() { return ComplexLiteral(BaseReal(0.0)); }
    static constexpr ComplexLiteral sqrt2() { return ComplexLiteral(std::numbers::sqrt2_v<BaseReal>); }
    static constexpr ComplexLiteral inv_sqrt2() { return ComplexLiteral(((BaseReal)1.0) / std::numbers::sqrt2_v<BaseReal>); }
    static constexpr ComplexLiteral i() { return ComplexLiteral(BaseReal(0.0), BaseReal(1.0)); }
    static constexpr ComplexLiteral minus_i() { return ComplexLiteral(BaseReal(0.0), BaseReal(-1.0)); }
    static constexpr ComplexLiteral i_inv_sqrt2() { return ComplexLiteral(BaseReal(0.0), ((BaseReal)1.0) / std::numbers::sqrt2_v<BaseReal>); }
    static constexpr ComplexLiteral minus_i_inv_sqrt2() { return ComplexLiteral(BaseReal(0.0), -((BaseReal)1.0) / std::numbers::sqrt2_v<BaseReal>); }

public: // -------------------------- INTERFACE --------------------------
    RealLiteral real() const { return value.real(); }
    RealLiteral imag() const { return value.imag(); }
    BaseComplex toStdComplex() const { return value; }
    bool hasReal() const { return std::abs(value.real()) >= ErrorTolerance; }
    bool hasImag() const { return std::abs(value.imag()) >= ErrorTolerance; }


public: // -------------------------- ARITHMETIC OPERATORS --------------------------
    // Overloaded operators for arithmetic
    ComplexLiteral operator+(const ComplexLiteral& other) const {return ComplexLiteral(value + other.value);}
    ComplexLiteral operator-(const ComplexLiteral& other) const {return ComplexLiteral(value - other.value);}
    ComplexLiteral operator*(const ComplexLiteral& other) const {return ComplexLiteral(value * other.value);}
    ComplexLiteral operator/(const ComplexLiteral& other) const {return ComplexLiteral(value / other.value);}

    // Overloaded operators for scalar RealLiteral
    ComplexLiteral operator+(RealLiteral r) const {return ComplexLiteral(value + r.value);}
    ComplexLiteral operator-(RealLiteral r) const {return ComplexLiteral(value - r.value);}
    ComplexLiteral operator*(RealLiteral r) const {return ComplexLiteral(value * r.value);}
    ComplexLiteral operator/(RealLiteral r) const {return ComplexLiteral(value / r.value);}

    // Overloaded operators for scalar BaseReal
    ComplexLiteral operator+(BaseReal r) const {return ComplexLiteral(value + r);}
    ComplexLiteral operator-(BaseReal r) const {return ComplexLiteral(value - r);}
    ComplexLiteral operator*(BaseReal r) const {return ComplexLiteral(value * r);}
    ComplexLiteral operator/(BaseReal r) const {return ComplexLiteral(value / r);}

    // Compound assignment operators
    ComplexLiteral& operator+=(const ComplexLiteral& other) { value += other.value; return *this; }
    ComplexLiteral& operator-=(const ComplexLiteral& other) { value -= other.value; return *this; }
    ComplexLiteral& operator*=(const ComplexLiteral& other) { value *= other.value; return *this; }
    ComplexLiteral& operator/=(const ComplexLiteral& other) { value /= other.value; return *this; }

    // Compound assignment operators for scalar RealLiteral
    ComplexLiteral& operator+=(RealLiteral r) { value += r.value; return *this; }
    ComplexLiteral& operator-=(RealLiteral r) { value -= r.value; return *this; }
    ComplexLiteral& operator*=(RealLiteral r) { value *= r.value; return *this; }
    ComplexLiteral& operator/=(RealLiteral r) { value /= r.value; return *this; }

    // Compound assignment operators for scalar BaseReal
    ComplexLiteral& operator+=(BaseReal r) { value += r; return *this; }
    ComplexLiteral& operator-=(BaseReal r) { value -= r; return *this; }
    ComplexLiteral& operator*=(BaseReal r) { value *= r; return *this; }
    ComplexLiteral& operator/=(BaseReal r) { value /= r; return *this; }

    // Right-operand friend operators for std::complex
    friend ComplexLiteral operator+(const BaseComplex& c, const ComplexLiteral& cl) { return ComplexLiteral(c) + cl; }
    friend ComplexLiteral operator-(const BaseComplex& c, const ComplexLiteral& cl) { return ComplexLiteral(c) - cl; }
    friend ComplexLiteral operator*(const BaseComplex& c, const ComplexLiteral& cl) { return ComplexLiteral(c) * cl; }
    friend ComplexLiteral operator/(const BaseComplex& c, const ComplexLiteral& cl) { return ComplexLiteral(c) / cl; }

    // Right-operand friend operators for scalar RealLiteral
    friend ComplexLiteral operator+(RealLiteral r, const ComplexLiteral& c) { return ComplexLiteral(r) + c; }
    friend ComplexLiteral operator-(RealLiteral r, const ComplexLiteral& c) { return ComplexLiteral(r) - c; }
    friend ComplexLiteral operator*(RealLiteral r, const ComplexLiteral& c) { return ComplexLiteral(r) * c; }
    friend ComplexLiteral operator/(RealLiteral r, const ComplexLiteral& c) { return ComplexLiteral(r) / c; }

    // Right-operand friend operators for scalar BaseReal
    friend ComplexLiteral operator+(double r, const ComplexLiteral& c) { return ComplexLiteral(r) + c; }
    friend ComplexLiteral operator-(double r, const ComplexLiteral& c) { return ComplexLiteral(r) - c; }
    friend ComplexLiteral operator*(double r, const ComplexLiteral& c) { return ComplexLiteral(r) * c; }
    friend ComplexLiteral operator/(double r, const ComplexLiteral& c) { return ComplexLiteral(r) / c; }


    // Unary operators
    ComplexLiteral operator-() const { return ComplexLiteral(-value); }
    ComplexLiteral operator+() const { return *this; }

public: // -------------------------- MATH FUNCTIONS --------------------------
    RealLiteral magnitude() const { return std::abs(value); }
    RealLiteral angle() const { return std::arg(value); }
    ComplexLiteral conjugate() const { return ComplexLiteral(std::conj(value)); }
    ComplexLiteral inverse() const { return ComplexLiteral(1.0 / value); }
    ComplexLiteral exp() const { return ComplexLiteral(std::exp(value)); }
    ComplexLiteral logn() const { return ComplexLiteral(std::log(value)); }
    RealLiteral norm() const { return std::norm(value); } // sqr(magnitude())
    ComplexLiteral sqrt() const { return ComplexLiteral(std::sqrt(value)); }
    ComplexLiteral sqr() const { return ComplexLiteral(value * value); }
    RealLiteral abs() const { return std::abs(value); } // sqrt(norm())= magnitude()
    bool NearZero() const { return std::abs(value) < ErrorTolerance; }
    bool Near(const ComplexLiteral& other) const { 
        return std::abs(value.real() - other.value.real()) < ErrorTolerance &&
                std::abs(value.imag() - other.value.imag()) < ErrorTolerance;
    }
    ComplexLiteral toleranceRound() const { 
        BaseComplex rounded(std::round(value.real() / ErrorTolerance) * ErrorTolerance,
                           std::round(value.imag() / ErrorTolerance) * ErrorTolerance);
        return ComplexLiteral(rounded); 
    }

public: // -------------------------- COMPARISON OPERATORS --------------------------
    bool operator==(const ComplexLiteral& other) const { 
        if (value == other.value) return true;
        BaseComplex diff = value - other.value;
        return std::abs(diff) < ErrorTolerance;
    }
    bool operator!=(const ComplexLiteral& other) const { return !(*this == other); }

public: // -------------------------- I/O --------------------------

    std::string toString() const {
        std::ostringstream os;
        bool has_real = value.real() != 0.0 && std::abs(value.real()) >= ErrorTolerance;
        bool has_imag = value.imag() != 0.0 && std::abs(value.imag()) >= ErrorTolerance;
        os.precision(decimals);
        if (has_real) {
            os << value.real();
        }
        if (has_imag) {
            if (value.imag() > 0 && has_real) {
                os << "+";
            }
            if (value.imag() < 0) {
                os << "-";
            }
            RealLiteral imag_abs = std::abs(value.imag());
            if (imag_abs != 1) {
                os << imag_abs;
            }
            os << "i";
            return os.str();
        }
        
        if (!has_real && !has_imag) {
            os << 0.0;
        }

        return os.str();
    }

    friend std::ostream & operator<<(std::ostream & os, const ComplexLiteral & c) { return (os << c.toString()); }
};





};

#endif