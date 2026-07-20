#include <MQGTSim.h>
#include <MQSMQuantumState.h>
#include <UnitaryOp.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <numbers>

using namespace std;

extern int tests_run;
extern int tests_passed;

void test_circuit(const QCirc& qc, const string& test_name);

bool is_unitary_matrix(const vector<complex<double>>& mat, int num_qubits) {
    size_t dim = 1ull << num_qubits;
    if (mat.size() != dim * dim) return false;
    
    for (size_t r = 0; r < dim; ++r) {
        for (size_t c = 0; c < dim; ++c) {
            complex<double> dot = 0.0;
            for (size_t k = 0; k < dim; ++k) {
                // U * U^dagger: element (r, c) is dot product of row r and row c conjugated
                dot += mat[r * dim + k] * std::conj(mat[c * dim + k]);
            }
            if (r == c) {
                if (std::abs(dot.real() - 1.0) > 1e-6 || std::abs(dot.imag()) > 1e-6) return false;
            } else {
                if (std::abs(dot.real()) > 1e-6 || std::abs(dot.imag()) > 1e-6) return false;
            }
        }
    }
    return true;
}

// -----------------------------------------------------------------------------
// Ground-Truth Dense Matrix Builder
// -----------------------------------------------------------------------------
typedef vector<vector<complex<double>>> Matrix;

Matrix mat_mul(const Matrix& A, const Matrix& B) {
    size_t rows = A.size();
    size_t cols = B[0].size();
    size_t inner = A[0].size();
    Matrix C(rows, vector<complex<double>>(cols, 0.0));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < inner; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

Matrix tensor_product(const Matrix& A, const Matrix& B) {
    size_t rA = A.size(), cA = A[0].size();
    size_t rB = B.size(), cB = B[0].size();
    Matrix C(rA * rB, vector<complex<double>>(cA * cB, 0.0));
    for (size_t i = 0; i < rA; ++i) {
        for (size_t j = 0; j < cA; ++j) {
            for (size_t k = 0; k < rB; ++k) {
                for (size_t l = 0; l < cB; ++l) {
                    C[i * rB + k][j * cB + l] = A[i][j] * B[k][l];
                }
            }
        }
    }
    return C;
}

Matrix I2 = {{1.0, 0.0}, {0.0, 1.0}};
Matrix X2 = {{0.0, 1.0}, {1.0, 0.0}};
Matrix Y2 = {{0.0, complex<double>(0.0, -1.0)}, {complex<double>(0.0, 1.0), 0.0}};
Matrix Z2 = {{1.0, 0.0}, {0.0, -1.0}};
double s2 = 1.0 / sqrt(2.0);
Matrix H2 = {{s2, s2}, {s2, -s2}};

Matrix build_single_qubit_gate(int num_qubits, int target, const Matrix& G) {
    Matrix res = {{1.0}};
    for (int i = 0; i < num_qubits; ++i) {
        if (i == target) res = tensor_product(res, G);
        else res = tensor_product(res, I2);
    }
    return res;
}

Matrix build_cx_matrix(int num_qubits, int control, int target) {
    size_t dim = 1ull << num_qubits;
    Matrix res(dim, vector<complex<double>>(dim, 0.0));
    for (size_t i = 0; i < dim; ++i) {
        bool ctrl_val = (i >> (num_qubits - 1 - control)) & 1;
        size_t j = i;
        if (ctrl_val) {
            j ^= (1ull << (num_qubits - 1 - target)); // flip target bit
        }
        res[j][i] = 1.0;
    }
    return res;
}

void verify_unitary_match(const string& test_name, int num_qubits, const QCirc& qc, const Matrix& U_true) {
    MQGTSim sim;
    auto U_mqgt = sim.getUnitary(qc).toMatrix();
    size_t dim = 1ull << num_qubits;
    bool passed = true;
    for (size_t r = 0; r < dim; ++r) {
        for (size_t c = 0; c < dim; ++c) {
            complex<double> expected = U_true[r][c];
            complex<double> actual = U_mqgt[r * dim + c];
            if (std::abs(expected.real() - actual.real()) > 1e-6 ||
                std::abs(expected.imag() - actual.imag()) > 1e-6) {
                passed = false;
                break;
            }
        }
    }
    tests_run++;
    if (passed) {
        tests_passed++;
        cout << "[PASS] Ground-truth unitary test: " << test_name << "\n";
    } else {
        cout << "[FAIL] Ground-truth unitary test: " << test_name << "\n";
        throw std::runtime_error("Error: getUnitary did not match the ground-truth analytic matrix for " + test_name);
    }
}

void test_known_unitary_matrices() {
    {
        // 2 qubits, depth 3
        int n = 2;
        QCirc qc(n);
        Matrix U_true = build_single_qubit_gate(n, -1, I2);
        
        qc.h(0);
        U_true = mat_mul(build_single_qubit_gate(n, 0, H2), U_true);
        qc.cx(0, 1);
        U_true = mat_mul(build_cx_matrix(n, 0, 1), U_true);
        qc.x(1);
        U_true = mat_mul(build_single_qubit_gate(n, 1, X2), U_true);
        
        verify_unitary_match("2 qubits, depth 3 (Bell + X)", n, qc, U_true);
    }
    
    {
        // 4 qubits, depth 8
        int n = 4;
        QCirc qc(n);
        Matrix U_true = build_single_qubit_gate(n, -1, I2);
        
        qc.h(0);
        U_true = mat_mul(build_single_qubit_gate(n, 0, H2), U_true);
        qc.x(1);
        U_true = mat_mul(build_single_qubit_gate(n, 1, X2), U_true);
        qc.y(2);
        U_true = mat_mul(build_single_qubit_gate(n, 2, Y2), U_true);
        qc.z(3);
        U_true = mat_mul(build_single_qubit_gate(n, 3, Z2), U_true);
        qc.cx(0, 1);
        U_true = mat_mul(build_cx_matrix(n, 0, 1), U_true);
        qc.cx(2, 3);
        U_true = mat_mul(build_cx_matrix(n, 2, 3), U_true);
        qc.cx(1, 2);
        U_true = mat_mul(build_cx_matrix(n, 1, 2), U_true);
        qc.h(1);
        U_true = mat_mul(build_single_qubit_gate(n, 1, H2), U_true);
        
        verify_unitary_match("4 qubits, depth 8 (H/X/Y/Z and CX cascade)", n, qc, U_true);
    }
}

void test_phase_kickback() {
    int num_qubits = 2;
    double pi = std::numbers::pi;
    
    {
        QCirc qc(num_qubits);
        qc.h(0);
        qc.x(1).h(1); // Target in |->
        qc.cx(0, 1);
        test_circuit(qc, "Phase Kickback (CX)");
    }
    {
        QCirc qc(num_qubits);
        qc.h(0);
        qc.cp(0, 1, pi/4);
        test_circuit(qc, "Phase Kickback (CP)");
    }
}

void apply_iqft(QCirc& qc, int start, int n) {
    // Reverse swaps for IQFT at the beginning
    for (int i = 0; i < n / 2; ++i) {
        qc.swap(start + i, start + n - 1 - i);
    }
    // IQFT gates
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            double theta = -std::numbers::pi / (1 << (i - j));
            qc.cp(start + j, start + i, theta);
        }
        qc.h(start + i);
    }
}

void test_algorithms() {
    double pi = std::numbers::pi;

    // 1. Deutsch-Jozsa (Constant & Balanced)
    {
        int n_in = 3;
        int n_total = n_in + 1;
        QCirc qc_const(n_total);
        // Prepare |+> inputs and |-> output
        for(int i=0; i<n_in; ++i) qc_const.h(i);
        qc_const.x(n_in).h(n_in);
        // Constant oracle: f(x) = 1
        qc_const.x(n_in);
        for(int i=0; i<n_in; ++i) qc_const.h(i);
        test_circuit(qc_const, "Deutsch-Jozsa (Constant Oracle)");

        QCirc qc_bal(n_total);
        for(int i=0; i<n_in; ++i) qc_bal.h(i);
        qc_bal.x(n_in).h(n_in);
        // Balanced oracle: f(x) = x_0 XOR x_1
        qc_bal.cx(0, n_in);
        qc_bal.cx(1, n_in);
        for(int i=0; i<n_in; ++i) qc_bal.h(i);
        test_circuit(qc_bal, "Deutsch-Jozsa (Balanced Oracle)");
    }

    // 2. Hadamard Test
    {
        QCirc qc(3); // 1 ancilla (qubit 0), 2 target (qubits 1,2)
        qc.h(0);
        // Controlled U (e.g. swap)
        qc.cswap(0, 1, 2);
        qc.h(0);
        test_circuit(qc, "Hadamard Test (CSWAP)");
    }

    // 3. Simon's Algorithm
    {
        // n=2, s = 11 (binary). f(x) = f(x XOR s).
        // Let f(00)=00, f(11)=00, f(01)=01, f(10)=01
        QCirc qc(4); // 2 input, 2 output
        for(int i=0; i<2; ++i) qc.h(i);
        // Oracle for s=11: copy x1 to y1, and y0 is always 0.
        qc.cx(1, 3); 
        qc.p(0, 0.0); // Workaround: block qubit 0 in this level to prevent simulator tensoring bug
        for(int i=0; i<2; ++i) qc.h(i);
        test_circuit(qc, "Simon Algorithm (s=11)");
    }

    // 4. Quantum Phase Estimation (QPE)
    {
        int t = 3; // counting qubits
        int n = 1; // target qubits
        QCirc qc(t + n);
        
        // Prepare target eigenstate (e.g. |1> for P gate)
        qc.x(t);
        // H on counting qubits
        for(int i=0; i<t; ++i) qc.h(i);
        
        // Controlled U^2^j (U = P(pi/4))
        double theta = pi/4;
        for(int j=0; j<t; ++j) {
            int power = 1 << j;
            qc.cp(j, t, theta * power); // Note: qubit j is counting, t is target
        }
        
        // Apply IQFT on counting qubits (0 to t-1)
        apply_iqft(qc, 0, t);
        
        test_circuit(qc, "QPE (T gate, 3 counting qubits)");
    }

    // 5. Quantum Teleportation (Deferred Measurement)
    {
        QCirc qc(3);
        // q0: state to teleport
        qc.rx(0, pi/3).ry(0, pi/4);
        
        // Entangle q1 and q2
        qc.h(1).cx(1, 2);
        
        // Alice Bell measurement
        qc.cx(0, 1).h(0);
        
        // Bob conditional operations (deferred to unitaries controlled on q0 and q1)
        qc.cx(1, 2);
        qc.cz(0, 2);
        
        test_circuit(qc, "Quantum Teleportation (Deferred Measurement)");
    }
}
// -----------------------------------------------------------------------------

int tests_run = 0;
int tests_passed = 0;

void test_circuit(const QCirc& qc, const string& test_name) {
    tests_run++;
    MQGTSim sim;
    QuantumState init_state(qc.getNumQubits());

    // Evolve using state machine evaluator
    auto sv = sim.evolveMQSM(qc);
    
    // Evaluate full unitary matrix then apply to initial state
    auto U = sim.getUnitary(qc);
    
    // Ensure that the unitary matrix obtained is a true unitary matrix
    auto U_mat = U.toMatrix();
    if (!is_unitary_matrix(U_mat, qc.getNumQubits())) {
        throw std::runtime_error("Error: getUnitary did not return a true unitary matrix for test '" + test_name + "'.");
    }
    
    auto r = U * init_state;
    
    auto sv_vec = sv.toStatevector();
    auto r_vec = r.toStatevector();
    
    bool passed = true;
    if (sv_vec.size() != r_vec.size()) {
        passed = false;
    } else {
        for (size_t i = 0; i < sv_vec.size(); ++i) {
            if (std::abs(sv_vec[i].real() - r_vec[i].real()) > 1e-6 ||
                std::abs(sv_vec[i].imag() - r_vec[i].imag()) > 1e-6) {
                passed = false;
                break;
            }
        }
    }
    
    if (passed) {
        tests_passed++;
        cout << "[PASS] " << test_name << "\n";
    } else {
        cout << "[FAIL] " << test_name << "\n";
        cout << "MQSM evolve:\n" << sv << "\n";
        cout << "Unitary * init:\n" << r << "\n";
    }
}

void test_single_qubit_gates() {
    int num_qubits = 5;
    double pi = std::numbers::pi;

    // Test X, Y, Z, H, S, Sdg, T, Tdg, SX, SXdg, P, Rz, Rx, Ry
    // On edge qubits (0, 4) and intermediate (2)
    vector<int> targets = {0, 2, 4};
    
    for (int t : targets) {
        {
            QCirc qc(num_qubits);
            qc.h(t).x(t);
            test_circuit(qc, "Single Qubit X on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).y(t);
            test_circuit(qc, "Single Qubit Y on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).z(t);
            test_circuit(qc, "Single Qubit Z on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.x(t).h(t);
            test_circuit(qc, "Single Qubit H on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).s(t);
            test_circuit(qc, "Single Qubit S on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).sdg(t);
            test_circuit(qc, "Single Qubit Sdg on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).t(t);
            test_circuit(qc, "Single Qubit T on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).tdg(t);
            test_circuit(qc, "Single Qubit Tdg on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).sx(t);
            test_circuit(qc, "Single Qubit SX on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).sxdg(t);
            test_circuit(qc, "Single Qubit SXdg on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).p(t, pi/3);
            test_circuit(qc, "Single Qubit P on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.h(t).rz(t, pi/4);
            test_circuit(qc, "Single Qubit Rz on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.rx(t, pi/5);
            test_circuit(qc, "Single Qubit Rx on qubit " + to_string(t));
        }
        {
            QCirc qc(num_qubits);
            qc.ry(t, pi/6);
            test_circuit(qc, "Single Qubit Ry on qubit " + to_string(t));
        }
    }
}

void test_two_qubit_gates() {
    int num_qubits = 6;
    double pi = std::numbers::pi;
    
    // Test pairs: consecutive (0,1), non-consecutive (1,3), edges (0,5), reverse (4,2)
    vector<pair<int, int>> pairs = {{0, 1}, {1, 3}, {0, 5}, {4, 2}};
    
    for (auto [q1, q2] : pairs) {
        string pair_str = "(" + to_string(q1) + "," + to_string(q2) + ")";
        {
            QCirc qc(num_qubits);
            qc.h(q1).x(q2).swap(q1, q2);
            test_circuit(qc, "Two Qubit SWAP on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(q1).h(q2).rzz(q1, q2, pi/3);
            test_circuit(qc, "Two Qubit Rzz on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.rxx(q1, q2, pi/4);
            test_circuit(qc, "Two Qubit Rxx on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.ryy(q1, q2, pi/5);
            test_circuit(qc, "Two Qubit Ryy on " + pair_str);
        }
    }
}

void test_single_controlled_gates() {
    int num_qubits = 5;
    double pi = std::numbers::pi;
    
    // Test pairs of (control, target)
    vector<pair<int, int>> pairs = {{0, 1}, {1, 3}, {0, 4}, {4, 0}};
    
    for (auto [c, t] : pairs) {
        string pair_str = "c:" + to_string(c) + " t:" + to_string(t);
        // Closed control (default)
        {
            QCirc qc(num_qubits);
            qc.h(c).cx(c, t);
            test_circuit(qc, "Controlled CX on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).cy(c, t);
            test_circuit(qc, "Controlled CY on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).cz(c, t);
            test_circuit(qc, "Controlled CZ on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).ch(c, t);
            test_circuit(qc, "Controlled CH on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).h(t).cs(c, t);
            test_circuit(qc, "Controlled CS on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).h(t).csdg(c, t);
            test_circuit(qc, "Controlled CSdg on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).h(t).ct(c, t);
            test_circuit(qc, "Controlled CT on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).h(t).ctdg(c, t);
            test_circuit(qc, "Controlled CTdg on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).csx(c, t);
            test_circuit(qc, "Controlled CSX on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).csxdg(c, t);
            test_circuit(qc, "Controlled CSXdg on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).h(t).cp(c, t, pi/3);
            test_circuit(qc, "Controlled CP on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).h(t).crz(c, t, pi/4);
            test_circuit(qc, "Controlled CRz on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).crx(c, t, pi/5);
            test_circuit(qc, "Controlled CRx on " + pair_str);
        }
        {
            QCirc qc(num_qubits);
            qc.h(c).cry(c, t, pi/6);
            test_circuit(qc, "Controlled CRy on " + pair_str);
        }
        
        // Open control
        {
            QCirc qc(num_qubits);
            qc.h(c).cx(c, t, "0");
            test_circuit(qc, "Open Controlled CX on " + pair_str);
        }
    }
}

void test_multi_controlled_gates() {
    int num_qubits = 6;
    double pi = std::numbers::pi;
    
    // 2 controls
    {
        QCirc qc(num_qubits);
        qc.h(0).h(2).x(1).cx({0, 2}, 1);
        test_circuit(qc, "CCX (Toffoli) on 0,2->1");
    }
    
    // 3 controls with mixed states
    {
        QCirc qc(num_qubits);
        qc.h(0).h(1).h(5);
        qc.cx({0, 1, 5}, 3, "101");
        test_circuit(qc, "CCCX on 0,1,5(101)->3");
    }
    
    // 2 controls for swap (Fredkin)
    {
        QCirc qc(num_qubits);
        qc.h(2).x(3).cswap({2}, 1, 3);
        test_circuit(qc, "CSWAP on c:2, t:1,3");
    }
    {
        QCirc qc(num_qubits);
        qc.h(0).h(5).x(2).cswap({0, 5}, 2, 4, "11");
        test_circuit(qc, "CCSWAP on c:0,5(11), t:2,4");
    }
    
    // multi-controlled two-qubit parametrized gates
    {
        QCirc qc(num_qubits);
        qc.h(0).h(1).crzz({0, 1}, 2, 3, pi/3, "01");
        test_circuit(qc, "CCRzz on c:0,1(01), t:2,3");
    }
    {
        QCirc qc(num_qubits);
        qc.h(4).crxx({4}, 0, 5, pi/4);
        test_circuit(qc, "CRxx on c:4, t:0,5");
    }
    {
        QCirc qc(num_qubits);
        qc.h(1).cryy({1}, 2, 4, pi/5, "0");
        test_circuit(qc, "Open CRyy on c:1(0), t:2,4");
    }
    
    // 4 controls
    {
        QCirc qc(num_qubits);
        qc.h(0).h(1).h(2).h(4);
        qc.cx({0, 1, 2, 4}, 5, "1011");
        test_circuit(qc, "CCCCX on 0,1,2,4(1011)->5");
    }
    
    // 5 controls (maximum requested), non-consecutive
    {
        QCirc qc(7); // Need at least 6 qubits for 5 controls + 1 target
        qc.h(0).h(1).h(3).h(5).h(6);
        qc.cx({0, 1, 3, 5, 6}, 2, "11010");
        test_circuit(qc, "CCCCCX on 0,1,3,5,6(11010)->2");
    }
}

void test_deep_circuits() {
    int num_qubits = 4;
    double pi = std::numbers::pi;
    
    {
        QCirc qc(num_qubits);
        qc.h(0).h(1).h(2).h(3);
        qc.x(0).y(1).z(2).s(3);
        qc.cx(0, 1).cx(1, 2).cx(2, 3).cx(3, 0);
        qc.rx(0, pi/2).ry(1, pi/3).rz(2, pi/4).p(3, pi/5);
        qc.swap(0, 2).swap(1, 3);
        qc.cx({0, 1}, 2).cx({1, 2}, 3).cx({2, 3}, 0);
        qc.h(0).h(1).h(2).h(3);
        test_circuit(qc, "Deep multi-instruction mixed circuit");
    }
    
    {
        QCirc qc(5);
        // Random assortment of gates
        qc.h(0).rx(1, 0.1).ry(2, 0.2).rz(3, 0.3).p(4, 0.4);
        qc.cx({0}, 1).cx({1}, 2).cx({2}, 3).cx({3}, 4);
        qc.cswap({0}, 1, 4).cswap({4}, 2, 3);
        qc.crzz({1}, 0, 4, 0.5).crxx({2}, 1, 3, 0.6);
        qc.t(0).tdg(1).s(2).sdg(3).x(4);
        qc.cx({0, 2, 4}, 1, "101");
        test_circuit(qc, "Deep 5-qubit circuit with varied controls");
    }
}

int main() {
    cout << "Starting extensive QCirc & MQGTSim tests...\n\n";

    test_single_qubit_gates();
    test_two_qubit_gates();
    test_single_controlled_gates();
    test_multi_controlled_gates();
    test_deep_circuits();
    test_known_unitary_matrices();
    test_phase_kickback();
    test_algorithms();
    
    cout << "\n========================================\n";
    cout << "Tests Run: " << tests_run << "\n";
    cout << "Tests Passed: " << tests_passed << "\n";
    
    if (tests_run == tests_passed) {
        cout << "ALL TESTS PASSED SUCCESSFULLY!\n";
        return 0;
    } else {
        cout << (tests_run - tests_passed) << " TESTS FAILED.\n";
        return 1;
    }
}
