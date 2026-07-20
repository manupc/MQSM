#include <MQGTSim.h>
#include <MQSMQuantumState.h>
#include <QCirc.h>
#include <UnitaryOp.h>
#include <utils/algorithms/Grover.h>

#include <iostream>
using namespace std;

void printExampleHeader(const string &title) {
  cout << "\n=================================================================="
          "==============\n";
  cout << " >>> EXAMPLE: " << title << "\n";
  cout << "===================================================================="
          "============\n";
}

void printSubSection(const string &title) {
  cout << "\n--- " << title << " ---\n";
}

int main() {

  // ---------------------------------------------------------
  printExampleHeader("Hadamard Gate (H)");
  Unitary h = Unitary::H();

  printSubSection("FSM Representation");
  cout << h.toFSMString() << endl;

  printSubSection("Extracted Unitary Matrix");
  cout << h << endl;

  // ---------------------------------------------------------
  printExampleHeader("Controlled-NOT Gate (CNOT)");
  Unitary cnot = Unitary::CX();

  printSubSection("FSM Representation");
  cout << cnot.toFSMString() << endl;

  printSubSection("Extracted Unitary Matrix");
  cout << cnot << endl;

  // ---------------------------------------------------------
  printExampleHeader("Evolution of |+10> with CNOT(0,2)");
  QuantumState ket = QuantumState::fromLabel("+10");

  printSubSection("Initial State: |+10>");
  cout << ket << endl;
  cout << "\nMQSM Representation:\n" << ket.toFSMString() << endl;

  QCirc qc(3);
  qc.cx(0, 2);

  MQGTSim sim;
  QuantumState new_ket = sim.evolveMQSM(qc, ket);

  printSubSection("Final State after CNOT(0,2)");
  cout << new_ket << endl;
  cout << "\nMQSM Representation:\n" << new_ket.toFSMString() << endl;

  // ---------------------------------------------------------
  printExampleHeader("Evolution of (|01+>+|100>)/sqrt{2} with SWAP(0,2)");
  QCirc qc2(3);
  qc2.h(0).x(1).h(2);
  qc2.cx(0, 1).ch(0, 2);

  ket = sim.evolveMQSM(qc2);

  printSubSection("Initial State: (|01+>+|100>)/sqrt{2}");
  cout << ket << endl;
  cout << "\nMQSM Representation:\n" << ket.toFSMString() << endl;

  qc2.swap(0, 2);
  new_ket = sim.evolveMQSM(qc2);

  printSubSection("Final State after SWAP(0,2)");
  cout << new_ket << endl;
  cout << "\nMQSM Representation:\n" << new_ket.toFSMString() << endl;

  // ---------------------------------------------------------
  printExampleHeader("Phase Kickback");
  ket = QuantumState::fromLabel("0-");

  printSubSection("Initial State: |0->");
  cout << ket << endl;
  cout << "\nMQSM Representation:\n" << ket.toFSMString() << endl;

  QCirc qc3(2);
  qc3.h(0).cx(0, 1).h(0);

  new_ket = sim.evolveMQSM(qc3, ket);

  printSubSection("Final State after Phase Kickback");
  cout << new_ket << endl;
  cout << "\nMQSM Representation:\n" << new_ket.toFSMString() << endl;

  // ---------------------------------------------------------
  printExampleHeader("Grover's Algorithm for a Simple Oracle");
  QCirc oracle(2);
  oracle.cz(0, 1);

  QCirc diffuser = GroverDiffuser(2);

  QCirc grover_example(2);
  grover_example.h(0).h(1);
  grover_example.append(oracle);
  grover_example.append(diffuser);

  new_ket = sim.evolveMQSM(grover_example);

  printSubSection("Final State (Oracle marks |11>)");
  cout << new_ket << endl;
  cout << "\nMQSM Representation:\n" << new_ket.toFSMString() << endl;

  cout << "\n=================================================================="
          "==============\n";
  cout << "                             END OF EXAMPLES\n";
  cout << "===================================================================="
          "============\n\n";

  return 0;
}