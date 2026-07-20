# MQSM

Source code for the paper: A theoretical framework for compact and efficient quantum state simulation using Finite State Machines

- *Python folder*: Basic/naive very short/compact implementation of MQGT and MQSM, with the examples of the paper.

**HOW TO USE IT: ** Just open the Jupyter Notebook and run.

- *C++20 Folder*: Optimized C++20 implementation using appropiate data structures to represent Finite State Machines, MQSMs and MQGTs. It also includes an interface to ease the creation of quantum circuits.

**HOW TO USE IT: ** It requires compiling the examples using a C++20 compatible g++ compiler. Makefile script may use linux-only bash commands for running. If you have a C++20 compiant g++ compiler, then download the sources, write "make" to build the examples, and make "run_examples" to run the program. Folder MQSMSim contains all library files (header-only implementation). The src folder contains the examples.cpp file with the examples of the paper.
