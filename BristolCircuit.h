/*
Copyright (c) 2017, The University of Bristol, Senate House, Tyndall Avenue,
Bristol, BS8 1TH, United Kingdom. Copyright (c) 2020, COSIC-KU Leuven,
Kasteelpark Arenberg 10, bus 2452, B-3001 Leuven-Heverlee, Belgium.

All rights reserved
*/
#ifndef _BristolCircuit
#define _BristolCircuit

/* This function defines a circuit
 *  A gate has one or two input wires (unless it is a MAND gate)
 *    - Only INV/NOT has one
 *  A gate has only one output wire (unless it is a MAND gate)
 *
 * The input file format is a little different from "classic"
 * Bristol format, we call the new format "Bristol Fashion"
 *
 * The new format is...
 *
 *   - A line defining the number of gates and then the number of wires in the
 *     circuit.
 *   - The number of input values niv  (e.g. if doing z=a+b mod p we have niv=3)
 *         - Then niv numbers defining the number of input wires per input value
 *           ni_1,..,ni_niv
 *   - The number of output values nov (e.g. if doing z=a+b mod p we have nov=1)
 *         - Then nov numbers defining the number of output wires per output
 * value, no_1,...,no_nov
 *   - The wire numbering is ordered so that the first i_0 wires correspond to
 * the first input value, the next i_1 wires correspond to the second input
 * value and so on.
 *   - With the last (o_0+...+o_{n-1}) wires corresponding to the outputs of
 *     the function, where n=no_1+...+no_nov
 *   - The gates are ordered topologically, so we can evaluate them in sequence.
 *   - Each gate is defined by
 *      -    Number input wires   (1 or 2, unless a MAND gate )
 *      -    Number output wires  (Always 1, unless a MAND gate)
 *      -    List of input wires
 *      -    List of output wires
 *      -    Gate operation (XOR, AND, INV, EQ, EQW or MAND).
 *     This is a bit redundant, as the first two entries can be inferred from
 *     the last, but we keep this for backwards compatibility reasons
 *   - So for example
 *          2 1 3 4 5 XOR
 *     corresponds to
 *          w_5 = XOR(w_3,w_4)
 *   - We also use
 *          1 1 0 3 EQ
 *          1 1 1 4 EQ
 *     to say that wire 3 is assigned the value 0 and wire 4 the value 1
 *   - And we use
 *          1 1 3 4 EQW
 *     to say wire 4 should equal wire 3
 *   - The MAND gate is a multiple AND gate it takes 2n input wires and
 *     produces n output wires. A gate of the form
 *          4 2 0 2 1 3 4 5 MAND
 *     executes the two MAND operations concurrently...
 *          2 1 0 1 4 AND
 *          2 1 2 3 5 AND
 *   - We call a circuit without MAND gates a `basic' Bristol Fashion circuit,
 *     and when we have MAND gates we say it is an `extended' Bristol Fashion
 *     circuit.
 */

#include "limbo_instances.h"
#include "gsl-lite.hpp"
#include "types.h"
#include <array>
#include <cstdint>
#include <vector>

using namespace std;

enum GateType { XOR, AND, INV, EQ, EQW, MAND };

// This function only works for type T != MAND
unsigned int cnt_numI(const GateType &T);

class Circuit {
  unsigned int nWires;

  // The number of input and output wires per variable
  vector<unsigned int> numI;
  vector<unsigned int> numO;

  // Each Gate is given as an entry in these arrays
  //   Type, Input Wires, Output Wires
  vector<GateType> GateT;
  vector<vector<unsigned int>> GateI; // nGates x v array
  vector<vector<unsigned int>> GateO; // nGates x v array

  // These maps are for AND gates only, not MAND gates
  vector<unsigned int> map;  // Map of n'th AND gate to the gate number
  vector<unsigned int> imap; // The inverse map

  unsigned int num_AND;       // Number of AND gates
  unsigned int total_num_AND; // Number of AND gates

  // Used for testing within the topological sort
  bool gate_is_ok(unsigned int i, const vector<bool> &used) const;

public:
  void recompute_map(); // Recomputes the mapping function

  void swap_gate(unsigned int i, unsigned int j); // Swaps two gates around

  // Applies a topological sort to the circuit
  // If test=true, just does a test
  void sort(bool flag = false);

  unsigned int get_nGates() const { return GateT.size(); }
  unsigned int get_nWires() const { return nWires; }

  // Returns the number of simple AND gates
  unsigned int num_AND_gates() const { return num_AND; }

  // Returns the total number of AND gates including MAND
  unsigned int total_num_AND_gates() const { return total_num_AND; }

  unsigned int num_inputs() const { return numI.size(); }
  unsigned int num_outputs() const { return numO.size(); }
  unsigned int total_num_out() const {
    unsigned int result = 0;
    for (size_t i = 0; i < numO.size(); i++) {
      result += numO[i];
    }
    return result;
  }

  // Return the number of wires in input/output variable i
  unsigned int num_iWires(unsigned int i) const { return numI[i]; }
  unsigned int num_oWires(unsigned int i) const { return numO[i]; }

  GateType get_GateType(unsigned int i) const { return GateT[i]; }

  unsigned int Gate_Wire_In(unsigned int i, unsigned int j) const {
    if ((i > GateI.size()) || (j >= GateI[i].size() && GateT[i] != MAND) ||
        (j == 1 && (GateT[i] == INV || GateT[i] == EQ || GateT[i] == EQW))) {
      throw std::runtime_error("Circuit error");
    }
    return GateI[i][j];
  }

  unsigned int Gate_Wire_Out(unsigned int i) const {
    if (i > GateO.size() || GateT[i] == MAND) {
      throw std::runtime_error("Circuit error");
    }
    return GateO[i][0];
  }

  unsigned int Gate_Wire_Out(unsigned int i, unsigned int j) const {
    if ((i > GateO.size()) || (GateT[i] != MAND) || (j >= GateO[i].size())) {
      throw std::runtime_error("Circuit error");
    }
    return GateO[i][j];
  }

  unsigned int MAND_Gate_Size(unsigned int i) const {
    if (GateT[i] != MAND) {
      throw std::runtime_error("Circuit error");
    }
    return GateO[i].size();
  }

  unsigned int get_nth_AND_Gate(unsigned int i) const { return map[i]; }

  unsigned int map_to_AND_Gate(unsigned int i) const {
    if (imap[i] >= GateT.size()) {
      throw std::runtime_error("Circuit error");
    }
    return imap[i];
  }

  // File IO
  friend ostream &operator<<(ostream &s, const Circuit &C);
  friend istream &operator>>(istream &s, Circuit &C);
  // Print gate i
  void output_gate(ostream &s, unsigned int i) const;

  // This function is for testing the circuits
  // outputs is resized by this function so can be blank on
  // entry
  bool base_circuit(const vector<uint8_t> &inputs,
                    vector<uint8_t> &outputs) const;

  std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, std::vector<uint8_t>>
  base_circuit_mult_output(const vector<uint8_t> &inputs,
                           vector<uint8_t> &outputs) const;

  void base_circuit_shares(
      const std::vector<std::array<uint64_t, PARTY64>> &witness_shares_in,
      const std::vector<std::array<uint64_t, PARTY64>> &output_mult_shares,
      std::vector<std::array<uint64_t, PARTY64>> &witness_shares_out,
      std::vector<std::array<uint64_t, PARTY64>> &right_input_mult_shares_out,
      std::vector<std::array<uint64_t, PARTY64>> &left_input_mult_shares_out)
      const;

  /* This function returns the AND gate depths for each gate
   *   Assumed Circuit is topologically sorted
   */
  vector<unsigned int> compute_depth() const;

  /* This merges all AND/MAND gates with a given depth
   *   Assumed Circuit is topologically sorted
   */
  void merge_AND_gates();
};

#endif
