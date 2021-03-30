/*
 *  This file is part of the optimized implementation of the BANQUET signature
 * scheme. See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef BANQUET_INSTANCES_H
#define BANQUET_INSTANCES_H

#include "config.h"
#include <cstdint>
#include <cstdlib>

struct limbo_instance_t {

  // limbo_aes_t aes_params;
  unsigned int nb_mult_gates;
  unsigned int nb_add_gates;
  unsigned int input_size;
  unsigned int output_size;

  uint32_t digest_size;     /* bytes */
  uint32_t seed_size;       /* bytes */
  uint32_t num_rounds;      // T
  uint32_t num_MPC_parties; // N
  uint32_t compression_factor;
  uint32_t lambda; // field expansion factor

  // limbo_params_t params;
};

const limbo_instance_t &limbo_instance_get(unsigned int i);

#endif
