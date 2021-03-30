#pragma once

#include "BristolCircuit.h"
#include "limbo_instances.h"
#include "config.h"
#include "types.h"
#include <array>
#include <cstdint>
#include <cstdlib>
#include <math.h>
#include <thread>
#include <vector>

limbo_proof_t limbo_prove(const limbo_instance_t &instance,
                              std::vector<uint8_t> &witness_in,
                              const Circuit &C);

std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
limbo_gen(const limbo_instance_t &instance,
            const Circuit &C);

bool limbo_verify(const limbo_instance_t &instance,
                    const std::vector<uint8_t> &witness_out,
                    const limbo_proof_t &signature,
                    const Circuit &C);

unsigned int compute_proof_size(limbo_proof_t proof,
                                limbo_instance_t instance);
