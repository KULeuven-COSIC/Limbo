/*
 * This file gives an example set of possible parameters for the Limbo scheme.
 * */

#include "limbo_instances.h"
#include <stdexcept>

static const limbo_instance_t instances[] = {
    {
        22573,  // nb_mult_gates
        110644, // nb_add_gates
        768,   // input_size
        256,   // output_size
        32,    // digest_size
        16,    // seed_size
        29,    // num_rounds
        64,    // num_MPC_parties
        16,     // compression_factor
        8      // lambda
    },
//Parameters below were defined to benchmark for different circuit size
//these are the parameters used in the Limbo paper.
#if SEC == 40
    // 1024 and gates
    {1024, 1024, 8, 8, 10, 5, 11, 16, 8, 8},  // 1 (16,64,8,11)
    {1024, 1024, 8, 8, 10, 5, 11, 16, 16, 8}, // 2 (16,64,16,11)
    {1024, 1024, 8, 8, 10, 5, 9, 32, 8, 8},   // 3 (32,64,8,9)
    {1024, 1024, 8, 8, 10, 5, 9, 32, 16, 8},  // 4 (32,64,8,9)
    {1024, 1024, 8, 8, 10, 5, 7, 64, 8, 8},   // 5 (64,64,8,7)
    {1024, 1024, 8, 8, 10, 5, 7, 64, 16, 8},  // 6 (64,64,16,7)
    {1024, 1024, 8, 8, 10, 5, 6, 128, 8, 8},  // 7 (128,64,8,6)
    {1024, 1024, 8, 8, 10, 5, 6, 128, 16, 8}, // 8 (128,64,16,6)

    // 16384 and gates
    {16384, 1024, 8, 8, 10, 5, 11, 16, 16, 8}, // 9 (16,64,16,11)
    {16384, 1024, 8, 8, 10, 5, 11, 16, 32, 8}, // 10 (16,64,32,11)
    {16384, 1024, 8, 8, 10, 5, 9, 32, 16, 8},  // 11 (32,64,16,9)
    {16384, 1024, 8, 8, 10, 5, 9, 32, 32, 8},  // 12 (32,64,32,9)
    {16384, 1024, 8, 8, 10, 5, 7, 64, 16, 8},  // 13 (64,64,16,7)
    {16384, 1024, 8, 8, 10, 5, 7, 64, 32, 8},  // 14 (64,64,32,7)
    {16384, 1024, 8, 8, 10, 5, 6, 128, 16, 8}, // 15 (128,64,16,6)
    {16384, 1024, 8, 8, 10, 5, 6, 128, 32, 8}, // 16 (128,64,32,6)

    // 65536 and gates
    {65536, 1024, 8, 8, 10, 5, 11, 16, 16, 8}, // 17 (16,64,16,11)
    {65536, 1024, 8, 8, 10, 5, 11, 16, 32, 8}, // 18 (16,64,32,11)
    {65536, 1024, 8, 8, 10, 5, 9, 32, 16, 8},  // 19 (32,64,16,9)
    {65536, 1024, 8, 8, 10, 5, 9, 32, 32, 8},  // 20 (32,64,32,9)
    {65536, 1024, 8, 8, 10, 5, 7, 64, 16, 8},  // 21 (64,64,16,7)
    {65536, 1024, 8, 8, 10, 5, 7, 64, 32, 8},  // 22 (64,64,32,7)
    {65536, 1024, 8, 8, 10, 5, 6, 128, 16, 8}, // 23 (128,64,16,6)
    {65536, 1024, 8, 8, 10, 5, 6, 128, 32, 8}, // 24 (128,64,32,6)

    // 1048576 and gates
    {1048576, 1024, 8, 8, 10, 5, 11, 16, 32, 8}, // 25 (16,64,32,11)
    {1048576, 1024, 8, 8, 10, 5, 11, 16, 64, 8}, // 26 (16,64,64,11)
    {1048576, 1024, 8, 8, 10, 5, 9, 32, 32, 8},  // 27 (32,64,32,9)
    {1048576, 1024, 8, 8, 10, 5, 9, 32, 64, 8},  // 28 (32,64,64,9)
    {1048576, 1024, 8, 8, 10, 5, 7, 64, 32, 8},  // 29 (64,64,32,7)
    {1048576, 1024, 8, 8, 10, 5, 7, 64, 64, 8},  // 30 (64,64,64,7)
    {1048576, 1024, 8, 8, 10, 5, 6, 128, 32, 8}, // 31 (128,40,32,6)
    {1048576, 1024, 8, 8, 10, 5, 6, 128, 64, 8}, // 32 (128,40,64,6)
#elif SEC == 128
    {1024, 1024, 8, 8, 32, 16, 40, 16, 8, 8},   // 1 (16,64,8,40)
    {1024, 1024, 8, 8, 32, 16, 38, 16, 16, 8},  // 2 (16,64,16,38)
    {1024, 1024, 8, 8, 32, 16, 34, 32, 8, 8},   // 3 (32,64,8,34)
    {1024, 1024, 8, 8, 32, 16, 32, 32, 16, 8},  // 4 (32,64,16,32)
    {1024, 1024, 8, 8, 32, 16, 30, 64, 8, 8},   // 5 (64,64,8,30)
    {1024, 1024, 8, 8, 32, 16, 28, 64, 16, 8},  // 6 (64,64,16,28)
    {1024, 1024, 8, 8, 32, 16, 27, 128, 8, 8},  // 7 (128,64,8,27)
    {1024, 1024, 8, 8, 32, 16, 25, 128, 16, 8}, // 8 (128,64,16,25)

    // 4096 and gates
    {4096, 1024, 8, 8, 32, 16, 42, 16, 8, 8},   // 9 (16,64,8,42)
    {4096, 1024, 8, 8, 32, 16, 40, 16, 16, 8},  // 10 (16,64,16,40)
    {4096, 1024, 8, 8, 32, 16, 36, 32, 8, 8},   // 11 (32,64,8,36)
    {4096, 1024, 8, 8, 32, 16, 34, 32, 16, 8},  // 12 (32,64,16,34)
    {4096, 1024, 8, 8, 32, 16, 32, 64, 8, 8},   // 13 (64,64,8,32)
    {4096, 1024, 8, 8, 32, 16, 30, 64, 16, 8},  // 14 (64,64,16,30)
    {4096, 1024, 8, 8, 32, 16, 29, 128, 8, 8},  // 15 (128,64,8,29)
    {4096, 1024, 8, 8, 32, 16, 27, 128, 16, 8}, // 16 (128,64,16,27)

    // 16384 and gates
    {16384, 1024, 8, 8, 32, 16, 40, 16, 16, 8},  // 17 (16,64,16,40)
    {16384, 1024, 8, 8, 32, 16, 38, 16, 32, 8},  // 18 (16,64,32,38)
    {16384, 1024, 8, 8, 32, 16, 34, 32, 16, 8},  // 19 (32,64,16,34)
    {16384, 1024, 8, 8, 32, 16, 32, 32, 32, 8},  // 20 (32,64,32,32)
    {16384, 1024, 8, 8, 32, 16, 30, 64, 16, 8},  // 21 (64,64,16,30)
    {16384, 1024, 8, 8, 32, 16, 28, 64, 32, 8},  // 22 (64,64,32,28)
    {16384, 1024, 8, 8, 32, 16, 27, 128, 16, 8}, // 23 (128,64,16,27)
    {16384, 1024, 8, 8, 32, 16, 25, 128, 32, 8}, // 24 (128,64,32,25)

    // 65536 and gates
    {65536, 1024, 8, 8, 32, 16, 42, 16, 16, 8},  // 25 (16,64,16,42)
    {65536, 1024, 8, 8, 32, 16, 40, 16, 32, 8},  // 26 (16,64,32,40)
    {65536, 1024, 8, 8, 32, 16, 36, 32, 16, 8},  // 27 (32,64,16,36)
    {65536, 1024, 8, 8, 32, 16, 34, 32, 32, 8},  // 28 (32,64,32,34)
    {65536, 1024, 8, 8, 32, 16, 32, 64, 16, 8},  // 29 (64,64,16,32)
    {65536, 1024, 8, 8, 32, 16, 30, 64, 32, 8},  // 30 (64,64,32,30)
    {65536, 1024, 8, 8, 32, 16, 29, 128, 16, 8}, // 31 (128,64,16,29)
    {65536, 1024, 8, 8, 32, 16, 27, 128, 32, 8}, // 32 (128,64,32,27)

    // 262144 and gates
    {262144, 1024, 8, 8, 32, 16, 43, 16, 16, 8},  // 33 (16,64,16,43)
    {262144, 1024, 8, 8, 32, 16, 41, 16, 32, 8},  // 34 (16,64,32,41)
    {262144, 1024, 8, 8, 32, 16, 37, 32, 16, 8},  // 35 (32,64,16,37)
    {262144, 1024, 8, 8, 32, 16, 35, 32, 32, 8},  // 36 (32,64,32,35)
    {262144, 1024, 8, 8, 32, 16, 33, 64, 16, 8},  // 37 (64,64,16,33)
    {262144, 1024, 8, 8, 32, 16, 31, 64, 32, 8},  // 38 (64,64,32,31)
    {262144, 1024, 8, 8, 32, 16, 30, 128, 16, 8}, // 39 (128,64,16,30)
    {262144, 1024, 8, 8, 32, 16, 28, 128, 32, 8}, // 40 (128,64,32,28)

    // 1048576 and gates
    {1048576, 1024, 8, 8, 32, 16, 43, 16, 32, 8},  // 41 (16,64,32,43)
    {1048576, 1024, 8, 8, 32, 16, 41, 16, 64, 8},  // 42 (16,64,64,41)
    {1048576, 1024, 8, 8, 32, 16, 37, 32, 32, 8},  // 43 (32,64,32,37)
    {1048576, 1024, 8, 8, 32, 16, 35, 32, 64, 8},  // 44 (32,64,64,35)
    {1048576, 1024, 8, 8, 32, 16, 33, 64, 32, 8},  // 45 (64,64,32,33)
    {1048576, 1024, 8, 8, 32, 16, 31, 64, 64, 8},  // 46 (64,64,64,31)
    {1048576, 1024, 8, 8, 32, 16, 30, 128, 32, 8}, // 47 (128,40,32,30)
    {1048576, 1024, 8, 8, 32, 16, 28, 128, 64, 8}, // 48 (128,40,64,28)
#endif
};

const limbo_instance_t &limbo_instance_get(unsigned int i) {
  return instances[i];
}