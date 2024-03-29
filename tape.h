#pragma once

#include "limbo.h"
extern "C" {
#include "kdf_shake.h"
}
#include "tree.h"
#include <cstdlib>

class RandomTape {
private:
  /* data */
  hash_context ctx;

public:
  RandomTape(const std::vector<uint8_t> &seed, const limbo_salt_t &salt,
             size_t rep_index, size_t party_index);
  ~RandomTape() = default;

  void squeeze_bytes(uint8_t *out, size_t len);
};
