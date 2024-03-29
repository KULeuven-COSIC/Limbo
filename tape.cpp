#include "tape.h"

RandomTape::RandomTape(const std::vector<uint8_t> &seed,
                       const limbo_salt_t &salt, size_t rep_index,
                       size_t party_index) {
  hash_init(&ctx, seed.size() * 2);
  hash_update(&ctx, seed.data(), seed.size());
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)rep_index);
  hash_update_uint16_le(&ctx, (uint16_t)party_index);
  hash_final(&ctx);
}

void RandomTape::squeeze_bytes(uint8_t *out, size_t len) {
  hash_squeeze(&ctx, out, len);
}
