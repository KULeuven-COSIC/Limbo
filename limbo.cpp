#include "limbo.h"

//#include "aes.h"
#include "field.h"
#include "tape.h"
#include "tree.h"
#include <algorithm>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cstring>
#include <iostream>

extern "C" {
#include "kdf_shake.h"
#include "randomness.h"
}

/* GLOBAL VARIABLES FOR TIMING */

#ifdef DETAILED_IN_CODE_TIMERS
extern std::vector<std::chrono::high_resolution_clock::time_point> start_timers;
extern std::vector<std::chrono::high_resolution_clock::time_point> stop_timers;
extern std::vector<std::chrono::duration<double>> time_spans;
#endif
/**/

namespace {
inline void hash_update_GF2E(hash_context *ctx,
                             const limbo_instance_t &instance,
                             const field::GF2E &element) {
  // 8 bytes is enough for supported field sizes
  std::array<uint8_t, 8> buffer;
  element.to_bytes(buffer.data());
  hash_update(ctx, buffer.data(), instance.lambda);
}

std::pair<limbo_salt_t, std::vector<std::vector<uint8_t>>>
generate_salt_and_seeds(const limbo_instance_t &instance,
                        const std::vector<uint8_t> &witness_in) {
  // salt, seed_1, ..., seed_r = H(instance||sk||pk||m)
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, witness_in.data(), witness_in.size());
  hash_final(&ctx);

  limbo_salt_t salt;
  hash_squeeze(&ctx, salt.data(), salt.size());
  std::vector<std::vector<uint8_t>> seeds;
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    std::vector<uint8_t> s(instance.seed_size);
    hash_squeeze(&ctx, s.data(), s.size());
    seeds.push_back(s);
  }
  return std::make_pair(salt, seeds);
}

std::vector<uint8_t> commit_to_party_seed(const limbo_instance_t &instance,
                                          const std::vector<uint8_t> &seed,
                                          const limbo_salt_t &salt,
                                          size_t rep_idx, size_t party_idx) {
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)rep_idx);
  hash_update_uint16_le(&ctx, (uint16_t)party_idx);
  hash_update(&ctx, seed.data(), seed.size());
  hash_final(&ctx);

  std::vector<uint8_t> commitment(instance.digest_size);
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<uint8_t> phase_1_commitment(
    const limbo_instance_t &instance, const limbo_salt_t &salt,
    const std::vector<uint8_t> &witness_out,
    const std::vector<std::vector<std::vector<uint8_t>>> &commitments,
    const std::vector<std::vector<uint8_t>> &witness_deltas,
    const std::vector<std::vector<uint8_t>> &mult_deltas,
    std::vector<std::vector<std::array<uint64_t, PARTY64>>>
        &output_broadcasts) {

  hash_context ctx;
  hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_1);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, witness_out.data(), witness_out.size());

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      hash_update(&ctx, commitments[repetition][party].data(),
                  commitments[repetition][party].size());
    }
    // dirty hack, need to check
    hash_update(&ctx, (uint8_t *)output_broadcasts[repetition].data(),
                PARTY64 * output_broadcasts[repetition].size() * 8);
    hash_update(&ctx, witness_deltas[repetition].data(),
                witness_deltas[repetition].size());
    hash_update(&ctx, mult_deltas[repetition].data(),
                mult_deltas[repetition].size());
  }
  hash_final(&ctx);

  std::vector<uint8_t> commitment(instance.digest_size);
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<field::GF2E> phase_1_expand(const limbo_instance_t &instance,
                                        const std::vector<uint8_t> &h_1,
                                        unsigned int nbChallenges) {
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, h_1.data(), h_1.size());
  hash_final(&ctx);

  std::vector<field::GF2E> r_ejs(nbChallenges);
  std::vector<uint8_t> lambda_sized_buffer(instance.lambda);
  r_ejs.reserve(instance.num_rounds);
  for (size_t e = 0; e < nbChallenges; e++) {
    hash_squeeze(&ctx, lambda_sized_buffer.data(), lambda_sized_buffer.size());
    r_ejs[e].from_bytes(lambda_sized_buffer.data());
  }
  return r_ejs;
}

std::vector<uint8_t> phase_2_commitment(
    const limbo_instance_t &instance, const limbo_salt_t &salt,
    const std::vector<uint8_t> &h_1,
    const std::vector<std::vector<uint8_t>> compression_c,
    size_t compression_nb,
    const std::vector<std::vector<std::vector<field::GF2E>>> &P_deltas) {

  hash_context ctx;
  hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_2);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, h_1.data(), h_1.size());

  for (size_t cn = 0; cn < compression_nb; cn++) {
    hash_update(&ctx, compression_c[cn].data(), compression_c[cn].size());
  }

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t k = 0; k < P_deltas[repetition][compression_nb].size(); k++) {
      hash_update_GF2E(&ctx, instance, P_deltas[repetition][compression_nb][k]);
    }
  }
  hash_final(&ctx);

  std::vector<uint8_t> commitment(instance.digest_size);
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<field::GF2E> phase_2_expand(
    const limbo_instance_t &instance, const std::vector<uint8_t> &h_2,
    const std::vector<field::GF2E> &forbidden_values, const size_t num_chall) {
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, h_2.data(), h_2.size());
  hash_final(&ctx);

  std::vector<uint8_t> lambda_sized_buffer(instance.lambda);
  std::vector<field::GF2E> R_es;
  for (size_t e = 0; e < num_chall; e++) { // instance.num_rounds; e++) {
    while (true) {
      hash_squeeze(&ctx, lambda_sized_buffer.data(),
                   lambda_sized_buffer.size());
      //  check that R is not in {0,...m2-1}
      field::GF2E candidate_R;
      candidate_R.from_bytes(lambda_sized_buffer.data());
      bool good = true;
      for (size_t k = 0; k < forbidden_values.size(); k++) {
        if (candidate_R == forbidden_values[k]) {
          good = false;
          break;
        }
      }
      if (good) {
        R_es.push_back(candidate_R);
        break;
      }
    }
  }
  return R_es;
}

std::vector<uint8_t>
phase_3_commitment(const limbo_instance_t &instance,
                   const limbo_salt_t &salt, const std::vector<uint8_t> &h_2,
                   const RepContainer<field::GF2E> &rep_shared_mult_left,
                   const RepContainer<field::GF2E> &rep_shared_mult_right,
                   const RepContainer<field::GF2E> &rep_shared_inner_prod,
                   size_t vector_size) {

  hash_context ctx;
  hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_3);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, h_2.data(), h_2.size());

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t j = 0; j < vector_size; j++) {
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        if (j == 0) {
          hash_update_GF2E(
              &ctx, instance,
              rep_shared_inner_prod.get(
                  repetition, party)[instance.compression_factor - 1]);
        }
        hash_update_GF2E(&ctx, instance,
                         rep_shared_mult_left.get(repetition, party)[j]);
        hash_update_GF2E(&ctx, instance,
                         rep_shared_mult_right.get(repetition, party)[j]);
      }
    }
  }
  hash_final(&ctx);

  std::vector<uint8_t> commitment(instance.digest_size);
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<uint16_t> phase_3_expand(const limbo_instance_t &instance,
                                     const std::vector<uint8_t> &h_3) {
  assert(instance.num_MPC_parties < (1ULL << 16));
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, h_3.data(), h_3.size());
  hash_final(&ctx);
  size_t num_squeeze_bytes = instance.num_MPC_parties > 256 ? 2 : 1;

  std::vector<uint16_t> opened_parties;
  uint16_t mask = (1ULL << ceil_log2(instance.num_MPC_parties)) - 1;
  for (size_t e = 0; e < instance.num_rounds; e++) {
    uint16_t party;
    do {
      hash_squeeze(&ctx, (uint8_t *)&party, num_squeeze_bytes);
      party = le16toh(party);
      party = party & mask;
    } while (party >= instance.num_MPC_parties);
    opened_parties.push_back(party);
  }
  return opened_parties;
}
} // namespace

std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
limbo_gen(const limbo_instance_t &instance,
            const Circuit &C)
{
  std::vector<uint8_t> witness_in(instance.input_size),
      witness_out(instance.output_size);

  // while (true) {
  rand_bytes(witness_in.data(), witness_in.size());
  uint8_t mask = 0b00000001;
  for (size_t i = 0; i < witness_in.size(); i++) {
    witness_in[i] &= mask;
  }
  if (not C.base_circuit(witness_in, witness_out))
  {
    throw std::runtime_error("invalied params");
  }

  std::pair<std::vector<uint8_t>, std::vector<uint8_t>> res =
      std::make_pair(witness_in, witness_out);
  return res;
}

limbo_proof_t limbo_prove(const limbo_instance_t &instance,
                              std::vector<uint8_t> &witness_in,
                              const Circuit &C)
{
  // init modulus of extension field F_{2^{8\lambda}}
  size_t maxp64 = instance.num_MPC_parties / 65;
  field::GF2E::init_extension_field(instance);

  // get mult gates inputs and outputs for circuit evaluation
  std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, std::vector<uint8_t>>
      mult_tuple;
  std::vector<uint8_t> witness_out;

#ifdef DETAILED_IN_CODE_TIMERS
  start_timers[0] = std::chrono::high_resolution_clock::now();
#endif

  mult_tuple = C.base_circuit_mult_output(witness_in, witness_out);

#ifdef DETAILED_IN_CODE_TIMERS
  stop_timers[0] = std::chrono::high_resolution_clock::now();
  time_spans[0] += std::chrono::duration_cast<std::chrono::duration<double>>(
      stop_timers[0] - start_timers[0]);

  start_timers[1] = std::chrono::high_resolution_clock::now();
#endif
  // generate salt and master seeds for each repetition
  auto [salt, master_seeds] = generate_salt_and_seeds(instance, witness_in);

  // buffer for squeezing field elements into
  std::vector<uint8_t> lambda_sized_buffer(instance.lambda);

  // do parallel repetitions
  // create seed trees and random tapes
  std::vector<SeedTree> seed_trees;
  std::vector<std::vector<RandomTape>> random_tapes;
  std::vector<std::vector<std::vector<uint8_t>>> party_seed_commitments;

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    // generate seed tree for the N parties
    seed_trees.push_back(SeedTree(master_seeds[repetition],
                                  instance.num_MPC_parties, salt, repetition));

    // commit to each party's seed;
    std::vector<std::vector<uint8_t>> current_party_seed_commitments;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      current_party_seed_commitments.push_back(commit_to_party_seed(
          instance, seed_trees[repetition].get_leaf(party).value(), salt,
          repetition, party));
    }
    party_seed_commitments.push_back(current_party_seed_commitments);

    // create random tape for each party
    std::vector<RandomTape> party_tapes;
    party_tapes.reserve(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      party_tapes.emplace_back(seed_trees[repetition].get_leaf(party).value(),
                               salt, repetition, party);
    }
    random_tapes.push_back(party_tapes);
  }

#ifdef DETAILED_IN_CODE_TIMERS
  stop_timers[1] = std::chrono::high_resolution_clock::now();
  time_spans[1] += std::chrono::duration_cast<std::chrono::duration<double>>(
      stop_timers[1] - start_timers[1]);
#endif
  /////////////////////////////////////////////////////////////////////////////
  // phase 1: commit to executions of the circuit
  /////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_shared_witness(
      instance.num_rounds);
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_output_broadcasts(
      instance.num_rounds);
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_shared_mult_left(
      instance.num_rounds);
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_shared_mult_right(
      instance.num_rounds);
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_shared_mult_out(
      instance.num_rounds);

  std::vector<std::vector<uint8_t>> rep_witness_deltas;
  std::vector<std::vector<uint8_t>> rep_mult_deltas;

// Loop over all repetitions (tau)
#ifdef IN_CODE_TIMERS
  std::chrono::high_resolution_clock::time_point start_time =
      std::chrono::high_resolution_clock::now();
#endif

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
#ifdef DETAILED_IN_CODE_TIMERS
    start_timers[2] = std::chrono::high_resolution_clock::now();
#endif
    rep_shared_witness[repetition].resize(instance.input_size);
    rep_output_broadcasts[repetition].resize(instance.output_size);
    rep_shared_mult_left[repetition].resize(instance.nb_mult_gates);
    rep_shared_mult_right[repetition].resize(instance.nb_mult_gates);
    rep_shared_mult_out[repetition].resize(instance.nb_mult_gates);

    // generate sharing of witness
    std::vector<uint8_t> witness_delta = witness_in;

    // Loops over all parties to sample the shares of witness
    std::vector<uint8_t> shared_witness(instance.input_size >> 3);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // auto shared_witness = rep_shared_witness.get(repetition, party);
      // sample as many random bits as needed
      random_tapes[repetition][party].squeeze_bytes(shared_witness.data(),
                                                    shared_witness.size());
      // set share of the player correctly
      for (size_t i = 0; i < shared_witness.size() - 1; i++) {
        for (size_t j = 0; j < 8; j++) {
          rep_shared_witness[repetition][i * 8 + j][party >> 6] ^=
              (((uint64_t)((shared_witness[i] >> j) & 1)) << (party % 64));
        }
      }
      for (size_t j = 0;
           j < std::max((long unsigned int)8,
                        instance.input_size - 8 * (shared_witness.size() - 1));
           j++) {
        rep_shared_witness[repetition][(shared_witness.size() - 1) * 8 + j]
                          [party >> 6] ^=
            (((uint64_t)((shared_witness[shared_witness.size() - 1] >> j) & 1))
             << (party));
      }
    }

    // Compute the witness_delta which will be xored to the first share
    uint8_t flip;
    for (size_t i = 0; i < instance.input_size; i++) {
      flip ^= flip;
      for (size_t p64 = 0; p64 <= maxp64; p64++) {
        flip ^=
            (__builtin_popcountll(rep_shared_witness[repetition][i][p64]) % 2);
      }
      witness_delta[i] ^= flip;
      // also fix share of the first party
      rep_shared_witness[repetition][i][0] ^= witness_delta[i];
    }

    rep_witness_deltas.push_back(witness_delta);

    // generate sharing of output of mult gates
    std::vector<uint8_t> mult_out_deltas = std::get<2>(mult_tuple);
    // Loops over all parties
    std::vector<uint8_t> shared_mult_out(instance.nb_mult_gates >> 3);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // fill share with random data
      random_tapes[repetition][party].squeeze_bytes(shared_mult_out.data(),
                                                    shared_mult_out.size());

      // set share of the player correctly
      for (size_t i = 0; i < shared_mult_out.size() - 1; i++) {
        for (size_t j = 0; j < 8;
             j++) // std::max((long unsigned int)8,instance.nb_mult_gates -
                  // 8*i); j++)
        {
          rep_shared_mult_out[repetition][i * 8 + j][party >> 6] ^=
              (((uint64_t)((shared_mult_out[i] >> j) & 1)) << (party % 64));
        }
      }
      for (size_t j = 0; j < std::max((long unsigned int)8,
                                      instance.nb_mult_gates -
                                          8 * (shared_mult_out.size() - 1));
           j++) {
        rep_shared_mult_out[repetition][(shared_mult_out.size() - 1) * 8 + j]
                           [party >> 6] ^=
            (((uint64_t)((shared_mult_out[shared_mult_out.size() - 1] >> j) &
                         1))
             << (party));
      }
    }
    // Compute the mult_out_delta which will be xored to the first share
    for (size_t i = 0; i < instance.nb_mult_gates; i++) {
      flip ^= flip;
      for (size_t p64 = 0; p64 <= maxp64; p64++) {
        flip ^=
            (__builtin_popcountll(rep_shared_mult_out[repetition][i][p64]) % 2);
      }
      mult_out_deltas[i] ^= flip;
      // also fix share of the first party
      rep_shared_mult_out[repetition][i][0] ^= mult_out_deltas[i];
    }

    rep_mult_deltas.push_back(mult_out_deltas);

#ifdef DETAILED_IN_CODE_TIMERS
    stop_timers[2] = std::chrono::high_resolution_clock::now();
    time_spans[2] += std::chrono::duration_cast<std::chrono::duration<double>>(
        stop_timers[2] - start_timers[2]);

    start_timers[3] = std::chrono::high_resolution_clock::now();
#endif
    C.base_circuit_shares(
        rep_shared_witness[repetition], rep_shared_mult_out[repetition],
        rep_output_broadcasts[repetition], rep_shared_mult_right[repetition],
        rep_shared_mult_left[repetition]);

#ifdef DETAILED_IN_CODE_TIMERS
    stop_timers[3] = std::chrono::high_resolution_clock::now();
    time_spans[3] += std::chrono::duration_cast<std::chrono::duration<double>>(
        stop_timers[3] - start_timers[3]);
#endif
  }

#ifdef IN_CODE_TIMERS
  std::chrono::high_resolution_clock::time_point tmp_time =
      std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(tmp_time -
                                                                start_time);
  std::cout << "PHASE 1 TOOK " << time_span.count() << std::endl;
#endif

  ///////////////////////////////////////////////////////////////////////////////
  //// phase 2: challenge the multiplications
  ///////////////////////////////////////////////////////////////////////////////

  // commit to salt, (all commitments of parties seeds, key_delta, t_delta)
  // for all repetitions
  std::vector<uint8_t> h_1 = phase_1_commitment(
      instance, salt, witness_out, party_seed_commitments, rep_witness_deltas,
      rep_mult_deltas, rep_output_broadcasts);

// expand challenge hash to 1 value R (still storing a vector in case we need
// more)
#ifdef DETAILED_IN_CODE_TIMERS
  start_timers[4] = std::chrono::high_resolution_clock::now();
#endif
#ifdef SAME_CHALLENGE
  std::vector<field::GF2E> Rs = phase_1_expand(instance, h_1, 1);
#else
  std::vector<field::GF2E> Rs =
      phase_1_expand(instance, h_1, instance.num_rounds);
#endif
  std::vector<std::vector<field::GF2E>> Rs_powers =
      field::compute_powers(Rs, instance.nb_mult_gates);

#ifdef DETAILED_IN_CODE_TIMERS
  stop_timers[4] = std::chrono::high_resolution_clock::now();
  time_spans[4] += std::chrono::duration_cast<std::chrono::duration<double>>(
      stop_timers[4] - start_timers[4]);
#endif

/////////////////////////////////////////////////////////////////////////////
// phase 2Bis: lift mult triples and make it an inner product
/////////////////////////////////////////////////////////////////////////////
#ifdef IN_CODE_TIMERS
  start_time = std::chrono::high_resolution_clock::now();
#endif
#ifdef DETAILED_IN_CODE_TIMERS
  start_timers[5] = std::chrono::high_resolution_clock::now();
#endif

  RepContainer<field::GF2E> rep_shared_mult_left_lifted(
      instance.num_rounds, instance.num_MPC_parties, instance.nb_mult_gates);
  RepContainer<field::GF2E> rep_shared_mult_right_lifted(
      instance.num_rounds, instance.num_MPC_parties, instance.nb_mult_gates);
  RepContainer<field::GF2E> rep_shared_inner_prod_lifted(
      instance.num_rounds, instance.num_MPC_parties,
      2 * instance.compression_factor + 1);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // auto share_mult_left = rep_shared_mult_left.get(repetition, party);
      auto share_mult_left_lifted =
          rep_shared_mult_left_lifted.get(repetition, party);

      // auto share_mult_right = rep_shared_mult_right.get(repetition, party);
      auto share_mult_right_lifted =
          rep_shared_mult_right_lifted.get(repetition, party);

      // auto share_mult_out = rep_shared_mult_out.get(repetition, party);
      auto share_inner_prod_lifted =
          rep_shared_inner_prod_lifted.get(repetition, party);
      size_t cp = party % 64;
      size_t ccp = party / 64;
      uint64_t mask = (1ULL << cp);
      for (size_t nbm = 0; nbm < instance.nb_mult_gates; nbm++) {
        share_mult_right_lifted[nbm] = field::GF2E(
            (rep_shared_mult_right[repetition][nbm][ccp] & (mask)) >> cp);
#ifdef SAME_CHALLENGE
        share_mult_left_lifted[nbm] = field::GF2E(
            ((rep_shared_mult_left[repetition][nbm][ccp] & (mask)) >> cp),
            Rs_powers[0][nbm].data);
        share_inner_prod_lifted[instance.compression_factor - 1] += field::GF2E(
            (rep_shared_mult_out[repetition][nbm][ccp] & (mask)) >> cp,
            Rs_powers[0][nbm].data);
#else
        share_mult_left_lifted[nbm] = field::GF2E(
            ((rep_shared_mult_left[repetition][nbm][ccp] & (mask)) >> cp),
            Rs_powers[repetition][nbm].data);
        share_inner_prod_lifted[instance.compression_factor - 1] += field::GF2E(
            (rep_shared_mult_out[repetition][nbm][ccp] & (mask)) >> cp,
            Rs_powers[repetition][nbm].data);
#endif
      }
    }
  }
#ifdef DETAILED_IN_CODE_TIMERS
  stop_timers[5] = std::chrono::high_resolution_clock::now();
  time_spans[5] += std::chrono::duration_cast<std::chrono::duration<double>>(
      stop_timers[5] - start_timers[5]);
#endif
#ifdef IN_CODE_TIMERS
  tmp_time = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
      tmp_time - start_time);
  std::cout << "PHASE 2BIS TOOK " << time_span.count() << std::endl;
#endif

/////////////////////////////////////////////////////////////////////////////
// Start compression phases !
/////////////////////////////////////////////////////////////////////////////
#ifdef IN_CODE_TIMERS
  start_time = std::chrono::high_resolution_clock::now();
#endif

  unsigned int nbCompressions =
      (unsigned int)(log(instance.nb_mult_gates) /
                     log(instance.compression_factor));

#ifdef DETAILED_IN_CODE_TIMERS
  start_timers[6] = std::chrono::high_resolution_clock::now();
#endif
  // Get k points for f,g
  std::vector<field::GF2E> eval_points =
      field::get_first_n_field_elements(instance.compression_factor);
  // Get 2*k - 1 points for h
  std::vector<field::GF2E> big_eval_points =
      field::get_first_n_field_elements(2 * instance.compression_factor - 1);
  // precompute lagrange pol
  std::vector<std::vector<field::GF2E>> lagrange_pol_eval =
      field::precompute_lagrange_polynomials(eval_points);
  // Evaluate those polynomoials at k+1,...,2k-1
  std::vector<std::vector<field::GF2E>> lagrange_pol_eval_inner_prod_s(
      instance.compression_factor);
  for (size_t point = instance.compression_factor;
       point < 2 * instance.compression_factor - 1; point++) {
    lagrange_pol_eval_inner_prod_s[point - instance.compression_factor].resize(
        instance.compression_factor);
    for (size_t slice = 0; slice < instance.compression_factor; slice++) {
      lagrange_pol_eval_inner_prod_s[point -
                                     instance.compression_factor][slice] =
          field::eval(lagrange_pol_eval[slice], big_eval_points[point]);
    }
  }

  // precompute lagrange pol
  std::vector<std::vector<field::GF2E>> big_lagrange_pol_eval =
      field::precompute_lagrange_polynomials(big_eval_points);

#ifdef DETAILED_IN_CODE_TIMERS
  stop_timers[6] = std::chrono::high_resolution_clock::now();
  time_spans[6] += std::chrono::duration_cast<std::chrono::duration<double>>(
      stop_timers[6] - start_timers[6]);

  start_timers[7] = std::chrono::high_resolution_clock::now();
#endif
  // vector to hold the injected inner products
  std::vector<std::vector<std::vector<field::GF2E>>> h_deltas(
      instance.num_rounds,
      std::vector<std::vector<field::GF2E>>(
          nbCompressions,
          std::vector<field::GF2E>(2 * instance.compression_factor - 1)));

  // vector to hold the commitments
  std::vector<std::vector<uint8_t>> compression_commitments(nbCompressions);

  // we expand to get the party challenge
  // And we open the last tuple
  std::vector<uint16_t> missing_parties;
  std::vector<uint8_t> h_3;
  std::vector<std::vector<field::GF2E>> f_final(instance.num_rounds),
      g_final(instance.num_rounds);
  std::vector<field::GF2E> h_final(instance.num_rounds);

  // additional point sampled at the last compression
  std::vector<std::vector<std::vector<field::GF2E>>> add_point_f(
      instance.num_rounds);
  std::vector<std::vector<std::vector<field::GF2E>>> add_point_g(
      instance.num_rounds);
  // hold the injected inner product for the last compression
  std::vector<std::vector<field::GF2E>> inner_product_output_last_compression(
      instance.num_rounds,
      std::vector<field::GF2E>(2 * instance.compression_factor + 1));

  std::vector<field::GF2E> f_i_poly(instance.compression_factor);
  std::vector<field::GF2E> g_i_poly(instance.compression_factor);
  std::vector<field::GF2E> h_poly(2 * instance.compression_factor - 1);

  unsigned int current_m = instance.nb_mult_gates;
  unsigned int vector_size = current_m;

#ifdef DETAILED_IN_CODE_TIMERS
  stop_timers[7] = std::chrono::high_resolution_clock::now();
  time_spans[7] += std::chrono::duration_cast<std::chrono::duration<double>>(
      stop_timers[7] - start_timers[7]);

#endif
  for (unsigned int compression = 1; compression <= nbCompressions;
       compression++) {

    current_m = vector_size;
    vector_size =
        (unsigned int)ceil(((double)current_m / instance.compression_factor));

#ifdef SAME_CHALLENGE
    std::vector<field::GF2E> left_input(instance.compression_factor *
                                        vector_size);
    std::vector<field::GF2E> right_input(instance.compression_factor *
                                         vector_size);
    std::vector<field::GF2E> inner_product_output(
        2 * instance.compression_factor - 1);
#endif

    // In the last compression round:
    // we add an additional point to f_j and g_j for zero knowledge
    if (compression == nbCompressions) {
#ifdef DETAILED_IN_CODE_TIMERS
      start_timers[8] = std::chrono::high_resolution_clock::now();
#endif
      // Get points 1,...,k+1 for f,g
      eval_points =
          field::get_first_n_field_elements(instance.compression_factor + 1);
      // Get 1,...,2*k+1 points for h
      big_eval_points = field::get_first_n_field_elements(
          2 * instance.compression_factor + 1);
      // precompute lagrange pol
      lagrange_pol_eval = field::precompute_lagrange_polynomials(eval_points);
      big_lagrange_pol_eval =
          field::precompute_lagrange_polynomials(big_eval_points);
      // Evaluate lagrange_pol_eval at k+2,...,2k+1 to interpolate f and g at
      // those points
      lagrange_pol_eval_inner_prod_s.resize(instance.compression_factor + 1);
      for (size_t point = instance.compression_factor + 1;
           point < 2 * instance.compression_factor + 1; point++) {
        lagrange_pol_eval_inner_prod_s[point - instance.compression_factor - 1]
            .resize(instance.compression_factor + 1);
        for (size_t slice = 0; slice < instance.compression_factor + 1;
             slice++) {
          lagrange_pol_eval_inner_prod_s[point - instance.compression_factor -
                                         1][slice] =
              field::eval(lagrange_pol_eval[slice], big_eval_points[point]);
        }
      }

#ifdef DETAILED_IN_CODE_TIMERS
      stop_timers[8] = std::chrono::high_resolution_clock::now();
      time_spans[8] +=
          std::chrono::duration_cast<std::chrono::duration<double>>(
              stop_timers[8] - start_timers[8]);
#endif
    }

    for (size_t repetition = 0; repetition < instance.num_rounds;
         repetition++) {

      // get the plain value of input/outputs
      auto share_mult_left_lifted =
          rep_shared_mult_left_lifted.get_repetition(repetition);
      auto share_mult_right_lifted =
          rep_shared_mult_right_lifted.get_repetition(repetition);
      auto share_inner_prod_lifted =
          rep_shared_inner_prod_lifted.get_repetition(repetition);
#ifndef SAME_CHALLENGE
      std::vector<field::GF2E> left_input(instance.compression_factor *
                                          vector_size);
      std::vector<field::GF2E> right_input(instance.compression_factor *
                                           vector_size);
      std::vector<field::GF2E> inner_product_output(
          2 * instance.compression_factor - 1);
#endif

#ifdef SAME_CHALLENGE
      if (repetition == 0)
#else
      if (true)
#endif
      {
#ifdef DETAILED_IN_CODE_TIMERS
        start_timers[9] = std::chrono::high_resolution_clock::now();
#endif
        for (size_t party = 0; party < instance.num_MPC_parties; party++) {
          share_mult_left_lifted[party] =
              share_mult_left_lifted[party].first(current_m);
          share_mult_right_lifted[party] =
              share_mult_right_lifted[party].first(current_m);

          for (size_t idx = 0; idx < instance.compression_factor * vector_size;
               idx++) {
            if (idx < current_m) {
              left_input[idx] += share_mult_left_lifted[party][idx];
              right_input[idx] += share_mult_right_lifted[party][idx];
            } else {
              left_input[idx] += 0;
              right_input[idx] += 0;
            }
          }
          // Result of the inner product is stored at
          // instance.compression_factor - 1 from previous compression round
          inner_product_output[instance.compression_factor - 1] +=
              share_inner_prod_lifted[party][instance.compression_factor - 1];
        }
#ifdef DETAILED_IN_CODE_TIMERS
        stop_timers[9] = std::chrono::high_resolution_clock::now();
        time_spans[9] +=
            std::chrono::duration_cast<std::chrono::duration<double>>(
                stop_timers[9] - start_timers[9]);
#endif
      }

#ifdef DETAILED_IN_CODE_TIMERS
      start_timers[10] = std::chrono::high_resolution_clock::now();
#endif
      // Inject inner product of all but the last slices
      for (size_t slice = 0; slice < instance.compression_factor - 1; slice++) {
#ifdef SAME_CHALLENGE
        if (repetition == 0)
#else
        if (true)
#endif
        {
          for (size_t slice_pos = 0; slice_pos < vector_size; slice_pos++) {
            inner_product_output[slice] +=
                left_input[vector_size * slice + slice_pos] *
                right_input[vector_size * slice + slice_pos];
          }
          inner_product_output[instance.compression_factor - 1] -=
              inner_product_output[slice];
        }
        h_deltas[repetition][compression - 1][slice] =
            inner_product_output[slice];
      }

      // Create sharing of the injected inner product
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        for (size_t slice = 0; slice < instance.compression_factor - 1;
             slice++) {
          random_tapes[repetition][party].squeeze_bytes(
              lambda_sized_buffer.data(), lambda_sized_buffer.size());
          share_inner_prod_lifted[party][slice].from_bytes(
              lambda_sized_buffer.data());
          h_deltas[repetition][compression - 1][slice] -=
              share_inner_prod_lifted[party][slice];
        }
      }

      // fix sharing of the first party
      for (size_t slice = 0; slice < instance.compression_factor - 1; slice++) {
        share_inner_prod_lifted[0][slice] +=
            h_deltas[repetition][compression - 1][slice];
      }

      // all parties fix their last inner product share
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        for (size_t slice = 0; slice < instance.compression_factor - 1;
             slice++) {
          share_inner_prod_lifted[party][instance.compression_factor - 1] -=
              share_inner_prod_lifted[party][slice];
        }
      }
#ifdef DETAILED_IN_CODE_TIMERS
      stop_timers[10] = std::chrono::high_resolution_clock::now();
      time_spans[10] +=
          std::chrono::duration_cast<std::chrono::duration<double>>(
              stop_timers[10] - start_timers[10]);
#endif

// Polynomials f_j and g_j are implicitly defined in left_input and right_input
// (evaluation domain) Polynomial h is partially defined in the first
// compression_factor coefficients of inner_product_output (evaluation domain)
// Need to fill the last compression_factor-1 coefficients using the lagrange
// polynomials and f_j,g_j
#ifdef SAME_CHALLENGE
      if (repetition == 0 && compression != nbCompressions)
#else
      if (compression != nbCompressions)
#endif
      {
#ifdef DETAILED_IN_CODE_TIMERS
        start_timers[11] = std::chrono::high_resolution_clock::now();
#endif

        for (size_t i = 0; i < vector_size; i++) {
          // First we want to interpolate the ith polynomials f,g of degree
          // $instance.compression_factor - 1$
          for (size_t point = 0; point < instance.compression_factor; point++) {
            f_i_poly[point] = left_input[vector_size * point + i];
            g_i_poly[point] = right_input[vector_size * point + i];
          }

          // Evaluate those polynomials at k+1,...,2k-1
          for (size_t point = instance.compression_factor;
               point < 2 * instance.compression_factor - 1; point++) {
            // inner_product_output[point] +=
            // (field::eval(f_i_poly,big_eval_points[point])) *
            // (field::eval(g_i_poly,big_eval_points[point]));
            inner_product_output[point] +=
                dot_product(
                    lagrange_pol_eval_inner_prod_s[point -
                                                   instance.compression_factor],
                    f_i_poly) *
                dot_product(
                    lagrange_pol_eval_inner_prod_s[point -
                                                   instance.compression_factor],
                    g_i_poly);
          }
        }
#ifdef DETAILED_IN_CODE_TIMERS
        stop_timers[11] = std::chrono::high_resolution_clock::now();
        time_spans[11] +=
            std::chrono::duration_cast<std::chrono::duration<double>>(
                stop_timers[11] - start_timers[11]);
#endif
      }

      if (compression == nbCompressions) {
#ifdef DETAILED_IN_CODE_TIMERS
        start_timers[12] = std::chrono::high_resolution_clock::now();
#endif
        f_i_poly.resize(instance.compression_factor + 1);
        g_i_poly.resize(instance.compression_factor + 1);

        // for each repetition, each party samples an additional share of point
        // for the f_j and g_j polynomials
        add_point_f[repetition].resize(vector_size);
        add_point_g[repetition].resize(vector_size);
        for (size_t i = 0; i < vector_size; i++) {
          add_point_f[repetition][i].resize(instance.num_MPC_parties + 1);
          add_point_g[repetition][i].resize(instance.num_MPC_parties + 1);
          // first points are defined by the inputs to the compression round
          // (f_i(1),...,f_i(k)) same across all repetitions
          for (size_t point = 0; point < instance.compression_factor; point++) {
            f_i_poly[point] = left_input[vector_size * point + i];
            g_i_poly[point] = right_input[vector_size * point + i];
          }
          // for each repetition Prover samples a random share for each parties
          // for each polynomials : f_i(k+1)
          for (size_t party = 0; party < instance.num_MPC_parties; party++) {
            random_tapes[repetition][party].squeeze_bytes(
                lambda_sized_buffer.data(), lambda_sized_buffer.size());
            add_point_f[repetition][i][party].from_bytes(
                lambda_sized_buffer.data());
            add_point_f[repetition][i][instance.num_MPC_parties] +=
                add_point_f[repetition][i][party];

            random_tapes[repetition][party].squeeze_bytes(
                lambda_sized_buffer.data(), lambda_sized_buffer.size());
            add_point_g[repetition][i][party].from_bytes(
                lambda_sized_buffer.data());
            add_point_g[repetition][i][instance.num_MPC_parties] +=
                add_point_g[repetition][i][party];
          }
          f_i_poly[instance.compression_factor] =
              (add_point_f[repetition][i][instance.num_MPC_parties]);
          g_i_poly[instance.compression_factor] =
              (add_point_g[repetition][i][instance.num_MPC_parties]);
          // set the first k point of inner product as usual
          for (size_t point = 0; point < instance.compression_factor; point++) {
            inner_product_output_last_compression[repetition][point] =
                inner_product_output[point];
          }
          // point k+1 is the inner prod of the additional points
          inner_product_output_last_compression[repetition]
                                               [instance.compression_factor] +=
              f_i_poly[instance.compression_factor] *
              g_i_poly[instance.compression_factor];
          // Evaluate those polynomials at k+1,...,2k+1
          for (size_t point = instance.compression_factor + 1;
               point < 2 * instance.compression_factor + 1; point++) {
            // inner_product_output[point] +=
            // (field::eval(f_i_poly,big_eval_points[point])) *
            // (field::eval(g_i_poly,big_eval_points[point]));
            inner_product_output_last_compression[repetition][point] +=
                dot_product(lagrange_pol_eval_inner_prod_s
                                [point - instance.compression_factor - 1],
                            f_i_poly) *
                dot_product(lagrange_pol_eval_inner_prod_s
                                [point - instance.compression_factor - 1],
                            g_i_poly);
          }
        }

#ifdef DETAILED_IN_CODE_TIMERS
        stop_timers[12] = std::chrono::high_resolution_clock::now();
        time_spans[12] +=
            std::chrono::duration_cast<std::chrono::duration<double>>(
                stop_timers[12] - start_timers[12]);
#endif
      }

#ifdef DETAILED_IN_CODE_TIMERS
      start_timers[13] = std::chrono::high_resolution_clock::now();
#endif
      // set the deltas
      size_t up_bound = 2 * instance.compression_factor - 1;
      if (compression == nbCompressions) {
        up_bound = 2 * instance.compression_factor + 1;
        h_deltas[repetition][compression - 1].resize(up_bound);
      }
      for (size_t point = instance.compression_factor; point < up_bound;
           point++) {
        if (compression == nbCompressions) {
          h_deltas[repetition][compression - 1][point] =
              inner_product_output_last_compression[repetition][point];
        } else {
          h_deltas[repetition][compression - 1][point] =
              inner_product_output[point];
        }
      }

      // Sample sharing
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        for (size_t slice = instance.compression_factor; slice < up_bound;
             slice++) {
          random_tapes[repetition][party].squeeze_bytes(
              lambda_sized_buffer.data(), lambda_sized_buffer.size());
          share_inner_prod_lifted[party][slice].from_bytes(
              lambda_sized_buffer.data());
          h_deltas[repetition][compression - 1][slice] -=
              share_inner_prod_lifted[party][slice];
        }
      }
      // fix sharing of the first party
      for (size_t slice = instance.compression_factor; slice < up_bound;
           slice++) {
        share_inner_prod_lifted[0][slice] +=
            h_deltas[repetition][compression - 1][slice];
      }

#ifdef DETAILED_IN_CODE_TIMERS
      stop_timers[13] = std::chrono::high_resolution_clock::now();
      time_spans[13] +=
          std::chrono::duration_cast<std::chrono::duration<double>>(
              stop_timers[13] - start_timers[13]);
#endif
    }

    // Commit to the injections of the last points of polynomial h
    compression_commitments[compression - 1] =
        phase_2_commitment(instance, salt, h_1, compression_commitments,
                           compression - 1, h_deltas);

// expand to get the challenge point s (expecting one for now)
#ifdef SAME_CHALLENGE
    std::vector<field::GF2E> s_challenges = phase_2_expand(
        instance, compression_commitments[compression - 1], eval_points, 1);
#else
    std::vector<field::GF2E> s_challenges =
        phase_2_expand(instance, compression_commitments[compression - 1],
                       eval_points, instance.num_rounds);
#endif

#ifdef SAME_CHALLENGE
#ifdef DETAILED_IN_CODE_TIMERS
    start_timers[14] = std::chrono::high_resolution_clock::now();
#endif
    std::vector<field::GF2E> precomputed_lagrange_coefs_at_s(
        instance.compression_factor);
    std::vector<field::GF2E> precomputed_lagrange_coefs_at_s_big(
        2 * instance.compression_factor - 1);
    if (compression == nbCompressions) {
      precomputed_lagrange_coefs_at_s.resize(instance.compression_factor + 1);
      precomputed_lagrange_coefs_at_s_big.resize(
          2 * instance.compression_factor + 1);
    }

    for (size_t point = 0; point < lagrange_pol_eval.size(); point++) {
      precomputed_lagrange_coefs_at_s[point] =
          field::eval(lagrange_pol_eval[point], s_challenges[0]);
    }

    for (size_t point = 0; point < big_lagrange_pol_eval.size(); point++) {
      precomputed_lagrange_coefs_at_s_big[point] =
          field::eval(big_lagrange_pol_eval[point], s_challenges[0]);
    }
#ifdef DETAILED_IN_CODE_TIMERS
    stop_timers[14] = std::chrono::high_resolution_clock::now();
    time_spans[14] += std::chrono::duration_cast<std::chrono::duration<double>>(
        stop_timers[14] - start_timers[14]);
#endif
#endif

    // Compute the shared inner_product tuple output, which is the shared
    // polynomials evaluated at s
    for (size_t repetition = 0; repetition < instance.num_rounds;
         repetition++) {
#ifndef SAME_CHALLENGE
#ifdef DETAILED_IN_CODE_TIMERS
      start_timers[14] = std::chrono::high_resolution_clock::now();
#endif
      std::vector<field::GF2E> precomputed_lagrange_coefs_at_s(
          instance.compression_factor);
      std::vector<field::GF2E> precomputed_lagrange_coefs_at_s_big(
          2 * instance.compression_factor - 1);
      if (compression == nbCompressions) {
        precomputed_lagrange_coefs_at_s.resize(instance.compression_factor + 1);
        precomputed_lagrange_coefs_at_s_big.resize(
            2 * instance.compression_factor + 1);
      }

      for (size_t point = 0; point < lagrange_pol_eval.size(); point++) {
        precomputed_lagrange_coefs_at_s[point] =
            field::eval(lagrange_pol_eval[point], s_challenges[repetition]);
      }

      for (size_t point = 0; point < big_lagrange_pol_eval.size(); point++) {
        precomputed_lagrange_coefs_at_s_big[point] =
            field::eval(big_lagrange_pol_eval[point], s_challenges[repetition]);
      }
#ifdef DETAILED_IN_CODE_TIMERS
      stop_timers[14] = std::chrono::high_resolution_clock::now();
      time_spans[14] +=
          std::chrono::duration_cast<std::chrono::duration<double>>(
              stop_timers[14] - start_timers[14]);
#endif
#endif

      size_t up_bound = 2 * instance.compression_factor - 1;
      if (compression == nbCompressions) {
        f_final[repetition].resize(vector_size);
        g_final[repetition].resize(vector_size);
        up_bound = 2 * instance.compression_factor + 1;
        h_poly.resize(up_bound);
      }

      auto share_mult_left_lifted =
          rep_shared_mult_left_lifted.get_repetition(repetition);
      auto share_mult_right_lifted =
          rep_shared_mult_right_lifted.get_repetition(repetition);
      auto share_inner_prod_lifted =
          rep_shared_inner_prod_lifted.get_repetition(repetition);

      // each party has to compute its share of the inner product tuple
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
#ifdef DETAILED_IN_CODE_TIMERS
        start_timers[15] = std::chrono::high_resolution_clock::now();
#endif

        share_mult_left_lifted[party] =
            share_mult_left_lifted[party].first(current_m);
        share_mult_right_lifted[party] =
            share_mult_right_lifted[party].first(current_m);
        share_inner_prod_lifted[party] =
            share_inner_prod_lifted[party].first(up_bound);

        // interpolate the ith polynomials f,g of degree
        // $instance.compression_factor - 1$
        for (size_t i = 0; i < vector_size; i++) {
          for (size_t point = 0; point < instance.compression_factor; point++) {
            if (vector_size * point + i < current_m) {
              f_i_poly[point] =
                  share_mult_left_lifted[party][vector_size * point + i];
              g_i_poly[point] =
                  share_mult_right_lifted[party][vector_size * point + i];
            } else {
              f_i_poly[point] = 0;
              g_i_poly[point] = 0;
            }
          }
          if (compression == nbCompressions) {
            f_i_poly.resize(instance.compression_factor + 1);
            g_i_poly.resize(instance.compression_factor + 1);
            f_i_poly[instance.compression_factor] =
                (add_point_f[repetition][i][party]);
            g_i_poly[instance.compression_factor] =
                (add_point_g[repetition][i][party]);
          }
          // evaluate at s
          share_mult_left_lifted[party][i] =
              dot_product(precomputed_lagrange_coefs_at_s, f_i_poly);
          share_mult_right_lifted[party][i] =
              dot_product(precomputed_lagrange_coefs_at_s, g_i_poly);
        }
        // interpolate h of degree $2*instance.compression_factor - 2$
        for (size_t point = 0; point < up_bound; point++) {
          h_poly[point] = share_inner_prod_lifted[party][point];
        }
        // evaluate at s
        share_inner_prod_lifted[party][instance.compression_factor - 1] =
            dot_product(precomputed_lagrange_coefs_at_s_big, h_poly);

        if (compression == nbCompressions) {
          for (size_t idx = 0; idx < vector_size; idx++) {
            f_final[repetition][idx] += share_mult_left_lifted[party][idx];
            g_final[repetition][idx] += share_mult_right_lifted[party][idx];
          }
          h_final[repetition] +=
              share_inner_prod_lifted[party][instance.compression_factor - 1];
        }

#ifdef DETAILED_IN_CODE_TIMERS
        stop_timers[15] = std::chrono::high_resolution_clock::now();
        time_spans[15] +=
            std::chrono::duration_cast<std::chrono::duration<double>>(
                stop_timers[15] - start_timers[15]);
#endif
      }
    }

    if (compression == nbCompressions) {
      h_3 = phase_3_commitment(
          instance, salt, compression_commitments[compression - 1],
          rep_shared_mult_left_lifted, rep_shared_mult_right_lifted,
          rep_shared_inner_prod_lifted, vector_size);
      missing_parties = phase_3_expand(instance, h_3);
    }
  }

#ifdef IN_CODE_TIMERS
  tmp_time = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
      tmp_time - start_time);
  std::cout << "PHASE 3 TOOK " << time_span.count() << std::endl;
#endif

  // The proof contains:
  //-seeds of all but one parties
  //-commitment to the last seed
  //-witness deltas
  //-mult deltas
  //-h deltas
  //-the last inner product at s

#ifdef SANITY_CHECK
  // Check if we do have an inner product in the end for each repetition
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    field::GF2E accumulator;
    for (size_t vec_s = 0; vec_s < f_final[repetition].size(); vec_s++) {
      accumulator += (f_final[repetition][vec_s] * g_final[repetition][vec_s]);
    }
    if (accumulator != h_final[repetition]) {
      std::cout << "ERROR IN CHECKING INNER PRODUCT AT REPETITION "
                << repetition << std::endl;
      std::cout << "SIZE WAS " << f_final[repetition].size() << std::endl;
      std::cout << accumulator.data << std::endl;
      std::cout << h_final[repetition].data << std::endl;
    }
  }
#endif

  std::vector<reveal_list_t> seeds;
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    seeds.push_back(
        seed_trees[repetition].reveal_all_but(missing_parties[repetition]));
  }

  std::vector<limbo_repetition_proof_t> proofs;
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    size_t missing_party = missing_parties[repetition];
    limbo_repetition_proof_t proof{
        seeds[repetition],
        party_seed_commitments[repetition][missing_party],
        rep_witness_deltas[repetition],
        rep_mult_deltas[repetition],
        h_deltas[repetition],
        f_final[repetition],
        g_final[repetition],
        h_final[repetition]};
    proofs.push_back(proof);
  }
  limbo_proof_t proof{salt, h_1, compression_commitments, h_3, proofs};
  return proof;
}

bool limbo_verify(const limbo_instance_t &instance,
                    const std::vector<uint8_t> &witness_out,
                    const limbo_proof_t &proof,
                    const Circuit &C)
{

  // init modulus of extension field F_{2^{8\lambda}}
  size_t maxp64 = instance.num_MPC_parties / 65;
  field::GF2E::init_extension_field(instance);

  // buffer for squeezing field elements into
  std::vector<uint8_t> lambda_sized_buffer(instance.lambda);

  // do parallel repetitions
  // create seed trees and random tapes
  std::vector<SeedTree> seed_trees;
  std::vector<std::vector<RandomTape>> random_tapes;
  std::vector<std::vector<std::vector<uint8_t>>> party_seed_commitments;

  // Recompute all compression_commitments
  unsigned int nbCompressions =
      (unsigned int)(log(instance.nb_mult_gates) /
                     log(instance.compression_factor));
  std::vector<std::vector<std::vector<field::GF2E>>> h_deltas(
      instance.num_rounds);
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    h_deltas[repetition] = proof.proofs[repetition].h_deltas;
  }
  std::vector<std::vector<uint8_t>> compression_commitments(nbCompressions);
  for (unsigned int compression = 1; compression <= nbCompressions;
       compression++) {
    compression_commitments[compression - 1] =
        phase_2_commitment(instance, proof.salt, proof.h_1,
                           compression_commitments, compression - 1, h_deltas);
  }

// Compute challenges based on hashes
// h1 expansion
#ifdef SAME_CHALLENGE
  std::vector<field::GF2E> Rs = phase_1_expand(instance, proof.h_1, 1);
#else
  std::vector<field::GF2E> Rs =
      phase_1_expand(instance, proof.h_1, instance.num_rounds);
#endif

  // compression rounds expansion
  std::vector<std::vector<field::GF2E>> s_challenges(nbCompressions);
  std::vector<field::GF2E> eval_points =
      field::get_first_n_field_elements(instance.compression_factor);
  for (unsigned int compression = 1; compression <= nbCompressions;
       compression++) {
#ifdef SAME_CHALLENGE
    s_challenges[compression - 1] = phase_2_expand(
        instance, compression_commitments[compression - 1], eval_points, 1);
#else
    s_challenges[compression - 1] =
        phase_2_expand(instance, compression_commitments[compression - 1],
                       eval_points, instance.num_rounds);
#endif
  }

  // last expansion
  std::vector<uint16_t> missing_parties = phase_3_expand(instance, proof.h_3);

  // rebuild SeedTrees
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const limbo_repetition_proof_t &lproof = proof.proofs[repetition];
    // regenerate generate seed tree for the N parties (except the missing
    // one)
    if (missing_parties[repetition] != lproof.reveallist.second)
      throw std::runtime_error(
          "modified signature between deserialization and verify");
    seed_trees.push_back(SeedTree(lproof.reveallist, instance.num_MPC_parties,
                                  proof.salt, repetition));
    // commit to each party's seed, fill up missing one with data from proof
    std::vector<std::vector<uint8_t>> current_party_seed_commitments;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
        current_party_seed_commitments.push_back(commit_to_party_seed(
            instance, seed_trees[repetition].get_leaf(party).value(),
            proof.salt, repetition, party));
      } else {
        current_party_seed_commitments.push_back(
            lproof.missing_seed_commitments);
      }
    }
    party_seed_commitments.push_back(current_party_seed_commitments);

    // create random tape for each party, dummy one for missing party
    std::vector<RandomTape> party_tapes;
    party_tapes.reserve(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
        party_tapes.emplace_back(seed_trees[repetition].get_leaf(party).value(),
                                 proof.salt, repetition, party);
      } else {
        party_tapes.emplace_back(std::vector<uint8_t>(instance.seed_size),
                                 proof.salt, repetition, party);
      }
    }
    random_tapes.push_back(party_tapes);
  }
  /////////////////////////////////////////////////////////////////////////////
  // phase 1: recompute commitments to executions of circuit
  /////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_shared_witness(
      instance.num_rounds);
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_output_broadcasts(
      instance.num_rounds);
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_shared_mult_left(
      instance.num_rounds);
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_shared_mult_right(
      instance.num_rounds);
  std::vector<std::vector<std::array<uint64_t, PARTY64>>> rep_shared_mult_out(
      instance.num_rounds);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    rep_shared_witness[repetition].resize(instance.input_size);
    rep_output_broadcasts[repetition].resize(instance.output_size);
    rep_shared_mult_left[repetition].resize(instance.nb_mult_gates);
    rep_shared_mult_right[repetition].resize(instance.nb_mult_gates);
    rep_shared_mult_out[repetition].resize(instance.nb_mult_gates);

    // Retrieve proof for this repetition
    const limbo_repetition_proof_t &lproof = proof.proofs[repetition];

    // generate sharing of witness
    std::vector<uint8_t> shared_witness(instance.input_size >> 3);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // auto shared_witness = rep_shared_witness.get(repetition, party);
      random_tapes[repetition][party].squeeze_bytes(shared_witness.data(),
                                                    shared_witness.size());

      // set share of the player correctly
      for (size_t i = 0; i < shared_witness.size(); i++) {
        for (size_t j = 0; j < 8; j++) {
          rep_shared_witness[repetition][i * 8 + j][party >> 6] ^=
              (((uint64_t)((shared_witness[i] >> j) & 1)) << (party % 64));
        }
      }
    }

    // fix first share
    for (size_t i = 0; i < instance.input_size; i++) {
      rep_shared_witness[repetition][i][0] ^= lproof.witness_deltas[i];
    }

    // generate sharing of output of mult gates
    std::vector<uint8_t> shared_mult_out(instance.nb_mult_gates >> 3);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // auto shared_mult_out = rep_shared_mult_out.get(repetition, party);
      random_tapes[repetition][party].squeeze_bytes(shared_mult_out.data(),
                                                    shared_mult_out.size());

      // set share of the player correctly
      for (size_t i = 0; i < shared_mult_out.size() - 1; i++) {
        for (size_t j = 0; j < 8;
             j++) // std::max((long unsigned int)8,instance.nb_mult_gates -
                  // 8*i); j++)
        {
          rep_shared_mult_out[repetition][i * 8 + j][party >> 6] ^=
              (((uint64_t)((shared_mult_out[i] >> j) & 1)) << (party % 64));
        }
      }
      for (size_t j = 0; j < std::max((long unsigned int)8,
                                      instance.nb_mult_gates -
                                          8 * (shared_mult_out.size() - 1));
           j++) {
        rep_shared_mult_out[repetition][(shared_mult_out.size() - 1) * 8 + j]
                           [party >> 6] ^=
            (((uint64_t)((shared_mult_out[shared_mult_out.size() - 1] >> j) &
                         1))
             << (party));
      }
    }
    // fix first share
    for (size_t i = 0; i < instance.nb_mult_gates; i++) {
      rep_shared_mult_out[repetition][i][0] ^= lproof.mult_deltas[i];
    }
    // get shares of mult gates inputs by executing circuit in MPC

    C.base_circuit_shares(
        rep_shared_witness[repetition], rep_shared_mult_out[repetition],
        rep_output_broadcasts[repetition], rep_shared_mult_right[repetition],
        rep_shared_mult_left[repetition]);

    // calculate missing output share from others shares and witness_out
    size_t ccp = missing_parties[repetition] / 64;
    size_t cp = missing_parties[repetition] % 64;
    uint64_t mask_missing = (1ULL << cp);
    uint8_t currsharing;
    for (size_t i = 0; i < instance.output_size; i++) {
      // shareofmissing = ((rep_output_broadcasts[repetition][i][ccp] &
      // mask_missing) >> cp);
      currsharing ^= currsharing;
      for (size_t p64 = 0; p64 <= maxp64; p64++) {
        currsharing ^=
            (__builtin_popcountll(rep_output_broadcasts[repetition][i][p64]) %
             2);
      }
      // if needed flip the share of the missing party
      if (currsharing != witness_out[i]) {
        rep_output_broadcasts[repetition][i][ccp] ^= mask_missing;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 2Bis: lift mult triples and make it an inner product
  /////////////////////////////////////////////////////////////////////////////
  RepContainer<field::GF2E> rep_shared_mult_left_lifted(
      instance.num_rounds, instance.num_MPC_parties, instance.nb_mult_gates);
  RepContainer<field::GF2E> rep_shared_mult_right_lifted(
      instance.num_rounds, instance.num_MPC_parties, instance.nb_mult_gates);
  RepContainer<field::GF2E> rep_shared_inner_prod_lifted(
      instance.num_rounds, instance.num_MPC_parties,
      2 * instance.compression_factor + 1);
  std::vector<std::vector<field::GF2E>> Rs_powers =
      field::compute_powers(Rs, instance.nb_mult_gates);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // auto share_mult_left = rep_shared_mult_left.get(repetition, party);
      auto share_mult_left_lifted =
          rep_shared_mult_left_lifted.get(repetition, party);

      // auto share_mult_right = rep_shared_mult_right.get(repetition, party);
      auto share_mult_right_lifted =
          rep_shared_mult_right_lifted.get(repetition, party);

      // auto share_mult_out = rep_shared_mult_out.get(repetition, party);
      auto share_inner_prod_lifted =
          rep_shared_inner_prod_lifted.get(repetition, party);
      size_t cp = party % 64;
      size_t ccp = party / 64;
      uint64_t mask = (1ULL << cp);
      for (size_t nbm = 0; nbm < instance.nb_mult_gates; nbm++) {
        share_mult_right_lifted[nbm] = field::GF2E(
            (rep_shared_mult_right[repetition][nbm][ccp] & (mask)) >> cp);
#ifdef SAME_CHALLENGE
        share_mult_left_lifted[nbm] = field::GF2E(
            ((rep_shared_mult_left[repetition][nbm][ccp] & (mask)) >> cp),
            Rs_powers[0][nbm].data);
        share_inner_prod_lifted[instance.compression_factor - 1] += field::GF2E(
            (rep_shared_mult_out[repetition][nbm][ccp] & (mask)) >> cp,
            Rs_powers[0][nbm].data);
#else
        share_mult_left_lifted[nbm] = field::GF2E(
            ((rep_shared_mult_left[repetition][nbm][ccp] & (mask)) >> cp),
            Rs_powers[repetition][nbm].data);
        share_inner_prod_lifted[instance.compression_factor - 1] += field::GF2E(
            (rep_shared_mult_out[repetition][nbm][ccp] & (mask)) >> cp,
            Rs_powers[repetition][nbm].data);
#endif
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Start compression phases !
  /////////////////////////////////////////////////////////////////////////////

  // precompute lagrange pol
  std::vector<std::vector<field::GF2E>> lagrange_pol_eval =
      field::precompute_lagrange_polynomials(eval_points);

  // Get 2*k - 1 points for h
  std::vector<field::GF2E> big_eval_points =
      field::get_first_n_field_elements(2 * instance.compression_factor - 1);

  // precompute lagrange pol
  std::vector<std::vector<field::GF2E>> big_lagrange_pol_eval =
      field::precompute_lagrange_polynomials(big_eval_points);

  // additional point sampled at the last compression
  std::vector<std::vector<std::vector<field::GF2E>>> add_point_f(
      instance.num_rounds);
  std::vector<std::vector<std::vector<field::GF2E>>> add_point_g(
      instance.num_rounds);
  std::vector<field::GF2E> f_i_poly(instance.compression_factor);
  std::vector<field::GF2E> g_i_poly(instance.compression_factor);

  unsigned int current_m = instance.nb_mult_gates;
  unsigned int vector_size = current_m;
  for (unsigned int compression = 1; compression <= nbCompressions;
       compression++) {
    current_m = vector_size;
    vector_size =
        (unsigned int)ceil(((double)current_m / instance.compression_factor));

    // std::cout << "Doing compression nb " << compression << " with vectors of
    // sizes " << vector_size << " requires padding " <<
    // (instance.compression_factor * vector_size) - current_m << std::endl;
    std::vector<field::GF2E> left_input(instance.compression_factor *
                                        vector_size);
    std::vector<field::GF2E> right_input(instance.compression_factor *
                                         vector_size);
    std::vector<field::GF2E> inner_product_output(
        2 * instance.compression_factor - 1);

    if (compression == nbCompressions) {
      // Get k+1 points for f,g
      eval_points =
          field::get_first_n_field_elements(instance.compression_factor + 1);
      // Get 2*k+1 points for h
      big_eval_points = field::get_first_n_field_elements(
          2 * instance.compression_factor + 1);
      // precompute lagrange pol
      lagrange_pol_eval = field::precompute_lagrange_polynomials(eval_points);
      big_lagrange_pol_eval =
          field::precompute_lagrange_polynomials(big_eval_points);
    }

    for (size_t repetition = 0; repetition < instance.num_rounds;
         repetition++) {
      auto share_mult_left_lifted =
          rep_shared_mult_left_lifted.get_repetition(repetition);
      auto share_mult_right_lifted =
          rep_shared_mult_right_lifted.get_repetition(repetition);
      auto share_inner_prod_lifted =
          rep_shared_inner_prod_lifted.get_repetition(repetition);

      // create sharing of the injected inner_product
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        for (size_t slice = 0; slice < instance.compression_factor - 1;
             slice++) {
          random_tapes[repetition][party].squeeze_bytes(
              lambda_sized_buffer.data(), lambda_sized_buffer.size());
          share_inner_prod_lifted[party][slice].from_bytes(
              lambda_sized_buffer.data());
        }
      }

      // fix sharing of the first party
      for (size_t slice = 0; slice < instance.compression_factor - 1; slice++) {
        share_inner_prod_lifted[0][slice] +=
            h_deltas[repetition][compression - 1][slice];
      }

      // all parties fix their last inner product share
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        for (size_t slice = 0; slice < instance.compression_factor - 1;
             slice++) {
          share_inner_prod_lifted[party][instance.compression_factor - 1] -=
              share_inner_prod_lifted[party][slice];
        }
      }

      size_t up_bound = 2 * instance.compression_factor - 1;
      // if last compression, sample also the additional point
      if (compression == nbCompressions) {
        up_bound = 2 * instance.compression_factor + 1;
        add_point_f[repetition].resize(vector_size);
        add_point_g[repetition].resize(vector_size);
        for (size_t i = 0; i < vector_size; i++) {
          add_point_f[repetition][i].resize(instance.num_MPC_parties + 1);
          add_point_g[repetition][i].resize(instance.num_MPC_parties + 1);
          for (size_t party = 0; party < instance.num_MPC_parties; party++) {
            random_tapes[repetition][party].squeeze_bytes(
                lambda_sized_buffer.data(), lambda_sized_buffer.size());
            add_point_f[repetition][i][party].from_bytes(
                lambda_sized_buffer.data());
            add_point_f[repetition][i][instance.num_MPC_parties] +=
                add_point_f[repetition][i][party];

            random_tapes[repetition][party].squeeze_bytes(
                lambda_sized_buffer.data(), lambda_sized_buffer.size());
            add_point_g[repetition][i][party].from_bytes(
                lambda_sized_buffer.data());
            add_point_g[repetition][i][instance.num_MPC_parties] +=
                add_point_g[repetition][i][party];
          }
        }
      }

      // Sample parties share for the last compression_factor-1 evaluation
      // points
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        for (size_t slice = instance.compression_factor; slice < up_bound;
             slice++) {
          random_tapes[repetition][party].squeeze_bytes(
              lambda_sized_buffer.data(), lambda_sized_buffer.size());
          share_inner_prod_lifted[party][slice].from_bytes(
              lambda_sized_buffer.data());
        }
      }
      // fix sharing of the first party
      for (size_t slice = instance.compression_factor; slice < up_bound;
           slice++) {
        share_inner_prod_lifted[0][slice] +=
            h_deltas[repetition][compression - 1][slice];
      }
    }
#ifdef SAME_CHALLENGE
    std::vector<field::GF2E> precomputed_lagrange_coefs_at_s(
        instance.compression_factor);
    std::vector<field::GF2E> precomputed_lagrange_coefs_at_s_big(
        2 * instance.compression_factor - 1);
    size_t up_bound = 2 * instance.compression_factor - 1;
    if (compression == nbCompressions) {
      precomputed_lagrange_coefs_at_s.resize(instance.compression_factor + 1);
      precomputed_lagrange_coefs_at_s_big.resize(
          2 * instance.compression_factor + 1);
      up_bound = 2 * instance.compression_factor + 1;
    }

    for (size_t point = 0; point < lagrange_pol_eval.size(); point++) {
      precomputed_lagrange_coefs_at_s[point] = field::eval(
          lagrange_pol_eval[point], s_challenges[compression - 1][0]);
    }

    for (size_t point = 0; point < big_lagrange_pol_eval.size(); point++) {
      precomputed_lagrange_coefs_at_s_big[point] = field::eval(
          big_lagrange_pol_eval[point], s_challenges[compression - 1][0]);
    }
#endif

    for (size_t repetition = 0; repetition < instance.num_rounds;
         repetition++) {
#ifndef SAME_CHALLENGE
      std::vector<field::GF2E> precomputed_lagrange_coefs_at_s(
          instance.compression_factor);
      std::vector<field::GF2E> precomputed_lagrange_coefs_at_s_big(
          2 * instance.compression_factor - 1);
      size_t up_bound = 2 * instance.compression_factor - 1;
      if (compression == nbCompressions) {
        precomputed_lagrange_coefs_at_s.resize(instance.compression_factor + 1);
        precomputed_lagrange_coefs_at_s_big.resize(
            2 * instance.compression_factor + 1);
        up_bound = 2 * instance.compression_factor + 1;
      }

      for (size_t point = 0; point < lagrange_pol_eval.size(); point++) {
        precomputed_lagrange_coefs_at_s[point] =
            field::eval(lagrange_pol_eval[point],
                        s_challenges[compression - 1][repetition]);
      }

      for (size_t point = 0; point < big_lagrange_pol_eval.size(); point++) {
        precomputed_lagrange_coefs_at_s_big[point] =
            field::eval(big_lagrange_pol_eval[point],
                        s_challenges[compression - 1][repetition]);
      }
#endif
      auto share_mult_left_lifted =
          rep_shared_mult_left_lifted.get_repetition(repetition);
      auto share_mult_right_lifted =
          rep_shared_mult_right_lifted.get_repetition(repetition);
      auto share_inner_prod_lifted =
          rep_shared_inner_prod_lifted.get_repetition(repetition);

      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        share_mult_left_lifted[party] =
            share_mult_left_lifted[party].first(current_m);
        share_mult_right_lifted[party] =
            share_mult_right_lifted[party].first(current_m);
        share_inner_prod_lifted[party] =
            share_inner_prod_lifted[party].first(up_bound);

        // interpolate the ith polynomials f,g of degree
        // $instance.compression_factor - 1$
        for (size_t i = 0; i < vector_size; i++) {
          for (size_t point = 0; point < instance.compression_factor; point++) {
            if (vector_size * point + i < current_m) {
              f_i_poly[point] =
                  share_mult_left_lifted[party][vector_size * point + i];
              g_i_poly[point] =
                  share_mult_right_lifted[party][vector_size * point + i];
            } else {
              f_i_poly[point] = 0;
              g_i_poly[point] = 0;
            }
          }
          if (compression == nbCompressions) {
            f_i_poly.resize(instance.compression_factor + 1);
            g_i_poly.resize(instance.compression_factor + 1);
            f_i_poly[instance.compression_factor] =
                (add_point_f[repetition][i][party]);
            g_i_poly[instance.compression_factor] =
                (add_point_g[repetition][i][party]);
          }

          // evaluate at s
          share_mult_left_lifted[party][i] =
              dot_product(precomputed_lagrange_coefs_at_s, f_i_poly);
          share_mult_right_lifted[party][i] =
              dot_product(precomputed_lagrange_coefs_at_s, g_i_poly);
        }
        // interpolate h of degree $2*instance.compression_factor - 2$
        std::vector<field::GF2E> h_poly(up_bound);
        for (size_t point = 0; point < up_bound; point++) {
          h_poly[point] = share_inner_prod_lifted[party][point];
        }
        // evaluate at s
        // share_inner_prod_lifted[party][instance.compression_factor - 1] =
        // field::eval(h_poly, s_challenges[compression - 1][0]);
        share_inner_prod_lifted[party][instance.compression_factor - 1] =
            dot_product(precomputed_lagrange_coefs_at_s_big, h_poly);
      }

      // at the end of the last compression compute the missing broadcast share
      if (compression == nbCompressions) {
        const limbo_repetition_proof_t &lproof = proof.proofs[repetition];
        std::copy(lproof.f_final.begin(), lproof.f_final.end(),
                  share_mult_left_lifted[missing_parties[repetition]].begin());
        std::copy(lproof.g_final.begin(), lproof.g_final.end(),
                  share_mult_right_lifted[missing_parties[repetition]].begin());
        share_inner_prod_lifted[missing_parties[repetition]]
                               [instance.compression_factor - 1] =
                                   lproof.h_final;

        for (size_t party = 0; party < instance.num_MPC_parties; party++) {
          if (party != missing_parties[repetition]) {
            for (size_t idx = 0; idx < vector_size; idx++) {
              if (idx < current_m) {
                share_mult_left_lifted[missing_parties[repetition]][idx] -=
                    share_mult_left_lifted[party][idx];
                share_mult_right_lifted[missing_parties[repetition]][idx] -=
                    share_mult_right_lifted[party][idx];
              }
            }
            share_inner_prod_lifted[missing_parties[repetition]]
                                   [instance.compression_factor - 1] -=
                share_inner_prod_lifted[party][instance.compression_factor - 1];
          }
        }
      }
    }
  }

  // Recompute h_1 and h_3
  std::vector<std::vector<uint8_t>> rep_witness_deltas;
  std::vector<std::vector<uint8_t>> rep_mult_deltas;
  for (const limbo_repetition_proof_t &lproof : proof.proofs) {
    rep_witness_deltas.push_back(lproof.witness_deltas);
    rep_mult_deltas.push_back(lproof.mult_deltas);
  }

  std::vector<uint8_t> h_1 = phase_1_commitment(
      instance, proof.salt, witness_out, party_seed_commitments,
      rep_witness_deltas, rep_mult_deltas, rep_output_broadcasts);

  std::vector<uint8_t> h_3 = phase_3_commitment(
      instance, proof.salt, compression_commitments[nbCompressions - 1],
      rep_shared_mult_left_lifted, rep_shared_mult_right_lifted,
      rep_shared_inner_prod_lifted, vector_size);

  // Check that these are ok
  if (memcmp(h_1.data(), proof.h_1.data(), h_1.size()) != 0) {
    std::cout << "ERROR IN COMPUTING h1 !" << std::endl;
    return false;
  }
  if (memcmp(h_3.data(), proof.h_3.data(), h_3.size()) != 0) {
    std::cout << "ERROR IN COMPUTING h3 !" << std::endl;
    return false;
  }

  // Check if we do have an inner product in the end for each repetition
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    field::GF2E accumulator;
    for (size_t vec_s = 0; vec_s < proof.proofs[repetition].f_final.size();
         vec_s++) {
      accumulator += (proof.proofs[repetition].f_final[vec_s] *
                      proof.proofs[repetition].g_final[vec_s]);
    }
    if (accumulator != proof.proofs[repetition].h_final) {
      std::cout << "ERROR IN CHECKING INNER PRODUCT AT REPETITION "
                << repetition << std::endl;
      std::cout << "SIZE WAS " << proof.proofs[repetition].f_final.size()
                << std::endl;
      return false;
    }
  }

  return true;
}

// compute size of the proof in bits
unsigned int compute_proof_size(limbo_proof_t proof,
                                limbo_instance_t instance) {
  unsigned int result;
  result = SALT_SIZE * 8;
  result += proof.h_1.size() * 8;
  for (unsigned int i = 0; i < proof.compression_commitments.size(); i++) {
    result += proof.compression_commitments[i].size() * 8;
  }
  result += proof.h_3.size() * 8;

  for (unsigned int i = 0; i < proof.proofs.size(); i++) {
    result += proof.proofs[i].reveallist.first.size() * instance.seed_size;
    result += proof.proofs[i].missing_seed_commitments.size() * 8;
    result += proof.proofs[i].witness_deltas.size();
    result += proof.proofs[i].mult_deltas.size();
    for (unsigned int j = 0; j < proof.proofs[i].h_deltas.size(); j++) {
      result += proof.proofs[i].h_deltas[j].size() * (8 * instance.lambda);
    }
    result += proof.proofs[i].f_final.size() * (8 * instance.lambda);
    result += proof.proofs[i].g_final.size() * (8 * instance.lambda);
    result += 8 * instance.lambda;
  }

  return result;
}

//  /////////////////////////////////////////////////////////////////////////////
//  // recompute shares of polynomials
//  /////////////////////////////////////////////////////////////////////////////
//  // a vector of the first m2+1 field elements for interpolation
//  std::vector<field::GF2E> x_values_for_interpolation_zero_to_m2 =
//      field::get_first_n_field_elements(instance.m2 + 1);
//
//  std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_m2 =
//      field::precompute_lagrange_polynomials(
//          x_values_for_interpolation_zero_to_m2);
//
//  std::vector<field::GF2E> x_values_for_interpolation_zero_to_2m2 =
//      field::get_first_n_field_elements(2 * instance.m2 + 1);
//
//  std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_2m2 =
//      field::precompute_lagrange_polynomials(
//          x_values_for_interpolation_zero_to_2m2);
//
//  std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> s_prime(
//      instance.num_rounds);
//  std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> t_prime(
//      instance.num_rounds);
//  std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> z_prime(
//      instance.num_rounds);
//
//  std::vector<std::vector<field::GF2E>> P_e(instance.num_rounds);
//  std::vector<std::vector<std::vector<field::GF2E>>> P_e_shares(
//      instance.num_rounds);
//
//  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++)
//  {
//    const limbo_repetition_proof_t &proof = signature.proofs[repetition];
//    // S_eji[repetition].resize(instance.num_MPC_parties);
//    // T_eji[repetition].resize(instance.num_MPC_parties);
//    s_prime[repetition].resize(instance.num_MPC_parties);
//    t_prime[repetition].resize(instance.num_MPC_parties);
//    z_prime[repetition].resize(instance.num_MPC_parties);
//    // P_ei[repetition].resize(instance.num_MPC_parties);
//    P_deltas[repetition].resize(instance.m2 + 1);
//
//    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
//      if (party != missing_parties[repetition]) {
//        // S_eji[repetition][party].resize(instance.m1);
//        // T_eji[repetition][party].resize(instance.m1);
//        s_prime[repetition][party].resize(instance.m1);
//        t_prime[repetition][party].resize(instance.m1);
//        z_prime[repetition][party].resize(instance.m1);
//        // lift shares from F_{2^8} to F_{2^{8\lambda}}
//        std::vector<std::reference_wrapper<const field::GF2E>> lifted_s;
//        std::vector<std::reference_wrapper<const field::GF2E>> lifted_t;
//        std::vector<std::reference_wrapper<const field::GF2E>> lifted_z;
//        lifted_s.reserve(instance.nb_mult_gates);
//        lifted_t.reserve(instance.nb_mult_gates);
//        lifted_z.reserve(instance.nb_mult_gates);
//        auto shared_s = rep_shared_mult_right.get(repetition, party);
//        auto shared_t = rep_shared_mult_left.get(repetition, party);
//        auto shared_z = rep_shared_mult_out.get(repetition, party);
//        for (size_t idx = 0; idx < instance.nb_mult_gates; idx++) {
//          lifted_s.push_back(field::lift_uint8_t(shared_s[idx]));
//          lifted_t.push_back(field::lift_uint8_t(shared_t[idx]));
//          lifted_z.push_back(field::lift_uint8_t(shared_z[idx]));
//        }
//
//        // rearrange shares
//        std::vector<field::GF2E> s_bar(instance.m2 + 1);
//        std::vector<field::GF2E> t_bar(instance.m2 + 1);
//        std::vector<field::GF2E> z_bar(instance.m2 + 1);
//
//        for (size_t j = 0; j < instance.m1; j++) {
//          for (size_t k = 0; k < instance.m2; k++) {
//            s_bar[k] = r_ejs[repetition][j] * lifted_s[j + instance.m1 * k];
//            t_bar[k] = lifted_t[j + instance.m1 * k];
//            z_bar[k] = lifted_z[j + instance.m1 * k];
//          }
//
//          // sample additional random points
//          random_tapes[repetition][party].squeeze_bytes(
//              lambda_sized_buffer.data(), lambda_sized_buffer.size());
//          s_bar[instance.m2].from_bytes(lambda_sized_buffer.data());
//
//          random_tapes[repetition][party].squeeze_bytes(
//              lambda_sized_buffer.data(), lambda_sized_buffer.size());
//          t_bar[instance.m2].from_bytes(lambda_sized_buffer.data());
//
//          // interpolate polynomials S_ej^i and T_ej^i
//          // S_eji[repetition][party][j] =
//          // utils::interpolate_with_precomputation(
//          // precomputation_for_zero_to_m2, s_bar);
//          // T_eji[repetition][party][j] =
//          // utils::interpolate_with_precomputation(
//          // precomputation_for_zero_to_m2, t_bar);
//          s_prime[repetition][party][j] = s_bar;
//          t_prime[repetition][party][j] = t_bar;
//          z_prime[repetition][party][j] = z_bar;
//        }
//      }
//    }
//
//    // compute sharing of P
//    std::vector<std::vector<field::GF2E>> &P_shares = P_e_shares[repetition];
//    P_shares.resize(instance.num_MPC_parties);
//    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
//      if (party != missing_parties[repetition]) {
//        // first m2 points: first party = sum of r_e,j, other parties = 0
//        P_shares[party].resize(2 * instance.m2 + 1);
//        //if (party == 0) {
//          field::GF2E sum_r;
//          //for (size_t j = 0; j < instance.m1; j++) {
//          //  sum_r += r_ejs[repetition][j] * z_prime[repetition][party][j];
//          //}
//          for (size_t k = 0; k < instance.m2; k++) {
//	     field::GF2E sum_r_mult;
//	     for(size_t j = 0; j < instance.m1; j++)
//	     {
//                sum_r_mult += r_ejs[repetition][j] *
//                z_prime[repetition][party][j][k];
//	     }
//            P_shares[party][k] = sum_r_mult;
//          }
//        //} else {
//        //  for (size_t k = 0; k < instance.m2; k++) {
//        //    P_shares[party][k] = field::GF2E(0);
//        //  }
//        //}
//
//        // second m2+1 points: sample from random tape
//        for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
//          random_tapes[repetition][party].squeeze_bytes(
//              lambda_sized_buffer.data(), lambda_sized_buffer.size());
//          P_shares[party][k].from_bytes(lambda_sized_buffer.data());
//        }
//      }
//    }
//    if (0 != missing_parties[repetition]) {
//      for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
//        // adjust first share with delta from signature
//        P_shares[0][k] += proof.P_delta[k - instance.m2];
//      }
//    }
//    // for (size_t party = 0; party < instance.num_MPC_parties; party++) {
//    //// iterpolate polynomial P_e^1 from 2m+1 points
//    // if (party != missing_parties[repetition]) {
//    // P_ei[repetition][party] = utils::interpolate_with_precomputation(
//    // precomputation_for_zero_to_2m2, P_shares[party]);
//    //}
//    //}
//  }
//
//  /////////////////////////////////////////////////////////////////////////////
//  // recompute views of polynomial checks
//  /////////////////////////////////////////////////////////////////////////////
//  std::vector<field::GF2E> c(instance.num_rounds);
//  std::vector<std::vector<field::GF2E>> c_shares(instance.num_rounds);
//  std::vector<std::vector<field::GF2E>> a(instance.num_rounds);
//  std::vector<std::vector<field::GF2E>> b(instance.num_rounds);
//  RepContainer<field::GF2E> a_shares(instance.num_rounds,
//                                     instance.num_MPC_parties, instance.m1);
//  RepContainer<field::GF2E> b_shares(instance.num_rounds,
//                                     instance.num_MPC_parties, instance.m1);
//
//  std::vector<field::GF2E> lagrange_polys_evaluated_at_Re_m2(instance.m2 + 1);
//  std::vector<field::GF2E> lagrange_polys_evaluated_at_Re_2m2(2 * instance.m2
//  +
//                                                              1);
//  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++)
//  {
//    const limbo_repetition_proof_t &proof = signature.proofs[repetition];
//    size_t missing_party = missing_parties[repetition];
//
//    for (size_t k = 0; k < instance.m2 + 1; k++) {
//      lagrange_polys_evaluated_at_Re_m2[k] =
//          field::eval(precomputation_for_zero_to_m2[k], R_es[repetition]);
//    }
//    for (size_t k = 0; k < 2 * instance.m2 + 1; k++) {
//      lagrange_polys_evaluated_at_Re_2m2[k] =
//          field::eval(precomputation_for_zero_to_2m2[k], R_es[repetition]);
//    }
//
//    c_shares[repetition].resize(instance.num_MPC_parties);
//    a[repetition].resize(instance.m1);
//    b[repetition].resize(instance.m1);
//    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
//      if (party != missing_party) {
//        auto a_shares_party = a_shares.get(repetition, party);
//        auto b_shares_party = b_shares.get(repetition, party);
//        for (size_t j = 0; j < instance.m1; j++) {
//          // compute a_ej^i and b_ej^i
//          //  a_shares[repetition][party][j] =
//          //  eval(S_eji[repetition][party][j], R_es[repetition]);
//          //  b_shares[repetition][party][j] =
//          //  eval(T_eji[repetition][party][j], R_es[repetition]);
//          a_shares_party[j] = dot_product(lagrange_polys_evaluated_at_Re_m2,
//                                          s_prime[repetition][party][j]);
//          b_shares_party[j] = dot_product(lagrange_polys_evaluated_at_Re_m2,
//                                          t_prime[repetition][party][j]);
//        }
//        // compute c_e^i
//        // c_shares[repetition][party] =
//        // eval(P_ei[repetition][party], R_es[repetition]);
//        c_shares[repetition][party] = dot_product(
//            lagrange_polys_evaluated_at_Re_2m2,
//            P_e_shares[repetition][party]);
//      }
//    }
//
//    // calculate missing shares
//    c[repetition] = proof.P_at_R;
//    c_shares[repetition][missing_party] = proof.P_at_R;
//    auto a_shares_missing = a_shares.get(repetition, missing_party);
//    auto b_shares_missing = b_shares.get(repetition, missing_party);
//    for (size_t j = 0; j < instance.m1; j++) {
//      a[repetition][j] = proof.S_j_at_R[j];
//      a_shares_missing[j] = proof.S_j_at_R[j];
//      b[repetition][j] = proof.T_j_at_R[j];
//      b_shares_missing[j] = proof.T_j_at_R[j];
//    }
//    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
//      if (party != missing_party) {
//        c_shares[repetition][missing_party] -= c_shares[repetition][party];
//        auto a_shares_party = a_shares.get(repetition, party);
//        auto b_shares_party = b_shares.get(repetition, party);
//        for (size_t j = 0; j < instance.m1; j++) {
//          a_shares_missing[j] -= a_shares_party[j];
//          b_shares_missing[j] -= b_shares_party[j];
//        }
//      }
//    }
//  }
//  /////////////////////////////////////////////////////////////////////////////
//  // recompute h_1 and h_3
//  /////////////////////////////////////////////////////////////////////////////
//  std::vector<std::vector<uint8_t>> sk_deltas;
//  std::vector<std::vector<uint8_t>> t_deltas;
//  for (const limbo_repetition_proof_t &proof : signature.proofs) {
//    sk_deltas.push_back(proof.sk_delta);
//    t_deltas.push_back(proof.t_delta);
//  }
//  std::vector<uint8_t> h_1 = phase_1_commitment(
//      instance, signature.salt, witness_out,
//      party_seed_commitments, sk_deltas, t_deltas, rep_output_broadcasts);
//
//  std::vector<uint8_t> h_3 = phase_3_commitment(
//      instance, signature.salt, h_2, c, c_shares, a, a_shares, b, b_shares);
//  // do checks
//  if (memcmp(h_1.data(), signature.h_1.data(), h_1.size()) != 0) {
//    return false;
//  }
//  if (memcmp(h_3.data(), signature.h_3.data(), h_3.size()) != 0) {
//    return false;
//  }
//
//  // check if P_e(R) = Sum_j S_e_j(R) * T_e_j(R) for all e
//  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++)
//  {
//    field::GF2E accum;
//    for (size_t j = 0; j < instance.m1; j++) {
//      accum += signature.proofs[repetition].S_j_at_R[j] *
//               signature.proofs[repetition].T_j_at_R[j];
//    }
//    if (accum != signature.proofs[repetition].P_at_R) {
//      return false;
//    }
//  }
//  return true;
//}
//
// std::vector<uint8_t>
// limbo_serialize_proof(const limbo_instance_t &instance,
//                            const limbo_proof_t &signature) {
//  std::vector<uint8_t> serialized;
//
//  serialized.insert(serialized.end(), signature.salt.begin(),
//                    signature.salt.end());
//  serialized.insert(serialized.end(), signature.h_1.begin(),
//                    signature.h_1.end());
//  serialized.insert(serialized.end(), signature.h_3.begin(),
//                    signature.h_3.end());
//
//  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++)
//  {
//    const limbo_repetition_proof_t &proof = signature.proofs[repetition];
//    for (const std::vector<uint8_t> &seed : proof.reveallist.first) {
//      serialized.insert(serialized.end(), seed.begin(), seed.end());
//    }
//    serialized.insert(serialized.end(), proof.C_e.begin(), proof.C_e.end());
//    serialized.insert(serialized.end(), proof.sk_delta.begin(),
//                      proof.sk_delta.end());
//    serialized.insert(serialized.end(), proof.t_delta.begin(),
//                      proof.t_delta.end());
//    for (size_t k = 0; k < instance.m2 + 1; k++) {
//      std::vector<uint8_t> buffer = proof.P_delta[k].to_bytes();
//      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
//    }
//    {
//      std::vector<uint8_t> buffer = proof.P_at_R.to_bytes();
//      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
//    }
//    for (size_t j = 0; j < instance.m1; j++) {
//      std::vector<uint8_t> buffer = proof.S_j_at_R[j].to_bytes();
//      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
//    }
//    for (size_t j = 0; j < instance.m1; j++) {
//      std::vector<uint8_t> buffer = proof.T_j_at_R[j].to_bytes();
//      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
//    }
//  }
//  return serialized;
//}
//
// limbo_proof_t
// limbo_deserialize_proof(const limbo_instance_t &instance,
//                              const std::vector<uint8_t> &serialized) {
//
//  size_t current_offset = 0;
//  limbo_salt_t salt;
//  memcpy(salt.data(), serialized.data() + current_offset, salt.size());
//  current_offset += salt.size();
//  std::vector<uint8_t> h_1(instance.digest_size), h_3(instance.digest_size);
//  memcpy(h_1.data(), serialized.data() + current_offset, h_1.size());
//  current_offset += h_1.size();
//  memcpy(h_3.data(), serialized.data() + current_offset, h_3.size());
//  current_offset += h_3.size();
//  std::vector<limbo_repetition_proof_t> proofs;
//
//  std::vector<uint16_t> missing_parties = phase_3_expand(instance, h_3);
//  size_t reveallist_size = ceil_log2(instance.num_MPC_parties);
//  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++)
//  {
//    reveal_list_t reveallist;
//    reveallist.first.reserve(reveallist_size);
//    reveallist.second = missing_parties[repetition];
//    for (size_t i = 0; i < reveallist_size; i++) {
//      std::vector<uint8_t> seed(instance.seed_size);
//      memcpy(seed.data(), serialized.data() + current_offset, seed.size());
//      current_offset += seed.size();
//      reveallist.first.push_back(seed);
//    }
//    std::vector<uint8_t> C_e(instance.digest_size);
//    memcpy(C_e.data(), serialized.data() + current_offset, C_e.size());
//    current_offset += C_e.size();
//
//    std::vector<uint8_t> sk_delta(instance.input_size);
//    memcpy(sk_delta.data(), serialized.data() + current_offset,
//           sk_delta.size());
//    current_offset += sk_delta.size();
//
//    std::vector<uint8_t> t_delta(instance.nb_mult_gates);
//    memcpy(t_delta.data(), serialized.data() + current_offset,
//    t_delta.size()); current_offset += t_delta.size();
//
//    field::GF2E tmp;
//    std::vector<field::GF2E> P_delta;
//    P_delta.reserve(instance.m2 + 1);
//    for (size_t k = 0; k < instance.m2 + 1; k++) {
//      std::vector<uint8_t> buffer(instance.lambda);
//      memcpy(buffer.data(), serialized.data() + current_offset,
//      buffer.size()); current_offset += buffer.size();
//      tmp.from_bytes(buffer.data());
//      P_delta.push_back(tmp);
//    }
//    field::GF2E P_at_R;
//    {
//      std::vector<uint8_t> buffer(instance.lambda);
//      memcpy(buffer.data(), serialized.data() + current_offset,
//      buffer.size()); current_offset += buffer.size();
//      P_at_R.from_bytes(buffer.data());
//    }
//    std::vector<field::GF2E> S_j_at_R;
//    S_j_at_R.reserve(instance.m1);
//    for (size_t j = 0; j < instance.m1; j++) {
//      std::vector<uint8_t> buffer(instance.lambda);
//      memcpy(buffer.data(), serialized.data() + current_offset,
//      buffer.size()); current_offset += buffer.size();
//      tmp.from_bytes(buffer.data());
//      S_j_at_R.push_back(tmp);
//    }
//    std::vector<field::GF2E> T_j_at_R;
//    T_j_at_R.reserve(instance.m1);
//    for (size_t j = 0; j < instance.m1; j++) {
//      std::vector<uint8_t> buffer(instance.lambda);
//      memcpy(buffer.data(), serialized.data() + current_offset,
//      buffer.size()); current_offset += buffer.size();
//      tmp.from_bytes(buffer.data());
//      T_j_at_R.push_back(tmp);
//    }
//    proofs.emplace_back(limbo_repetition_proof_t{reveallist, C_e, sk_delta,
//                                                   t_delta, P_delta, P_at_R,
//                                                   S_j_at_R, T_j_at_R});
//  }
//  assert(current_offset == serialized.size());
//  limbo_proof_t signature{salt, h_1, h_3, proofs};
//  return signature;
//}
