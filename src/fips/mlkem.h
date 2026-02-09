/* Copyright 2026 (c) Yao Zeran, Zhang Chenzhi
 * 
 * This file defines mlkem core functions. */

#ifndef MLKEM_H
#define MLKEM_H

#include <string>
#include <iostream>

#include "poly.h"
#include "rand.h"

namespace mlkem
{

static inline void pack_ek(uint8_t ek[poly_vec_len + seed_len], poly_vec* ekpv, const uint8_t seed[seed_len]) {
  poly_vec_to_bytes(ekpv, ek);
  memcpy(ek + poly_vec_len, seed, seed_len);
}

static inline void unpack_ek(poly_vec* ekpv, uint8_t seed[seed_len], const uint8_t ek[poly_vec_len + seed_len]) {
  bytes_to_poly_vec(ek, ekpv);
  memcpy(seed, ek + poly_vec_len, seed_len);
}

static inline void pack_dk(uint8_t dk[poly_vec_len], const poly_vec* dkpv) { 
  poly_vec_to_bytes(dkpv, dk); 
}

static inline void unpack_dk(poly_vec* dkpv, const uint8_t dk[poly_vec_len]) { 
  bytes_to_poly_vec(dk, dkpv); 
}

/* A cipher is comprised of a poly vec u:
 *   u = A^Ty + e_1 
 * and a poly v:
 *   v = t^Ty + e_2 + mu */

static inline void pack_cipher(uint8_t c[compressed_poly_vec_len + compressed_poly_len], 
    const poly_vec* u, const poly* v) {
  poly_vec_compress(c, u);
  poly_compress(c + compressed_poly_vec_len, v);
}

static inline void unpack_cipher(poly_vec* u, poly* v, 
  const uint8_t c[compressed_poly_vec_len + compressed_poly_len]) {
  poly_vec_decompress(u, c);
  poly_decompress(v, c + compressed_poly_vec_len);
}

/* Helper: generate key pairs given two random seeds
 * 
 *   ek: the encapsulation key (pub key) that sender will use to encap 
 *   dk: the decapsulation key (pri key) that receiver will use to decap */
static inline void indcpa_key_gen(uint8_t ek[poly_vec_len + seed_len], 
    uint8_t dk[poly_vec_len], const uint8_t seeds[seed_len]) {

  uint8_t buf[2 * seed_len]; /* used to store seeds */
  const uint8_t* pub_seed = buf;
  const uint8_t* noise_seed = buf + seed_len;

  memcpy(buf, seeds, seed_len);
  buf[seed_len] = mlkem_k;
  sha::sha3_512(buf, buf, seed_len + 1); /* expand 32+1 bytes to two pseudorandom 32-byte seeds */
  
  uint8_t nonce = 0;
  poly_vec a[mlkem_k], e, ekpv, dkpv, b;

  /* the process of calculating 
   *   t = As + e 
   * where A is a random matrix, t is ek, e is noise, s is dk 
   * specified as noisy linear system in ntt domain as in fips203 algo 13 */
  gen_matrix(a, pub_seed, 0);
  for (int i = 0; i < mlkem_k; ++i) { gen_noise_poly_eta1(&dkpv.vec[i], noise_seed, nonce++); }
  for (int i = 0; i < mlkem_k; ++i) { gen_noise_poly_eta1(&e.vec[i], noise_seed, nonce++); }
  poly_vec_ntt(&dkpv);
  poly_vec_ntt(&e);
  for (int i = 0; i < mlkem_k; ++i) {
    poly_vec_basemul(&ekpv.vec[i], &a[i], &dkpv);
    poly_tomont(&ekpv.vec[i]);
  }
  poly_vec_add(&ekpv, &ekpv, &e);
  poly_vec_reduce(&ekpv);

  pack_ek(ek, &ekpv, pub_seed);
  pack_dk(dk, &dkpv);
}

/* Key generation scheme, specified as algo 16 and 19
 *
 *   ek: output encapsulation key 
 *   dk: output decapsulation key */
inline void key_gen(uint8_t* ek, uint8_t* dk) {
  uint8_t seeds[2 * seed_len];
  gen_rand_bytes(seeds, 2 * seed_len);

  indcpa_key_gen(ek, dk, seeds);

  memcpy(dk + poly_vec_len, ek, poly_vec_len + seed_len); /* store ek(ek and pub seed to gen matrix) after dk */
  sha::sha3_256(dk + poly_vec_len + poly_vec_len + seed_len, ek, poly_vec_len + seed_len); /* also store ek hash */
  memcpy(dk + 2 * poly_vec_len + seed_len + sha::hash256_len, seeds + seed_len, seed_len); /* rejection value z */
}

/* Helper: encapsulate a msg using the encapsulation key and random seed
 * */
static inline void indcpa_enc(uint8_t c[compressed_poly_vec_len + compressed_poly_len],
  const uint8_t m[msg_len], const uint8_t ek[poly_vec_len + seed_len], const uint8_t seed[seed_len]) {
  
  uint8_t seed_received[seed_len];
  uint8_t nonce = 0;
  
  poly_vec at[mlkem_k], y, u, e1, ekpv;
  poly v, e2, mu;

  unpack_ek(&ekpv, seed_received, ek);
  msg_to_poly(&mu, m); /* convert msg to poly form */
  
  /* the encap process of calculating
   *   u = A^{T}y + e1, v = t^{T}y + e2 + u */
  gen_matrix(at, seed_received, 1); /* regen mat a used in key gen */
  for (int i = 0; i < mlkem_k; ++i) {
    gen_noise_poly_eta1(&y.vec[i], seed, nonce++);
  }
  for (int i = 0; i < mlkem_k; ++i) {
    gen_noise_poly_eta2(&e1.vec[i], seed, nonce++);
  }
  gen_noise_poly_eta2(&e2, seed, nonce++);

  poly_vec_ntt(&y);
  for (int i = 0; i < mlkem_k; ++i) {
    poly_vec_basemul(&u.vec[i], &at[i], &y); /* u = A^{T}y + e1 */
  }
  poly_vec_basemul(&v, &ekpv, &y); /* v = t^{T}y + e2 + mu */

  poly_vec_invntt(&u);
  poly_invntt(&v);
  
  poly_vec_add(&u, &u, &e1);
  poly_add(&v, &v, &e2);
  poly_add(&v, &v, &mu);
  poly_vec_reduce(&u);
  poly_reduce(&v);

  pack_cipher(c, &u, &v); /* byte encode */  
}

/* Encapsulate, specified as algo 17 and 20
 * 
 *   cipher: cipher text
 *   sk: the shared secret key used for future symmetric encryption 
 *   ek: the encryption key received */
inline void encap(uint8_t* cipher, uint8_t* sk, const uint8_t* ek) {
  uint8_t seed[seed_len];
  gen_rand_bytes(seed, seed_len);

  uint8_t buf[seed_len + sha::hash256_len];
  uint8_t kr[64]; /* hash of key and randomness */

  memcpy(buf, seed, seed_len);

  sha::sha3_256(buf + seed_len, ek, poly_vec_len + seed_len);
  sha::sha3_512(kr, buf, seed_len + sha::hash256_len);

  indcpa_enc(cipher, buf, ek, kr + seed_len);

  memcpy(sk, kr, seed_len);
}

/* Helper: decapsulate 
 * 
 *   m: output message 
 *   cipher: cipher received 
 *   dk: decapsulation key */
static inline void indcpa_dec(uint8_t m[msg_len], 
  const uint8_t cipher[compressed_poly_vec_len + compressed_poly_len], uint8_t dk[poly_vec_len]) {
  poly_vec u, dkpv;
  poly v, w;

  unpack_cipher(&u, &v, cipher);
  unpack_dk(&dkpv, dk); // dkpv is now in the ntt domain

  /* process of calculating
   *   w = v' - invntt(s^T(u')) 
   * notice the u and v will be different from u and v in encap as compression 
   * and decompression alter their values */
  poly_vec_ntt(&u); // u to ntt domain

  poly_vec_basemul(&w, &dkpv, &u);
  poly_invntt(&w);
  poly_sub(&w, &v, &w);
  poly_reduce(&w);
  
  poly_to_msg(m, &w);
}

/* Decapsulation, specified as algo 18 and 21
 * 
 *   ss: output shared secret key */
inline void decap(uint8_t* ss, uint8_t* cipher, uint8_t* dk) {
  int fail;
  /* cipher text decrypted, contains shared key and random seed for fo transform*/
  uint8_t buf[msg_len + sha::hash256_len];
  /* key and random seed used to re-encrypt */
  uint8_t kr[msg_len + seed_len];
  /* used to store newly generated cipher */
  uint8_t cmp[compressed_poly_vec_len + compressed_poly_len];
  const uint8_t* ek = dk + poly_vec_len;

  indcpa_dec(buf, cipher, dk);

  memcpy(buf + msg_len, dk + 2 * poly_vec_len + seed_len, sha::hash256_len); /* extract hash of encap key */
  sha::sha3_512(kr, buf, msg_len + sha::hash256_len);

  indcpa_enc(cmp, buf, ek, kr + msg_len); /* recalculate encrypted msg */

  fail = ccmp(cipher, cmp, compressed_poly_vec_len + compressed_poly_len);

  shake256_rkprf(ss, dk + 2 * poly_vec_len + seed_len + sha::hash256_len, cipher); /* compute rejection key */

  cmov(ss, kr, seed_len, !fail); /* constant copy kr */
}

} /* namespace mlkem */

#endif /* MLKEM_H */