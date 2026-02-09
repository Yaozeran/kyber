#ifndef WEAK_IMPL_H
#define WEAK_IMPL_H

#include <cstddef>
#include <array>
#include <cstdio>

#include "fips/include.h"

namespace weak_impl
{

constexpr size_t mlkem_q = 3329;

constexpr size_t mlkem_k = 2;
constexpr size_t mlkem_n = 16;

constexpr size_t mlkem_eta1 = 2;
constexpr size_t mlkem_eta2 = 2;

constexpr size_t poly_len = mlkem_n * 16 / 8;
constexpr size_t poly_vec_len = mlkem_k * poly_len;

constexpr size_t seed_len = 32;

typedef struct {
  int16_t coeffs[mlkem_n];
} poly;

typedef struct {
  poly vec[mlkem_k];
} poly_vec;

inline void poly_reduce(poly* p) { for (int i = 0; i < mlkem_n; ++i) { p->coeffs[i] %= mlkem_q; } }

inline void poly_add(poly* p, const poly* a, const poly* b) {
  for (int i = 0; i < mlkem_n; ++i) { p->coeffs[i] = a->coeffs[i] + b->coeffs[i]; }
  poly_reduce(p);
}

inline void poly_mult(poly* p, const poly* a, const poly* b) {
  for (int i = 0; i < mlkem_n; ++i) {
    for (int j = 0; j < mlkem_n; ++j) {
      int idx = i + j;
      int32_t prod = (int32_t)a->coeffs[i] * (int32_t)b->coeffs[j];
      if (idx < mlkem_n) {
        p->coeffs[idx] += prod;
      } else { // x^n = -1
        p->coeffs[idx - mlkem_n] -= prod;
      }
    }
  }
  poly_reduce(p);
}

inline void poly_vec_reduce(poly_vec* pv) { for (int i = 0; i < mlkem_k; ++i) { poly_reduce(&pv->vec[i]); } }

inline void poly_vec_add(poly_vec* pv, const poly_vec* a, const poly_vec* b) {
  for (int i = 0; i < mlkem_k; ++i) { poly_add(&pv->vec[i], &a->vec[i], &b->vec[i]); }
  poly_vec_reduce(pv);
}

inline void poly_vec_mult(poly* p, const poly_vec* a, const poly_vec* b) {
  poly t{0};
  for (int v = 0; v < mlkem_k; ++v) {
    for (int i = 0; i < mlkem_n; ++i) { t.coeffs[i] = 0; }
    poly_mult(&t, &a->vec[v], &b->vec[v]);
    poly_add(p, p, &t);
  }
}

inline void poly_to_bytes(const poly* p, uint8_t out[poly_len]) {
  // Serialize each coefficient as 16 bits (little-endian), canonicalize to [0, mlkem_q-1]
  for (size_t i = 0; i < mlkem_n; ++i) {
    int32_t v = p->coeffs[i] % (int32_t)mlkem_q;
    if (v < 0) { v += (int32_t)mlkem_q; }
    out[2*i + 0] = (uint8_t)(v & 0xff);
    out[2*i + 1] = (uint8_t)((v >> 8) & 0xff);
  }
}

inline void poly_vec_to_bytes(const poly_vec* pv, uint8_t out[poly_vec_len]) {
  // Concatenate serialized polys
  for (size_t i = 0; i < mlkem_k; ++i) {
    poly_to_bytes(&pv->vec[i], out + i * poly_len);
  }
}

static uint32_t load32_le(const uint8_t bytes[4]) {
  return (uint32_t)bytes[0] | (uint32_t)bytes[1] << 8 | (uint32_t)bytes[2] << 16 | (uint32_t)bytes[3] << 24;
}

static uint32_t load24_le(const uint8_t bytes[3]) {
  return (uint32_t)bytes[0] | (uint32_t)bytes[1] << 8 | (uint32_t)bytes[2] << 16;
}

/* Central binomial distribution (eta = 2) to generate small noise poly
 * 
 *   for eta = 2, 
 *   for each coefficient, we need 4 bits as 
 *     coeff = sum_{i=1}^eta(a_i) - sum_{i=(eta+1)}^{2*eta}(a_i) */
static inline void cbd2(poly* p, const uint8_t rand_bytes[4*mlkem_n/8]) {
  uint32_t t, d;
  int16_t a, b;
  for (int i = 0; i < mlkem_n / 8; ++i) {
    t = load32_le(rand_bytes + 4 * i);
    d = t & 0x55555555;
    d += (t >> 1) & 0x55555555;
    for (int j = 0; j < 8; ++j) {
      a = (d >> (4 * j + 0)) & 0x3;
      b = (d >> (4 * j + 2)) & 0x3;
      p->coeffs[8*i+j] = a - b;
    }
  }
}

/* Central binomial distribution (eta = 3) to generate small noise poly
 *
 *   for eta = 3, 
 *   for each coeff, we have
 *     coeff = (x_1 + x_2 + x_3) - (x_4 + x_5 + x_6) 
 *   so we need 6 * n / 8 bytes in total */
static inline void cbd3(poly* p, const uint8_t rand_bytes[6*mlkem_n/8]) {
  uint32_t t, d;
  int16_t a, b;
  for(int i = 0; i < mlkem_n / 4; ++i) {
    t = load24_le(rand_bytes + 3 * i);
    d = t & 0x00249249;
    d += (t >> 1) & 0x00249249;
    d += (t >> 2) & 0x00249249;
    for(int j = 0; j < 4; ++j) {
      a = (d >> (6 * j + 0)) & 0x7;
      b = (d >> (6 * j + 3)) & 0x7;
      p->coeffs[4*i+j] = a - b;
    }
  }
}

/* Central binomial distribution to generate small noise polynomial e
 * */
static inline void cbd_eta1(poly* p, const uint8_t rand_bytes[mlkem_eta1*mlkem_n/4]) {
  if (mlkem_eta1 == 2) {
    cbd2(p, rand_bytes);
  } else if (mlkem_eta1 == 3) {
    cbd3(p, rand_bytes);
  }
}

/* Central binomial distribution to generate small ephemeral noise polynomial 
 * */
static inline void cbd_eta2(poly* p, const uint8_t rand_bytes[mlkem_eta2*mlkem_n/4]) { cbd2(p, rand_bytes); }

inline void gen_rand_poly_eta1(poly* p, const uint8_t seed[seed_len], uint8_t nonce) {
  uint8_t buf[mlkem_eta1 * mlkem_n / 4];
  mlkem::shake256_prf(buf, sizeof(buf), seed, nonce);
  cbd_eta1(p, buf);
}

static inline unsigned int rej_sample_uniform(int16_t* ptr, unsigned int len, const uint8_t* buf, unsigned int buflen) {
  unsigned int cnt = 0, pos = 0;
  uint16_t v0, v1;
  while (cnt < len && pos + 3 <= buflen) {
    v0 = ((buf[pos + 0] >> 0) | ((uint16_t)buf[pos + 1] << 8)) & 0xfff;
    v1 = ((buf[pos + 1] >> 4) | ((uint16_t)buf[pos + 2] << 4)) & 0xfff;
    pos += 3;
    if (v0 < mlkem_q) { ptr[cnt++] = v0; }
    if (cnt < len && v1 < mlkem_q) { ptr[cnt++] = v1; }
  }
  return cnt;
}

static constexpr size_t mat_nblocks = ( 12 * mlkem_n / 8 * (1 << 12) / mlkem_q + sha::shake128_rate ) / sha::shake128_rate;

static inline void gen_matrix(poly_vec* a, const uint8_t seed[seed_len], int transposed) {
  unsigned int cnt, buflen;
  uint8_t buf[mat_nblocks * sha::shake128_rate];
  sha::keccak_ctx ctx;
  for (int i = 0; i < mlkem_k; ++i) {
    for (int j = 0; j < mlkem_k; ++j) {
      if (transposed) {
        mlkem::shake128_absorb(&ctx, seed, i, j);
      } else {
        mlkem::shake128_absorb(&ctx, seed, j, i);
      }
      sha::shake128_squeeze_blocks(buf, mat_nblocks, &ctx);
      buflen = mat_nblocks * sha::shake128_rate;
      cnt = rej_sample_uniform(a[i].vec[j].coeffs, mlkem_n, buf, buflen);
      while (cnt < mlkem_n) {
        sha::shake128_squeeze_blocks(buf, 1, &ctx);
        buflen = sha::shake128_rate;
        cnt += rej_sample_uniform(a[i].vec[j].coeffs + cnt, mlkem_n - cnt, buf, buflen);
      }
    }
  }
}

inline void key_gen(uint8_t* ek, uint8_t* dk, uint8_t* a, uint8_t* e) {
  // generate two random seeds
  std::array<uint8_t, 2 * seed_len> seeds;
  mlkem::gen_rand_bytes(seeds.data(), 2 * seed_len);
  // generate random components
  uint8_t nonce = 0;
  poly_vec mat[mlkem_k], ekpv{0}, dkpv, epv;
  gen_matrix(mat, seeds.data(), 0);
  for (int i = 0; i < mlkem_k; ++i) { gen_rand_poly_eta1(&dkpv.vec[i], seeds.data() + seed_len, nonce++); }
  for (int i = 0; i < mlkem_k; ++i) { gen_rand_poly_eta1(&epv.vec[i], seeds.data() + seed_len, nonce++); }
  for (int i = 0; i < mlkem_k; ++i) {
    poly_vec_mult(&ekpv.vec[i], &mat[i], &dkpv);
  }
  poly_vec_add(&ekpv, &ekpv, &epv);
  // serialize to output stream
  poly_vec_to_bytes(&ekpv, ek);
  poly_vec_to_bytes(&dkpv, dk);
  for (int i = 0; i < mlkem_k; ++i) {
    poly_vec_to_bytes(&mat[i], a);
    a += poly_vec_len;
  }
  poly_vec_to_bytes(&epv, e);
  memcpy(ek + poly_vec_len, seeds.data(), seed_len);
}

inline void encap() {
  
}

inline void decap() {

}

} /* namespace weak_impl */

#endif /* WEAK_H */