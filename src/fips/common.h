/* Copyright 2026, Yao Zeran, Zhang Chenzhi 
 * 
 * The <common.h> file defines constants used in mlkem and const time helper
 *   functions. */

#ifndef COMMON_H
#define COMMON_H

#include <cstdint>

/* random seed byte len */

namespace sha 
{
  constexpr int seed_len = 32;

  constexpr int hash128_len = 16;
   
  constexpr int hash256_len = 32;

  constexpr int hash512_len = 64;

} /* namespace sha */

namespace mlkem
{

constexpr int mlkem_config_set = 512;

/* mlkem operates in ring R_q: Z_q^n*/
constexpr int mlkem_q = 3329;

/* when using montgomery reduction, need q^{-1} mod 2^16 */
constexpr int mlkem_inverse_q = -3327;

/* number of coefficients in a polynomial */
constexpr int mlkem_n = 256;

/* number of poly in a vector of polys */
constexpr int mlkem_k = (mlkem_config_set == 512) ? 2 :
                        (mlkem_config_set == 768) ? 3 : 
                        (mlkem_config_set == 1024) ? 4 : -1;

constexpr int mlkem_eta1 = (mlkem_config_set == 512) ? 3 : 
                           (mlkem_config_set == 768) ? 2 : 
                           (mlkem_config_set == 1024) ? 2 : -1;

constexpr int mlkem_eta2 = 2;

/* polynomial byte len: 12 bits/coefficient, 256 * 12 / 8 = 384 */
constexpr int poly_len = 384;

/* polynomial vector byte len */
constexpr int poly_vec_len = mlkem_k * poly_len;

constexpr int compressed_poly_len = (mlkem_config_set == 512) ? 128 :
                                    (mlkem_config_set == 768) ? 128 :
                                    (mlkem_config_set == 1024) ? 160 : -1;

constexpr int compressed_poly_vec_len = (mlkem_config_set == 512) ? mlkem_k * 320 : 
                                        (mlkem_config_set == 768) ? mlkem_k * 320 :
                                        (mlkem_config_set == 1024) ? mlkem_k * 352 : -1;

/* other lens */

constexpr int seed_len = 32;

constexpr int msg_len = 32;

constexpr int ek_len = poly_vec_len + seed_len;
constexpr int dk_len = poly_vec_len + poly_vec_len + seed_len + sha::hash256_len + seed_len;

} /* namepsace mlkem */

namespace mlkem
{

/* Check if two given bytes array are the same in const time
 * */ 
inline int ccmp(const uint8_t* a, const uint8_t* b, size_t len) {
  uint8_t r = 0;
  for (int i = 0; i < len; ++i) {
    r |= a[i] ^ b[i];
  }
  return (-(uint64_t)r) >> 63;
  // return (((uint32_t)r | -(uint32_t)r) >> 31);
}

/* Copy content in 'in' to 'out' in const time 
 * */
inline void cmov(uint8_t* out, const uint8_t* in, size_t len, uint8_t b) {
  b = -b;
  for (int i = 0; i < len; ++i) {
    out[i] ^= b & (out[i] ^ in[i]);
  }
}

inline void cmov_int16(int16_t* out, int16_t in, uint16_t b) {
  b = -b;
  *out ^= b & ((*out) ^ in);
}

} /* namepsace mlkem */


#endif /* COMMON_H */