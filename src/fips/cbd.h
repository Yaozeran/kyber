/* Copyright 2026 (c), Yao Zeran, Zhang Chenzhi 
 * 
 * The <cbd.h> file implements usage of central binomial distribution function 
 *   to generate small noise polynomial e. */

#ifndef CBD_H
#define CBD_H

#include "common.h"
#include "poly.h"

namespace mlkem
{

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
inline void cbd_eta1(poly* p, const uint8_t rand_bytes[mlkem_eta1*mlkem_n/4]) {
  if (mlkem_eta1 == 2) {
    cbd2(p, rand_bytes);
  } else if (mlkem_eta1 == 3) {
    cbd3(p, rand_bytes);
  }
}

/* Central binomial distribution to generate small ephemeral noise polynomial 
 * */
inline void cbd_eta2(poly* p, const uint8_t rand_bytes[mlkem_eta2*mlkem_n/4]) { cbd2(p, rand_bytes); }

} /* namespace mlkem */

#endif /* CBD_H */