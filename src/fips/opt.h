/* Copyright 2026, Yao Zeran, Zhang Chenzhi 
 *
 * The <opt.h> file defines ntt and barrett reduction related optmization logic. */

#ifndef OPT_H
#define OPT_H

#include "common.h"

namespace mlkem
{

/* Calc t = a * 2^{-16} mod q, where t in (-q, q)
 * */
inline int16_t montgomery_reduce(int32_t a) {
  int16_t t;
  // find t = m s.t. a + m * q = 0 mod 2^16
  t = (int16_t)a * mlkem_inverse_q;
  t = (a - (int32_t)t * mlkem_q) >> 16;
  return t;
}

/* Calc a mod q, works by estimate floor(a/q),
 * to calc, need to find a - q * floor(a/q) normally, by replace a/q by a * (v / 2^k) where
 *   v /approx (2^k / q)
 * */
inline int16_t barrett_reduce(int16_t a) {
  int16_t t;
  // v = floor(2^26/q + 1/2)
  const int16_t v = ((1 << 26) + mlkem_q / 2) / mlkem_q;
  t = ((int32_t)v * a + (1 << 25)) >> 26;
  t *= mlkem_q;
  return a - t;
}

static inline int16_t fqmul(int16_t a, int16_t b) {
  return montgomery_reduce((int32_t)a * b);
}

/* zetas_impl * 2^16 mod 3329 = zeta_fips203 */
const int16_t zetas[128] = {
  -1044,  -758,  -359, -1517,  1493,  1422,   287,   202,
   -171,   622,  1577,   182,   962, -1202, -1474,  1468,
    573, -1325,   264,   383,  -829,  1458, -1602,  -130,
   -681,  1017,   732,   608, -1542,   411,  -205, -1571,
   1223,   652,  -552,  1015, -1293,  1491,  -282, -1544,
    516,    -8,  -320,  -666, -1618, -1162,   126,  1469,
   -853,   -90,  -271,   830,   107, -1421,  -247,  -951,
   -398,   961, -1508,  -725,   448, -1065,   677, -1275,
  -1103,   430,   555,   843, -1251,   871,  1550,   105,
    422,   587,   177,  -235,  -291,  -460,  1574,  1653,
   -246,   778,  1159,  -147,  -777,  1483,  -602,  1119,
  -1590,   644,  -872,   349,   418,   329,  -156,   -75,
    817,  1097,   603,   610,  1322, -1285, -1465,   384,
  -1215,  -136,  1218, -1335,  -874,   220, -1187, -1659,
  -1185, -1530, -1278,   794, -1510,  -854,  -870,   478,
   -108,  -308,   996,   991,   958, -1460,  1522,  1628
};

inline void ntt(int16_t v[256]) {
  unsigned int j, len, start, k = 1;
  int16_t zeta, t;
  for (len = 128; len >= 2; len >>= 1) {
    for (start = 0; start < 256; start = j + len) {
      zeta = zetas[k++];
      for (j = start; j < start + len; ++j) {
        t = fqmul(zeta, v[j + len]);
        v[j + len] = v[j] - t;
        v[j] = v[j] + t;
      }
    }
  }
}

inline void invntt(int16_t v[256]) {
  unsigned int start, len, j, k = 127;
  int16_t zeta, t;
  const int16_t f = 1441; /* mont^2/128 */
  for (len = 2; len <= 128; len <<= 1) {
    for (start = 0; start < 256; start = j + len) {
      zeta = zetas[k--];
      for (j = start; j < start + len; ++j) {
        t = v[j];
        v[j] = barrett_reduce(t + v[j + len]);
        v[j + len] = v[j + len] - t;
        v[j + len] = fqmul(zeta, v[j + len]);
      }
    }
  }
  for (j = 0; j < 256; ++j) {
    v[j] = fqmul(v[j], f);
  }
}

/* Multiplication in the ring Z_q[x]/(x^2-zeta)
 *   where a(x) = a_0 + a_1x, b(x) = b_0 + b_1(x) 
 *   and a(x)b(x) = (a_0b_0) + (a_0b_1 + b_0a_1)x + (a_1b_1)x^2, x^2 = zeta
 */
inline void basemul(int16_t p[2], const int16_t a[2], const int16_t b[2], int16_t zeta) {
  p[0]  = fqmul(a[1], b[1]);
  p[0]  = fqmul(p[0], zeta);
  p[0] += fqmul(a[0], b[0]);
  p[1]  = fqmul(a[0], b[1]);
  p[1] += fqmul(a[1], b[0]);
}

} /* namespace mlkem */

#endif /* OPT_H */