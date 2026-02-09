/* Copyright 2026 (c), Yao Zeran, Zhang Chenzhi 
 * 
 * The <rand.h> file implements get entropy src from os system call, and how to 
 *   gen random matrix and polynomials */

#ifndef RAND_H
#define RAND_H

#ifdef _WIN32
#include <windows.h>
#include <wincrypt.h>
#else
#include <fcntl.h>
#include <errno.h>
#endif

#if defined(__linux__) || defined(__APPLE__)
#include <unistd.h>
#include <sys/random.h>
#endif

#include <cstdint>
#include <cstdlib>

#include "poly.h"
#include "shake.h"
#include "cbd.h"

namespace mlkem
{

/* Generate random bytes from os entropy src
 * */
inline void gen_rand_bytes(uint8_t* out, size_t outlen) {
#ifdef _WIN32
  // todo
#elif defined(__linux__) || defined(__APPLE__)
  static int fd = -1;
  ssize_t ret;

  while(fd == -1) {
    fd = open("/dev/urandom", O_RDONLY);
    if(fd == -1 && errno == EINTR)
      continue;
    else if(fd == -1)
      abort();
  }

  while(outlen > 0) {
    ret = read(fd, out, outlen);
    if(ret == -1 && errno == EINTR)
      continue;
    else if(ret == -1)
      abort();

    out += ret;
    outlen -= ret;
  }
#endif
}

static constexpr size_t mat_nblocks = ( 12 * mlkem_n / 8 * (1 << 12) / mlkem_q + sha::shake128_rate ) / sha::shake128_rate;

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

inline void gen_matrix(poly_vec* a, const uint8_t seed[seed_len], int transposed) {
  unsigned int cnt, buflen;
  uint8_t buf[mat_nblocks * sha::shake128_rate];
  sha::keccak_ctx ctx;
  for (int i = 0; i < mlkem_k; ++i) {
    for (int j = 0; j < mlkem_k; ++j) {
      if (transposed) {
        shake128_absorb(&ctx, seed, i, j);
      } else {
        shake128_absorb(&ctx, seed, j, i);
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

inline void gen_noise_poly_eta1(poly* p, const uint8_t seed[seed_len], uint8_t nonce) {
  uint8_t buf[mlkem_eta1 * mlkem_n / 4];
  shake256_prf(buf, sizeof(buf), seed, nonce);
  cbd_eta1(p, buf);
}

inline void gen_noise_poly_eta2(poly* p, const uint8_t seed[seed_len], uint8_t nonce) {
  uint8_t buf[mlkem_eta2 * mlkem_n / 4];
  shake256_prf(buf, sizeof(buf), seed, nonce);
  cbd_eta2(p, buf);
}

} /* namespace mlkem */

#endif /* RAND_H */