/* Copyright 2026 (c), Yao Zeran, Zhang Chenzhi */

#include <iostream>
#include <format>

#include "../src/fips/include.h"

static inline void print_bytes(const uint8_t* in, int len) {
  for (int i = 0; i < len; ++i) {
    std::cout << std::format("{:02x}", *(in + i));
  }
  std::cout << std::endl;
}

static inline void print_numbers(const uint8_t* in, int len) {
  for (int i = 0; i < len; ++i) {
    std::cout << std::format("{:02d} ", *(in + i));
  }
  std::cout << std::endl;
}

uint8_t drb[mlkem::seed_len]{
  0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
  0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F, 0x10,
  0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18,
  0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F, 0x20
};

int main(int argc, char* argv[]) {

  /* test mlkem: key gen, encap, decap */
  uint8_t public_key[mlkem::poly_vec_len + mlkem::seed_len];
  uint8_t private_key[2 * mlkem::poly_vec_len + 2 * mlkem::seed_len + sha::hash256_len];

  std::cout << std::endl;
  std::cout << std::endl;
  mlkem::key_gen(public_key, private_key);
  std::cout << "Public key generated: " << std::endl;
  print_bytes(public_key, mlkem::poly_vec_len + mlkem::seed_len);
  std::cout << std::endl;
  std::cout << "Private key generated: "<< std::endl;
  print_bytes(private_key, mlkem::poly_vec_len);
  std::cout << std::endl;

  uint8_t alice_secret_key[32];
  uint8_t cipher[mlkem::compressed_poly_len + mlkem::compressed_poly_vec_len];
  mlkem::encap(cipher, alice_secret_key, public_key);

  uint8_t bob_secret_key[32];
  mlkem::decap(bob_secret_key, cipher, private_key);
  std::cout << std::endl;
  std::cout << "Alice secret key" << std::endl;
  print_bytes(alice_secret_key, 32);
  std::cout << std::endl;
  std::cout << "Bob secret key" << std::endl;
  print_bytes(bob_secret_key, 32);
  std::cout << std::endl;

  return 0;
}