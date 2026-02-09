#include "../src/weak_impl.h"


int main(int argc, char* argv[]) {

  uint8_t encap_key[weak_impl::poly_vec_len + weak_impl::seed_len];
  uint8_t decap_key[weak_impl::poly_vec_len];
  
  uint8_t rand_noise_e[weak_impl::poly_vec_len];

  uint8_t matrix_a[weak_impl::poly_vec_len * weak_impl::mlkem_k];


  weak_impl::key_gen(encap_key, decap_key, matrix_a, rand_noise_e);

  printf(" ---- the encapsulation key ---- \n");
  for (int i = 0; i < weak_impl::poly_vec_len; ++i) { printf("%02x ", encap_key[i]); }
  printf("\n");
  printf(" ---- the encapsulation key ---- \n");
  int16_t* encap_coeffs = (int16_t*)encap_key;
  for (int i = 0; i < weak_impl::poly_vec_len / 2; ++i) {
      printf("%d ", encap_coeffs[i]);
  }
  printf("\n");

  printf(" ----  the random matrix a  ---- \n");
  for (int i = 0; i < weak_impl::poly_vec_len * weak_impl::mlkem_k; ++i) {
    printf("%02x ", matrix_a[i]); 
  }
  printf("\n");
  printf(" ----  the random matrix a  ---- \n");
  int16_t* mat_coeffs = (int16_t*)matrix_a;
  for (int i = 0; i < weak_impl::poly_vec_len * weak_impl::mlkem_k / 2; ++i) {
      printf("%d ", mat_coeffs[i]);
  }
  printf("\n");

  printf(" ---- the decapsulation key ---- \n");
  for (int i = 0; i < weak_impl::poly_vec_len; ++i) { printf("%02x ", decap_key[i]); }
  printf("\n");
  printf(" ---- the decapsulation key ---- \n");
  int16_t* decap_coeffs = (int16_t*)decap_key;
  for (int i = 0; i < weak_impl::poly_vec_len / 2; ++i) {
      printf("%d ", decap_coeffs[i]);
  }
  printf("\n");

  printf(" ---- the random noise poly ---- \n");
  for (int i = 0; i < weak_impl::poly_vec_len; ++i) { printf("%02x ", rand_noise_e[i]); }
  printf("\n");
  printf(" ---- the random noise poly ---- \n");
  int16_t* rand_noise_e_coeffs = (int16_t*)rand_noise_e;
  for (int i = 0; i < weak_impl::poly_vec_len / 2; ++i) {
      printf("%d ", rand_noise_e_coeffs[i]);
  }
  printf("\n");
};