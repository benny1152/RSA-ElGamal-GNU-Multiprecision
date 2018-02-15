/* Copyright (C) 2017 Daniel Page <csdsp@bristol.ac.uk>
 *
 * Use of this source code is restricted per the CC BY-NC-ND license, a copy of
 * which can be found via http://creativecommons.org (and should be included as
 * LICENSE.txt within the associated archive or repository).
 */

#include "modmul.h"
#include "assert.h"


// Used Wikipedia article extensively (as well as lecture notes)
void montgomery_mult(mpz_t r, mpz_t x, mpz_t y, mp_limb_t omega, mpz_t N) {
    mp_limb_t r_0, y_i, x_0, u;
    // mpz_t t; // Using temp t therefore r and x can be passed as same variable (adkin to GMP functions)
    // mpz_init(t);
    // mpz_set_ui(t, 0);

    for (mp_size_t i = 0; i < mpz_size(N); i++) {
      x_0 = mpz_getlimbn(x, 0); // 0-th limb x
      y_i = mpz_getlimbn(y, i); // i-th limb y
      r_0 = mpz_getlimbn(r, 0); // 0-th limb r
      u = (r_0 + (y_i * x_0)) * omega; // u = (r_0 + y_i.x_0). omega (mod b)

      mpz_addmul_ui(r, x, y_i); // r = r + y_i.x
      mpz_addmul_ui(r, N, u);   // r = r + u.N
      mpz_tdiv_q_2exp(r, r, mp_bits_per_limb);  // r = r / b = (r + y_i.x + u.N) / b
    }

    // if (t > N) then t = t - N
    if(mpz_cmp(r,N) >= 0) {
      mpz_sub(r,r,N);
    }

    // mpz_clear(t);
}

void precompute_omega(mp_limb_t omega, mpz_t N) {
  omega = 1;
  mp_limb_t b = mpz_getlimbn(N, 0);

  for (mp_size_t i = 1; i <= mp_bits_per_limb; i++) {
    omega = omega * b;
  }

  omega = -omega;
}


/* Perform stage 1:
 *
 * - read each 3-tuple of N, e and m from stdin,
 * - compute the RSA encryption c, then
 * - write the ciphertext c to stdout.
 */

void stage1(int count) {
  mpz_t N, e, m, c;
  mpz_inits(N, e, m, c, NULL);

  for (int i = 0; i<(count/3); i++) {
    gmp_scanf("%Zx\n", N);
    gmp_scanf("%Zx\n", e);
    gmp_scanf("%Zx\n", m);

    assert(mpz_cmp(m, N) < 0);
    assert(mpz_cmp_ui(e, 1) > 0);
    mpz_powm(c, m, e, N);
    gmp_printf("%Zx\n", c);
  }
  mpz_clears(N, e, m, c, NULL);

}


/* Perform stage 2:
 *
 * - read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
 * - compute the RSA decryption m, then
 * - write the plaintext m to stdout.
 */

// RSA Decryption using CRT, sourced method from https://www.di-mgt.com.au/crt_rsa.html
// and Crypto A slides.
void stage2(int count) {
  mpz_t N, d, p, q, d_p, d_q, i_p, i_q, c, m, m_1, m_2;
  mpz_inits(N, d, p, q, d_p, d_q, i_p, i_q, c, m, m_1, m_2, NULL);
  for (int i = 0; i < (count/9); i++) {
    gmp_scanf("%Zx\n", N);
    gmp_scanf("%Zx\n", d);
    gmp_scanf("%Zx\n", p);
    gmp_scanf("%Zx\n", q);
    gmp_scanf("%Zx\n", d_p);
    gmp_scanf("%Zx\n", d_q);
    gmp_scanf("%Zx\n", i_p);
    gmp_scanf("%Zx\n", i_q);
    gmp_scanf("%Zx\n", c);

    mpz_powm(m_1, c, d_p, p);
    mpz_powm(m_2, c, d_q, q);
    mpz_sub(m_1, m_1, m_2);
    mpz_mul(m_1, m_1, i_q);
    mpz_mod(m_1, m_1, p);
    mpz_mul(m_1, m_1, q);
    mpz_add(m_1, m_1, m_2);

    gmp_printf("%Zx\n", m_1);
  }
  mpz_clears(N, d, p, q, d_p, d_q, i_p, i_q, c, m, m_1, m_2, NULL);

}

/* Perform stage 3:
 *
 * - read each 5-tuple of p, q, g, h and m from stdin,
 * - compute the ElGamal encryption c = (c_1,c_2), then
 * - write the ciphertext c to stdout.
 */

#define key_length 256
mpz_t random_numbers[50];
int numbers_index = 0;

// Using /dev/urandom to generate 50 random ephemeral key values, used
// https://stackoverflow.com/a/2572373 as a basis
void generate_numbers() {
  FILE *file;
  uint8_t ints[key_length/8];
  for (int i = 0; i < 50; i++) {
    mpz_init(random_numbers[i]);
    file = fopen("/dev/urandom", "r");
    fread(&ints, 1, sizeof ints, file);
    mpz_import(random_numbers[i], sizeof ints, 1, sizeof ints[0], 0, 0, ints);
  }
  fclose(file);
}

void clear_numbers() {
  for (int i = 0; i < 50; i++) {
    mpz_clear(random_numbers[i]);
  }
}

// void get_numbers_index() {
//   numbers_index++;
//   return numbers_index;
// }

void pseudorandom_generator(mpz_t k, mpz_t q) {
  mpz_t seed;
  mpz_init(seed);
  gmp_randstate_t state;
  gmp_randinit_mt(state); // 'Mersene Twister' algorithm random state creation (gmp docs)

  mpz_swap(seed, random_numbers[numbers_index]); // efficient number swapping
  numbers_index++;
  gmp_randseed(state, seed);

  mpz_urandomm(k, state, q);

  mpz_clear(seed);
  gmp_randclear(state);
}

void stage3(int count) {
  mpz_t p, q, g, h, m, c1, c2, k;
  mpz_inits(p, q, g, h, m, c1, c2, k, NULL);

  generate_numbers();
  // Uncomment below if a fixed ephemeral key of 1 is wanted
  // mpz_set_str(k, "1", 10);
  for (int i = 0; i < (count/5); i++) {
    gmp_scanf("%Zx\n", p);
    gmp_scanf("%Zx\n", q);
    gmp_scanf("%Zx\n", g);
    gmp_scanf("%Zx\n", h);
    gmp_scanf("%Zx\n", m);
    pseudorandom_generator(k, q);

    mpz_powm_sec(c1, g, k, p);
    mpz_powm_sec(c2, h, k, p);
    mpz_mul(c2, m, c2);
    mpz_mod(c2, c2, p);

    gmp_printf("%Zx\n", c1);
    gmp_printf("%Zx\n", c2);
  }
  mpz_clears(p, q, g, h, m, c1, c2, k, NULL);
  clear_numbers();
}

/* Perform stage 4:
 *
 * - read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
 * - compute the ElGamal decryption m, then
 * - write the plaintext m to stdout.
 */

void stage4(int count) {
  mpz_t p, q, g, x, c1, c2, m;
  mpz_inits(p, q, g, x, c1, c2, m, NULL);
  for (int i = 0; i < (count/5); i ++) {
    gmp_scanf("%Zx\n", p);
    gmp_scanf("%Zx\n", q);
    gmp_scanf("%Zx\n", g);
    gmp_scanf("%Zx\n", x);
    gmp_scanf("%Zx\n", c1);
    gmp_scanf("%Zx\n", c2);

    mpz_powm_sec(c1, c1, x, p);
    mpz_invert(c1, c1, p);
    mpz_mul(m, c2, c1);
    mpz_mod(m, m, p);
    gmp_printf("%Zx\n", m);
  }
  mpz_clears(p, q, g, x, c1, c2, m, NULL);
}

// Sums the number of lines in a file, made with the help of https://stackoverflow.com/a/12733630
int lines(char *filename) {
  FILE *file = fopen(filename, "r");
  int line_count = 0;
  int ch = 0;
  while (!feof(file)) {
    ch = fgetc(file);
    if (ch == '\n') {
        line_count++;
    }
  }
  fclose(file);
  return line_count;
}

/* The main function acts as a driver for the assignment by simply invoking the
 * correct function for the requested stage.
 */
int main( int argc, char* argv[] ) {

  if( 2 != argc ) {
    abort();
  }

  // Get number of lines in input file
  char stage[5];
  strcpy(stage, argv[1]);
  char *filename = strcat(stage, ".input");
  int count = lines(filename);

  if     ( !strcmp( argv[ 1 ], "stage1") ) {
    stage1(count);
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2(count);
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3(count);
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    stage4(count);
  }
  else {
    abort();
  }

  return 0;
}
