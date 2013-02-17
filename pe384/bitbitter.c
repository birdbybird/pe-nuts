#include <stdio.h>
#include <stdint.h>

#include "bitbitter.h"

/* count leading zeroes */
int clz(uint32_t x) {
  int zeros = 0;

  if (x == 0) { 
    zeros = (sizeof x) * 8;
    return zeros;
  }
  if ((x >> 16) & 0x0000FFFF)
    x = x >> 16;
  else 
    zeros += 16;
  if ((x >> 8) & 0x000000FF)
    x = x >> 8;
  else 
    zeros += 8;
  if ((x >> 4) & 0x0000000F)
    x = x >> 4;
  else zeros += 4;
  if ((x >> 2) & 0x00000003)
    x = x >> 2;
  else 
    zeros += 2;
  //if (x == 1)
  // zeros++;
  zeros += 1 - (x >> 1); /* +1 if (x == 1) */ 
  return zeros;
}

/* count leading zeros, again */
int nlz(uint32_t x) {
  int n = 0;

  if (x == 0) return 32;
  if ((x >> 16) == 0) { n += 16; x = x << 16; }
  if ((x >> 24) == 0) { n += 8; x = x << 8; }
  if ((x >> 28) == 0) { n += 4; x = x << 4; }
  if ((x >> 30) == 0) { n += 2; x = x << 2; }
  n = n + 1 - (x >> 31); // if shifted x still has a 0 leading bit, add 1
  return n;
}

int rounddown_pow2(uint32_t x) {
  return (32 - nlz(x) - 1);
}

int is_pow2(uint32_t x) {
  return ((x & (x - 1)) == 0 ? 1 : 0);
}

int is_even(uint32_t x) {
  return ((x & (-x)) == 1 ? 0 : 1);
}

#if TEST_MODE
int main(int argc, char **argv) {
  int i; 

  for (i = 0; i < 100; i++) {
    fprintf(stdout, "%d, %1$#x, nlz: %d, %d\n", i, nlz(i), clz(i));
  }
  return EXIT_SUCCESS;
}
#endif
