#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "bitbitter.h"

/* count 1-bits in x */
int c1s(uint32_t x) {
  x = (x & 0x55555555) + ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
  x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
  x = (x & 0x0000FFFF) + ((x >>16) & 0x0000FFFF);
  return x;
}

int gen_a(uint32_t i) {
  int n11s = i & (i << 1);
  return c1s(n11s);
}

int gen_b(uint32_t i) {
  int n11s = i & (i << 1); /* i's adjacent pairs of 11s as 1-bits in n11s */
  /* then count the 1-bits; b = (-1)^c1s() */
  if (is_even(c1s(n11s)))
    return 1;
  else 
    return -1;
}

/* generate a sequence s of partial sums of Rudin-Shapiro sequence b */
uint32_t gen_s(uint32_t count, int all) {
  uint32_t s = 1; /* a(0) = 0, b(0) = (-1)^0 = 1 => s(0) = 1 */
  uint32_t i = 0;

  if (count == 0 || all == 1)
    printf("%u\n", s);
  for (i = 1; i <= count; i++) {
    s += gen_b(i);
    if (i == count || all == 1)
      printf("%u\n", s);
  }
  return s;
}

int main(int argc, char *argv[]) {
  uint32_t n = 0;

  /* generate s elements up to n */
  if (argc == 3 && (strcmp("-c", argv[1]) == 0)) {
    if (sscanf(argv[2],"%u", &n) == 1)
      (void) gen_s(n, 1);
  }
  /* generate nth element of s */
  else if (argc == 3 && (strcmp("-n", argv[1]) == 0)) {
    if (sscanf(argv[2],"%u", &n) == 1)
      (void) gen_s(n, 0);
  }
#if 0
  /* generate 1st-order nval of s == n */
  else if (argc == 3 && (strcmp("-s", argv[1]) == 0)) 
    if (sscanf(argv[2],"%u", &n) == 1)
      gen_s_nval(n, 1);
#endif
  else {
      printf("no-go\n");
  }
  return 0;
}

