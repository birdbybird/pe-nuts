/* PE 384
 *
 * calculate Sum_t{ G(F(t), F(t-1)) }, 2 <= t <= 45,
 *
 * where G(x, i) -- ith-order n-value in the sequence s_n, s.t. s(G(x, i)) = x
 * and
 * s_n -- sequence of partial sums of Rudin-Shapiro sequence
 *
 * F(0) = F(1) = 1;
 * F(2) = 2;
 *
 * argument seq.: 1  1  2  3  5  8  13  21  34  55..
 *                      ^
 * 	            start t=2
 *                   F(2)
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h> /* explicit integer types */
#include <stdbool.h>

#include <assert.h>

#include "bitbitter.h"

#define LOGMIN 1
#define LOG 0

enum { Ik1, Ik2, IkLast };
enum { s1, s2, s3, sLast };

typedef struct {
  uint32_t val;
  uint32_t nvals_count[IkLast];
} DetSum; /* determinant subsum */

typedef struct {
  uint32_t val;
#if TEST_MODE
  uint32_t nval_count;
#endif
  uint32_t rdown_pow2, rup_pow2;
  uint8_t flpow2, clpow2;
  uint8_t k1_pow4, k2_pow4;
  uint8_t nl_pow2, nr_pow2;
// void (*pop_s) (SnElem *, uint32_t sval); 
// void (*set_bounds) (SnElem *); 
} SnElem;

typedef struct {
  SnElem sn_elem;
  uint64_t nl, nr;
} SnElemExt;

void set_s_bounds(SnElem *s);
void set_s_bounds_ext(SnElemExt *s);
void pop_s(SnElem *s, uint32_t sval); /* populate an element */
uint64_t calculate_ith_nval_wake1a(uint32_t s, uint32_t target_nval_count);
uint64_t calculate_ith_nval_take1a(SnElem *s, uint32_t target_nval_count,
                                  uint8_t tg_nl_pow2, uint8_t tg_nr_pow2,
                                  uint32_t *add_nnvals_ptr
                                  /*, int *error_code */);
uint64_t calculate_ith_nval_take2(uint32_t sval, uint32_t target_nval_count);
uint64_t Gval(SnElem *s, uint32_t target_nval_count);
int calculate_Ik_nnvals_take2(uint32_t s,
                              uint8_t tg_nl_pow2, unsigned char tg_nr_pow2,
                              uint32_t *nval_num_fix_ptr);
int calculate_Ik_nnvals(SnElem *s,
                        uint8_t tg_nl_pow2, uint8_t tg_nr_pow2,
                        uint32_t *nval_count_fix_ptr);
uint64_t sum_Gvals_for_fibs(int *error); 
void print_fib_info(int t, SnElemExt *s, uint64_t tg_nval);

/* function implementations */
void set_s_bounds(SnElem *s) {
  s->flpow2 = rounddown_pow2(s->val); 
  s->clpow2 = s->flpow2 + 1;

  s->k1_pow4 = s->flpow2 - 1; 
  s->k2_pow4 = s->flpow2;

  s->rdown_pow2 = 1 << s->flpow2; 
  s->rup_pow2 = 1 << s->clpow2;

  s->nl_pow2 = 0;
  if (is_pow2(s->val) && s->val != 1) {
    s->nl_pow2 = (s->flpow2 - 1) << 1;
    s->nr_pow2 = s->flpow2 << 1;
  }
  else {
    if (s->flpow2 != 0) 
      s->nl_pow2 = (s->flpow2 << 1) - 1;
    s->nr_pow2 = (s->flpow2 << 1) + 1;
  }
}

void set_s_bounds_ext(SnElemExt *s) {
  set_s_bounds(&(s->sn_elem));
  s->nl = (uint64_t)1 << s->sn_elem.nl_pow2;
  s->nr = (uint64_t)1 << s->sn_elem.nr_pow2;
}

void pop_s(SnElem *s, uint32_t sval) {
  s->val = sval;
  // inline set_s_bounds(s);
  s->flpow2 = rounddown_pow2(s->val); 
  s->clpow2 = s->flpow2 + 1;

  s->k1_pow4 = s->flpow2 - 1; 
  s->k2_pow4 = s->flpow2;

  s->rdown_pow2 = 1 << s->flpow2; 
  s->rup_pow2 = 1 << s->clpow2;

  s->nl_pow2 = 0;
  if (is_pow2(s->val) && s->val != 1) {
    s->nl_pow2 = (s->flpow2 - 1) << 1;
    s->nr_pow2 = s->flpow2 << 1;
  }
  else {
    if (s->flpow2 != 0) // s != 1
      s->nl_pow2 = (s->flpow2 << 1) - 1;
    s->nr_pow2 = (s->flpow2 << 1) + 1; // tg_nl_pow2 + 2
  }
}

int calculate_Ik_nnvals_take2(uint32_t s,
                              uint8_t tg_nl_pow2, uint8_t tg_nr_pow2,
                              uint32_t *nval_num_fix_ptr) {
  uint32_t det_s1 = 0;
  uint32_t det_s2 = 0;

  uint32_t s_rdown_pow2, s_rup_pow2;

  uint8_t s_flpow2 = 0;
  uint8_t s_clpow2 = 0;

  uint8_t s_nl_pow2 = 0;
  uint8_t s_nr_pow2 = 0;
  uint8_t s_nmid_pow2 = 0;

  s_flpow2 = rounddown_pow2(s);  // k1_pow4 = s_flpow2 - 1;
  s_clpow2 = s_flpow2 + 1;       // k2_pow4 = s_flpow2;

  s_rdown_pow2 = 1 << s_flpow2;
  s_rup_pow2 = 1 << s_clpow2;

  assert(nval_num_fix_ptr);

  if (is_pow2(s) && s != 1) {
    s_nl_pow2 = (s_flpow2 - 1) << 1;
    s_nr_pow2 = s_flpow2 << 1;
    s_nmid_pow2 = (s_flpow2 << 1) - 1;
  }
  if (s_flpow2 != 0) // s != 1
    s_nl_pow2 = (s_flpow2 << 1) - 1;
  s_nr_pow2 = (s_flpow2 << 1) + 1; // tg_nl_pow2 + 2
  s_nmid_pow2 = s_flpow2 << 1;

  det_s1 = s - s_rdown_pow2;
  det_s2 = s_rup_pow2 - s;
#if LOG
  printf("\n%s called with s: %lu\n\t"
         "tg_nl_pow2: %u, tg_nr_pow2: %u\n\t"
         "s_nl_pow2: %u, s_nmid_pow2: %u, s_nr_pow2: %u\n",
         __func__,
         s,
         tg_nl_pow2, tg_nr_pow2,
         s_nl_pow2, s_nmid_pow2, s_nr_pow2);
#endif
  if (tg_nl_pow2 <= s_nl_pow2 && s_nr_pow2 <= tg_nr_pow2) {
    // valid s interval is inside target interval -- e.g initial case of calling for F(n)
#if LOG
    printf("-> standard valid interval\n");
#endif
    *nval_num_fix_ptr = 0;
    return s;
  }
  else if (s_nr_pow2 < tg_nl_pow2 || tg_nr_pow2 < s_nl_pow2) {
#if LOG
    printf("-> no valid nvals in the target invterval\n");
#endif
    return 0; // error code: #define INVALID_TARGET_INTERVAL -2
  }
  else if (s_nl_pow2 < tg_nl_pow2 && s_nr_pow2 <= tg_nr_pow2) {
    // assert: tg_nl_pow2 == s_nmid_pow2
    if (tg_nl_pow2 != s_nmid_pow2) {
#if LOG
      printf("-> cannot resolve nnvals in the target interval, cond. 1\n");
#endif
      return -3;
    }
    // => focus on s-values in Ik2, # of nvals = s - det_s2; // == 2 * s - s_rup_pow2
#if LOG
    printf("-> nvals resolved in Ik2 for s\n");
#endif
    *nval_num_fix_ptr = det_s2;
    return (s - det_s2);
  }
  else if (tg_nl_pow2 <= s_nl_pow2 && tg_nr_pow2 < s_nr_pow2) {
    // assert: tg_nr_pow2 == s_nmid_pow2
    if (tg_nr_pow2 != s_nmid_pow2) {
#if LOG
      printf("-> cannot resolve nnvals in the target interval, cond. 2\n");
#endif
      return -3;
    }
    // => work with s-values only in Ik1; # of nvals = det_s2
#if LOG
    printf("-> nvals resolved in Ik1 for s\n");
#endif
    *nval_num_fix_ptr = 0;
    return det_s2;
  }
}

uint64_t calculate_ith_nval_wake1a(uint32_t sval, uint32_t target_nval_count) {
  SnElem s = {0};
  uint32_t add_nnvals_dummy = 0;

  pop_s(&s, sval);
#if LOGMIN
  printf("\nF_cur: %u, tg_nval_count: %u\n", s.val, target_nval_count);
#endif
  return calculate_ith_nval_take1a(&s, target_nval_count, s.nl_pow2, s.nr_pow2, &add_nnvals_dummy);
}

/* ith-order nval (ith-nval) is calculated relative to the interval specified; 
 * topmost call for s = F(k) calculates ith-nval in the standard boundaries for s
 * (in the two adjacent intervals Ik1, Ik2) */
uint64_t calculate_ith_nval_take1a(SnElem *s, uint32_t target_nval_count,
                                   uint8_t tg_nl_pow2, uint8_t tg_nr_pow2,
                                   uint32_t *add_nnvals_ptr
                                   /*, int *error_code */) {
  uint32_t dets[sLast] = {0};
  SnElem subsum = {0};
  uint64_t target_nval = 0;

  uint8_t s_nmid_pow2 = 0;
  uint32_t add_nnvals = 0;
  uint8_t beta = 0;

  assert(add_nnvals_ptr != NULL);
  //assert(s->val >= target_nval_count || target_nval_count <= 2);

  /* base case */
  if (s->val == 1) {
    if (target_nval_count > 1)
      *add_nnvals_ptr = 1;
    return 0;
  }
  /* special processing for powers of 2 */
  if (is_pow2(s->val)) {
#if LOGMIN
    printf("s: %u is a power of 2\n", s->val);
#endif
    /* beta coef. directly from the solution to the recurrence equation */
    beta = (is_even(target_nval_count) ? 3 : 1); 
    /* interval is correct by default, by (18), (*) */
    target_nval = beta + 4 * calculate_ith_nval_wake1a(s->val >> 1, (target_nval_count >> 1) + 1);
#if LOGMIN
    printf("s: %10u, s_target_nval: %llu\n", 
           s->val, target_nval);
#endif
    return target_nval;
  }

  s_nmid_pow2 = s->flpow2 << 1;
  dets[s1] = s->val - s->rdown_pow2;
  dets[s2] = s->rup_pow2 - s->val; // confirmed: hypothesis: dets[s2] equals the number of nvals for s in Ik1
#if LOGMIN
  printf("take1a(): s: %u, i: %u, dets[s1]: %u, dets[s2]: %u\n", 
         s->val, target_nval_count, dets[s1], dets[s2]);
#endif

  if (tg_nl_pow2 <= s->nl_pow2 && s->nr_pow2 <= tg_nr_pow2) {
    // valid s interval is inside target interval -- e.g initial case of calling for F(k)
    // proceed to phase 2
    ;
  }
  else if (s->nr_pow2 < tg_nl_pow2 || tg_nr_pow2 < s->nl_pow2) {
    *add_nnvals_ptr = 0;
    printf("now here, INVALID_TARGET_INTERVAL !\n");
    return 0;
  }
  else if (s->nl_pow2 < tg_nl_pow2 && s->nr_pow2 <= tg_nr_pow2) {
    // assert(tg_nl_pow2 == s_nmid_pow2) =>
    // focus on s-values in Ik2, # of nvals = s - dets[s2]; (== 2 * s - s_rup_pow2)
    if (s->val - dets[s2] < target_nval_count) {
      *add_nnvals_ptr = s->val - dets[s2];
      return 0;
    }
    // Ik2_only = 1;
    target_nval_count += dets[s2]; // proceed for s with updated target_nval_count
  }
  else if (tg_nl_pow2 <= s->nl_pow2 && tg_nr_pow2 < s->nr_pow2) {
    // assert(tg_nr_pow2 == s_nmid_pow2) =>
    // work with s-values only in Ik1; # of nvals = dets[s2]
    if (dets[s2] < target_nval_count) {
      *add_nnvals_ptr = dets[s2];
      return 0;
    }
    // Ik1_only = 1;
    // target_nval_count is not changed, proceed for s
  }

  /* phase 2 */
  if (target_nval_count > dets[s2]) {
    // (15), target_nval falls in Ik2 interval for s
    if (dets[s1] < target_nval_count - dets[s2])
      add_nnvals = dets[s1];
    else {
      pop_s(&subsum, dets[s1]);
      target_nval = calculate_ith_nval_take1a(&subsum, target_nval_count - dets[s2],
                                              // 0, 2 * k2_pow4 - 1
                                              0, (s->k2_pow4 << 1) - 1,
                                              &add_nnvals);
    }
    if (add_nnvals != 0) {
      // by (16):
      dets[s3] = 3 * s->rdown_pow2  - s->val;
      // note: not using add_nnvals since dets[s3] must generate enough nvals
      // by properties of s_n sequence
      pop_s(&subsum, dets[s3]);
      target_nval = calculate_ith_nval_take1a(&subsum,
                                              target_nval_count - dets[s2] - add_nnvals,
                                              (s->k2_pow4 << 1) - 1, s->k2_pow4 << 1,
                                              &add_nnvals);
    }
    target_nval += (uint64_t)1 << (s->k2_pow4 << 1); // s(2^(2*k2_pow4)), 1st half of Ik2
  }
  else { // (dets[s2] >= target_nval_count)
    // searched s_target_nval is in Ik1 interval for s
    if (dets[s1] < target_nval_count)
      add_nnvals = dets[s1];
    else {
      pop_s(&subsum, dets[s1]);
      target_nval = calculate_ith_nval_take1a(&subsum, target_nval_count,
                                              0, s->k1_pow4 << 1, // 0, 2 * k1_pow4
                                              &add_nnvals);
    }
    if (add_nnvals != 0) { // onto dets[s2]:
      // s_target_nval has to be derived from (dets[s2], target_nval_count - dets[s1]_nvals_Ik1)
      // in Ik
      pop_s(&subsum, dets[s2]);
      target_nval = calculate_ith_nval_take1a(&subsum,
                                              target_nval_count - add_nnvals,
                                              s->k1_pow4 << 1, (s->k1_pow4 << 1) + 1,
                                              &add_nnvals);
    }
    target_nval += ((uint64_t)1 << (s->k1_pow4 << 1)) << 1; // 2nd half of Ik1
  }
#if LOGMIN
  printf("\ts: %10u, s_target_nval: %llu\n", s->val, target_nval);
#endif
  return target_nval;
}

enum { 
  INVALID_TARGET_INTERVAL = -2,
  NO_VALUES_IN_TARGET_INTERVAL = -3
};

// TODO: check return type 
int calculate_Ik_nnvals(SnElem *s,
                        uint8_t tg_nl_pow2, uint8_t tg_nr_pow2,
                        uint32_t *nval_count_fix_ptr) {
  uint32_t det_s1 = 0;
  uint32_t det_s2 = 0;

  uint8_t s_nl_pow2 = 0;
  uint8_t s_nr_pow2 = 0;
  uint8_t s_nmid_pow2 = 0;

  assert(nval_count_fix_ptr);

  // target interval:
  // -- assert overlapping with Ik1, Ik2
  // -- assert target interval is entirely in (1) Ik1 or (2) Ik2 or in (3) their union

  // given: target pow2 interval : [ tg_nl_pow2, tg_nr_pow2 ]
  // s valid interval : [ 2 * 4^(s_flpow2 - 1), 2 * 4^(s_flpow2) ] or
  //                    [ 2 *  s_flpow2 - 1; 2 * s_flpow2 + 1 ]
  if (is_pow2(s->val) && s->val != 1) {
    s_nl_pow2 = (s->flpow2 - 1) << 1;
    s_nr_pow2 = s->flpow2 << 1;
    s_nmid_pow2 = (s->flpow2 << 1) - 1;
  }
  if (s->flpow2 != 0) // s != 1
    s_nl_pow2 = (s->flpow2 << 1) - 1;
  s_nr_pow2 = (s->flpow2 << 1) + 1; // tg_nl_pow2 + 2
  s_nmid_pow2 = s->flpow2 << 1;

  det_s1 = s->val - s->rdown_pow2;
  det_s2 = s->rup_pow2 - s->val; // hypthesis: det_s2 equals the number of nvals for s in Ik1
#if LOG
  printf("\n%s called with s: %lu\n\t"
         "tg_nl_pow2: %u, tg_nr_pow2: %u\n\t"
         "s_nl_pow2: %u, s_nmid_pow2: %u, s_nr_pow2: %u\n",
         __func__,
         s->val,
         tg_nl_pow2, tg_nr_pow2,
         s_nl_pow2, s_nmid_pow2, s_nr_pow2);
#endif
  if (tg_nl_pow2 <= s_nl_pow2 && s_nr_pow2 <= tg_nr_pow2) {
    /* valid s interval is inside target interval -- e.g initial case of calling for F(n) */
#if LOG
  printf("-> standard valid interval\n");
#endif
    *nval_count_fix_ptr = 0;
    return s->val;
  }
  else if (s_nr_pow2 < tg_nl_pow2 || tg_nr_pow2 < s_nl_pow2) {
#if LOG
    printf("-> no valid nvals in the target invterval\n");
#endif
    return 0; // return INVALID_TARGET_INTERVAL;
  }
  else if (s_nl_pow2 < tg_nl_pow2 && s_nr_pow2 <= tg_nr_pow2) {
    if (tg_nl_pow2 != s_nmid_pow2) {
#if LOG
      printf("-> cannot resolve nnvals in the target interval, cond. 1\n");
#endif
      return NO_VALUES_IN_TARGET_INTERVAL;
    }
    /* => focus on s-values in Ik2, # of nvals = s - det_s2; // == 2 * s - s_rup_pow2 */
#if LOG
    printf("-> nvals resolved in Ik2 for s\n");
#endif
    *nval_count_fix_ptr = det_s2;
    return (s->val - det_s2);
  }
  else if (tg_nl_pow2 <= s_nl_pow2 && tg_nr_pow2 < s_nr_pow2) {
    if (tg_nr_pow2 != s_nmid_pow2) {
#if LOG
      printf("-> cannot resolve nnvals in the target interval, cond. 2\n");
#endif
      return NO_VALUES_IN_TARGET_INTERVAL;
    }
    /* => work with s-values only in Ik1; # of nvals = det_s2 */
#if LOG
    printf("-> nvals resolved in Ik1 for s\n");
#endif
    *nval_count_fix_ptr = 0;
    return det_s2;
  }
}

// TODO: use pop_s()
/* calculate_ith_nval_take2(), s via pointer */
uint64_t Gval(SnElem *s, uint32_t target_nval_count) {
  DetSum dets[sLast] = {{0}};
  SnElem subsum = {0};
  uint32_t nval_count_fix = 0;
  uint64_t target_nval = 0;
  uint8_t beta = 0;

  /* base case i.e. target_nval == 0 */
  if (s->val == 1) {
    return 0;
  }
  /* special processing for powers of 2 */
  if (is_pow2(s->val)) {
#if LOG
    printf("s: %lu is a power of 2\n", s->val);
#endif
    /* beta coef. directly from the solution to the recurrence equation */
    beta = (is_even(target_nval_count) ? 3 : 1); 
    /* interval is correct by default, see (18), (*) */
    subsum.val = s->val >> 1;
    set_s_bounds(&subsum);
    target_nval = beta + 4 * Gval(&subsum, (target_nval_count >> 1) + 1);
#if LOGMIN
    printf("s: %10u, s_target_nval: %llu\n", s->val, target_nval);
#endif
    return target_nval;
  }

  dets[s1].val = s->val - s->rdown_pow2;
  dets[s2].val = s->rup_pow2 - s->val; /* det_s2 equals the number of nvals for s in Ik1 */
#if LOG 
  printf("Gval(): s: %u, i: %lu; dets[s1].val: %lu, dets[s2]: %lu\n",
         s->val, target_nval_count,
         dets[s1].val, dets[s2].val);
#endif
  if (target_nval_count > dets[s2].val) { 
    /* target_nval falls into Ik2 interval for s */
    subsum.val = dets[s1].val;
    set_s_bounds(&subsum);
    dets[s1].nvals_count[Ik2] = calculate_Ik_nnvals(&subsum,
                                                    // 0, 2 * s->k2_pow4 - 1
                                                    0, (s->k2_pow4 << 1) - 1,
                                                    &nval_count_fix);
    if (dets[s1].nvals_count[Ik2] >= target_nval_count - dets[s2].val) {
      target_nval = Gval(&subsum,
                         target_nval_count - dets[s2].val + nval_count_fix);
    }
    else {
      /* by (16) */
      dets[s3].val = 3 * s->rdown_pow2 - s->val;
      subsum.val = dets[s3].val;
      set_s_bounds(&subsum);
      dets[s3].nvals_count[Ik2] = calculate_Ik_nnvals(&subsum,
                                                     // 2 * s->k2_pow4 - 1, 2 * s->k2_pow4
                                                     (s->k2_pow4 << 1) - 1, s->k2_pow4 << 1,
                                                     &nval_count_fix);
      target_nval = Gval(&subsum,
                         target_nval_count - dets[s2].val - dets[s1].nvals_count[Ik2]
                         + nval_count_fix);
    }
    target_nval += (uint64_t)1 << (s->k2_pow4 << 1); /* s(2^(2*s->k2_pow4)), 1st half of Ik2 */
  }
  else { /* (det_s2 >= target_nval_count); (17) applies */
    /* target_nval is in Ik1 interval for s */
    subsum.val = dets[s1].val;
    set_s_bounds(&subsum);
    dets[s1].nvals_count[Ik1] = calculate_Ik_nnvals(&subsum,
                                                    // 0, 2 * s->k1_pow4
                                                    0, s->k1_pow4 << 1,
                                                    &nval_count_fix);
    if (dets[s1].nvals_count[Ik1] >= target_nval_count) {  /* _unlikely ? */
      target_nval = Gval(&subsum,
                         target_nval_count + nval_count_fix);
    }
    else { /* onto det. sum == det_s2: apply (18) */
      subsum.val = dets[s2].val; 
      set_s_bounds(&subsum);
      dets[s2].nvals_count[Ik1] = calculate_Ik_nnvals(&subsum,
                                                      // 2 * s->k1_pow4, 2 * s->k1_pow4 + 1
                                                      s->k1_pow4 << 1, (s->k1_pow4 << 1) + 1,
                                                      &nval_count_fix);
      target_nval = Gval(&subsum,
                         target_nval_count - dets[s1].nvals_count[Ik1]
                         + nval_count_fix);
    }
    target_nval += ((uint64_t)1 << (s->k1_pow4 << 1)) << 1; /* 2nd half of Ik1 */
  }
#if LOGMIN
  printf("s: %10u, s_target_nval: %llu\n", s->val, target_nval);
#endif
  return target_nval;
}

/* Find ith-order nval, s.t. s_n(nval) == s
 * (nval is always found in Ik1 or Ik2 for a given s val);
 * version which uses calculate_Ik_nnvals() */
uint64_t calculate_ith_nval_take2(uint32_t sval, uint32_t target_nval_count) {
  DetSum dets[sLast] = {{0}};
  SnElem s = { .val = sval };
  uint32_t nval_count_fix = 0;
  uint64_t target_nval = 0;
  uint8_t beta = 0;

  /* base case i.e. target_nval == 0; */
      if (s.val == 1) {
    return 0;
  }

  set_s_bounds(&s);

  if (is_pow2(s.val)) {
#if LOG
    printf("s: %lu is a power of 2\n", s.val);
#endif
    /* beta coef. directly from the solution to the recurrence equation */
    beta = (is_even(target_nval_count) ? 3 : 1); 
    /* interval is correct by default, see (18), (*) */
    target_nval = beta + 4 * calculate_ith_nval_take2(s.val >> 1, (target_nval_count >> 1) + 1);
#if LOGMIN
    printf("s: %10u, s_target_nval: %llu\n", s.val, target_nval);
#endif
    return target_nval;
  }

  dets[s1].val = s.val - s.rdown_pow2;
  dets[s2].val = s.rup_pow2 - s.val; /* hypothesis: det_s2 equals the number of nvals for s in Ik1 */

  if (target_nval_count > dets[s2].val) { /* target_nval falls into Ik2 interval for s */
    /* (15) applies (det_s1 reappears in it, so we do not calculate;
     * continue recursively for appropriate sum det_s1 in this interval */
    dets[s1].nvals_count[Ik2] = calculate_Ik_nnvals_take2(dets[s1].val,
                                                    // 0, 2 * s.k2_pow4 - 1
                                                    0, (s.k2_pow4 << 1) - 1,
                                                    &nval_count_fix); 
    if (dets[s1].nvals_count[Ik2] >= target_nval_count - dets[s2].val) {
      target_nval = calculate_ith_nval_take2(dets[s1].val,
                                             target_nval_count - dets[s2].val + nval_count_fix);
    }
    else {
      dets[s3].val = 3 * s.rdown_pow2 - s.val;
      dets[s3].nvals_count[Ik2] = calculate_Ik_nnvals_take2(dets[s3].val,
                                                     // 2 * s.k2_pow4 - 1, 2 * s.k2_pow4
                                                     (s.k2_pow4 << 1) - 1, s.k2_pow4 << 1,
                                                     &nval_count_fix);
      target_nval = calculate_ith_nval_take2(dets[s3].val,
                                             target_nval_count - dets[s2].val - dets[s1].nvals_count[Ik2]
                                             + nval_count_fix);
    }
    target_nval += (uint64_t)1 << (s.k2_pow4 << 1); /* s(2^(2*s.k2_pow4)), 1st half of Ik2 */
  }
  else { /* (det_s2 >= target_nval_count); apply (17) */
    /* target_nval is in Ik1 interval for s */
    dets[s1].nvals_count[Ik1] = calculate_Ik_nnvals_take2(dets[s1].val,
                                                          // 0, 2 * s.k1_pow4
                                                          0, s.k1_pow4 << 1,
                                                          &nval_count_fix);
    if (dets[s1].nvals_count[Ik1] >= target_nval_count) {  /* _unlikely ? */
      target_nval = calculate_ith_nval_take2(dets[s1].val,
                                             target_nval_count + nval_count_fix);
    }
    else { /* onto det. sum == det_s2: apply (18) */
      dets[s2].nvals_count[Ik1] = calculate_Ik_nnvals_take2(dets[s2].val,
                                                            // 2 * s.k1_pow4, 2 * s.k1_pow4 + 1
                                                            s.k1_pow4 << 1, (s.k1_pow4 << 1) + 1,
                                                            &nval_count_fix);
      target_nval = calculate_ith_nval_take2(dets[s2].val,
                                             target_nval_count - dets[s1].nvals_count[Ik1]
                                             + nval_count_fix);
    }
    target_nval += ((uint64_t)1 << (s.k1_pow4 << 1)) << 1; /* 2nd half of Ik1 */
  }
#if LOGMIN
  printf("s: %10u, s_target_nval: %llu\n", s.val, target_nval);
#endif
  return target_nval;
}

void print_fib_info(int t, SnElemExt *fib, uint64_t tg_nval) {
  printf("t: %2d: F_cur: %10u, clpow2: %u,\n\t"
          "nl_pow2: %u  nr_pow2: %u, nl: %llu  nr: %llu\n\t"
          "tg_nval: %llu\n",
          t, fib->sn_elem.val, fib->sn_elem.clpow2,
          fib->sn_elem.nl_pow2, fib->sn_elem.nr_pow2, 
          fib->nl, fib->nr,
          tg_nval);
  bool tg_inrange = (fib->nl <= tg_nval && tg_nval <= fib->nr);
  printf("\t->%s in range\n\n", tg_inrange ? "" : " NOT");
}

#define COLOR_RESET       "\033[0m"
#define COLOR_BOLDGREEN   "\033[1m\033[32m" /* Bold Green */

uint64_t sum_Gvals_for_fibs(int *error) {
  /* our fibs:
   * about 3 fibs per Ik'th (first occurrence), interval: Ik, I(k+1)
   * about 45/3 = 15 intervals Ik, max_k = 15 (== max_k_pow4), max_k_pow2 = 30
   * => 32 bit numbers shall be enough to store fibs nvals */
  uint32_t F_prev = 1;
  uint32_t F_cur = 2;
#if LOGMIN
  SnElemExt fib = {{0}};
  uint64_t tg_nval = 0;
#else
  SnElem fib = {0};
#endif
  int t;
  uint64_t nvals_grand_sum = 0;

  for (t = 2; t <= 45; t++) {
#if LOGMIN
    fib.sn_elem.val = F_cur;
    set_s_bounds_ext(&fib);
    // TODO: error codes by Gval()
    //tg_nval = Gval(&(fib.sn_elem), F_prev);
    //tg_nval = calculate_ith_nval_take2(F_cur, F_prev);
    tg_nval = calculate_ith_nval_wake1a(F_cur, F_prev);
#if TEST_MODE
    assert(F_cur == gen_s(tg_nval, 0));
#endif
    print_fib_info(t, &fib, tg_nval);
    nvals_grand_sum += tg_nval;
#else
    fib.val = F_cur;
    set_s_bounds(&fib);
    //nvals_grand_sum += Gval(&fib, F_prev);
    //nvals_grand_sum += calculate_ith_nval_take2(F_cur, F_prev);
    nvals_grand_sum += calculate_ith_nval_wake1a(F_cur, F_prev);
#endif
    F_cur += F_prev;
#if LOGMIN
    F_prev = fib.sn_elem.val;
#else
    F_prev = fib.val;
#endif
  }
  printf(COLOR_BOLDGREEN "\nANSWER: nvals_grand_sum: " COLOR_RESET "%llu\n", nvals_grand_sum);

  return nvals_grand_sum; 
}

#define PICK_FUNC(func, take) PICK_FUNC_(func, take)
#define PICK_FUNC_(func, take) func ## _ ## take

#include <string.h> /* strcmp */

int main(int argc, char *argv[]) {
  int error = 0;

  if (argc == 2 && (strcmp("-t", argv[1]) == 0)) {
    printf("gotta test 'em all\n");
#if TEST_MODE
    test_harness(PICK_FUNC(sum_Gvals_for_fibs, take2), &error);
#endif
  }
  else if (argc == 2 && (strcmp("-s", argv[1]) == 0)) {
    printf("test test: we'd like to call.. %p\n", PICK_FUNC_(calculate_ith_nval, take2));
  }
  else if (argc == 1) {
    (void) sum_Gvals_for_fibs(&error);
  }
  else {
    printf("intention unclear ..halt\n");
  }
  return error;
}

