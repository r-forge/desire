/*
 * truncnorm.c - Implementation of truncated normal distribution
 *
 * Authors:
 *  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
 *  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
 *  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#define UNPACK_REAL_VECTOR(ARGS, SXP, DBL) \
  ARGS = CDR(ARGS); \
  SXP = CAR(ARGS); \
  PROTECT(SXP = coerceVector(SXP, REALSXP)); \
  DBL = REAL(SXP);

#define UNPACK_REAL_VALUE(ARGS, SXP, DBL) \
  ARGS = CDR(ARGS); \
  SXP = CAR(ARGS); \
  PROTECT(SXP = coerceVector(SXP, REALSXP)); \
  DBL = *REAL(SXP);

#define UNPACK_INTEGER_VALUE(ARGS, SXP, INT) \
  ARGS = CDR(ARGS); \
  SXP = CAR(ARGS); \
  PROTECT(SXP = coerceVector(SXP, INTSXP)); \
  INT = *INTEGER(SXP);

#define ALLOC_REAL_VECTOR(SIZE, SXP, DBL) \
  PROTECT(SXP = allocVector(REALSXP, SIZE)); \
  DBL = REAL(SXP);


SEXP dtruncnorm(SEXP args) {
  R_len_t i, n;
  SEXP s_x, s_a, s_b, s_mean, s_sd, s_ret;
  double *x, a, b, mean, sd, c1, c2, *ret;

  UNPACK_REAL_VECTOR(args, s_x, x);
  UNPACK_REAL_VALUE(args, s_a, a);
  UNPACK_REAL_VALUE(args, s_b, b);
  UNPACK_REAL_VALUE(args, s_mean, mean);
  UNPACK_REAL_VALUE(args, s_sd, sd);

  n = length(s_x);
  ALLOC_REAL_VECTOR(n, s_ret, ret);

  c2 = sd * pnorm(b, mean, sd, FALSE, FALSE) - pnorm(a, mean, sd, FALSE, FALSE);
  for (i = 0; i < n; ++i) {
    double t = x[i];
    if (a <= t && t <= b) { /* In range: */
      c1 = dnorm(t, mean, sd, FALSE);
      ret[i] = - c1 / c2;
    } else { /* Truncated: */
      ret[i] = 0.0;
    }
  }
  UNPROTECT(6);
  return s_ret;
}


SEXP ptruncnorm(SEXP args) {
  R_len_t i, n;
  SEXP s_q, s_a, s_b, s_mean, s_sd, s_ret;
  double *q, a, b, mean, sd, cx, ca, cb, *ret;
  
  UNPACK_REAL_VECTOR(args, s_q, q);
  UNPACK_REAL_VALUE(args, s_a, a);
  UNPACK_REAL_VALUE(args, s_b, b);
  UNPACK_REAL_VALUE(args, s_mean, mean);
  UNPACK_REAL_VALUE(args, s_sd, sd);

  n = length(s_q);
  ALLOC_REAL_VECTOR(n, s_ret, ret);

  ca = pnorm(a, mean, sd, FALSE, FALSE);
  cb = pnorm(b, mean, sd, FALSE, FALSE);
  
  for (i = 0; i < n; ++i) {
    if (q[i] < a) {
      ret[i] = 0.0;
    } else if (q[i] > b) {
      ret[i] = 1.0;
    } else {
      cx = pnorm(q[i], mean, sd, FALSE, FALSE);
      ret[i] = (cx - ca) / (cb - ca);
    }
  }
  UNPROTECT(6);
  return s_ret;
}


SEXP rtruncnorm(SEXP args) {
 R_len_t i;
 int n;
 SEXP s_n, s_a, s_b, s_mean, s_sd, s_ret;
 double a, b, mean, sd, *ret, tmp;

 UNPACK_INTEGER_VALUE(args, s_n, n);
 UNPACK_REAL_VALUE(args, s_a, a);
 UNPACK_REAL_VALUE(args, s_b, b);
 UNPACK_REAL_VALUE(args, s_mean, mean);
 UNPACK_REAL_VALUE(args, s_sd, sd);
 
 ALLOC_REAL_VECTOR(n, s_ret, ret);
 GetRNGstate();
 for (i = 0; i < n; ++i) {
   while (1) {
     tmp = rnorm(mean, sd);
     if (a <= tmp && tmp <= b)
       break;
   }
   ret[i] = tmp;
 }
 PutRNGstate();
 UNPROTECT(6);
 return s_ret;
}


SEXP etruncnorm(SEXP args) {
  R_len_t i, na, nb;
  SEXP s_a, s_b, s_mean, s_sd, s_ret;
  double *a, *b, mean, sd, *ret;

  UNPACK_REAL_VECTOR(args, s_a, a);
  UNPACK_REAL_VECTOR(args, s_b, b);
  UNPACK_REAL_VALUE(args, s_mean, mean);
  UNPACK_REAL_VALUE(args, s_sd, sd);
  
  na = length(s_a);
  nb = length(s_b);
  
  if (na != nb) 
    error("Length of a and b differ (%i != %i)", na, nb);

  ALLOC_REAL_VECTOR(na, s_ret, ret);
  
  for (i = 0; i < na; ++i) {
    double ca = dnorm(a[i], mean, sd, FALSE);
    double cb = dnorm(b[i], mean, sd, FALSE);
    double Ca = pnorm(a[i], mean, sd, FALSE, FALSE);
    double Cb = pnorm(b[i], mean, sd, FALSE, FALSE);

    ret[i] = mean + sd * ((ca - cb) / (Cb - Ca));
  } 
  UNPROTECT(5);
  return s_ret;
}
