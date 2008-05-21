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
  double *x, a, b, mean, sd, *ret;

  UNPACK_REAL_VECTOR(args, s_x, x);
  UNPACK_REAL_VALUE(args, s_a, a);
  UNPACK_REAL_VALUE(args, s_b, b);
  UNPACK_REAL_VALUE(args, s_mean, mean);
  UNPACK_REAL_VALUE(args, s_sd, sd);

  n = length(s_x);
  ALLOC_REAL_VECTOR(n, s_ret, ret);
  const double c1 = pnorm((a-mean)/sd, 0.0, 1.0, TRUE, FALSE);
  const double c2 = pnorm((b-mean)/sd, 0.0, 1.0, TRUE, FALSE);
  const double c3 = sd * (c2 - c1);
  for (i = 0; i < n; ++i) {
    const double cx = x[i];
    if (a <= cx && cx <= b) { /* In range: */
      const double c4 = dnorm((cx - mean)/sd, 0.0, 1.0, FALSE);
      ret[i] = c4 / c3;
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
  double *q, a, b, mean, sd, ca, cb, *ret;
  
  UNPACK_REAL_VECTOR(args, s_q, q);
  UNPACK_REAL_VALUE(args, s_a, a);
  UNPACK_REAL_VALUE(args, s_b, b);
  UNPACK_REAL_VALUE(args, s_mean, mean);
  UNPACK_REAL_VALUE(args, s_sd, sd);

  n = length(s_q);
  ALLOC_REAL_VECTOR(n, s_ret, ret);

  ca = pnorm(a, mean, sd, TRUE, FALSE);
  cb = pnorm(b, mean, sd, TRUE, FALSE);
  
  for (i = 0; i < n; ++i) {
    if (q[i] < a) {
      ret[i] = 0.0;
    } else if (q[i] > b) {
      ret[i] = 1.0;
    } else {
      const double cx = pnorm(q[i], mean, sd, TRUE, FALSE);
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
  R_len_t i, n_a, n_b, n_mean, n_sd;
  SEXP s_a, s_b, s_mean, s_sd, s_ret;
  double *a, *b, *mean, *sd, *ret;

  UNPACK_REAL_VECTOR(args, s_a, a);
  UNPACK_REAL_VECTOR(args, s_b, b);
  UNPACK_REAL_VECTOR(args, s_mean, mean);
  UNPACK_REAL_VECTOR(args, s_sd, sd);
  
  n_a = length(s_a);
  n_b = length(s_b);
  n_mean = length(s_mean);
  n_sd = length(s_sd);  
  if (n_a != n_b || n_b != n_mean || n_mean != n_sd) 
    error("Length of a, b, mean or sd differ.");

  ALLOC_REAL_VECTOR(n_a, s_ret, ret);
  
  for (i = 0; i < n_a; ++i) {
    const double ca = dnorm(a[i], mean[i], sd[i], FALSE);
    const double cb = dnorm(b[i], mean[i], sd[i], FALSE);
    const double Ca = pnorm(a[i], mean[i], sd[i], TRUE, FALSE);
    const double Cb = pnorm(b[i], mean[i], sd[i], TRUE, FALSE);
    
    ret[i] = mean[i] + sd[i] * ((ca - cb) / (Cb - Ca));
  } 
  UNPROTECT(5);
  return s_ret;
}


SEXP vtruncnorm(SEXP args) {
  R_len_t i, n_a, n_b, n_mean, n_sd;
  SEXP s_a, s_b, s_mean, s_sd, s_ret;
  double *a, *b, *mean, *sd, *ret;
  
  UNPACK_REAL_VECTOR(args, s_a, a);
  UNPACK_REAL_VECTOR(args, s_b, b);
  UNPACK_REAL_VECTOR(args, s_mean, mean);
  UNPACK_REAL_VECTOR(args, s_sd, sd);
  
  n_a = length(s_a);
  n_b = length(s_b);
  n_mean = length(s_mean);
  n_sd = length(s_sd);
  if (n_a != n_b || n_b != n_mean || n_mean != n_sd) 
    error("Length of a, b, mean or sd differ.");

  ALLOC_REAL_VECTOR(n_a, s_ret, ret);
  
  for (i = 0; i < n_a; ++i) {
    const double am = (a[i] - mean[i])/sd[i];
    const double bm = (b[i] - mean[i])/sd[i];
    const double ca = dnorm(am, 0.0, 1.0, FALSE);
    const double cb = dnorm(bm, 0.0, 1.0, FALSE);
    const double Ca = pnorm(am, 0.0, 1.0, TRUE, FALSE);
    const double Cb = pnorm(bm, 0.0, 1.0, TRUE, FALSE);
    const double v  = sd[i] * sd[i];
    
    const double d = Cb - Ca;
    const double m1 = (ca - cb)/d;
    const double m2 = (am*ca - bm*cb)/d;
    ret[i] = (1.0 + m2 + m1*m1)*v;
  } 
  UNPROTECT(5);
  return s_ret;
}
