/*
 * loglognorm.c - Implementation of double log normal functions
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


SEXP dloglognorm(SEXP args) {
  R_len_t i, n;
  SEXP s_x, s_mean, s_sd, s_ret;
  double *x, mean, sd, *ret, c1, c2, c3, c4;
  
  UNPACK_REAL_VECTOR(args, s_x, x);
  UNPACK_REAL_VALUE(args, s_mean, mean);
  UNPACK_REAL_VALUE(args, s_sd, sd);
  
  n = length(s_x);
  ALLOC_REAL_VECTOR(n, s_ret, ret);

  c1 = - M_1_SQRT_2PI / sd;
  c2 = -1.0/(2.0 * sd * sd);
  for (i = 0; i < n; ++i) {
    if (0 <= x[i] && x[i] <= 1) {
      c3 = c1 / (log(x[i]) * x[i]);
      c4 = c2 * pow(log(-log(x[i])) - mean, 2);
      ret[i] = c3 * exp(c4);
    } else {
      ret[i] = 0.0;
    }
  }
  UNPROTECT(4);
  return s_ret;
}


SEXP ploglognorm(SEXP args) {
  R_len_t i, n;
  SEXP s_q, s_mean, s_sd, s_ret;
  double *q, mean, sd, *ret;
  
  UNPACK_REAL_VECTOR(args, s_q, q);
  UNPACK_REAL_VALUE(args, s_mean, mean);
  UNPACK_REAL_VALUE(args, s_sd, sd);
  
  n = length(s_q);
  ALLOC_REAL_VECTOR(n, s_ret, ret);

  for (i = 0; i < n; ++i) {
    if (q[i] < 0.0) {
      ret[i] = 0.0;
    } else if (q[i] > 1.0) {
      ret[i] = 1.0;
    } else {
      /* FIXME: Better use log=TRUE? */
      ret[i] = 1.0 - pnorm((log(-log(q[i])) - mean)/sd, 0.0, 1.0, FALSE, FALSE);    
    }
  }
  UNPROTECT(4);
  return s_ret;
}


SEXP qloglognorm(SEXP args) {
  R_len_t i, n;
  SEXP s_p, s_mean, s_sd, s_ret;
  double *p, mean, sd, *ret;
  
  UNPACK_REAL_VECTOR(args, s_p, p);
  UNPACK_REAL_VALUE(args, s_mean, mean);
  UNPACK_REAL_VALUE(args, s_sd, sd);
  
  n = length(s_p);
  ALLOC_REAL_VECTOR(n, s_ret, ret);
  
  for (i = 0; i < n; ++i) {
    if (0 <= p[i] && p[i] <= 1.0) {
      ret[i] = exp(-exp(sd * qnorm(1-p[i], 0.0, 1.0, FALSE, FALSE) + mean));
    } else {
      ret[i] = R_NaN;
    }      
  }
  UNPROTECT(4);
  return s_ret;
}


typedef struct {
  double mean, sd, r;
} loglognorm_param;

static void loglognorm_intgr(double *x, int n, void *ex) {
  int i;
  loglognorm_param *lp = (loglognorm_param *)ex;
  const double mean = lp->mean;
  const double sd = lp->sd;
  const double r = lp->r;

  /* Taken from Trautmann (2004) p. 54 */
  for (i = 0; i < n; ++i) {
    x[i] = exp(-r*exp(mean + sd * x[i])) * M_1_SQRT_2PI * exp(-0.5 * pow(x[i], 2.0));
  }
}

SEXP mloglognorm(SEXP args) {
  R_len_t i, n_mean, n_sd, n_moment;
  SEXP s_mean, s_sd, s_moment, s_ret;
  double *mean, *sd, *moment, *ret, tmp;
  
  UNPACK_REAL_VECTOR(args, s_mean, mean);
  UNPACK_REAL_VECTOR(args, s_sd, sd);
  UNPACK_REAL_VECTOR(args, s_moment, moment);
  n_mean = length(s_mean);
  n_sd = length(s_sd);
  n_moment = length(s_moment);
  if (n_mean != n_sd)
    error("Length of mean, sd and moment differ (%i != %i)", n_mean, n_sd);

  ALLOC_REAL_VECTOR(n_mean, s_ret, ret);

  /* Parameters for Rdqagi() */
  int limit = 100;
  int lenw = 4 * limit;
  int *iwork = (int *) R_alloc(limit, sizeof(int));
  double *work = (double *) R_alloc(lenw,  sizeof(double));
  loglognorm_param lp;
  double bound = 0; /* not used */
  int inf = 2; /* 2 = -1nf - inf */
  double epsabs = 0.0000001;
  double epsrel = epsabs;
  double result, abserr;
  int neval, ier, last;
      
  for (i = 0; i < n_mean; ++i) {
    lp.mean = mean[i];
    lp.sd = sd[i];
    lp.r = moment[i];
    Rdqagi(loglognorm_intgr, (void *)&lp, &bound, &inf,
	   &epsabs, &epsrel, 
	   &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    /* FIXME: Possibly check agains lower bound given in
     * Trautmann (2004):
     *
     *   E(X^r) >= exp(-r) * Phi(-mean/sd)
     */
    if (ier >= 1) { /* Failure */
      ret[i] = R_NaN;
    } else { /* No error */
      ret[i] = result;
    }
  }
  UNPROTECT(4); /* s_mean, s_sd, s_moment, s_ret */
  return s_ret;
}

