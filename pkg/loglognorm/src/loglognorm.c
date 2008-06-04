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

#ifndef MAX
#define MAX(A, B) ((A > B)?(A):(B))
#endif

#define UNPACK_REAL_VECTOR(S, D, N)	\
  double *D = REAL(S);                  \
  R_len_t N = length(S);

#define ALLOC_REAL_VECTOR(SXP, DBL, SIZE)    \
  SEXP SXP;				     \
  PROTECT(SXP = allocVector(REALSXP, SIZE)); \
  double *DBL = REAL(SXP);


SEXP dloglognorm(SEXP s_x, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_x   , x   , n_x);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  
  n = MAX(MAX(n_x, n_mean), n_sd);
  ALLOC_REAL_VECTOR(s_ret, ret, n);

  for (i = 0; i < n; ++i) {
    const double cx = x[i % n_x];
    if (0 <= cx && cx <= 1) {
      const double cmean = mean[i % n_mean];
      const double csd = sd[i % n_sd];

      const double c1 = -1.0/(2.0 * csd * csd);
      const double c2 = log(cx);
      const double c3 = -M_1_SQRT_2PI / (csd * c2 * cx);
      const double c4 = c1 * pow(log(-c2) - cmean, 2);
      ret[i] = c3 * exp(c4);
    } else {
      ret[i] = 0.0;
    }
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}


SEXP ploglognorm(SEXP s_q, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_q   , q   , n_q);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  
  n = MAX(MAX(n_q, n_mean), n_sd);
  ALLOC_REAL_VECTOR(s_ret, ret, n);

  for (i = 0; i < n; ++i) {
    const double cq = q[i % n_q];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    if (cq < 0.0) { 
      ret[i] = 0.0;
    } else if (cq > 1.0) {
      ret[i] = 1.0;
    } else { /* q \in [0, 1] */
      /* Directly return the upper tail: */  
      ret[i] = pnorm(log(-log(cq)), cmean, csd, FALSE, FALSE);      
    }
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}


SEXP qloglognorm(SEXP s_p, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_p   , p   , n_p);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  
  n = MAX(MAX(n_p, n_mean), n_sd);
  ALLOC_REAL_VECTOR(s_ret, ret, n);
  
  for (i = 0; i < n; ++i) {
    const double cp = p[i % n_p];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    if (0 <= cp && cp <= 1.0) { /* p \in [0, 1] */
      ret[i] = exp(-exp(qnorm(1-cp, cmean, csd, TRUE, FALSE)));
    } else { 
      ret[i] = 0.0;
    }      
  }
  UNPROTECT(1); /* s_ret */
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

SEXP mloglognorm(SEXP s_moment, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_moment, moment, n_moment);
  UNPACK_REAL_VECTOR(s_mean  , mean  , n_mean);
  UNPACK_REAL_VECTOR(s_sd    , sd    , n_sd);
  
  n = MAX(MAX(n_moment, n_mean), n_sd);
  ALLOC_REAL_VECTOR(s_ret, ret, n);
  
  /* Parameters for Rdqagi: */
  int limit = 100;
  int lenw = 4 * limit;
  int *iwork = (int *) R_alloc(limit, sizeof(int));
  double *work = (double *) R_alloc(lenw,  sizeof(double));
  loglognorm_param lp;
  double bound = 0; /* not used */
  int inf = 2; /* 2 = -Inf - inf */
  double epsabs = 0.00000001;
  double epsrel = epsabs;
  double result, abserr;
  int neval, ier, last;
      
  for (i = 0; i < n_mean; ++i) {
    lp.mean = mean[i % n_mean];
    lp.sd = sd[i % n_sd];
    lp.r = moment[i % n_moment];
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
  UNPROTECT(1); /* s_ret */
  return s_ret;
}

