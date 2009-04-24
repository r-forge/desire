/*
 * loglognorm.c - Implementation of double log normal functions
 *
 * Authors:
 *  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
 *  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
 *  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#ifndef MAX
#define MAX(A, B) ((A > B)?(A):(B))
#endif

#define UNPACK_REAL_VECTOR(S, D, N)		\
  double *D = REAL(S);				\
  const R_len_t N = length(S);

/*
 * derringer_suich(x, y, d, beta)
 *
 * Assumes 
 *  + length(y) = length(d) = length(beta) - 1 
 *  + y is sorted
 */
SEXP ds_eval(SEXP s_x, SEXP s_y, SEXP s_d, SEXP s_beta) {
  R_len_t i, j;
  double *x, *y, *d, *beta, *ret;
  SEXP s_ret;

  /* Unpack arguments */  
  const R_len_t k = length(s_x);
  const R_len_t n = length(s_y);
  x = REAL(s_x);
  y = REAL(s_y);
  d = REAL(s_d);
  beta = REAL(s_beta);

  /* Allocate return vector */
  PROTECT(s_ret = allocVector(REALSXP, k));
  ret = REAL(s_ret);

  /* Declare constants */
  const double ymin = y[0];
  const double ymax = y[n-1];
  /* Iterate over input */
  for (i = 0; i < k; ++i) {
    if (x[i] < ymin || x[i] > ymax) {
      ret[i] = 0.0;
    } else {
      /* Find j so that y[j] is smaller than x[i] */
      for (j=0; y[j] < x[i]; ++j) {};
      const double ym = y[j];
      const double yp = y[j-1];
      if (ym == R_PosInf || yp == R_NegInf) { /* LB Type edge cases */
	ret[i] = 1.0;
      } else {
	const double dm = d[j];
	const double dp = d[j-1];
	const double b = beta[j-1];
	if (dm <= dp) {
	  const double c1 = (x[i] - ym)/(yp - ym);
	  ret[i] = dm + (dp - dm)*pow(c1, b);
	} else {
	  const double c1 = (x[i] - yp)/(ym - yp);
	  ret[i] = dp + (dm - dp)*pow(c1, b);
	}
      }
    }
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}


/* 
 * dsLTU11:
 */

SEXP ddsLTU11(SEXP s_x, SEXP s_l, SEXP s_t, SEXP s_u, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  SEXP s_res;
  /* Unpack arguments */
  UNPACK_REAL_VECTOR(s_x   , x   , n_x);
  UNPACK_REAL_VECTOR(s_l   , l   , n_l);
  UNPACK_REAL_VECTOR(s_t   , t   , n_t);
  UNPACK_REAL_VECTOR(s_u   , u   , n_u);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  /* Maximum argument length == result size */
  n = MAX(MAX(n_mean, n_sd), MAX(MAX(n_x, n_l), MAX(n_t, n_u)));

  /* Allocate result vector */
  PROTECT(s_res = allocVector(REALSXP, n));
  double *res = REAL(s_res);
  
  for (i = 0; i < n; ++i) {
    const double cx = x[i % n_x];
    if (cx < 0.0 || cx >= 1.0) { /* x \nin [0, 1) */
      res[i] = 0.0;
    } else {
      const double cmean = mean[i % n_mean];
      const double csd = sd[i % n_sd];
      const double cl = l[i % n_l];
      const double ct = t[i % n_t];
      const double cu = u[i % n_u];
      if (cx == 0.0) { /* x == 0 */
	const double c1 = pnorm((cl - cmean)/csd, 0.0, 1.0, TRUE, FALSE);
	const double c2 = pnorm((cu - cmean)/csd, 0.0, 1.0, TRUE, FALSE);
	res[i] = c1 + 1 - c2;
      } else { /* x \in (0, 1) */
	const double c1 = dnorm((cl + cx * (ct - cl) - cmean)/csd, 0.0, 1.0, FALSE);
	const double c2 = dnorm((cu - cx * (cu - ct) - cmean)/csd, 0.0, 1.0, FALSE);
	res[i] = (ct - cl)/csd * c1 + (cu - ct)/csd * c2;
      }
    }
  }
  UNPROTECT(1); /* s_res */
  return s_res;
}


SEXP pdsLTU11(SEXP s_q, SEXP s_l, SEXP s_t, SEXP s_u, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  SEXP s_res;
  /* Unpack arguments */
  UNPACK_REAL_VECTOR(s_q   , q   , n_q);
  UNPACK_REAL_VECTOR(s_l   , l   , n_l);
  UNPACK_REAL_VECTOR(s_t   , t   , n_t);
  UNPACK_REAL_VECTOR(s_u   , u   , n_u);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  /* Maximum argument length == result size */
  n = MAX(MAX(n_mean, n_sd), MAX(MAX(n_q, n_l), MAX(n_t, n_u)));
  
  /* Allocate result vector */
  PROTECT(s_res = allocVector(REALSXP, n));
  double *res = REAL(s_res);
  
  for(i = 0; i < n; ++i) {
    const double cq = q[i % n_q];
    if (cq < 0.0) {
      res[i] = 0.0;
    } else if (cq > 1.0) {
      res[i] = 1.0;
    } else {
      const double cl = l[i % n_l];
      const double ct = t[i % n_t];
      const double cu = u[i % n_u];
      const double cmean = mean[i % n_mean];
      const double csd = sd[i % n_sd];
      
      const double c1 = cl + cq * (ct - cl) - cmean;
      const double c2 = pnorm(c1/csd, 0.0, 1.0, TRUE, FALSE);
      const double c3 = cu - cq * (cu - ct) - cmean;
      const double c4 = pnorm(c3/csd, 0.0, 1.0, TRUE, FALSE);
      res[i] = c2 + 1 - c4;
    }
  }
  UNPROTECT(1); /* s_res */
  return s_res;
}


SEXP edsLTU11(SEXP s_l, SEXP s_t, SEXP s_u, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  SEXP s_res;
  /* Unpack arguments */
  UNPACK_REAL_VECTOR(s_l   , l   , n_l);
  UNPACK_REAL_VECTOR(s_t   , t   , n_t);
  UNPACK_REAL_VECTOR(s_u   , u   , n_u);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  /* Maximum argument length == result size */
  n = MAX(MAX(n_mean, n_sd), MAX(n_l, MAX(n_t, n_u)));
  
  PROTECT(s_res = allocVector(REALSXP, n));
  double *res = REAL(s_res);

  for(i = 0; i < n; ++i) {
    const double cl = l[i % n_l];
    const double ct = t[i % n_t];
    const double cu = u[i % n_u];
    const double cy = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    
    const double yml = cy - cl;
    const double tml = ct - cl;
    const double umy = cu - cy;
    const double umt = cu - ct;
 
    const double ncl = (cl - cy)/csd;
    const double nct = (ct - cy)/csd;
    const double ncu = (cu - cy)/csd;
    
    const double pl = pnorm(ncl, 0.0, 1.0, TRUE, FALSE);
    const double dl = dnorm(ncl, 0.0, 1.0, FALSE);
    const double pt = pnorm(nct, 0.0, 1.0, TRUE, FALSE);
    const double dt = dnorm(nct, 0.0, 1.0, FALSE);
    const double pu = pnorm(ncu, 0.0, 1.0, TRUE, FALSE);
    const double du = dnorm(ncu, 0.0, 1.0, FALSE);

    const double Gl = pt - pl;
    const double gl = dl - dt;
    const double Gr = pu - pt;
    const double gr = dt - du;

    res[i] = (Gl*yml + gl*csd)/tml  + (Gr*umy - gr*csd)/umt;
  }
  UNPROTECT(1); /* s_res */
  return s_res;
}


/* 
 * dsLTI11:
 */
SEXP ddsLTI11(SEXP s_x, SEXP s_l, SEXP s_t, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  SEXP s_res;
  /* Unpack arguments */
  UNPACK_REAL_VECTOR(s_x   , x   , n_x);
  UNPACK_REAL_VECTOR(s_l   , l   , n_l);
  UNPACK_REAL_VECTOR(s_t   , t   , n_t);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  /* Maximum argument length == result size */
  n = MAX(MAX(n_mean, n_sd), MAX(MAX(n_x, n_l),  n_t));

  /* Allocate result vector */
  PROTECT(s_res = allocVector(REALSXP, n));
  double *res = REAL(s_res);
  
  for (i = 0; i < n; ++i) {
    const double cx = x[i % n_x];
    if (cx < 0.0 || cx >= 1.0) { /* x \nin [0, 1) */
      res[i] = 0.0;
    } else {
      const double cmean = mean[i % n_mean];
      const double csd = sd[i % n_sd];
      const double cl = l[i % n_l];
      const double ct = t[i % n_t];
      if (cx == 0.0) { /* x == 0 */
	res[i] = pnorm((cl - cmean)/csd, 0.0, 1.0, TRUE, FALSE);
      } else if (cx == 1.0) { /* x == 0 */
	res[i] = 1.0 - pnorm((ct - cmean)/csd, 0.0, 1.0, TRUE, FALSE);
      } else { /* x \in (0, 1) */
	const double c1 = dnorm((cl + cx * (ct - cl) - cmean)/csd, 0.0, 1.0, FALSE);
	res[i] = (ct - cl)/csd * c1;
      }
    }
  }
  UNPROTECT(1); /* s_res */
  return s_res;
}


SEXP pdsLTI11(SEXP s_q, SEXP s_l, SEXP s_t, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  SEXP s_res;
  /* Unpack arguments */
  UNPACK_REAL_VECTOR(s_q   , q   , n_q);
  UNPACK_REAL_VECTOR(s_l   , l   , n_l);
  UNPACK_REAL_VECTOR(s_t   , t   , n_t);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  /* Maximum argument length == result size */
  n = MAX(MAX(n_mean, n_sd), MAX(MAX(n_q, n_l), n_t));

  /* Allocate result vector */
  PROTECT(s_res = allocVector(REALSXP, n));
  double *res = REAL(s_res);
  
  for (i = 0; i < n; ++i) {
    const double cq = q[i % n_q];
    if (cq < 0.0) {
      res[i] = 0.0;
    } else if (cq >= 1.0) {
      res[i] = 1.0;
    } else {
      const double cmean = mean[i % n_mean];
      const double csd = sd[i % n_sd];
      const double cl = l[i % n_l];
      const double ct = t[i % n_t];
      res[i] = pnorm((cl + cq*(ct - cl) - cmean)/csd, 0.0, 1.0, TRUE, FALSE);
    }
  }
  UNPROTECT(1); /* s_res */
  return s_res;
}


SEXP edsLTI11(SEXP s_l, SEXP s_t, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  SEXP s_res;
  /* Unpack arguments */
  UNPACK_REAL_VECTOR(s_l   , l   , n_l);
  UNPACK_REAL_VECTOR(s_t   , t   , n_t);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  /* Maximum argument length == result size */
  n = MAX(MAX(n_mean, n_sd), MAX(n_l, n_t));

  /* Allocate result vector */
  PROTECT(s_res = allocVector(REALSXP, n));
  double *res = REAL(s_res);
  
  for (i = 0; i < n; ++i) { /* see Steuer p. 74 */
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    const double cl = l[i % n_l];
    const double ct = t[i % n_t];
    
    const double dml = cmean - cl;
    const double dtl = ct - cl;
    
    const double qlm = (cl - cmean)/csd;
    const double qtm = (ct - cmean)/csd;
    
    const double ptm = pnorm(qtm, 0.0, 1.0, TRUE, FALSE);
    const double dtm = dnorm(qtm, 0.0, 1.0, FALSE);

    const double Gl = ptm - pnorm(qlm, 0.0, 1.0, TRUE, FALSE);
    const double gl = dnorm(qlm, 0.0, 1.0, FALSE) - dtm;

    res[i] = (Gl * dml + gl*csd)/dtl + 1.0 - ptm;
  }
  UNPROTECT(1); /* s_res */
  return s_res;
}

/*
 * dsA1
 *
 * Assumptions:
 *  length(y) == length(d)
 */

SEXP edsA1(SEXP s_y, SEXP s_d, SEXP s_mean, SEXP s_sd) {
  R_len_t i, k, n;
  SEXP s_res;
  /* Unpack arguments */
  UNPACK_REAL_VECTOR(s_y   , y   , n_y);
  UNPACK_REAL_VECTOR(s_d   , d   , n_d);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  /* Maximum argument length == result size */
  n = MAX(n_mean, n_sd);
  
  /* Allocate result vector */
  PROTECT(s_res = allocVector(REALSXP, n));
  double *res = REAL(s_res);
  
  for (i = 0; i < n; ++i) { /* see Steuer p. 88 */
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    res[i] = 0.0;
    for (k = 1; k < n_y; ++k) {
      const double yk  = y[k];
      const double ykm = y[k-1];

      const double dk  = d[k];
      const double dkm = d[k-1];

      const double pyk  = pnorm(yk, cmean, csd, TRUE, FALSE);
      const double pykm = pnorm(ykm, cmean, csd, TRUE, FALSE);
      const double pyd  = pyk - pykm;
      if (dkm == dk) {
	res[i] += dk * pyd;
      } else if (pyd != 0){
	/* We skip parts where pyd is numerically zero because
	 * they would lead to a division by zero and therefor a 
	 * result of NaN. Since their weight is almost zero, this
	 * does not change the result by much...
	 */
	const double dyk  = dnorm(yk, cmean, csd, FALSE);
	const double dykm = dnorm(ykm, cmean, csd, FALSE);
	const double frac = csd * csd * (dykm - dyk)/pyd;
	if (dkm > dk) {
	  const double c1 = (dk - dkm)/(yk - ykm);
	  res[i] += (dkm + c1 * (cmean + frac - ykm)) * pyd;
	} else { /* dkm > dk */
	  const double c1 = (dkm - dk)/(yk - ykm);
	  res[i] += (dk + c1 * (cmean + frac + yk)) * pyd;
	} 
      }
    }
  }
  UNPROTECT(1); /* s_res */
  return s_res;
}
