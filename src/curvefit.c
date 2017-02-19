
#include "rjmcmc/rjmcmc_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rjmcmc/curvefit.h"

#include "rjmcmc/rjmcmc_util.h"
#include "rjmcmc/rjmcmc_debug.h"

#define M_PI 3.14159265358979323846

struct _curvefit_result {
  
  int maxorder;

  double *alpha;
  double *beta;

  double **L;

  /* Contains Cholesky Decomposition of Covariance Matrix */
  double **Z;

  /* Contains Inverse Covariance Matrix (sigma^-1 or Cm^-1) */
  double **S;
  double **Si;

  /* Contains Mean best fit */
  double *mu;

  double *x;
  double *b;

};

static int
compute_hankel(const point1d_t *points,
	       int n,
	       double lambda,
	       int order,
	       curvefit_result_t *cf);

static void
set_I(double **m,
      int n);

static int
compute_hankel_cholesky(const double *alpha,
			int order,
			double **L);

static int
compute_forward_substitution(double **L,
			     int order,
			     double *b,
			     double *x);

static int
compute_backward_substitution(double **L,
			      int order,
			      double *b,
			      double *x);


static int
compute_inverse(double **L,
		int order,
		double **Z,
		double *b,
		double *x);

static int
compute_square(double **L,
	       int order,
	       double **L2);

	       
static int
compute_cholesky(double **A,
		 int order,
		 double **L);

curvefit_result_t *
curvefit_create(int maxorder)
{
  curvefit_result_t *r;

  r = (curvefit_result_t*)malloc(sizeof(curvefit_result_t));
  if (r == NULL) {
    return NULL;
  }

  r->maxorder = maxorder;

  r->alpha = rjmcmc_create_array_1d(2*(maxorder + 1));
  if (r->alpha == NULL) {
    return NULL;
  }

  r->beta = rjmcmc_create_array_1d(2*(maxorder + 1));
  if (r->beta == NULL) {
    return NULL;
  }

  r->L = rjmcmc_create_array_2d(maxorder + 1, maxorder + 1);
  if (r->L == NULL) {
    return NULL;
  }

  r->Z = rjmcmc_create_array_2d(maxorder + 1, maxorder + 1);
  if (r->Z == NULL) {
    return NULL;
  }

  r->S = rjmcmc_create_array_2d(maxorder + 1, maxorder + 1);
  if (r->S == NULL) {
    return NULL;
  }

  r->Si = rjmcmc_create_array_2d(maxorder + 1, maxorder + 1);
  if (r->Si == NULL) {
    return NULL;
  }

  r->mu = rjmcmc_create_array_1d(maxorder + 1);
  if (r->mu == NULL) {
    return NULL;
  }

  r->x = rjmcmc_create_array_1d(2*(maxorder + 1));
  if (r->x == NULL) {
    return NULL;
  }

  r->b = rjmcmc_create_array_1d(2*(maxorder + 1));
  if (r->b == NULL) {
    return NULL;
  }

  return r;
}

void
curvefit_destroy(curvefit_result_t *cf)
{
  if (cf != NULL) {
    rjmcmc_destroy_array_1d(cf->alpha);
    rjmcmc_destroy_array_1d(cf->beta);
    rjmcmc_destroy_array_2d(cf->maxorder + 1, cf->L);
    rjmcmc_destroy_array_2d(cf->maxorder + 1, cf->Z);
    rjmcmc_destroy_array_2d(cf->maxorder + 1, cf->S);
    rjmcmc_destroy_array_2d(cf->maxorder + 1, cf->Si);

    rjmcmc_destroy_array_1d(cf->mu);

    rjmcmc_destroy_array_1d(cf->x);
    rjmcmc_destroy_array_1d(cf->b);

    free(cf);
  }
}

int
curvefit_compute(const dataset1d_t *d,
		 int di,
		 int dj,
		 int order,
		 curvefit_result_t *cf)
{
  int n;

  if (cf == NULL) {
    rjmcmc_error("curvefit_compute: result is not allocated\n");
    return -1;
  }

  if (order > cf->maxorder) {
    rjmcmc_error("curvefit_compute: requested order is too large (%d > %d)\n",
		 order, cf->maxorder);
    return -1;
  }
  
  if (dj <= di) {
    rjmcmc_error("curvefit_compute: invalid range\n");
    return -1;
  }

  n = dj - di + 1;
  if (n < order) {
    /* Insufficient no. of points to solve */
    rjmcmc_error("curvefit_compute: insufficient points for order\n");
    return -1;
  }
  
  compute_hankel(d->points + di,
		 n,
		 1.0,
		 order,
		 cf);

  if (compute_hankel_cholesky(cf->alpha,
			      order, 
			      cf->L) < 0) {
    /* Singular/Near singular */
    rjmcmc_error("curvefit_compute: "
		 "failed to compute hankel/cholesky for order %d\n",
		 order);
    return -1;
  }

  compute_square(cf->L,
		 order,
		 cf->Si);

  compute_forward_substitution(cf->L, order, cf->beta, cf->x);
  compute_backward_substitution(cf->L, order, cf->x, cf->mu);

  compute_inverse(cf->L,
		  order,
		  cf->Z,
		  cf->b,
		  cf->x);

  compute_square(cf->Z,
		 order,
		 cf->S);

  if (compute_cholesky(cf->S,
		       order,
		       cf->Z) < 0) {
    /* rjmcmc_error( */
    /* 	    "curvefit_compute: failed to compute cholesky\n"); */
    return -1;
  }

  return 0;
}

int
curvefit_compute_lambda(const dataset1d_t *d,
			double lambda,
			int di,
			int dj,
			int order,
			curvefit_result_t *cf)
{
  int n;

  if (cf == NULL) { 
    rjmcmc_error("curvefit_compute_lambda: result not allocated\n");
    return -1;
  }

  if (order > cf->maxorder) { 
    rjmcmc_error("curvefit_compute_lambda: "
		 "requested order is too large (%d > %d)\n",
		 order, cf->maxorder);
    return -1;
  }
  
  if (dj <= di) {
    rjmcmc_error("curvefit_compute_lambda: invalid range (%d %d)\n", di, dj);
    return -1;
  }

  n = dj - di + 1;
  if (n < order) {
    /* Insufficient no. of points to solve */
    rjmcmc_error("curvefit_compute_lambda: insufficient points\n");
    return -1;
  }
  
  compute_hankel(d->points + di,
		 n,
		 lambda,
		 order,
		 cf);

  if (compute_hankel_cholesky(cf->alpha,
			      order, 
			      cf->L) < 0) {
    /* rjmcmc_error("curvefit_compute_lambda: failed to compute hankel cholesky\n (%d %d %d)\n", di, dj, order); */
    return -1;
  }

  compute_forward_substitution(cf->L, order, cf->beta, cf->x);
  compute_backward_substitution(cf->L, order, cf->x, cf->mu);

  compute_inverse(cf->L,
		  order,
		  cf->Z,
		  cf->b,
		  cf->x);

  compute_square(cf->Z,
		 order,
		 cf->S);

  if (compute_cholesky(cf->S,
		       order,
		       cf->Z) < 0) {
    /* rjmcmc_error("curvefit_compute_lambda: failed to compute cholesky\n"); */
    return -1;
  }

  return 0;
}

int 
curvefit_sample(curvefit_result_t *cf,
		rjmcmc_normal_rand_t normal,
		double *coeff,
		int coeff_len,
		double *prob)
{
  int i;
  int j;

  double sigmah;
  double sigma2;
  double sigma2t;

  /* printf("N: "); */
  for (i = 0; i < coeff_len; i ++) {
    cf->x[i] = normal();
    /* printf("%g ", cf->x[i]); */
  }
  /* printf("\n"); */

  /*
   * Work out the coefficients 
   */
  /* printf("b: "); */
  for (i = 0; i < coeff_len; i ++) {
    coeff[i] = cf->mu[i];

    for (j = 0; j < coeff_len; j ++) {
      coeff[i] += cf->x[j] * cf->Z[i][j];
    }

    /* b is the x - mu vector */
    cf->b[i] = coeff[i] - cf->mu[i];
    /* printf("%g ", cf->b[i]); */
  }
  /* printf("\n"); */

  /* printf("--- coeff ---\n"); */
  /* for (i = 0; i < coeff_len; i ++) { */
  /*   printf("%f\n", coeff[i]); */
  /* } */
  /* printf("-------------\n"); */
				    
  /*
   * Work out (x - mu)T sigma^-1 (x - mu)
   */
  sigma2 = 0.0;
  /* printf("S:\n"); */
  for (i = 0; i < coeff_len; i ++) {
    sigma2t = 0.0;
    for (j = 0; j < coeff_len; j ++) {
      sigma2t += cf->Si[i][j]*cf->b[j];
      /* printf("  %g", cf->Si[i][j]); */
    }
    /* printf("\n"); */
    sigma2 += cf->b[i] * sigma2t;
  }

  /* printf("Si:\n"); */
  for (i = 0; i < coeff_len; i ++) {
    for (j = 0; j < coeff_len; j ++) {
      /* printf("  %g", cf->Si[i][j]); */
    }
    /* printf("\n"); */
  }      

  /*
   * Work out det of cholesky
   */
  sigmah = 1.0;
  /* printf("Z: "); */
  for (i = 0; i < coeff_len; i ++) {
    /* printf("%g ", cf->Z[i][i]); */
    sigmah *= cf->Z[i][i];
  }
  /* printf("\n"); */

  if (sigmah < 0.0) {
    /* Not positive definite: probably shouldn't happen */
    rjmcmc_error("curvefit_sample: det less than zero\n");
    return -1;
  }

  /* printf("p_one_d_exp: %f\n", sigma2); */
  /* printf("tr(Z): %g (%g)\n", sigmah, 1.0/sigmah); */
  /* printf("cl: %d\n", coeff_len); */
  
  /* printf("s: %g %g\n", sigma2, sigmah); */

  *prob = 
    exp(-0.5 * sigma2)/
    (pow(2.0 * M_PI, (double)(coeff_len)/2.0) * (1.0/sigmah));
  return 0;
}

int
curvefit_sample_mean(curvefit_result_t *cf,
		     double *coeff,
		     int coeff_len)
{
  int i;

  for (i = 0; i < coeff_len; i ++) {
    coeff[i] = cf->mu[i];
  }

  return 0;
}

int
curvefit_sample_sigma(curvefit_result_t *cf,
		      double *sigma,
		      int sigma_len)
{
  int i;

  for (i = 0; i < sigma_len; i ++) {
    sigma[i] = sqrt(cf->S[i][i]);
  }

  return 0;
}

int
curvefit_sample_detCm(curvefit_result_t *cf,
		      double *detCm,
		      int order)
{
  int i;
  double sigmah;

  sigmah = 1.0;
  for (i = 0; i <= order; i ++) {
    sigmah *= cf->Z[i][i];
  }

  *detCm = (sigmah * sigmah);
  return 0;
}

static int
compute_hankel(const point1d_t *points,
	       int n,
	       double lambda,
	       int order,
	       curvefit_result_t *cf)
{
  int m;
  int i;
  int j;

  double ai;
  double bi;

  double l2;

  m = 2 * (order + 1);

  for (j = 0; j < m; j ++) {
    cf->alpha[j] = 0.0;
    cf->beta[j] = 0.0;
  }

  l2 = lambda*lambda;
  
  for (i = 0; i < n; i ++) {

    ai = 1.0/(l2 * points[i].n * points[i].n);
    bi = ai * points[i].y;

    cf->alpha[0] += ai;
    cf->beta[0] += bi;

    for (j = 1; j < m; j ++) {
      ai *= points[i].x;
      bi *= points[i].x;

      cf->alpha[j] += ai;
      cf->beta[j] += bi;
    }
  }

  return 0;
}

static void
set_I(double **m,
      int n)
{
  int i;
  int j;

  for (i = 0; i < n; i ++) {
    for (j = 0; j < n; j ++) {
      m[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
}

static int
compute_hankel_cholesky(const double *alpha,
			int order,
			double **L)
{
  int i;
  int j;
  int k;
  int m;

  m = order + 1;

  set_I(L, m);

  for (j = 0; j < m; j ++) {

    L[j][j] = alpha[2*j];
    for (k = 0; k < j; k ++) {
      L[j][j] -= L[j][k]*L[j][k];
    }
    if (L[j][j] <= 0.0) {
      /* rjmcmc_error("compute_hankel_cholesky: %g %d %g\n", L[j][j], j, */
      /* 	      alpha[2*j]); */
      return -1;
    }
    L[j][j] = sqrt(L[j][j]);
      
    for (i = j + 1; i < m; i ++) {
      L[i][j] = alpha[i + j];
      for (k = 0; k < j; k ++) {
	L[i][j] -= L[i][k]*L[j][k];
      }
      L[i][j] /= L[j][j];
    }
  }

  return 0;
}

static int
compute_forward_substitution(double **L,
			     int order,
			     double *b,
			     double *x)
{
  int m;
  int j;
  int i;

  m = order + 1;

  for (j = 0; j < m; j ++) {

    x[j] = b[j];

    for (i = 0; i < j; i ++) {
      x[j] -= L[j][i]*x[i];
    }

    x[j] /= L[j][j];

  }

  return 0;
}

static int
compute_backward_substitution(double **L,
			      int order,
			      double *b,
			      double *x)
{
  int m;
  int j;
  int i;

  m = order + 1;

  for (j = m - 1; j >= 0; j --) {

    x[j] = b[j];

    for (i = m - 1; i > j; i --) {
      x[j] -= L[i][j]*x[i];
    }

    x[j] /= L[j][j];

  }

  return 0;
}

static int
compute_inverse(double **L,
		int order,
		double **Z,
		double *b,
		double *x)
{
  int m;
  
  int i;
  int j;

  m = order + 1;

  for (j = 0; j < m; j ++) {
    for (i = 0; i < m; i ++) {
      b[i] = (j == i) ? 1.0 : 0;
    }

    compute_forward_substitution(L, order, b, x);

    for (i = 0; i < m; i ++) {
      Z[j][i] = x[i];
    }
  }
      
  return 0;
}

static int
compute_square(double **L,
	       int order,
	       double **L2)
{
  int m;
  int i;
  int j;
  int k;
  double t;


  m = order + 1;
  for (i = 0; i < m; i ++) {
    for (j = 0; j < m; j ++) {

      t = 0.0;

      for (k = 0; k < m; k ++) {
	t += L[i][k] * L[j][k];
      }

      L2[i][j] = t;
    }
  }

  return 0;
}

static int
compute_cholesky(double **A,
		 int order,
		 double **L)
{
  int m;
  int i;
  int j;
  int k;

  m = order + 1;

  set_I(L, m);

  for (j = 0; j < m; j ++) {

    L[j][j] = A[j][j];
    for (k = 0; k < j; k ++) {
      L[j][j] -= L[j][k]*L[j][k];
    }
    if (L[j][j] <= 0.0) {
      /* rjmcmc_error( */
      /* 	      "compute_cholesky: not positive definite %g\n", */
      /* 	      L[j][j]); */
      return -1;
    }
    L[j][j] = sqrt(L[j][j]);
      
    for (i = j + 1; i < m; i ++) {
      L[i][j] = A[i][j];
      for (k = 0; k < j; k ++) {
	L[i][j] -= L[i][k]*L[j][k];
      }
      L[i][j] /= L[j][j];
    }
  }

  return 0;
  
}

static double fx(double *a, int order, double x)
{
  double xp;
  double y;
  int i;

  xp = 1.0;
  y = 0.0;
  for (i = 0; i <= order; i ++) {
    y += xp * a[i];
    xp *= x;
  }

  return y;
}

int
curvefit_compute_mean_misfit(curvefit_result_t *cf,
			     const dataset1d_t *data,
			     int di,
			     int dj,
			     double lambda,
			     int order,
			     double *mean,
			     double *sigma,
			     double *mean_misfit,
			     double *detCm)
{
  int i;
  double dy;
  double n;
  double misfit;

  if (curvefit_compute_lambda(data,
			      lambda,
			      di,
			      dj,
			      order,
			      cf) < 0) {
    /* rjmcmc_error( */
    /* 	    "curvefit_compute_mean_misfit: " */
    /* 	    "failed to compute (%d %d, %d (%g))\n", */
    /* 	    di, */
    /* 	    dj, */
    /* 	    order, */
    /* 	    lambda); */
    return -1;
  }
  
  if (curvefit_sample_mean(cf,
			   mean,
			   order + 1) < 0) {
    return -1;
  }

  if (curvefit_sample_sigma(cf,
			    sigma,
			    order + 1) < 0) {
    return -1;
  }
  
  if (curvefit_sample_detCm(cf,
			    detCm,
			    order) < 0) {
    return -1;
  }
  
  /* Compute mean misfit */
  misfit = 0.0;
  
  for (i = di; i <= dj; i ++) {
    dy = fx(mean, order, data->points[i].x) - 
      data->points[i].y;
    n = data->points[i].n;
    misfit += dy*dy/(2.0 * n * n);
  }
  
  *mean_misfit = misfit;
  return 0;
}

int
curvefit_evaluate_pk(curvefit_result_t *cf,
		     const dataset1d_t *data,
		     int di, 
		     int dj,
		     int max_order,
		     const double *fixed_prior,
		     double auto_z,
		     double *mean_misfit,
		     double *detCm,
		     double *autoprior,
		     double **S,
		     double *pk,
		     double *kcdf,
		     double **mean,
		     double **sigma)
{
  int k;
  int j;

  for (k = 0; k <= max_order; k ++) {
    mean_misfit[k] = 0.0;
    detCm[k] = 0.0;

    if (curvefit_compute_mean_misfit(cf,
				     data,
				     di,
				     dj,
				     1.0,
				     k,
				     mean[k],
				     sigma[k],
				     mean_misfit + k,
				     detCm + k) < 0) {
      /* rjmcmc_error( */
      /* 	      "curvefit_evaluate_pk: " */
      /* 	      "failed to compute mean misfit (%d) f\n", */
      /* 	      k); */

      /* It is safe to continue here as we use detCm[k] == 0.0 as a test
	 to indicate if an order is unavailable */

      /* return -1; */
    }
  }
  
  if (fixed_prior != NULL) {
    /*
     * If using a fixed prior, we overwrite the auto prior before 
     * calculating the probabilities.
     */
    for (k = 0; k <= max_order; k ++) {
      if (k == 0) {
	autoprior[k] = fixed_prior[k];
      } else {
	autoprior[k] = autoprior[k - 1] * fixed_prior[k];
      }
    }
  } else {
    /*
     * Compute the auto prior from the std dev of the coefficients.
     */
    for (k = 0; k <= max_order; k ++) {
      autoprior[k] = 1.0;
      for (j = 0; j <= k; j ++) {
	autoprior[k] *= (sigma[k][j] * 2.0 * auto_z);
      }
    }

    /* rjmcmc_debug("az = %g\n", auto_z); */
    /* for (k = 0; k <= max_order; k ++) { */
    /*   rjmcmc_debug("%g\n", autoprior[k]); */
    /* } */
    /* for (k = 0; k <= max_order; k ++) { */
    /*   for (j = 0; j <= k; j ++) { */
    /* 	rjmcmc_debug("%g ", sigma[k][j]); */
    /*   } */
    /*   rjmcmc_debug("\n"); */
    /* } */
  }

  for (k = 0; k <= max_order; k ++) {
    
    S[k][k] = 1.0;
    
    if (detCm[k] > 0.0) {
      
      for (j = k + 1; j <= max_order; j ++) {
	
	if (detCm[j] > 0.0) {
	  S[k][j] = exp(mean_misfit[j] - mean_misfit[k]) *
	    sqrt(pow(2.0 * M_PI, k - j) * detCm[k]/detCm[j]);
	  
	  S[k][j] *= autoprior[j]/autoprior[k];
	  
	  S[j][k] = 1.0/S[k][j];
	  
	} else {
	  S[k][j] = 0.0;
	  S[j][k] = 0.0;
	}
	
      }
    } else {
      
      for (j = 0; j <= max_order; j ++) {
	S[k][j] = 0.0;
	S[j][k] = 0.0;
      }
    }
  }

  for (k = 0; k <= max_order; k ++) {
    if (k == 0) {
      kcdf[k] = 0.0;
    } else {
      kcdf[k] = kcdf[k - 1];
    }

    pk[k] = 0.0;
    for (j = 0; j <= max_order; j ++) {
      pk[k] += S[j][k];
    }
    
    if (pk[k] > 0.0) {
      pk[k] = 1.0/pk[k];
    }

    kcdf[k] += pk[k];
  }
    
  return 0;
}

double **cf_L(curvefit_result_t *cf)
{
  return cf->L;
}

double **cf_Z(curvefit_result_t *cf)
{
  return cf->Z;
}

double **cf_S(curvefit_result_t *cf)
{
  return cf->S;
}

double **cf_Si(curvefit_result_t *cf)
{
  return cf->Si;
}

double *cf_mu(curvefit_result_t *cf)
{
  return cf->mu;
}

