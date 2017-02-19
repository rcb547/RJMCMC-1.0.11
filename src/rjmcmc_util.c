
#include "rjmcmc/rjmcmc_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <rjmcmc/rjmcmc_debug.h>

#if defined(HAVE_STDINT_H)
#include <stdint.h>
#endif

#include "rjmcmc/rjmcmc_util.h"

double *rjmcmc_create_array_1d(int w)
{
  double *r;
  if (w == 0) {
    return NULL;
  }

  r = (double*)malloc(sizeof(double) * w);
  memset(r, 0, sizeof(double) * w);

  return r;
}

double **rjmcmc_create_array_2d(int w, int h)
{
  double **r;
  int i;

  if (w == 0 || h == 0) {
    return NULL;
  }

  r = (double**)malloc(sizeof(double*) * w);
  for (i = 0; i < w; i ++) {
    r[i] = rjmcmc_create_array_1d(h);
    if (r[i] == NULL) {
      return NULL;
    }
  }

  return r;
}

double ***rjmcmc_create_array_3d(int w, int h, int d)
{
  double ***r;
  int i;

  if (w == 0 || h == 0 || d == 0) {
    return NULL;
  }

  r = (double***)malloc(sizeof(double**) * w);
  for (i = 0; i < w; i ++) {
    r[i] = rjmcmc_create_array_2d(h, d);
    if (r[i] == NULL) {
      return NULL;
    }
  }

  return r;
}

void rjmcmc_destroy_array_1d(double *a)
{
  free(a);
}

void rjmcmc_destroy_array_2d(int w, double **a)
{
  int i;
  if (a != NULL) {
    for (i = 0; i < w; i ++) {
      rjmcmc_destroy_array_1d(a[i]);
    }
    free(a);
  }
}

void rjmcmc_destroy_array_3d(int w, int h, double ***a)
{
  int i;
  if (a != NULL) {
    for (i = 0; i < w; i ++) {
      rjmcmc_destroy_array_2d(h, a[i]);
    }
    free(a);
  }
}


int *rjmcmc_create_int_array_1d(int w)
{
  int *r;

  if (w == 0) {
    return NULL;
  }

  r = (int*)malloc(sizeof(int) * w);
  memset(r, 0, sizeof(int) * w);

  return r;
}

int **rjmcmc_create_int_array_2d(int w, int h)
{
  int **r;
  int i;

  if (w == 0 || h == 0) {
    return NULL;
  }

  r = (int**)malloc(sizeof(int*) * w);
  for (i = 0; i < w; i ++) {
    r[i] = rjmcmc_create_int_array_1d(h);
    if (r[i] == NULL) {
      return NULL;
    }
  }

  return r;
}

int ***rjmcmc_create_int_array_3d(int w, int h, int d)
{
  int ***r;
  int i;

  if (w == 0 || h == 0 || d == 0) {
    return NULL;
  }

  r = (int***)malloc(sizeof(int**) * w);
  for (i = 0; i < w; i ++) {
    r[i] = rjmcmc_create_int_array_2d(h, d);
    if (r[i] == NULL) {
      return NULL;
    }
  }

  return r;
}

int ****rjmcmc_create_int_array_4d(int w, int h, int d, int e)
{
  int ****r;
  int i;

  if (w == 0 || h == 0 || d == 0 || e == 0) {
    return NULL;
  }

  r = (int****)malloc(sizeof(int**) * w);
  for (i = 0; i < w; i ++) {
    r[i] = rjmcmc_create_int_array_3d(h, d, e);
    if (r[i] == NULL) {
      return NULL;
    }
  }

  return r;
}

void rjmcmc_destroy_int_array_1d(int *a)
{
  free(a);
}

void rjmcmc_destroy_int_array_2d(int w, int **a)
{
  int i;
  if (a != NULL) {
    for (i = 0; i < w; i ++) {
      rjmcmc_destroy_int_array_1d(a[i]);
    }
    free(a);
  }
}

void rjmcmc_destroy_int_array_3d(int w, int h, int ***a)
{
  int i;
  if (a != NULL) {
    for (i = 0; i < w; i ++) {
      rjmcmc_destroy_int_array_2d(h, a[i]);
    }
    free(a);
  }
}

void rjmcmc_destroy_int_array_4d(int w, int h, int d, int ****a)
{
  int i;
  if (a != NULL) {
    for (i = 0; i < w; i ++) {
      rjmcmc_destroy_int_array_3d(h, d, a[i]);
    }
    free(a);
  }
}

int rjmcmc_map_to_index(double v, double vmin, double vmax, int n)
{
  int i;
  
  i = (int)((double)n * (v - vmin)/(vmax - vmin));

  if (i < 0) {
    return 0;
  }

  if (i > (n - 1)) {
    return n - 1;
  }

  return i;
}

double rjmcmc_mean(const double *v, int n)
{
  return rjmcmc_mean_skip(v, 0, n);
}

double rjmcmc_mean_skip(const double *v, int skip, int n)
{
  int i;
  double s;

  s = 0.0;
  for (i = skip; i < n; i ++) {
    s += v[i];
  }

  return s/(double)n;
}

int rjmcmc_mean_variance(const double *v, 
			 int n, 
			 double *_mean, 
			 double *_variance)
{
  int t;
  double mean;
  double m;
  double delta;
  int i;

  if (n < 2) {
    return -1;
  }

  t = 0;
  mean = 0.0;
  m = 0.0;

  for (i = 0; i < n; i ++) {
    
    t ++;
    delta = v[i] - mean;
    mean += delta/(double)t;
    m += delta * (v[i] - mean);

  }

  *_mean = mean;
  *_variance = m/(double)(t - 1);

  return 0;
}

int
rjmcmc_vector_to_histogram(int s, int n, const double *v,
			   int hn, double vmin, double vmax, int *hist)
{
  int i;

  for (i = 0; i < hn; i ++) {
    hist[i] = 0;
  }

  for (i = s; i < n; i ++) {
    hist[rjmcmc_map_to_index(v[i], vmin, vmax, hn)] ++;
  }

  return 0;
}


double rjmcmc_median_from_histogram(int *hist, double vmin, double vmax, int n)
{
  int i;
  int j;
  int ci;
  int cj;

  i = 0;
  j = n - 1;
  ci = 0;
  cj = 0;

  while (i != j) {
    if (ci < cj) {
      ci += hist[i];
      i ++;
    } else {
      cj += hist[j];
      j --;
    }
  }

  return (double)i/(double)n * (vmax - vmin) + vmin;
}

double rjmcmc_head_from_histogram(int *hist, double vmin, double vmax, int n, int drop)
{
  int i;
  int ci;

  i = 0; 
  ci = 0;
  while(i < n && ci < drop) {
    if (hist[i] + ci >= drop) {
      break;
    }

    ci += hist[i];
    i ++;
  }

  return (double)i/(double)n * (vmax - vmin) + vmin;
}

double rjmcmc_tail_from_histogram(int *hist, double vmin, double vmax, int n, int drop)
{
  int i;
  int ci;

  i = n - 1; 
  ci = 0;
  while(i > 0 && ci < drop) {
    if (hist[i] + ci >= drop) {
      break;
    }

    ci += hist[i];
    i --;
  }

  return (double)i/(double)n * (vmax - vmin) + vmin;
}

double rjmcmc_mode_from_histogram(int *hist, double vmin, double vmax, int n)
{
  int i;
  int m;
  int mi;

  m = 0;
  mi = -1;

  for (i = 0; i < n; i ++) {
    if (hist[i] > m) {
      m = hist[i];
      mi = i;
    }
  }
  
  if (mi < 0) {
    return 0.0;
  }

  return (double)mi/(double)n * (vmax - vmin) + vmin;
}

int rjmcmc_save_vector(const char *filename,
		       const double *v,
		       int n)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    rjmcmc_error("rjmcmc_save_vector: failed to open file for writing\n");
    return -1;
  }

  for (i = 0; i < n; i ++) {
    fprintf(fp, "%f\n", v[i]);
  }

  fclose(fp);
  return 0;
}


int rjmcmc_save_coords(const char *filename,
		       const double *x,
		       const double *y,
		       int n)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    rjmcmc_error("rjmcmc_save_coords: failed to open file for writing\n");
    return -1;
  }

  for (i = 0; i < n; i ++) {
    fprintf(fp, "%f %f\n", x[i], y[i]);
  }

  fclose(fp);
  return 0;
}

int rjmcmc_save_vector_as_histogram(const char *filename,
				    double minv,
				    double maxv,
				    int bins,
				    const double *v,
				    int n)
{
  double *x;
  double *hist;
  int i;
  int j;
  int r;

  x = rjmcmc_create_array_1d(bins);
  if (x == NULL) {
    return -1;
  }

  hist = rjmcmc_create_array_1d(bins);
  if (hist == NULL) {
    return -1;
  }

  for (i = 0; i < n; i ++) {
    j = (int)((v[i] - minv)/(maxv - minv) * (double)bins);
    if (j >= 0 && j < bins) {
      hist[j] ++;
    }
  }

  for (i = 0; i < bins; i ++) {
    x[i] = minv + (maxv - minv) * ((double)i + 0.5)/(double)bins;
  }

  r = rjmcmc_save_coords(filename, x, hist, bins);
  rjmcmc_destroy_array_1d(x);
  rjmcmc_destroy_array_1d(hist);
  
  return r;
}

int rjmcmc_save_int_vector(const char *filename,
			   const int *v,
			   int n)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    rjmcmc_error("rjmcmc_save_int_vector: failed to open file for writing\n");
    return -1;
  }

  for (i = 0; i < n; i ++) {
    fprintf(fp, "%d\n", v[i]);
  }

  fclose(fp);
  return 0;
}

int rjmcmc_save_int_vector_as_histogram(const char *filename,
					int minv,
					int maxv,
					const int *v,
					int n)
{
  int *hist;
  int i;
  int j;
  int l;

  int r;

  l = maxv - minv + 1;
  hist = rjmcmc_create_int_array_1d(l);
  if (hist == NULL) {
    return -1;
  }

  for (i = 0; i < n; i ++) {
    j = v[i] - minv;
    if (j >= 0 && j < l) {
      hist[j] ++;
    }
  }

  r = rjmcmc_save_int_vector(filename, hist, l);
  rjmcmc_destroy_int_array_1d(hist);
  
  return r;
}

int rjmcmc_save_int_coords(const char *filename,
			   const double *x,
			   const int *iy,
			   int n)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    rjmcmc_error("rjmcmc_save_coords: failed to open file for writing\n");
    return -1;
  }

  for (i = 0; i < n; i ++) {
    fprintf(fp, "%f %d\n", x[i], iy[i]);
  }

  fclose(fp);
  return 0;
}

int rjmcmc_save_matrix(const char *filename,
		       const double **m,
		       int c,
		       int r)
{
  int i;
  int j;
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    rjmcmc_error("rjmcmc_save_matrix: failed to open file for writing\n");
    return -1;
  }

  for (j = 0; j < r; j ++) {
    for (i = 0; i < c; i ++) {
      fprintf(fp, "%f ", m[i][j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}

int rjmcmc_save_int_matrix(const char *filename,
			   const int **m,
			   int c,
			   int r)
{
  int i;
  int j;
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    rjmcmc_error("rjmcmc_save_int_matrix: failed to open file for writing\n");
    return -1;
  }

  for (j = 0; j < r; j ++) {
    for (i = 0; i < c; i ++) {
      fprintf(fp, "%d ", m[i][j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}

#define SQRT2PI 2.5066282746310002
double rjmcmc_gaussian_probability(double phi,
				   double sigma)
{
  return exp(-(phi*phi)/(2.0 * sigma*sigma))/(SQRT2PI * sigma);
}

double rjmcmc_log_gaussian_probability(double phi,
				       double sigma)
{
  return -log(sigma * SQRT2PI) - (phi*phi)/(2.0 * sigma*sigma);
}

double rjmcmc_gaussian_interval_probability(double phi,
					    double sigma,
					    double deltasigma)
{
  double d = deltasigma * sigma;

#if 0
printf("erf: %g %g %g\n",
	 (0.5 * (1.0 + erf((phi + d)/(sqrt(2.0) * sigma)))),

	 (0.5 * (1.0 + erf((phi - d)/(sqrt(2.0) * sigma)))),
	 (0.5 * (1.0 + erf((phi + d)/(sqrt(2.0) * sigma)))) -
	 (0.5 * (1.0 + erf((phi - d)/(sqrt(2.0) * sigma)))));

  return 
    (0.5 * (1.0 + erf((phi + d)/(sqrt(2.0) * sigma)))) -
    (0.5 * (1.0 + erf((phi - d)/(sqrt(2.0) * sigma))));
#endif
  return 0.0;
}

double rjmcmc_polynomial_value(const double *coeff,
			       int ncoeff,
			       double x)
{
  int i;
  double xs = coeff[0];
  double xp = x;

  for (i = 1; i < ncoeff; i ++) {
    xs += coeff[i] * xp;
    xp *= x;
  }

  return xs;
}

int rjmcmc_factorial(int n)
{
  int t = 1;

  while (n > 0) {
    t *= n;
    n --;
  }

  return t;
}

int rjmcmc_binomial(int _n, int _k)
{
  uint64_t r = 1, d = _n - _k;
  uint64_t n = _n;
  uint64_t k = _k;

  if (d > k) { 
    k = d; d = n - k; 
  }
 
  while (n > k) {
    if (r >= (UINT64_MAX/n)) {
      return INT_MAX;
    }
    r *= n--;
 
    while (d > 1 && !(r % d)) {
      r /= d--;
    }
  }

  return (int)r;
}

int
rjmcmc_fill_coord_vector(double xmin,
			 double xmax,
			 int xsamples,
			 double *x)
{
  int i;

  for (i = 0; i < xsamples; i ++) {
    x[i] = xmin + 
      (xmax - xmin) * ((double)i + 0.5)/(double)(xsamples);
  }

  return xsamples;
}
