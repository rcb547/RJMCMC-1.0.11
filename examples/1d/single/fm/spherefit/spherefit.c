
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <rjmcmc/forwardmodel.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/resultsetfm.h>

#include <rjmcmc/rjmcmc_util.h>

static double unif(void)
{
  return rjmcmc_uniform();
}

static double norm(void)
{
  return rjmcmc_normal();
}

struct spherefit_data {
  int npoints;
  double *x;
  double *y;
  double *z;

  double fixed_sigma;
};

typedef enum {
  P_XC = 0,
  P_YC,
  P_ZC,
  P_R,
  P_TOTAL
} parameter_t;

static double spherefit_misfit(void *user_arg,
			       int nvalues,
			       const double *values);

static int write_points(const double *x,
			const double *y,
			const double *z,
			int n,
			const char *filename);

int main(int argc, char *argv[]) 
{
  struct spherefit_data data;

  int burnin = 1000;
  int total = 20000;
  int samples = 100;
  double confidence_interval = 0.95;
  int requested_results = RESULTSETFM_MEAN;

  forwardmodelparameter_t parameters[P_TOTAL];

  int nproc;
  const double *v;
  const int *iv;
  int i;

  double *xcoords;
  double *y;
  int xcl;

  resultsetfm_t *results;

  double real_xc;
  double real_yc;
  double real_zc;
  double real_r;
  double real_rsigma;

  double theta;
  double phi;
  double gamma;

  int npoints;

  double r;


  /*
   * Create the synthetic data
   */
  real_xc = -5.0;
  real_yc = 2.0;
  real_zc = 1.0;
  real_r = 7.5;
  real_rsigma = 1.0;

  npoints = 100;

  data.x = rjmcmc_create_array_1d(npoints);
  data.y = rjmcmc_create_array_1d(npoints);
  data.z = rjmcmc_create_array_1d(npoints);
  if (data.x == NULL ||
      data.y == NULL || 
      data.z == NULL) {
    fprintf(stderr, "error: failed to allocate memory for data\n");
    return -1;
  }

  for (i = 0; i < npoints; i ++) {

    theta = rjmcmc_random_choose_double(-M_PI, M_PI, unif);
    phi = rjmcmc_random_choose_double(-M_PI, M_PI, unif);
    gamma = rjmcmc_random_choose_double(-M_PI, M_PI, unif);

    r = real_r + real_rsigma * norm();

    data.x[i] = real_xc + r * (cos(gamma)*cos(theta) - sin(gamma)*sin(phi)*sin(theta));
    data.y[i] = real_yc - r * sin(gamma) * cos(phi);
    data.z[i] = real_zc + r*(cos(gamma)*sin(theta) + sin(gamma)*sin(phi)*cos(theta));

  }

  data.npoints = npoints;
  data.fixed_sigma = real_rsigma;

  if (write_points(data.x, 
		   data.y,
		   data.z,
		   data.npoints,
		   "points.txt") < 0) {
    fprintf(stderr, "error: failed to save points data\n");
    return -1;
  }
  
  /*
   * Initialize the search space for the parameters
   */
  
  parameters[P_XC].fmin = -10.0;
  parameters[P_XC].fmax = 10.0;
  parameters[P_XC].fstd_value = 1.0;
  parameters[P_XC].fstd_bd = 0.0;

  parameters[P_YC].fmin = -10.0;
  parameters[P_YC].fmax = 10.0;
  parameters[P_YC].fstd_value = 1.0;
  parameters[P_YC].fstd_bd = 0.0;

  parameters[P_ZC].fmin = -10.0;
  parameters[P_ZC].fmax = 10.0;
  parameters[P_ZC].fstd_value = 1.0;
  parameters[P_ZC].fstd_bd = 0.0;

  parameters[P_R].fmin = 1.0;
  parameters[P_R].fmax = 20.0;
  parameters[P_R].fstd_value = 1.0;
  parameters[P_R].fstd_bd = 0.0;


  /*
   * Run the forward model
   */

  results = single_forwardmodel(burnin,
				total,
				unif, /*rjmcmc_uniform,*/
				norm, /*rjmcmc_normal,*/
				P_TOTAL,
				parameters,
				spherefit_misfit,
				&data,
				samples,
				confidence_interval,
				requested_results);

  if (results == NULL) {
    fprintf(stderr, 
	    "error: failed to run functionfit\n");
    return -1;
  }

  v = resultsetfm_get_misfit(results);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get misfit data\n");
    return -1;
  }
  if (rjmcmc_save_vector("spherefit.misfit", v, total) < 0) {
    fprintf(stderr, "error: failed to save misfit data\n");
    return -1;
  }

  iv = resultsetfm_get_propose(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get propose counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  iv = resultsetfm_get_accept(results, &nproc);
  if (iv == NULL) {
    fprintf(stderr, "error: failed to get accept counts\n");
    return -1;
  }
  for (i = 0; i < nproc; i ++) {
    printf("%6d ", iv[i]);
  }
  printf("\n");

  printf("x centre: %f %f\n", 
	 real_xc, 
	 resultsetfm_get_parameter_mean(results, P_XC));
  printf("y centre: %f %f\n", 
	 real_yc, 
	 resultsetfm_get_parameter_mean(results, P_YC));
  printf("z centre: %f %f\n", 
	 real_zc, 
	 resultsetfm_get_parameter_mean(results, P_ZC));
  printf("radius  : %f %f\n", 
	 real_r, 
	 resultsetfm_get_parameter_mean(results, P_R));

  /*
   * X coordinate
   */
  v = resultsetfm_get_parameter_history(results, P_XC);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get parameter x\n");
    return -1;
  }
  if (rjmcmc_save_vector("spherefit.x", v, total) < 0) {
    fprintf(stderr, "error: failed to save x paramter\n");
    return -1;
  }

  /*
   * Y coordinate
   */
  v = resultsetfm_get_parameter_history(results, P_YC);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get parameter y\n");
    return -1;
  }
  if (rjmcmc_save_vector("spherefit.y", v, total) < 0) {
    fprintf(stderr, "error: failed to save y paramter\n");
    return -1;
  }

  /*
   * Z coordinate
   */
  v = resultsetfm_get_parameter_history(results, P_ZC);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get parameter z\n");
    return -1;
  }
  if (rjmcmc_save_vector("spherefit.z", v, total) < 0) {
    fprintf(stderr, "error: failed to save z paramter\n");
    return -1;
  }

  /*
   * Radius
   */
  v = resultsetfm_get_parameter_history(results, P_R);
  if (v == NULL) {
    fprintf(stderr, "error: failed to get parameter radius\n");
    return -1;
  }
  if (rjmcmc_save_vector("spherefit.r", v, total) < 0) {
    fprintf(stderr, "error: failed to save radius paramter\n");
    return -1;
  }

  resultsetfm_destroy(results);

  return 0;
}

static double sqr(double x)
{
  return x*x;
}

static double spherefit_misfit(void *user_arg,
			       int nvalues,
			       const double *values)
{
  struct spherefit_data *data = (struct spherefit_data *)user_arg;

  double xc = values[P_XC];
  double yc = values[P_YC];
  double zc = values[P_ZC];
  double r = values[P_R];

  double dr;

  int i;

  double sum;

  sum = 0.0;

  for (i = 0; i < data->npoints; i ++) {

    dr = r - sqrt(sqr(data->x[i] - xc) + 
		  sqr(data->y[i] - yc) + 
		  sqr(data->z[i] - zc));

    sum += sqr(dr);
  }

  return sum/(sqr(data->fixed_sigma) * 2.0);
}

static int write_points(const double *x,
			const double *y,
			const double *z,
			int n,
			const char *filename)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    return -1;
  }

  for (i = 0; i < n; i ++) {
    fprintf(fp, "%f %f %f\n", x[i], y[i], z[i]);
  }

  fclose(fp);
  return 0;
}

