
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/forwardmodel_util.h>

int approx_equal(double a, double b)
{
  static const double epsilon = 1e-5;

  if (fabs(a - b) > epsilon) {
    printf("approx_equal: %f != %f (e = %g)\n", a, b, fabs(a - b));
    return 0;
  } 

  return -1;
}

int within_1percent(double expected, double value)
{
  double p = fabs(expected - value)/expected;

  return p < 0.01;
}

int init_util_suite(void)
{
  return 0;
}

int clean_util_suite(void)
{
  return 0;
}

void test_util_polynomial(void);
void test_util_sigmapower(void);
void test_util_gaussian(void);
void test_fill_coord_vector(void);

#define ADD_TEST(name, function) \
  if (CU_add_test(pSuite, name, function) == NULL) {\
  CU_cleanup_registry(); \
  return CU_get_error(); \
  }

int main(int argc, char *argv[])
{
  CU_pSuite pSuite = NULL;

  if (CU_initialize_registry() != CUE_SUCCESS) {
    return CU_get_error();
  }

  pSuite = CU_add_suite("Util Suite", init_util_suite, clean_util_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Polynomial", test_util_polynomial);
  ADD_TEST("Forwardmodel Sigma Power Error", test_util_sigmapower);
  ADD_TEST("Gaussian Probability", test_util_gaussian);
  ADD_TEST("Fill Coord Vector", test_fill_coord_vector);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

double one(void)
{
  return 1.0;
}

void test_util_polynomial(void)
{
  double poly[2] = {5.0, 1.0};

  CU_ASSERT(approx_equal(rjmcmc_polynomial_value(poly, 2, 0.0), 5.0));
  CU_ASSERT(approx_equal(rjmcmc_polynomial_value(poly, 2, 1.0), 6.0));
  CU_ASSERT(approx_equal(rjmcmc_polynomial_value(poly, 2, -1.0), 4.0));
}

void test_util_sigmapower(void)
{
  int n = 10;
  double sigma = 2.0;
  double r = 0.5;
  int i;

  double *phi;
  double error;

  double phi_ex[20] = {
    31,26,10,13,9.2,8.5,8.2,42,7.3,6.6,6.2,5.6,2.9,4,4.6,2.6,2.6,1.7,2.2,2.2
  };

  phi = rjmcmc_create_array_1d(n);
  CU_ASSERT_FATAL(phi != NULL);

  for (i = 0; i < n; i ++) {
    phi[i] = 1.0;
  }

  error = forwardmodel_misfit_sigma_power(phi,
					  n,
					  sigma,
					  r);

  CU_ASSERT(error == 1.0);

  for (i = 0; i < n; i ++) {
    phi[i] = (double)(i + 1);
  }

  error = forwardmodel_misfit_sigma_power(phi,
					  n,
					  sigma,
					  r);

  CU_ASSERT(error == 42.0);
  
  rjmcmc_destroy_array_1d(phi);
  
  n = 20;

  error = forwardmodel_misfit_sigma_power(phi_ex,
					  n,
					  1.003220,
					  0.796472);

  CU_ASSERT(within_1percent(6661.9, error));
  

}

void test_util_gaussian(void)
{
  CU_ASSERT(within_1percent(exp(-0.5)/sqrt(2.0 * M_PI),
			    rjmcmc_gaussian_probability(1.0, 1.0)));

  CU_ASSERT(within_1percent(log(rjmcmc_gaussian_probability(1.0, 1.0)),
			    rjmcmc_log_gaussian_probability(1.0, 1.0)));
}

void test_fill_coord_vector(void)
{
  #define COORDS_SIZE 20
  double coords[COORDS_SIZE];
  double xmin;
  double xmax;
  double halfstep;
  int i;

  xmin = -3.0;
  xmax = 7.0;
  
  rjmcmc_fill_coord_vector(xmin, xmax, COORDS_SIZE, coords);

  /* printf("\n"); */
  /* for (i = 0; i < COORDS_SIZE; i ++) { */
  /*   printf("%2d %10.5f\n", i, coords[i]); */
  /* } */
  /* printf("\n"); */

  halfstep = 0.5*(xmax - xmin)/(double)(COORDS_SIZE);

  CU_ASSERT(approx_equal(xmin + halfstep, coords[0]));
  CU_ASSERT(approx_equal(xmax - halfstep, coords[COORDS_SIZE - 1]));
}
