
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/dataset1d.h>

int approx_equal(double a, double b)
{
  static const double epsilon = 1e-5;

  if (fabs(a - b) > epsilon) {
    printf("approx_equal: %f != %f (e = %g)\n", a, b, fabs(a - b));
    return 0;
  } 

  return -1;
}

int init_dataset1d_suite(void)
{
  return 0;
}

int clean_dataset1d_suite(void)
{
  return 0;
}

void test_dataset1d_mean_variance(void);

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

  pSuite = CU_add_suite("Dataset1d Suite", init_dataset1d_suite, clean_dataset1d_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Mean Variance", test_dataset1d_mean_variance);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

double one(void)
{
  return 1.0;
}

void test_dataset1d_mean_variance(void)
{
  double x[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 
		  6.0, 7.0, 8.0, 9.0, 10.0};
  double y[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 
		  6.0, 7.0, 8.0, 9.0, 10.0};
  double n[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 
		  1.0, 1.0, 1.0, 1.0, 1.0};
  int npoints = 10;


  dataset1d_t *p;

  double mean;
  double variance;
  

  p = dataset1d_create_from_array(x, y, n, npoints);
  CU_ASSERT_FATAL(p != NULL);

  CU_ASSERT_FATAL(dataset1d_mean_variance(p,
					  0,
					  npoints - 1,
					  &mean,
					  &variance) >= 0);

  /*
   * Values obtained from R:
   *> mean(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
   *[1] 5.5
   *> var(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
   *[1] 9.166667
   *
   */

  CU_ASSERT(approx_equal(mean, 5.5));
  CU_ASSERT(approx_equal(variance, 9.166667));

  dataset1d_destroy(p);
}


