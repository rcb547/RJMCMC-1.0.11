
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/dataset1d.h>
#include <rjmcmc/part1d_regression_rj.h>

int approx_equal(double a, double b)
{
  static const double epsilon = 1e-5;

  if (fabs(a - b) > epsilon) {
    printf("approx_equal: %f != %f (e = %g)\n", a, b, fabs(a - b));
    return 0;
  } 

  return -1;
}

static dataset1d_t *data = NULL;

int init_part1d_regression_suite(void)
{
  int i;

  #define DATASIZE 20
  double datapoints[DATASIZE][3] = {
    {0.00000000e+00, -7.93237964e-01, 5.00000000e-02},
    {5.26315789e-02, -7.41269756e-01, 5.00000000e-02},
    {1.05263158e-01, -6.87997715e-01, 5.00000000e-02},
    {1.57894737e-01, -6.32798220e-01, 5.00000000e-02},
    {2.10526316e-01, -5.75047653e-01, 5.00000000e-02},
    {2.63157895e-01, -5.14122396e-01, 5.00000000e-02},
    {3.15789474e-01, -4.49398828e-01, 5.00000000e-02},
    {3.68421053e-01, -3.80253332e-01, 5.00000000e-02},
    {4.21052632e-01, -3.06062289e-01, 5.00000000e-02},
    {4.73684211e-01, -2.26202078e-01, 5.00000000e-02},
    {5.26315789e-01, -1.40049083e-01, 5.00000000e-02},
    {5.78947368e-01, -4.69796827e-02, 5.00000000e-02},
    {6.31578947e-01,  5.36297405e-02, 5.00000000e-02},
    {6.84210526e-01,  1.62402806e-01, 5.00000000e-02},
    {7.36842105e-01,  2.79963132e-01, 5.00000000e-02},
    {7.89473684e-01,  4.06934339e-01, 5.00000000e-02},
    {8.42105263e-01,  5.43940044e-01, 5.00000000e-02},
    {8.94736842e-01,  6.91603867e-01, 5.00000000e-02},
    {9.47368421e-01,  8.50549427e-01, 5.00000000e-02},
    {1.00000000e+00,  1.02140034e+00, 5.00000000e-02}
  };

  data = dataset1d_create(DATASIZE);
  if (data == NULL) {
    return -1;
  }

  data->xmin = 0.0;
  data->xmax = 1.0;

  data->ymin = -2.0;
  data->ymax = 2.0;

  for (i = 0; i < DATASIZE; i ++) {
    data->points[i].x = datapoints[i][0];
    data->points[i].y = datapoints[i][1];
    data->points[i].n = datapoints[i][2];
  }

  return 0;
}

int clean_part1d_regression_suite(void)
{
  dataset1d_destroy(data);

  return 0;
}

void test_part1d_regression_initialize(void);

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

  pSuite = CU_add_suite("Partition 1D Regression Suite", init_part1d_regression_suite, clean_part1d_regression_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Initialize", test_part1d_regression_initialize);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

double one(void)
{
  return 1.0;
}

void test_part1d_regression_initialize(void)
{
  part1d_regression_rj_t *p;

  p = part1d_regression_rj_create(2, 5, 5, 1, 0.0, 1.0, 0.1);
  CU_ASSERT_FATAL(p != NULL);

  CU_ASSERT_FATAL(part1d_regression_rj_initialize(p,
						  (const dataset1d_t**)&data,
						  1,
						  rjmcmc_uniform,
						  rjmcmc_normal) == 0);
}

