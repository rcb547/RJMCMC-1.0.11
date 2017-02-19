
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/rjmcmc_config.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

int init_random_suite(void)
{
  return 0;
}

int clean_random_suite(void)
{
  return 0;
}

void test_random_uniform(void);
void test_random_normal(void);
void test_random_choose_int(void);

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

  pSuite = CU_add_suite("Random Suite", 
			init_random_suite, 
			clean_random_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Uniform", test_random_uniform);
  ADD_TEST("Normal", test_random_normal);
  ADD_TEST("Choose Int", test_random_choose_int);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

void test_random_uniform(void)
{
  int n;
  double *v;
  int i;

  double mean;
  double variance;

  int histsize;
  int *hist;

  double histmin;
  double histmax;

  FILE *fp;

  n = 10000000;
  v = rjmcmc_create_array_1d(n);
  CU_ASSERT_FATAL(v != NULL);

  for (i = 0; i < n; i ++) {
    //    v[i] = rjmcmc_uniform();
    v[i] = rjmcmc_uniform();
    CU_ASSERT(v[i] >= 0.0 && v[i] <= 1.0);
  }

  /* 
   * Write out a histogram for visual inspection
   */
  histsize = 100;
  hist = rjmcmc_create_int_array_1d(histsize);
  CU_ASSERT_FATAL(hist != NULL);

  histmin = 0.0;
  histmax = 1.0;
  for (i = 0; i < n; i ++) {
    hist[rjmcmc_map_to_index(v[i], histmin, histmax, histsize)] ++;
  }

  fp = fopen("uniform_histogram.txt", "w");
  CU_ASSERT_FATAL(fp != NULL);
  
  for (i = 0; i < histsize; i ++) {

    fprintf(fp, "%f %d\n", 
	    ((histmax - histmin) * ((double)i + 0.5))/(double)histsize + histmin, 
	    hist[i]);

  }

  fclose(fp);

  fp = fopen("uniform_data.txt", "w");
  CU_ASSERT_FATAL(fp != NULL);
  
  for (i = 0; i < n; i ++) {
    fprintf(fp, "%g\n", v[i]);
  }

  fclose(fp);

}

void test_random_normal(void)
{
  int n;
  double *v;
  int i;

  double mean;
  double variance;

  int histsize;
  int *hist;

  double histmin;
  double histmax;

  FILE *fp;

  n = 10000000;
  v = rjmcmc_create_array_1d(n);
  CU_ASSERT_FATAL(v != NULL);

  for (i = 0; i < n; i ++) {
    v[i] = rjmcmc_normal();
  }

  CU_ASSERT_FATAL(rjmcmc_mean_variance(v, n, &mean, &variance) >= 0);

  CU_ASSERT(fabs(mean) < 0.01);
  CU_ASSERT(fabs(sqrt(variance) - 1.0) < 0.01);


  /* printf("%f %f\n", mean, variance); */

  /* 
   * Write out a histogram for visual inspection
   */
  histsize = 1024;
  hist = rjmcmc_create_int_array_1d(histsize);
  CU_ASSERT_FATAL(hist != NULL);

  histmin = -6.0;
  histmax = 6.0;
  for (i = 0; i < n; i ++) {
    hist[rjmcmc_map_to_index(v[i], histmin, histmax, histsize)] ++;
  }

  fp = fopen("normal_histogram.txt", "w");
  CU_ASSERT_FATAL(fp != NULL);
  
  for (i = 0; i < histsize; i ++) {

    fprintf(fp, "%f %d\n", 
	    ((histmax - histmin) * ((double)i + 0.5))/(double)histsize + histmin, 
	    hist[i]);

  }

  fclose(fp);

  fp = fopen("normal_data.txt", "w");
  CU_ASSERT_FATAL(fp != NULL);
  
  for (i = 0; i < n; i ++) {
    fprintf(fp, "%g\n", v[i]);
  }

  fclose(fp);

}

void test_random_choose_int(void)
{
  int vmin;
  int vmax;
  int c;
  int cmin;
  int cmax;
  int i;
  int t;

  vmin = 13;
  vmax = 27;
  c = 0;
  cmin = 0;
  cmax = 0;
  for (i = 0; i < 1000; i ++) {

    t = rjmcmc_random_choose_int(vmin, vmax, rjmcmc_uniform);
    if (t < vmin || t > vmax) {
      c ++;
    }

    if (t == vmin) {
      cmin ++;
    }
    if (t == vmax) {
      cmax ++;
    }
  }

  CU_ASSERT(c == 0);
  CU_ASSERT(cmin > 0);
  CU_ASSERT(cmax > 0);
}

