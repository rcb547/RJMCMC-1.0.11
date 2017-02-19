
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/resultset2dfm.h>

int init_dataset1d_suite(void)
{
  return 0;
}

int clean_dataset1d_suite(void)
{
  return 0;
}

void test_save_load_1(void);

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

  ADD_TEST("Save/Load 01", test_save_load_1);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

void test_save_load_1(void)
{
  const int PI = 1;
  const int NPROPOSE = 31;
  const int AI = 2;
  const int NACCEPT = 17;

  resultset2dfm_t *src;
  resultset2dfm_t *dst;

  int burnin;
  int total;
  int thin;
  int nglobalparameters = 2;
  forwardmodelparameter_t global_parameters[2];
  int nlocalparameters = 3;
  forwardmodelparameter_t local_parameters[3];
  int nhierarchical;
  int xsamples;
  int ysamples;
  int zsamples;
  int max_partitions;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int nprocesses;
  double credible_interval;
  int results;

  /* Loaded parameters for comparison */
  int l_results;
  int l_burnin;
  int l_total;
  int l_thin;
  int l_xsamples;
  int l_ysamples;
  int l_zsamples;
  int l_nprocesses;
  int l_maxpartitions;
  double l_xmin;
  double l_xmax;
  double l_ymin;
  double l_ymax;
  double l_credibleinterval;
  int l_nhierarchical;
  int l_nglobal;
  int l_nlocal;

  forwardmodelparameter_t *dst_global;
  forwardmodelparameter_t *dst_local;

  int i;
  int dstnprocesses;
  const int *p;

  /*
   * Set some reasonable values.
   */
  burnin = 1000;
  total = 10000;
  thin = 3;
  nhierarchical = 1;
  xsamples = 100;
  ysamples = 150;
  zsamples = 200;
  max_partitions = 13;
  
  xmin = -1.5;
  xmax = 2.5;
  ymin = -3.5;
  ymax = 4.5;
  nprocesses = 5;
  credible_interval = 0.625;
  results = RESULTSET2DFM_MEAN;

  /*
   * Create some tests parameters
   */
  for (i = 0; i < nglobalparameters; i ++) {
    global_parameters[i].fmin = (double)i*10.0 + 1.0;
    global_parameters[i].fmax = (double)i*10.0 + 2.0;
    global_parameters[i].fstd_value = (double)i*10.0 + 3.0;
    global_parameters[i].fstd_bd = (double)i*10.0 + 4.0;
  }

  for (i = 0; i < nlocalparameters; i ++) {
    local_parameters[i].fmin = (double)(i + nglobalparameters)*10.0 + 1.0;
    local_parameters[i].fmax = (double)(i + nglobalparameters)*10.0 + 2.0;
    local_parameters[i].fstd_value = (double)(i + nglobalparameters)*10.0 + 3.0;
    local_parameters[i].fstd_bd = (double)(i + nglobalparameters)*10.0 + 4.0;
  }
  
  /*
   * Create the resultset
   */

  src = resultset2dfm_create(burnin,
			     total,
			     thin,
			     nglobalparameters,
			     global_parameters,
			     nlocalparameters,
			     local_parameters,
			     nhierarchical,
			     xsamples,
			     ysamples,
			     zsamples,
			     max_partitions,
			     xmin,
			     xmax,
			     ymin,
			     ymax,
			     nprocesses,
			     credible_interval,
			     results);
  CU_ASSERT_FATAL(src != NULL);

  /*
   * Fill in some data 
   */
  for (i = 0; i < NPROPOSE; i ++) {
    resultset2dfm_propose(src, PI);
  }

  for (i = 0; i < NACCEPT; i ++) {
    resultset2dfm_accept(src, AI);
  }

  /*
   * Save the results
   */

  CU_ASSERT_FATAL(resultset2dfm_save(src, 
				     "resultset2dfm_saveload_1.dat",
				     RESULTSET2DFM_BINARY) == 0);

  /*
   * Load the results.
   */
  dst = resultset2dfm_load("resultset2dfm_saveload_1.dat",
			   RESULTSET2DFM_BINARY,
			   &l_results,
			   &l_burnin,
			   &l_total,
			   &l_thin,
			   &l_xsamples,
			   &l_ysamples,
			   &l_zsamples,
			   &l_nprocesses,
			   &l_maxpartitions,
			   &l_xmin,
			   &l_xmax,
			   &l_ymin,
			   &l_ymax,
			   &l_credibleinterval,
			   &l_nhierarchical,
			   &l_nglobal,
			   &l_nlocal,
			   &dst_global,
			   &dst_local);
  CU_ASSERT_FATAL(dst != NULL);

  CU_ASSERT(results == l_results);
  CU_ASSERT(burnin == l_burnin);
  CU_ASSERT(total == l_total);
  CU_ASSERT(thin == l_thin);
  CU_ASSERT(xsamples == l_xsamples);
  CU_ASSERT(ysamples == l_ysamples);
  CU_ASSERT(zsamples == l_zsamples);
  CU_ASSERT(nprocesses == l_nprocesses);
  CU_ASSERT(max_partitions == l_maxpartitions);
  CU_ASSERT(xmin == l_xmin);
  CU_ASSERT(xmax == l_xmax);
  CU_ASSERT(ymin == l_ymin);
  CU_ASSERT(ymax == l_ymax);
  CU_ASSERT(credible_interval == l_credibleinterval);
  CU_ASSERT(nhierarchical == l_nhierarchical);
  CU_ASSERT(nglobalparameters == l_nglobal);
  CU_ASSERT(nlocalparameters == l_nlocal);

  /*
   * Verify the results match the source
   */
  CU_ASSERT(resultset2dfm_get_nparameters(dst) ==
	    resultset2dfm_get_nparameters(src));

  p = resultset2dfm_get_propose(dst, &dstnprocesses);
  CU_ASSERT(dstnprocesses == nprocesses);
  CU_ASSERT_FATAL(p != NULL);
  for (i = 0; i < dstnprocesses; i ++) {
    if (i == PI) {
      CU_ASSERT(p[i] == NPROPOSE);
    } else {
      CU_ASSERT(p[i] == 0);
    }
  }

  p = resultset2dfm_get_accept(dst, &dstnprocesses);
  CU_ASSERT(dstnprocesses == nprocesses);
  CU_ASSERT_FATAL(p != NULL);
  for (i = 0; i < dstnprocesses; i ++) {
    if (i == AI) {
      CU_ASSERT(p[i] == NACCEPT);
    } else {
      CU_ASSERT(p[i] == 0);
    }
  }

  for (i = 0; i < nglobalparameters; i ++) {
    CU_ASSERT(global_parameters[i].fmin = dst_global[i].fmin);
    CU_ASSERT(global_parameters[i].fmax = dst_global[i].fmax);
    CU_ASSERT(global_parameters[i].fstd_value = dst_global[i].fstd_value);
    CU_ASSERT(global_parameters[i].fstd_bd = dst_global[i].fstd_bd);
  }

  for (i = 0; i < nlocalparameters; i ++) {
    CU_ASSERT(local_parameters[i].fmin = dst_local[i].fmin);
    CU_ASSERT(local_parameters[i].fmax = dst_local[i].fmax);
    CU_ASSERT(local_parameters[i].fstd_value = dst_local[i].fstd_value);
    CU_ASSERT(local_parameters[i].fstd_bd = dst_local[i].fstd_bd);
  }

  forwardmodelparameter_destroy(dst_global);
  forwardmodelparameter_destroy(dst_local);

  resultset2dfm_destroy(src);
  resultset2dfm_destroy(dst);
}
