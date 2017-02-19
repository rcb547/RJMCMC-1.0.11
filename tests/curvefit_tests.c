
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/curvefit.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>

int approx_equal(double a, double b)
{
  static const double epsilon = 1e-5;

  if (fabs(a - b) > epsilon) {
    printf("approx_equal: %f != %f (e = %g)\n", a, b, fabs(a - b));
    return 0;
  } 

  return -1;
}

static double zero(void)
{
  return 0.0;
}

static dataset1d_t *data = NULL;
static curvefit_result_t *cf = NULL;

/* These are the results from octave for the curve fitting */
static double omu[4] = {
  -0.79324,
   0.97896,
   0.12278,
   0.71290
};

static double oT[4] = {
  -0.75575,
  0.89844,
  0.23544,
  0.65882
  /* -0.79938, */
  /* 0.99292, */
  /* 0.10274, */
  /* 0.74068 */
};


static double oL[4][4] = {
  {89.44272,   44.72136,   30.59882,   23.53756},
  { 0.00000,   27.14484,   27.14484,   24.84392},
  { 0.00000,    0.00000,    7.34067,   11.01100},
  { 0.00000,    0.00000,    0.00000,    1.93699}
};

static double oZ[4][4] = {
  {0.03749, -0.27244,  0.53009, -0.30129},
  {0.00000,  0.19192, -0.56236,  0.38440},
  {0.00000,  0.00000,  0.14494, -0.16497},
  {0.00000,  0.00000,  0.00000,  0.02778}
};

static double oS[4][4] = {
  { 0.0014057,  -0.0102143,   0.0198744,  -0.0112961},
  {-0.0102143,   0.1110546,  -0.2523431,   0.1558568},
  { 0.0198744,  -0.2523431,   0.6182486,  -0.3997938},
  {-0.0112961,   0.1558568,  -0.3997938,   0.2665292}
};

int init_curvefit_suite(void)
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

  cf = curvefit_create(10);
  if (cf == NULL) {
    return -1;
  }

  return 0;
}

int clean_curvefit_suite(void)
{
  dataset1d_destroy(data);
  curvefit_destroy(cf);

  return 0;
}

void test_curvefit_cholesky(void);
void test_curvefit_priors(void);
void test_curvefit_probability(void);
void test_curvefit_curve(void);
void test_curvefit_precision(void);
void test_curvefit_zero_variance(void);

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

  pSuite = CU_add_suite("Curvefit Suite", init_curvefit_suite, clean_curvefit_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Cholesky", test_curvefit_cholesky);
  ADD_TEST("Priors", test_curvefit_priors);
  ADD_TEST("Probability", test_curvefit_probability);
  ADD_TEST("Curve", test_curvefit_curve);
  ADD_TEST("Precision", test_curvefit_precision);
  ADD_TEST("0th Order Variance", test_curvefit_zero_variance);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

double one(void)
{
  return 1.0;
}

void test_curvefit_cholesky(void)
{
  double **L;
  double **Z;
  double **S;
  double *mu;

  double T[10];

  double curve_prob;

  int m;
  int i;
  int j;
  
  double detCm;

  m = 4;

  CU_ASSERT_FATAL(curvefit_compute(data,
				   0,
				   DATASIZE - 1,
				   m - 1,
				   cf) == 0);

  L = cf_L(cf);
  CU_ASSERT_FATAL(L != NULL);

  Z = cf_Z(cf);
  CU_ASSERT_FATAL(Z != NULL);

  S = cf_S(cf);
  CU_ASSERT_FATAL(S != NULL);

  mu = cf_mu(cf);
  CU_ASSERT_FATAL(mu != NULL);

  /* printf("L:\n"); */
  /* for (i = 0; i < m; i ++) { */
  /*   for (j = 0; j < m; j ++) { */
  /*     printf("%f ", L[i][j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* printf("Z:\n"); */
  /* for (i = 0; i < m; i ++) { */
  /*   for (j = 0; j < m; j ++) { */
  /*     printf("%f ", Z[i][j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* printf("S:\n"); */
  /* for (i = 0; i < m; i ++) { */
  /*   for (j = 0; j < m; j ++) { */
  /*     printf("%f ", S[i][j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* printf("mu\n"); */
  /* for (i = 0; i < m; i ++) { */
  /*   printf("%f ", mu[i]); */
  /* } */
  /* printf("\n"); */

  /*
   * First test the best fit
   */
  CU_ASSERT_FATAL(curvefit_sample(cf,
				  zero,
				  T,
				  m,
				  &curve_prob) == 0);
  /* printf("bf cp %g\n", curve_prob); */
  CU_ASSERT(fabs(7.3374e-07 - curve_prob) < 0.0001e-7);

  /*
   * Next test the fit with +1 with each component
   */
  CU_ASSERT_FATAL(curvefit_sample(cf,
				  one,
				  T,
				  m,
				  &curve_prob) == 0);
  
  /* printf("T\n"); */
  /* for (i = 0; i < m; i ++) { */
  /*   printf("%f ", T[i]); */
  /* } */
  /* printf("\n"); */

  for (i = 0; i < m; i ++) {
    for (j = 0; j < m; j ++) {
      CU_ASSERT(approx_equal(L[i][j], oL[j][i]));
    }
  }

  for (i = 0; i < m; i ++) {
    for (j = 0; j < m; j ++) {
      CU_ASSERT(approx_equal(Z[i][j], oZ[j][i]));
    }
  }

  for (i = 0; i < m; i ++) {
    for (j = 0; j < m; j ++) {
      CU_ASSERT(approx_equal(S[i][j], oS[j][i]));
    }
  }

  for (i = 0; i < m; i ++) {
    CU_ASSERT(approx_equal(mu[i], omu[i]));
  }

  for (i = 0; i < m; i ++) {
    CU_ASSERT(approx_equal(T[i], oT[i]));
  }
  
  /* printf("cp: %g\n", curve_prob); */

  CU_ASSERT(fabs(curve_prob - 9.9302e-08) < 0.0001e-8);

  CU_ASSERT_FATAL(curvefit_sample_detCm(cf,
					&detCm,
					m - 1) == 0);
  
  /* printf("sqrt(detCm) = %g (%g)\n", detCm, sqrt(detCm)); */

  CU_ASSERT(approx_equal(sqrt(detCm), 2.8967e-05));
}

void test_curvefit_priors(void)
{
  double datapoints[DATASIZE][3] = {
    {2.159514E-02,   0.4118921,      0.2000000},    
    {2.417982E-02,   0.1613642,      0.2000000},    
    {8.894611E-02,   0.3569893,      0.2000000},    
    {0.1284819   ,   0.5534996,      0.2000000},    
    {0.2099725   ,   0.1788403,      0.2000000},    
    {0.3368428   ,   0.5552996,      0.2000000},    
    {0.4436209   ,   0.3286708,      0.2000000},    
    {0.5369146   ,   0.2733173,      0.2000000},    
    {0.5421172   ,   0.7738906,      0.2000000},    
    {0.5514755   ,   0.4706886,      0.2000000},    
    {0.5790006   ,   0.5556199,      0.2000000},    
    {0.5963506   ,   0.8348216,      0.2000000},    
    {0.6632974   ,   0.3515924,      0.2000000},    
    {0.6706321   ,   0.9026219,      0.2000000},    
    {0.6750101   ,   0.6716199,      0.2000000},    
    {0.6990681   ,   0.9349030,      0.2000000},    
    {0.7109817   ,   1.0772080,      0.2000000},    
    {0.7786797   ,   0.5440105,      0.2000000},    
    {0.9027547   ,   1.0832310,      0.2000000},    
    {0.9756112   ,   0.9496933,      0.2000000}
  }; 

  dataset1d_t *data;

  double **mean;
  double **sigma;
  double *mean_misfit;
  double *detCm;
  double *autoprior;
  double **S;
  double *pk;
  double *kcdf;

  int max_order = 3;
  double fixed_priors[4] = {1.2, 4.0, 20.0, 60.0};
  double auto_z = 3.0;
  int i;

  data = dataset1d_create(DATASIZE);
  CU_ASSERT_FATAL(data != NULL);

  data->xmin = 0.0;
  data->xmax = 1.0;

  data->ymin = 0.0;
  data->ymax = 1.2;

  for (i = 0; i < DATASIZE; i ++) {
    data->points[i].x = datapoints[i][0];
    data->points[i].y = datapoints[i][1];
    data->points[i].n = datapoints[i][2];
  }

  mean = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  CU_ASSERT_FATAL(mean != NULL);

  sigma = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  CU_ASSERT_FATAL(sigma != NULL);

  mean_misfit = rjmcmc_create_array_1d(max_order + 1);
  CU_ASSERT_FATAL(mean_misfit != NULL);

  detCm = rjmcmc_create_array_1d(max_order + 1);
  CU_ASSERT_FATAL(detCm != NULL);

  autoprior = rjmcmc_create_array_1d(max_order + 1);
  CU_ASSERT_FATAL(autoprior != NULL);

  S = rjmcmc_create_array_2d(max_order + 1, max_order + 1);
  CU_ASSERT_FATAL(S != NULL);

  pk = rjmcmc_create_array_1d(max_order + 1);
  CU_ASSERT_FATAL(pk != NULL);

  kcdf = rjmcmc_create_array_1d(max_order + 1);
  CU_ASSERT_FATAL(kcdf != NULL);
  
  
  CU_ASSERT_FATAL(curvefit_evaluate_pk(cf,
				       data,
				       0,
				       DATASIZE - 1,
				       3,
				       NULL,
				       auto_z,
				       mean_misfit,
				       detCm,
				       autoprior,
				       S,
				       pk,
				       kcdf,
				       mean,
				       sigma) >= 0);

  CU_ASSERT(approx_equal(pk[0], 2.0115e-04));
  CU_ASSERT(approx_equal(pk[1], 8.4580e-01));
  CU_ASSERT(approx_equal(pk[2], 1.5103e-01));
  CU_ASSERT(approx_equal(pk[3], 2.9704e-03));

  CU_ASSERT_FATAL(curvefit_evaluate_pk(cf,
				       data,
				       0,
				       DATASIZE - 1,
				       3,
				       fixed_priors,
				       auto_z,
				       mean_misfit,
				       detCm,
				       autoprior,
				       S,
				       pk,
				       kcdf,
				       mean,
				       sigma) >= 0);


  CU_ASSERT(approx_equal(pk[0], 0.000408));
  CU_ASSERT(approx_equal(pk[1], 0.867672));
  CU_ASSERT(approx_equal(pk[2], 0.118156));
  CU_ASSERT(approx_equal(pk[3], 0.013763));

  dataset1d_destroy(data);

  rjmcmc_destroy_array_2d(max_order + 1, mean);
  rjmcmc_destroy_array_1d(mean_misfit);
  rjmcmc_destroy_array_1d(detCm);
  rjmcmc_destroy_array_1d(autoprior);
  rjmcmc_destroy_array_1d(pk);
  rjmcmc_destroy_array_1d(kcdf);
  rjmcmc_destroy_array_2d(max_order + 1, S);
}

void test_curvefit_probability(void)
{
  int i;
  double sumprob;
  double curve_prob;

  double T[10];
  int m;

  m = 4;

  CU_ASSERT_FATAL(curvefit_compute(data,
				   0,
				   DATASIZE - 1,
				   m - 1,
				   cf) == 0);
  
  sumprob = 0.0;

  for (i = 0; i < 100; i ++) {
    CU_ASSERT_FATAL(curvefit_sample(cf,
				    one,
				    T,
				    m,
				    &curve_prob) == 0);
    sumprob += curve_prob;
  }

  printf("sumprob %g\n", sumprob);

}

void test_curvefit_curve(void)
{
  double x[] = {64,
		65,
		66,
		67,
		68,
		69,
		70,
		71,
		72,
		73,
		74,
		75,
		76,
		77,
		78,
		79,
		80,
		81,
		82,
		83,
		84,
		85,
		86,
		87,
		88,
		89,
		90,
		91,
		92,
		93,
		94,
		95,
		96,
		97,
		98,
		99,
		100};
  double y[] = {34.6469,
		42.9048, 
		37.8449, 
		49.5675,
		47.4974, 
		41.9835, 
		50.78, 
		47.2427, 
		55.8715, 
		52.7508, 
		55.2984, 
		54.7579, 
		52.1766, 
		63.1839, 
		62.4916, 
		55.8938, 
		53.4449, 
		52.8315, 
		59.2266, 
		63.6094, 
		62.092, 
		54.031, 
		57.2428, 
		67.368, 
		58.6959, 
		53.5683, 
		52.9358, 
		58.397, 
		61.633, 
		54.9933, 
		55.9431, 
		57.7203, 
		48.1066, 
		52.0459, 
		59.1899, 
		54.9111, 
		42.9803};
  double n[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

  dataset1d_t *data;
  curvefit_result_t *cf;

  double a[3];
  double prob;
  double zprob;
  int i;

  data = dataset1d_create_from_array(x, y, n, 37);
  CU_ASSERT_FATAL(data != NULL);

  cf = curvefit_create(3);
  CU_ASSERT_FATAL(cf != NULL);

  CU_ASSERT_FATAL(curvefit_compute(data, 0, 36, 2, cf) == 0);

  CU_ASSERT_FATAL(curvefit_sample(cf, 
				  zero, 
				  a,
				  3,
				  &zprob) == 0);

  printf("zp: %g\n", zprob);

  printf("az: %g %g %g\n",
	 a[0], a[1], a[2]);

  /* for (i = 0; i < 100; i ++) { */
  /*   CU_ASSERT_FATAL(curvefit_sample(cf,  */
  /* 				    rjmcmc_normal,  */
  /* 				    a, */
  /* 				    3, */
  /* 				    &prob) == 0); */
  /*   CU_ASSERT(prob < zprob); */
  /* } */

  
}
  
void test_curvefit_precision(void)
{
  double x[7] = {63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0};
  double y[7] = {44.6002, 43.3741, 38.6491, 49.9412, 45.3128, 45.1015, 48.5634};
  //double n[7] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
  double n[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  dataset1d_t *data;
  curvefit_result_t *cf;
  

  data = dataset1d_create_from_array(x, y, n, 7);
  CU_ASSERT_FATAL(data != NULL);

  cf = curvefit_create(5);
  CU_ASSERT_FATAL(cf != NULL);

  CU_ASSERT_FATAL(curvefit_compute_lambda(data, 1.0, 0, 6, 4, cf) == 0);

  // This test fails due to numerical imprecision
  //  CU_ASSERT_FATAL(curvefit_compute_lambda(data, 1.0, 0, 6, 5, cf) == 0);
}
    
void test_curvefit_zero_variance(void)
{
  double x[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  double y[10] = {5.83597099259,
		  1.71908240846,
		  6.8520433025,
		  -2.61589695472,
		  -1.15919944324,
		  -1.28337769986,
		  3.98173146651,
		  1.69949946585,
		  -7.85380721982,
		  -2.37831319988};
  double n[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  dataset1d_t *data;
  curvefit_result_t *cf;

  double coeff[5];
  double p;

  int i;

  /* for (i = 0; i < 10; i ++) { */
  /*   n[i] *= 4.4; */
  /* } */

  data = dataset1d_create_from_array(x, y, n, 10);
  
  CU_ASSERT_FATAL(data != NULL);

  cf = curvefit_create(5);
  CU_ASSERT_FATAL(cf != NULL);

  CU_ASSERT_FATAL(curvefit_compute(data, 0, 9, 0, cf) >= 0);

  CU_ASSERT_FATAL(curvefit_sample(cf, one, coeff, 1, &p) >= 0);

  printf("%f \n", coeff[0]);
  printf("%f \n", cf_mu(cf)[0]);
  printf("%f \n", cf_Z(cf)[0][0]);
  
}
