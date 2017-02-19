
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/rjmcmc_config.h>
#include <rjmcmc/rjmcmc_random.h>
#include <rjmcmc/rjmcmc_util.h>
#include <rjmcmc/part1d_forwardmodel.h>

static double zero(void) {
  return 0.0;
}

static double half(void) {
  return 0.5;
}

static double quarter(void) {
  return 0.25;
}

int init_random_suite(void)
{
  return 0;
}

int clean_random_suite(void)
{
  return 0;
}

void test_simple_natural(void);
void test_complex_natural(void);

void test_simple_zero(void);
void test_complex_zero(void);

void test_birth_natural(void);
void test_death_natural(void);

void test_birth_zero(void);
void test_death_zero(void);
void test_death_zero_complex(void);


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

  pSuite = CU_add_suite("Part1d Forward Model Suite", 
			init_random_suite, 
			clean_random_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Simple Natural", test_simple_natural);
  ADD_TEST("Complex Natural", test_complex_natural);

  ADD_TEST("Simple Zero", test_simple_zero);
  ADD_TEST("Complex Zero", test_complex_zero);

  ADD_TEST("Natural Birth", test_birth_natural);
  ADD_TEST("Natural Death", test_death_natural);

  ADD_TEST("Zero Birth", test_birth_zero);
  ADD_TEST("Zero Death", test_death_zero);
  ADD_TEST("Zero Death Complex", test_death_zero_complex);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

void test_simple_natural(void)
{
  part1d_forwardmodel_t *p;
  double **local_parameters;

  double t;

  local_parameters = rjmcmc_create_array_2d(2, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 1.0;

  p = part1d_forwardmodel_create(PART1D_FM_NATURAL,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       2,
						       NULL,
						       local_parameters) >= 0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 1.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.5);

}

void test_complex_natural(void)
{
  part1d_forwardmodel_t *p;
  double partitions[1] = {0.5};
  double **local_parameters;

  double t;

  local_parameters = rjmcmc_create_array_2d(3, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 0.0;
  local_parameters[2][0] = 1.0;

  p = part1d_forwardmodel_create(PART1D_FM_NATURAL,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       3,
						       partitions,
						       local_parameters) >= 0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 1.0);

}

void test_simple_zero(void)
{
  part1d_forwardmodel_t *p;
  double partitions[2] = {0.0, 1.0};
  double **local_parameters;

  double t;

  local_parameters = rjmcmc_create_array_2d(2, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 1.0;

  p = part1d_forwardmodel_create(PART1D_FM_ZERO,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       2,
						       partitions,
						       local_parameters) >= 0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

}

void test_complex_zero(void)
{
  part1d_forwardmodel_t *p;
  double partitions[1] = {0.5};
  double **local_parameters;

  double t;

  local_parameters = rjmcmc_create_array_2d(3, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 0.0;
  local_parameters[2][0] = 1.0;

  p = part1d_forwardmodel_create(PART1D_FM_ZERO,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       3,
						       partitions,
						       local_parameters) >= 0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 1.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(p, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 1.0);

}

void test_birth_natural(void)
{
  part1d_forwardmodel_t *p, *q;
  double partitions[2] = {0.0, 1.0};
  double **local_parameters;

  double t;

  forwardmodelparameter_t flp;

  local_parameters = rjmcmc_create_array_2d(2, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 1.0;

  p = part1d_forwardmodel_create(PART1D_FM_NATURAL,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  q = part1d_forwardmodel_create(PART1D_FM_NATURAL,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(q != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       2,
						       partitions,
						       local_parameters) >= 0);

  part1d_forwardmodel_clone(p, q);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 1.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.5);

  flp.fmin = -0.25;
  flp.fmax = 1.25;
  flp.fstd_value = 0.25;
  flp.fstd_bd = 0.25;

  CU_ASSERT(part1d_forwardmodel_propose_birth(p, 
					      q,
					      0,
					      NULL,
					      1,
					      &flp,
					      half,
					      quarter,
					      &t) == 1);

}

void test_death_natural(void)
{
  part1d_forwardmodel_t *p, *q;
  double partitions[1] = {0.5};
  double **local_parameters;

  double t;

  forwardmodelparameter_t flp;

  local_parameters = rjmcmc_create_array_2d(3, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 0.125;
  local_parameters[2][0] = 0.0;

  p = part1d_forwardmodel_create(PART1D_FM_NATURAL,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  q = part1d_forwardmodel_create(PART1D_FM_NATURAL,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(q != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       3,
						       partitions,
						       local_parameters) >= 0);
  part1d_forwardmodel_clone(p, q);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.125);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  flp.fmin = -0.25;
  flp.fmax = 1.25;
  flp.fstd_value = 0.25;
  flp.fstd_bd = 0.25;

  CU_ASSERT(part1d_forwardmodel_propose_death(p, 
					      q,
					      0,
					      NULL,
					      1,
					      &flp,
					      half,
					      quarter,
					      &t) == 1);
}

void test_birth_zero(void)
{
  part1d_forwardmodel_t *p, *q;
  double partitions[2] = {0.0, 1.0};
  double **local_parameters;

  double t;

  forwardmodelparameter_t flp;

  local_parameters = rjmcmc_create_array_2d(2, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 1.0;

  p = part1d_forwardmodel_create(PART1D_FM_ZERO,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  q = part1d_forwardmodel_create(PART1D_FM_ZERO,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(q != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       2,
						       partitions,
						       local_parameters) >= 0);

  part1d_forwardmodel_clone(p, q);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  flp.fmin = -0.25;
  flp.fmax = 1.25;
  flp.fstd_value = 0.25;
  flp.fstd_bd = 0.25;

  CU_ASSERT(part1d_forwardmodel_propose_birth(p, 
					      q,
					      0,
					      NULL,
					      1,
					      &flp,
					      half,
					      quarter,
					      &t) == 1);


  CU_ASSERT(rjmcmc_gaussian_probability(0.0625, flp.fstd_bd) == t);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0625);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.0625);


}

void test_death_zero(void)
{
  part1d_forwardmodel_t *p, *q;
  double partitions[1] = {0.5};
  double **local_parameters;

  double t;

  forwardmodelparameter_t flp;

  local_parameters = rjmcmc_create_array_2d(3, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 0.125;
  local_parameters[2][0] = 0.0625;

  p = part1d_forwardmodel_create(PART1D_FM_ZERO,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  q = part1d_forwardmodel_create(PART1D_FM_ZERO,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(q != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       3,
						       partitions,
						       local_parameters) >= 0);
  part1d_forwardmodel_clone(p, q);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0625);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.0625);

  flp.fmin = -0.25;
  flp.fmax = 1.25;
  flp.fstd_value = 0.25;
  flp.fstd_bd = 0.25;

  CU_ASSERT(part1d_forwardmodel_propose_death(p, 
					      q,
					      0,
					      NULL,
					      1,
					      &flp,
					      half,
					      quarter,
					      &t) == 1);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

}

void test_death_zero_complex(void)
{
  part1d_forwardmodel_t *p, *q;
  double partitions[3] = {0.5, 0.25, 0.75};
  double **local_parameters;

  double t;

  forwardmodelparameter_t flp;

  local_parameters = rjmcmc_create_array_2d(5, 1);
  CU_ASSERT_FATAL(local_parameters != NULL);

  local_parameters[0][0] = 0.0;
  local_parameters[1][0] = 0.125;
  local_parameters[2][0] = 0.0625;
  local_parameters[3][0] = 0.0;
  local_parameters[4][0] = 0.125;

  p = part1d_forwardmodel_create(PART1D_FM_ZERO,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(p != NULL);

  q = part1d_forwardmodel_create(PART1D_FM_ZERO,
				 2,
				 10,
				 0.0,
				 1.0,
				 0.1,
				 0,
				 1,
				 0);
  CU_ASSERT_FATAL(q != NULL);

  CU_ASSERT_FATAL(part1d_forwardmodel_initialize_fixed(p,
						       NULL,
						       NULL,
						       5,
						       partitions,
						       local_parameters) >= 0);
  part1d_forwardmodel_clone(p, q);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.125);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.0625);

  flp.fmin = -0.25;
  flp.fmax = 1.25;
  flp.fstd_value = 0.25;
  flp.fstd_bd = 0.25;

  CU_ASSERT(part1d_forwardmodel_propose_death(p, 
					      q,
					      0,
					      NULL,
					      1,
					      &flp,
					      zero,
					      quarter,
					      &t) == 1);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 1.0, 1, &t) >= 0);
  CU_ASSERT(t == 0.125);

  CU_ASSERT_FATAL(part1d_forwardmodel_value_at(q, 0.5, 1, &t) >= 0);
  CU_ASSERT(t == 0.0);

}
