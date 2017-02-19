
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/position_map2d.h>
#include <rjmcmc/part2d_forwardmodel.h>
#include <rjmcmc/rjmcmc_util.h>

int init_part2d_forwardmodel_suite(void)
{
  return 0;
}

int clean_part2d_forwardmodel_suite(void)
{
  return 0;
}

void test_save_load_1(void);
void test_birth(void);
void test_birth_internal(void);

void test_death(void);
void test_death_internal(void);

void test_move_internal(void);

#define ADD_TEST(name, function) \
  if (CU_add_test(pSuite, name, function) == NULL) {\
  CU_cleanup_registry(); \
  return CU_get_error(); \
  }


int main(int argc, char *argv[])
{
  CU_pSuite pSuite = NULL;

  position_map2d_set_type(POSITION_MAP2D_LINEAR);
  //position_map2d_set_type(POSITION_MAP2D_DELAUNAY);

  if (CU_initialize_registry() != CUE_SUCCESS) {
    return CU_get_error();
  }

  pSuite = CU_add_suite("Part2d Forwardmodel Suite", 
			init_part2d_forwardmodel_suite, 
			clean_part2d_forwardmodel_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Save/Load 01", test_save_load_1);
  ADD_TEST("Birth", test_birth);
  ADD_TEST("Birth Internal", test_birth_internal);
  ADD_TEST("Death", test_death);
  ADD_TEST("Death Internal", test_death_internal);
  ADD_TEST("Move Internal", test_move_internal);


  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

void test_save_load_1(void)
{
  static const int NPARTITIONS = 5;
  static const int NSPOTCHECKS = 100;

  part2d_forwardmodel_t *src;
  part2d_forwardmodel_t *prop;
  part2d_forwardmodel_t *t;

  part2d_forwardmodel_t *dst;

  int nglobalparameters = 2;
  int nlocalparameters = 3;
  int nhierarchicalparameters = 1;
  forwardmodelparameter_t global_parameters[2];
  forwardmodelparameter_t local_parameters[3];
  forwardmodelparameter_t hierarchical_parameters[1];
  int max_partitions;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double pdx;
  double pdy;
  int includecorners;

  bbox2d_t bound;

  const double *a;
  const double *b;

  double src_localparameters[3];
  double dst_localparameters[3];

  double prob;

  double sx, sy;
  double dx, dy;

  int i;
  int j;

  /*
   * Set some reasonable values.
   */
  max_partitions = 13;
  
  xmin = -1.5;
  xmax = 2.5;
  ymin = -3.5;
  ymax = 4.5;

  pdx = 0.75;
  pdy = 1.25;

  includecorners = 0;

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
    local_parameters[i].fmax = (double)(i + nglobalparameters)*10.0 + 20.0;
    local_parameters[i].fstd_value = (double)(i + nglobalparameters)/10.0 + 0.03;
    local_parameters[i].fstd_bd = (double)(i + nglobalparameters)/10.0 + 0.04;
  }

  for (i = 0; i < nhierarchicalparameters; i ++) {
    hierarchical_parameters[i].fmin = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 1.0;
    hierarchical_parameters[i].fmax = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 2.0;
    hierarchical_parameters[i].fstd_value = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 3.0;
    hierarchical_parameters[i].fstd_bd = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 4.0;
  }

  /*
   * Create the initial model
   */

  src = part2d_forwardmodel_create(PART2D_FM_ZERO,
				   1,
				   max_partitions,
				   xmin, 
				   xmax,
				   ymin,
				   ymax,
				   pdx,
				   pdy,
				   nglobalparameters,
				   nlocalparameters,
				   nhierarchicalparameters,
				   includecorners);

  prop = part2d_forwardmodel_create(PART2D_FM_ZERO,
				   1,
				   max_partitions,
				   xmin, 
				   xmax,
				   ymin,
				   ymax,
				   pdx,
				   pdy,
				   nglobalparameters,
				   nlocalparameters,
				   nhierarchicalparameters,
				   includecorners);

  CU_ASSERT_FATAL(src != NULL);
  CU_ASSERT_FATAL(prop != NULL);

  CU_ASSERT_FATAL(part2d_forwardmodel_initialize(src,
						 global_parameters,
						 nglobalparameters,
						 local_parameters,
						 nlocalparameters,
						 hierarchical_parameters,
						 nhierarchicalparameters,
						 0,
						 rjmcmc_uniform,
						 rjmcmc_normal) == 0);

  /*
   * Create a few points 
   */
  for (i = 0; i < (NPARTITIONS - 1); i ++) {
    CU_ASSERT(part2d_forwardmodel_propose_birth(src,
						prop,
						nglobalparameters,
						global_parameters,
						nlocalparameters,
						local_parameters,
						rjmcmc_uniform,
						rjmcmc_normal,
						&bound,
						&prob) != 1);
    t = src;
    src = prop;
    prop = t;
  }
  CU_ASSERT_FATAL(part2d_forwardmodel_partitions(src) == NPARTITIONS);
  for (i = 0; i < NPARTITIONS; i ++) {
    CU_ASSERT_FATAL(part2d_forwardmodel_partition_centre(src, i, &sx, &sy) == 0);
  }
    

  part2d_forwardmodel_destroy(prop);
				      
  /*
   * Save the results
   */

  CU_ASSERT_FATAL(part2d_forwardmodel_save(src, 
					   "part2d_forwardmodel_saveload_1.dat") == 0);

  /*
   * Load the results.
   */
  dst = part2d_forwardmodel_load("part2d_forwardmodel_saveload_1.dat");
  CU_ASSERT_FATAL(dst != NULL);

  /*
   * Verify the results match the source
   */
  CU_ASSERT(part2d_forwardmodel_type(src) ==
	    part2d_forwardmodel_type(dst));
  CU_ASSERT(part2d_forwardmodel_min_partitions(src) ==
	    part2d_forwardmodel_min_partitions(dst));
  CU_ASSERT(part2d_forwardmodel_max_partitions(src) ==
	    part2d_forwardmodel_max_partitions(dst));

  CU_ASSERT(part2d_forwardmodel_pdx(src) ==
	    part2d_forwardmodel_pdx(dst));
  CU_ASSERT(part2d_forwardmodel_pdy(src) ==
	    part2d_forwardmodel_pdy(dst));

  /* Global Parameters */
  a = part2d_forwardmodel_global_parameters(src);
  b = part2d_forwardmodel_global_parameters(dst);
  CU_ASSERT_FATAL(a != NULL);
  CU_ASSERT_FATAL(b != NULL);
  for (i = 0; i < nglobalparameters; i ++) {
    CU_ASSERT(a[i] == b[i]);
  }

  /* Hierarchical Parameters */
  a = part2d_forwardmodel_global_parameters(src);
  b = part2d_forwardmodel_global_parameters(dst);
  CU_ASSERT_FATAL(a != NULL);
  CU_ASSERT_FATAL(b != NULL);
  for (i = 0; i < nglobalparameters; i ++) {
    CU_ASSERT(a[i] == b[i]);
  }

  /* Partitions */
  CU_ASSERT(part2d_forwardmodel_partitions(src) == 
	    part2d_forwardmodel_partitions(dst));

  for (i = 0; i < NPARTITIONS; i ++) {
    CU_ASSERT_FATAL(part2d_forwardmodel_partition_centre(src, i + 4, &sx, &sy) == 0);
    CU_ASSERT_FATAL(part2d_forwardmodel_partition_centre(dst, i + 4, &dx, &dy) == 0);

    CU_ASSERT(sx == dx);
    CU_ASSERT(sy == dy);

  }

  for (i = 0; i < NSPOTCHECKS; i ++) {

    sx = rjmcmc_random_choose_double(xmin, xmax, rjmcmc_uniform);
    sy = rjmcmc_random_choose_double(ymin, ymax, rjmcmc_uniform);

    CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src,
  						 sx,
  						 sy,
  						 nlocalparameters,
  						 src_localparameters) == 0);
 
    CU_ASSERT_FATAL(part2d_forwardmodel_value_at(dst,
  						 sx,
  						 sy,
  						 nlocalparameters,
  						 dst_localparameters) == 0);

    for (j = 0; j < nlocalparameters; j ++) {
      CU_ASSERT(src_localparameters[j] == dst_localparameters[j]);
    }
  }
   
  part2d_forwardmodel_destroy(src);
  part2d_forwardmodel_destroy(dst);
}

static double half(void)
{
  return 0.5;
}

static double one(void)
{
  return 1.0;
}

void test_birth(void)
{
  static const int NPARTITIONS = 5;
  static const int NSPOTCHECKS = 100;

  part2d_forwardmodel_t *src;
  part2d_forwardmodel_t *prop;
  part2d_forwardmodel_t *t;

  part2d_forwardmodel_t *dst;

  int nglobalparameters = 2;
  int nlocalparameters = 3;
  int nhierarchicalparameters = 1;
  forwardmodelparameter_t global_parameters[2];
  forwardmodelparameter_t local_parameters[3];
  forwardmodelparameter_t hierarchical_parameters[1];
  int max_partitions;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double pdx;
  double pdy;
  int includecorners;

  bbox2d_t bound;

  const double *a;
  const double *b;

  double src_localparameters[3];
  double dst_localparameters[3];

  double prob;

  double sx, sy;
  double dx, dy;

  int i;
  int j;

  double expected_prob;

  /*
   * Set some reasonable values.
   */
  max_partitions = 13;
  
  xmin = -1.5;
  xmax = 2.5;
  ymin = -3.5;
  ymax = 4.5;

  pdx = 0.75;
  pdy = 1.25;

  includecorners = 0;

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
    local_parameters[i].fmax = (double)(i + nglobalparameters)*10.0 + 20.0;
    local_parameters[i].fstd_value = (double)(i + nglobalparameters)*0.01 + 0.03;
    local_parameters[i].fstd_bd = (double)(i + nglobalparameters)*0.01 + 0.04;
  }

  for (i = 0; i < nhierarchicalparameters; i ++) {
    hierarchical_parameters[i].fmin = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 1.0;
    hierarchical_parameters[i].fmax = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 2.0;
    hierarchical_parameters[i].fstd_value = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 3.0;
    hierarchical_parameters[i].fstd_bd = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 4.0;
  }

  /*
   * Create the initial model
   */

  src = part2d_forwardmodel_create(PART2D_FM_ZERO,
				   1,
				   max_partitions,
				   xmin, 
				   xmax,
				   ymin,
				   ymax,
				   pdx,
				   pdy,
				   nglobalparameters,
				   nlocalparameters,
				   nhierarchicalparameters,
				   includecorners);

  prop = part2d_forwardmodel_create(PART2D_FM_ZERO,
				    1,
				    max_partitions,
				    xmin, 
				    xmax,
				    ymin,
				    ymax,
				    pdx,
				    pdy,
				    nglobalparameters,
				    nlocalparameters,
				    nhierarchicalparameters,
				    includecorners);

  sx = (xmax + xmin)/2.0;
  sy = (ymax + ymin)/2.0;

  CU_ASSERT_FATAL(src != NULL);
  CU_ASSERT_FATAL(prop != NULL);

  CU_ASSERT_FATAL(part2d_forwardmodel_initialize(src,
						 global_parameters,
						 nglobalparameters,
						 local_parameters,
						 nlocalparameters,
						 hierarchical_parameters,
						 nhierarchicalparameters,
						 0,
						 rjmcmc_uniform,
						 rjmcmc_normal) == 0);


  /*
   * Evaulate the model at the centre
   */
  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src,
					       sx,
					       sy,
					       nlocalparameters,
					       src_localparameters) == 0);


  /*
   * Propose a birth with fixed random functions. This should create
   * a new point in the centre, one standard deviation above the existing
   * level.
   */
  CU_ASSERT(part2d_forwardmodel_propose_birth(src,
					      prop,
					      nglobalparameters,
					      global_parameters,
					      nlocalparameters,
					      local_parameters,
					      half,
					      one,
					      &bound,
					      &prob) != 0);

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(prop,
					       sx,
					       sy,
					       nlocalparameters,
					       dst_localparameters) == 0);

  CU_ASSERT_FATAL((part2d_forwardmodel_partitions(src) + 1) ==
		  part2d_forwardmodel_partitions(prop));

  expected_prob = 1.0;

  for (i = 0; i < nlocalparameters; i ++) {
    
    CU_ASSERT(dst_localparameters[i] > src_localparameters[i]);
    CU_ASSERT(fabs((dst_localparameters[i] - src_localparameters[i]) - 
		   local_parameters[i].fstd_bd) < 1.0e-3);

    expected_prob *= rjmcmc_gaussian_probability(local_parameters[i].fstd_bd,
						 local_parameters[i].fstd_bd);
  }

  CU_ASSERT(fabs(expected_prob - prob) < 1.0e-6);

}

static double quarter(void)
{
  return 0.25;
}

static double negativehalf(void)
{
  return -0.5;
}

static double zero(void)
{
  return 0.0;
}


void test_death(void) 
{
  part2d_forwardmodel_t *src;
  part2d_forwardmodel_t *prop;
  part2d_forwardmodel_t *t;

  part2d_forwardmodel_t *dst;

  int nglobalparameters = 2;
  int nlocalparameters = 3;
  int nhierarchicalparameters = 1;
  forwardmodelparameter_t global_parameters[2];
  forwardmodelparameter_t local_parameters[3];
  forwardmodelparameter_t hierarchical_parameters[1];
  int max_partitions;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double pdx;
  double pdy;
  int includecorners;

  bbox2d_t bound;

  double src_c_localparameters[3];
  double src_bl_localparameters[3];

  double c_localparameters[3];
  double bl_localparameters[3];

  double test_localparameters[3];

  double prob;

  double sx, sy;
  double dx, dy;

  double px, py;

  int i;

  /*
   * Set some reasonable values.
   */
  max_partitions = 13;
  
  xmin = -1.5;
  xmax = 2.5;
  ymin = -3.5;
  ymax = 4.5;

  pdx = 0.75;
  pdy = 1.25;

  includecorners = 0;

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
    local_parameters[i].fmax = (double)(i + nglobalparameters)*10.0 + 20.0;
    local_parameters[i].fstd_value = (double)(i + nglobalparameters)*0.01 + 0.03;
    local_parameters[i].fstd_bd = (double)(i + nglobalparameters)*0.01 + 0.04;
  }

  for (i = 0; i < nhierarchicalparameters; i ++) {
    hierarchical_parameters[i].fmin = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 1.0;
    hierarchical_parameters[i].fmax = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 2.0;
    hierarchical_parameters[i].fstd_value = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 3.0;
    hierarchical_parameters[i].fstd_bd = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 4.0;
  }

  /*
   * Create the initial model
   */

  src = part2d_forwardmodel_create(PART2D_FM_ZERO,
				   1,
				   max_partitions,
				   xmin, 
				   xmax,
				   ymin,
				   ymax,
				   pdx,
				   pdy,
				   nglobalparameters,
				   nlocalparameters,
				   nhierarchicalparameters,
				   includecorners);

  prop = part2d_forwardmodel_create(PART2D_FM_ZERO,
				    1,
				    max_partitions,
				    xmin, 
				    xmax,
				    ymin,
				    ymax,
				    pdx,
				    pdy,
				    nglobalparameters,
				    nlocalparameters,
				    nhierarchicalparameters,
				    includecorners);


  sx = (xmax + xmin)/2.0;
  sy = (ymax + ymin)/2.0;

  dx = xmin + (xmax - xmin)/4.0;
  dy = ymin + (ymax - ymin)/4.0;

  CU_ASSERT_FATAL(src != NULL);
  CU_ASSERT_FATAL(prop != NULL);

  CU_ASSERT_FATAL(part2d_forwardmodel_initialize(src,
						 global_parameters,
						 nglobalparameters,
						 local_parameters,
						 nlocalparameters,
						 hierarchical_parameters,
						 nhierarchicalparameters,
						 0,
						 rjmcmc_uniform,
						 rjmcmc_normal) == 0);

  /*
   * Evaulate the model at the centre before adding points
   */
  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src,
					       sx,
					       sy,
					       nlocalparameters,
					       src_c_localparameters) == 0);
  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src,
					       sx,
					       sy,
					       nlocalparameters,
					       src_bl_localparameters) == 0);

  /*
   * Propose a birth with fixed random functions. This should create
   * a new point in the centre, one standard deviation above the existing
   * level.
   */
  CU_ASSERT(part2d_forwardmodel_propose_birth(src,
					      prop,
					      nglobalparameters,
					      global_parameters,
					      nlocalparameters,
					      local_parameters,
					      half,
					      one,
					      &bound,
					      &prob) != 0);

  CU_ASSERT_FATAL((part2d_forwardmodel_partitions(src) + 1) ==
		  part2d_forwardmodel_partitions(prop));

  /* Swap models around */
  t = src;
  src = prop;
  prop = t;

  /*
   * Propose a second birth with fixed random functions. This should create
   * a new point in bottom left , half a standard deviation below the existing
   * level.
   */
  CU_ASSERT(part2d_forwardmodel_propose_birth(src,
					      prop,
					      nglobalparameters,
					      global_parameters,
					      nlocalparameters,
					      local_parameters,
					      quarter,
					      negativehalf,
					      &bound,
					      &prob) != 0);

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(prop,
					       sx,
					       sy,
					       nlocalparameters,
					       c_localparameters) == 0);

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(prop,
					       dx,
					       dy,
					       nlocalparameters,
					       bl_localparameters) == 0);

  CU_ASSERT_FATAL((part2d_forwardmodel_partitions(src) + 1) ==
		  part2d_forwardmodel_partitions(prop));

  t = src;
  src = prop;
  prop = t;

  CU_ASSERT(part2d_forwardmodel_propose_death(src,
					      prop,
					      nglobalparameters,
					      global_parameters,
					      nlocalparameters,
					      local_parameters,
					      half,
					      rjmcmc_normal,
					      &bound,
					      &prob) != 0);

  CU_ASSERT_FATAL((part2d_forwardmodel_partitions(src) - 1) ==
		  part2d_forwardmodel_partitions(prop));

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(prop, 
					       sx, sy,
					       nlocalparameters,
					       test_localparameters) == 0);

  for (i = 0; i < nlocalparameters; i ++) {
    CU_ASSERT(src_c_localparameters[i] == test_localparameters[i]);
  }

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(prop, 
					       dx, dy,
					       nlocalparameters,
					       test_localparameters) == 0);
  for (i = 0; i < nlocalparameters; i ++) {
    CU_ASSERT(bl_localparameters[i] == test_localparameters[i]);
  }



}

void test_birth_internal(void)
{
  part2d_forwardmodel_t *src;

  int nglobalparameters = 2;
  int nlocalparameters = 3;
  int nhierarchicalparameters = 1;
  forwardmodelparameter_t global_parameters[2];
  forwardmodelparameter_t local_parameters[3];
  forwardmodelparameter_t hierarchical_parameters[1];
  int max_partitions;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double pdx;
  double pdy;
  int includecorners;

  bbox2d_t bound;

  double src_localparameters[3];
  double dst_localparameters[3];

  double parameters[3];

  double sx, sy;

  int i;

  /*
   * Set some reasonable values.
   */
  max_partitions = 13;
  
  xmin = -1.5;
  xmax = 2.5;
  ymin = -3.5;
  ymax = 4.5;

  pdx = 0.75;
  pdy = 1.25;

  includecorners = 0;

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
    local_parameters[i].fmax = (double)(i + nglobalparameters)*10.0 + 20.0;
    local_parameters[i].fstd_value = (double)(i + nglobalparameters)*0.01 + 0.03;
    local_parameters[i].fstd_bd = (double)(i + nglobalparameters)*0.01 + 0.04;

    /* Set new point parameters */
    parameters[i] = local_parameters[i].fmin;
  }

  for (i = 0; i < nhierarchicalparameters; i ++) {
    hierarchical_parameters[i].fmin = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 1.0;
    hierarchical_parameters[i].fmax = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 2.0;
    hierarchical_parameters[i].fstd_value = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 3.0;
    hierarchical_parameters[i].fstd_bd = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 4.0;
  }

  /*
   * Create the initial model
   */

  src = part2d_forwardmodel_create(PART2D_FM_ZERO,
				   1,
				   max_partitions,
				   xmin, 
				   xmax,
				   ymin,
				   ymax,
				   pdx,
				   pdy,
				   nglobalparameters,
				   nlocalparameters,
				   nhierarchicalparameters,
				   includecorners);

  sx = (xmax + xmin)/2.0;
  sy = (ymax + ymin)/2.0;

  CU_ASSERT_FATAL(src != NULL);

  CU_ASSERT_FATAL(part2d_forwardmodel_initialize(src,
						 global_parameters,
						 nglobalparameters,
						 local_parameters,
						 nlocalparameters,
						 hierarchical_parameters,
						 nhierarchicalparameters,
						 0,
						 rjmcmc_uniform,
						 rjmcmc_normal) == 0);


  /*
   * Evaulate the model at the centre
   */
  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src,
					       sx,
					       sy,
					       nlocalparameters,
					       src_localparameters) == 0);


  /*
   * Propose a birth with fixed random functions. This should create
   * a new point in the centre, one standard deviation above the existing
   * level.
   */
  CU_ASSERT(part2d_forwardmodel_addpoint(src,
					 sx, 
					 sy,
					 nlocalparameters,
					 parameters,
					 &bound) == 0);


  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx, 
					       sy,
					       nlocalparameters,
					       dst_localparameters) == 0);

  
  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters[i]);
  }

}

void test_death_internal(void)
{
  part2d_forwardmodel_t *src;

  int nglobalparameters = 2;
  int nlocalparameters = 3;
  int nhierarchicalparameters = 1;
  forwardmodelparameter_t global_parameters[2];
  forwardmodelparameter_t local_parameters[3];
  forwardmodelparameter_t hierarchical_parameters[1];
  int max_partitions;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double pdx;
  double pdy;
  int includecorners;

  bbox2d_t bound;

  double src_localparameters[3];
  double dst_localparameters[3];

  double parameters[3];
  double parameters_s[3];

  double sx, sy;

  int i;

  /*
   * Set some reasonable values.
   */
  max_partitions = 13;
  
  xmin = -1.5;
  xmax = 2.5;
  ymin = -3.5;
  ymax = 4.5;

  pdx = 0.75;
  pdy = 1.25;

  includecorners = 0;

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
    local_parameters[i].fmax = (double)(i + nglobalparameters)*10.0 + 20.0;
    local_parameters[i].fstd_value = (double)(i + nglobalparameters)*0.01 + 0.03;
    local_parameters[i].fstd_bd = (double)(i + nglobalparameters)*0.01 + 0.04;

    /* Set new point parameters */
    parameters[i] = local_parameters[i].fmin;
    parameters_s[i] = local_parameters[i].fmax;
  }

  for (i = 0; i < nhierarchicalparameters; i ++) {
    hierarchical_parameters[i].fmin = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 1.0;
    hierarchical_parameters[i].fmax = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 2.0;
    hierarchical_parameters[i].fstd_value = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 3.0;
    hierarchical_parameters[i].fstd_bd = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 4.0;
  }

  /*
   * Create the initial model
   */

  src = part2d_forwardmodel_create(PART2D_FM_ZERO,
				   1,
				   max_partitions,
				   xmin, 
				   xmax,
				   ymin,
				   ymax,
				   pdx,
				   pdy,
				   nglobalparameters,
				   nlocalparameters,
				   nhierarchicalparameters,
				   includecorners);

  sx = (xmax + xmin)/2.0;
  sy = (ymax + ymin)/2.0;

  CU_ASSERT_FATAL(src != NULL);

  CU_ASSERT_FATAL(part2d_forwardmodel_initialize(src,
						 global_parameters,
						 nglobalparameters,
						 local_parameters,
						 nlocalparameters,
						 hierarchical_parameters,
						 nhierarchicalparameters,
						 0,
						 rjmcmc_uniform,
						 rjmcmc_normal) == 0);


  /*
   * Evaulate the model at the centre
   */
  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src,
					       sx,
					       sy,
					       nlocalparameters,
					       src_localparameters) == 0);


  /*
   * Propose a birth with fixed random functions. This should create
   * a new point in the centre, one standard deviation above the existing
   * level.
   */
  CU_ASSERT(part2d_forwardmodel_addpoint(src,
					 sx, 
					 sy,
					 nlocalparameters,
					 parameters,
					 &bound) == 0);


  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx, 
					       sy,
					       nlocalparameters,
					       dst_localparameters) == 0);

  
  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters[i]);
  }

  /*
   * Add another point near the centre
   */
  CU_ASSERT(part2d_forwardmodel_addpoint(src, 
					 sx + 0.01,
					 sy + 0.01,
					 nlocalparameters,
					 parameters_s,
					 &bound) == 0);

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx + 0.01, 
					       sy + 0.01,
					       nlocalparameters,
					       dst_localparameters) == 0);

  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters_s[i]);
  }

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx, 
					       sy,
					       nlocalparameters,
					       dst_localparameters) == 0);

  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters[i]);
  }

  /*
   * Remove the first added point
   */
  CU_ASSERT(part2d_forwardmodel_delpoint(src,
					 1,
					 &bound) == 0);


  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx, 
					       sy,
					       nlocalparameters,
					       dst_localparameters) == 0);

  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters_s[i]);
  }


}

void test_move_internal(void)
{
  part2d_forwardmodel_t *src;

  int nglobalparameters = 2;
  int nlocalparameters = 3;
  int nhierarchicalparameters = 1;
  forwardmodelparameter_t global_parameters[2];
  forwardmodelparameter_t local_parameters[3];
  forwardmodelparameter_t hierarchical_parameters[1];
  int max_partitions;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double pdx;
  double pdy;
  int includecorners;

  bbox2d_t bound;

  double src_localparameters[3];
  double dst_localparameters[3];

  double parameters[3];
  double parameters_s[3];

  double sx, sy;

  int i;

  /*
   * Set some reasonable values.
   */
  max_partitions = 13;
  
  xmin = -1.5;
  xmax = 2.5;
  ymin = -3.5;
  ymax = 4.5;

  pdx = 0.75;
  pdy = 1.25;

  includecorners = 0;

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
    local_parameters[i].fmax = (double)(i + nglobalparameters)*10.0 + 20.0;
    local_parameters[i].fstd_value = (double)(i + nglobalparameters)*0.01 + 0.03;
    local_parameters[i].fstd_bd = (double)(i + nglobalparameters)*0.01 + 0.04;

    /* Set new point parameters */
    parameters[i] = local_parameters[i].fmin;
    parameters_s[i] = local_parameters[i].fmax;
  }

  for (i = 0; i < nhierarchicalparameters; i ++) {
    hierarchical_parameters[i].fmin = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 1.0;
    hierarchical_parameters[i].fmax = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 2.0;
    hierarchical_parameters[i].fstd_value = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 3.0;
    hierarchical_parameters[i].fstd_bd = (double)(i + nglobalparameters + nlocalparameters)*10.0 + 4.0;
  }

  /*
   * Create the initial model
   */

  src = part2d_forwardmodel_create(PART2D_FM_ZERO,
				   1,
				   max_partitions,
				   xmin, 
				   xmax,
				   ymin,
				   ymax,
				   pdx,
				   pdy,
				   nglobalparameters,
				   nlocalparameters,
				   nhierarchicalparameters,
				   includecorners);

  sx = (xmax + xmin)/2.0;
  sy = (ymax + ymin)/2.0;

  CU_ASSERT_FATAL(src != NULL);

  CU_ASSERT_FATAL(part2d_forwardmodel_initialize(src,
						 global_parameters,
						 nglobalparameters,
						 local_parameters,
						 nlocalparameters,
						 hierarchical_parameters,
						 nhierarchicalparameters,
						 0,
						 rjmcmc_uniform,
						 rjmcmc_normal) == 0);


  /*
   * Evaulate the model at the centre
   */
  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src,
					       sx,
					       sy,
					       nlocalparameters,
					       src_localparameters) == 0);


  /*
   * Propose a birth with fixed random functions. This should create
   * a new point in the centre, one standard deviation above the existing
   * level.
   */
  CU_ASSERT(part2d_forwardmodel_addpoint(src,
					 sx, 
					 sy,
					 nlocalparameters,
					 parameters,
					 &bound) == 0);


  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx, 
					       sy,
					       nlocalparameters,
					       dst_localparameters) == 0);

  
  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters[i]);
  }

  /*
   * Add another point near the centre
   */
  CU_ASSERT(part2d_forwardmodel_addpoint(src, 
					 sx + 0.01,
					 sy + 0.01,
					 nlocalparameters,
					 parameters_s,
					 &bound) == 0);

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx + 0.01, 
					       sy + 0.01,
					       nlocalparameters,
					       dst_localparameters) == 0);

  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters_s[i]);
  }

  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx, 
					       sy,
					       nlocalparameters,
					       dst_localparameters) == 0);

  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters[i]);
  }

  /*
   * Move the first added point
   */
  CU_ASSERT(part2d_forwardmodel_movepoint(src,
					  1,
					  sx + 0.5,
					  sy + 0.5,
					  &bound) == 0);


  CU_ASSERT_FATAL(part2d_forwardmodel_value_at(src, 
					       sx, 
					       sy,
					       nlocalparameters,
					       dst_localparameters) == 0);

  for (i = 0; i < nlocalparameters; i ++) { 
    CU_ASSERT(dst_localparameters[i] == parameters_s[i]);
  }


}
