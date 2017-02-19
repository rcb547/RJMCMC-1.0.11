
#include <stdio.h>
#include <stdlib.h>

#include <CUnit/Basic.h>

#include <rjmcmc/quadtree.h>

#define DEFAULT_LEVELS 3

int init_util_suite(void)
{
  return 0;
}

int clean_util_suite(void)
{
  return 0;
}

void test_quadtree_create(void);
void test_quadtree_random_add(void);
void test_quadtree_random_delete(void);

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

  pSuite = CU_add_suite("Quadtree Suite", init_util_suite, clean_util_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Quadtree Create", test_quadtree_create);
  ADD_TEST("Quadtree Random Add", test_quadtree_random_add);
  ADD_TEST("Quadtree Random Delete", test_quadtree_random_delete);


  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

void test_quadtree_create(void)
{
  quadtree_t *d;

  d = quadtree_create(10, 
		      DEFAULT_LEVELS,
		      -1.0,
		      1.0,
		      -1.0,
		      1.0);
  CU_ASSERT_FATAL(d != NULL);

  quadtree_destroy(d);

  d = quadtree_create(10, 
		      DEFAULT_LEVELS,
		      -2.0,
		      2.0,
		      -1.0,
		      1.0);
  CU_ASSERT(d != NULL);

}

void test_quadtree_basic_add(void)
{
  quadtree_t *d;
  bbox2d_t bound;

  d = quadtree_create(10, 
		      DEFAULT_LEVELS,
		      -1.0,
		      1.0,
		      -1.0,
		      1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(quadtree_add(d, -0.5, 0.5, &bound) >= 0);

}

void test_quadtree_random_add(void)
{
  #define NPOINTS 100

  quadtree_t *d;
  int i;
  int j;
  double x[NPOINTS];
  double y[NPOINTS];
  bbox2d_t bound;

  d = quadtree_create(100, 
		      DEFAULT_LEVELS,
		      -1.0,
		      1.0,
		      -1.0,
		      1.0);
  CU_ASSERT_FATAL(d != NULL);
  
  /*
   * Add all the random points
   */
  for (i = 0; i < NPOINTS; i ++) {

    x[i] = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;
    y[i] = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;

    CU_ASSERT(quadtree_add(d, x[i], y[i], &bound) >= 0);
  }

  CU_ASSERT(quadtree_npoints(d) == (NPOINTS + 4));

  /*
   * Check that the nearest point is the correct one
   */
  for (i = 0; i < NPOINTS; i ++) {

    j = quadtree_nearest(d, 0, x[i] + 0.0001, y[i] - 0.0001);

    CU_ASSERT((j - 4) == i);
  }

  quadtree_destroy(d);
}


void test_quadtree_random_delete(void)
{
  #define NPOINTS 100
  #define DELETEPOINTS 10

  quadtree_t *d;
  int i;
  int j;
  double x[NPOINTS];
  double y[NPOINTS];
  bbox2d_t bound;

  d = quadtree_create(100, 
		      DEFAULT_LEVELS,
		      -1.0,
		      1.0,
		      -1.0,
		      1.0);
  CU_ASSERT_FATAL(d != NULL);
  
  /*
   * Add all the random points
   */
  for (i = 0; i < NPOINTS; i ++) {

    x[i] = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;
    y[i] = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;

    CU_ASSERT(quadtree_add(d, x[i], y[i], &bound) >= 0);
  }

  CU_ASSERT(quadtree_npoints(d) == (NPOINTS + 4));

  /*
   * Check that the nearest point is the correct one
   */
  for (i = 0; i < NPOINTS; i ++) {

    j = quadtree_nearest(d, 0, x[i] + 0.001, y[i] - 0.001);

    CU_ASSERT((j - 4) == i);
  }

  /*
   * Delete the first n points
   */
  for (i = 0; i < DELETEPOINTS; i ++) {
    CU_ASSERT(quadtree_delete(d, 4, &bound) >= 0);
  }

  CU_ASSERT(quadtree_npoints(d) == (NPOINTS - DELETEPOINTS + 4));

  /*
   * Check that the nearest point is the correct one
   */
  for (i = DELETEPOINTS; i < NPOINTS; i ++) {

    j = quadtree_nearest(d, 0, x[i], y[i]);

    CU_ASSERT((j - 4) == (i - DELETEPOINTS));
  }

  quadtree_destroy(d);
}
