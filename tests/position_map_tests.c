
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/position_map1d.h>

int approx_equal(double a, double b)
{
  static const double epsilon = 1e-5;

  if (fabs(a - b) > epsilon) {
    printf("approx_equal: %f != %f (e = %g)\n", a, b, fabs(a - b));
    return 0;
  } 

  return -1;
}

int init_position_map1d_suite(void)
{
  return 0;
}

int clean_position_map1d_suite(void)
{
  return 0;
}

void test_position_map1d_insert(void);
void test_position_map1d_small_move(void);
void test_position_map1d_intervals(void);
void test_position_map1d_fill(void);
void test_position_map1d_delete(void);


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

  pSuite = CU_add_suite("Position_Map1d Suite", init_position_map1d_suite, clean_position_map1d_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Insert", test_position_map1d_insert);
  ADD_TEST("Small Move", test_position_map1d_small_move);
  ADD_TEST("Intervals", test_position_map1d_intervals);
  ADD_TEST("Fill", test_position_map1d_fill);
  ADD_TEST("Delete", test_position_map1d_delete);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

double one(void)
{
  return 1.0;
}

void test_position_map1d_insert(void)
{
  double x[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  int npoints = 5;

  int i;

  position_map1d_t *p;
  double nx;

  p = position_map1d_create(2*npoints, 0.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_insert(p, x[i], i + 2) >= 0);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_position_of_index(p, i + 2) == x[i]);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map1d_nearest(p, x[i] + 0.025, &nx) >= 0);
    CU_ASSERT(nx == x[i]);
  }

  CU_ASSERT(position_map1d_predecessor_of_index(p, 2) == 0);
  CU_ASSERT(position_map1d_predecessor_of_index(p, 3) == 5);
  CU_ASSERT(position_map1d_predecessor_of_index(p, 4) == 6);
  CU_ASSERT(position_map1d_predecessor_of_index(p, 5) == 4);
  CU_ASSERT(position_map1d_predecessor_of_index(p, 6) == 2);

  CU_ASSERT(position_map1d_successor_of_index(p, 2) == 6);
  CU_ASSERT(position_map1d_successor_of_index(p, 3) == 1);
  CU_ASSERT(position_map1d_successor_of_index(p, 4) == 5);
  CU_ASSERT(position_map1d_successor_of_index(p, 5) == 3);
  CU_ASSERT(position_map1d_successor_of_index(p, 6) == 4);
}


void test_position_map1d_small_move(void)
{
  double x[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  int npoints = 5;

  int i;

  position_map1d_t *p;
  double nx;

  p = position_map1d_create(2*npoints, 0.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_insert(p, x[i], i + 2) >= 0);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_position_of_index(p, i + 2) == x[i]);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map1d_nearest(p, x[i] + 0.025, &nx) >= 0);
    CU_ASSERT(nx == x[i]);
  }

  /* CU_ASSERT_FATAL(position_map1d_small_move(p, x[2], x[2] + 0.1) < 0); */
  /* CU_ASSERT_FATAL(position_map1d_small_move(p, x[2], x[2] - 0.2) < 0); */

  /* CU_ASSERT_FATAL(position_map1d_small_move(p, x[2], x[2] + 0.11) < 0); */
  /* CU_ASSERT_FATAL(position_map1d_small_move(p, x[2], x[2] - 0.21) < 0); */

  /* Move the 0.5 point close to 0.6 */
  CU_ASSERT_FATAL(position_map1d_small_move(p, x[2], x[2] + 0.09) >= 0);
}

struct cb_data {
  int ii;
  
  double xmin[10];
  double xmax[10];
  int iy[10];
};
  
static int icb(void *user_arg, 
	       double xmin,
	       double xmax,
	       int iy,
	       int riy)
{
  struct cb_data *c = (struct cb_data*)user_arg;

  c->xmin[c->ii] = xmin;
  c->xmax[c->ii] = xmax;
  c->iy[c->ii] = iy;

  c->ii ++;

  return 0;
}

void test_position_map1d_intervals(void)
{
  double x[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  int npoints = 5;

  struct cb_data c;

  int i;

  position_map1d_t *p;
  double nx;

  p = position_map1d_create(2*npoints, 0.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_insert(p, x[i], i + 2) >= 0);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_position_of_index(p, i + 2) == x[i]);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map1d_nearest(p, x[i] + 0.025, &nx) >= 0);
    CU_ASSERT(nx == x[i]);
  }

  c.ii = 0;
  CU_ASSERT_FATAL(position_map1d_traverse_intervals(p, icb, &c) >= 0);

  CU_ASSERT(c.ii == 6);

  /* for (i = 0; i < 5; i ++) { */
  /*   printf("%4d %12.6f %12.6f\n", c.iy[i], c.xmin[i], c.xmax[i]); */
  /* } */

  CU_ASSERT(c.iy[0] == 0);
  CU_ASSERT(c.xmin[0] == 0.0);
  CU_ASSERT(c.xmax[0] == 0.1);

  CU_ASSERT(c.iy[1] == 2);
  CU_ASSERT(c.xmin[1] == 0.1);
  CU_ASSERT(c.xmax[1] == 0.3);

  CU_ASSERT(c.iy[2] == 6);
  CU_ASSERT(c.xmin[2] == 0.3);
  CU_ASSERT(c.xmax[2] == 0.5);

  CU_ASSERT(c.iy[3] == 4);
  CU_ASSERT(c.xmin[3] == 0.5);
  CU_ASSERT(c.xmax[3] == 0.6);

  CU_ASSERT(c.iy[4] == 5);
  CU_ASSERT(c.xmin[4] == 0.6);
  CU_ASSERT(c.xmax[4] == 0.9);

  CU_ASSERT(c.iy[5] == 3);
  CU_ASSERT(c.xmin[5] == 0.9);
  CU_ASSERT(c.xmax[5] == 1.0);
}

void test_position_map1d_fill(void)
{
  double x[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  
  int npoints = 5;

  struct cb_data c;

  int i;

  position_map1d_t *p;
  double nx;

  double xf[10];
  int fillpoints;

  p = position_map1d_create(2*npoints, 0.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_insert(p, x[i], i + 2) >= 0);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_position_of_index(p, i + 2) == x[i]);
  }

  fillpoints = 10;
  CU_ASSERT_FATAL(position_map1d_fill_list(p, xf, &fillpoints) >= 0);

  CU_ASSERT(fillpoints == 7);

  CU_ASSERT(xf[0] == 0.0);
  CU_ASSERT(xf[1] == 0.1);
  CU_ASSERT(xf[2] == 0.3);
  CU_ASSERT(xf[3] == 0.5);
  CU_ASSERT(xf[4] == 0.6);
  CU_ASSERT(xf[5] == 0.9);
  CU_ASSERT(xf[6] == 1.0);
}

void test_position_map1d_delete(void)
{
  double x[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  
  int npoints = 5;

  struct cb_data c;

  int i;

  position_map1d_t *p;
  double nx;

  double xf[10];
  int fillpoints;

  p = position_map1d_create(2*npoints, 0.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_insert(p, x[i], i + 2) >= 0);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map1d_position_of_index(p, i + 2) == x[i]);
  }

  CU_ASSERT_FATAL(position_map1d_delete(p, 0.5, 4) >= 0);

  fillpoints = 10;
  CU_ASSERT_FATAL(position_map1d_fill_list(p, xf, &fillpoints) >= 0);

  CU_ASSERT(fillpoints == 6);

  CU_ASSERT(xf[0] == 0.0);
  CU_ASSERT(xf[1] == 0.1);
  CU_ASSERT(xf[2] == 0.3);
  CU_ASSERT(xf[3] == 0.6);
  CU_ASSERT(xf[4] == 0.9);
  CU_ASSERT(xf[5] == 1.0);

  CU_ASSERT(position_map1d_position_of_index(p, 2) == 0.1);
  CU_ASSERT(position_map1d_position_of_index(p, 3) == 0.9);
  CU_ASSERT(position_map1d_position_of_index(p, 4) == 0.6);
  CU_ASSERT(position_map1d_position_of_index(p, 5) == 0.3);

}
