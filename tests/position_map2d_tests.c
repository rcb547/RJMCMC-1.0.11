
#include <stdio.h>

#include <CUnit/Basic.h>

#include <rjmcmc/position_map2d.h>
#include <rjmcmc/rjmcmc_random.h>

int approx_equal(double a, double b)
{
  static const double epsilon = 1e-5;

  if (fabs(a - b) > epsilon) {
    printf("approx_equal: %f != %f (e = %g)\n", a, b, fabs(a - b));
    return 0;
  } 

  return -1;
}

int init_position_map_suite(void)
{
  return 0;
}

int clean_position_map_suite(void)
{
  return 0;
}

void test_position_map_insert(void);
void test_position_map_delete(void);
void test_position_map_move(void);

void test_position_map_nearest(void);

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

  pSuite = CU_add_suite("Position Map 2D Suite", init_position_map_suite, clean_position_map_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Insert", test_position_map_insert);
  ADD_TEST("Delete", test_position_map_delete);
  ADD_TEST("Move", test_position_map_move);
  ADD_TEST("Nearest", test_position_map_nearest);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

double one(void)
{
  return 1.0;
}

void test_position_map_insert(void)
{
  double x[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  double y[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  int npoints = 5;
  bbox2d_t bound;

  int i;

  position_map2d_t *p;
  double nx;
  double ny;
  int ni;

  p = position_map2d_create(2*npoints, -1.0, 1.0, -1.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map2d_insert(p, x[i], y[i], &bound) >= 0);
  }

  CU_ASSERT_FATAL(position_map2d_validate(p) == 0);

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map2d_position_of_index(p, i + 4, &nx, &ny) >= 0);
    CU_ASSERT(nx == x[i]);
    CU_ASSERT(ny == y[i]);
  }

  for (i = 0; i < npoints; i ++) {
    ni = position_map2d_nearest(p, 
				x[i] + 0.025, 
				y[i] + 0.025, 
				&nx, &ny,
				0);

    /* printf("%d %d : %f %f = %f %f\n", */
    /* 	   i, ni, */
    /* 	   nx, ny, */
    /* 	   x[i], y[i]); */

    CU_ASSERT_FATAL(ni >= 0);
    CU_ASSERT(nx == x[i]);
    CU_ASSERT(ny == y[i]);
  }
}

void test_position_map_delete(void)
{
  double x[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  double y[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  int npoints = 5;

  int i;

  position_map2d_t *p;
  bbox2d_t bound;
  double nx;
  double ny;

  int del_i = 2;

  p = position_map2d_create(2*npoints, -1.0, 1.0, -1.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map2d_insert(p, x[i], y[i], &bound) >= 0);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map2d_position_of_index(p, i + 4, &nx, &ny) >= 0);
    CU_ASSERT(nx == x[i]);
    CU_ASSERT(ny == y[i]);
  }

  CU_ASSERT_FATAL(position_map2d_delete(p,
					del_i + 4,
					&bound) >= 0);

  CU_ASSERT_FATAL(position_map2d_position_of_index(p, del_i + 4, &nx, &ny) >= 0);
  CU_ASSERT(nx != x[del_i]);
  CU_ASSERT(ny != y[del_i]);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map2d_nearest(p, x[i] + 0.025, y[i] + 0.025, 
					   &nx, &ny, 0) >= 0);

    if (i == del_i) {
      CU_ASSERT(nx == x[3]);
      CU_ASSERT(ny == y[3]);
    } else {
      CU_ASSERT(nx == x[i]);
      CU_ASSERT(ny == y[i]);
    }
  }
}

void test_position_map_move(void)
{
  double x[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  double y[5] = {0.1, 0.9, 0.5, 0.6, 0.3};
  int npoints = 5;

  int i;

  position_map2d_t *p;
  bbox2d_t bound;
  double nx;
  double ny;

  int move_i = 2;
  double mx = 0.7;
  double my = 0.7;

  p = position_map2d_create(2*npoints, -1.0, 1.0, -1.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);
  
  for (i = 0; i < npoints; i ++) {
    CU_ASSERT(position_map2d_insert(p, x[i], y[i], &bound) >= 0);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map2d_position_of_index(p, i + 4, &nx, &ny) >= 0);
    CU_ASSERT(nx == x[i]);
    CU_ASSERT(ny == y[i]);
  }

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map2d_nearest(p, x[i] + 0.025, y[i] + 0.025, 
					   &nx, &ny, 0) >= 0);
    CU_ASSERT(nx == x[i]);
    CU_ASSERT(ny == y[i]);
  }

  CU_ASSERT_FATAL(position_map2d_move(p,
				      move_i + 4,
				      mx, my,
				      &bound) >= 0);

  CU_ASSERT_FATAL(position_map2d_nearest(p, mx + 0.025, my + 0.025,
					 &nx, &ny, 0) >= 0);
  CU_ASSERT(nx == mx);
  CU_ASSERT(ny == my);

  for (i = 0; i < npoints; i ++) {
    CU_ASSERT_FATAL(position_map2d_nearest(p, x[i] + 0.025, y[i] + 0.025, 
					   &nx, &ny, 0) >= 0);
    if (i == move_i) {
      CU_ASSERT(nx == x[3]);
      CU_ASSERT(ny == y[3]);
    } else {
      CU_ASSERT(nx == x[i]);
      CU_ASSERT(ny == y[i]);
    }
  }
}

double dist2(double x1, double y1, double x2, double y2)
{
  double dx;
  double dy;

  dx = x2 - x1;
  dy = y2 - y1;

  return dx*dx + dy*dy;
}

int linear_nearest(const double *px,
		   const double *py,
		   int npoints,
		   double x,
		   double y)
{
  int di;
  double d;
  double mind;
  int i;

  di = 0;
  mind = dist2(px[0], py[0], x, y);

  for (i = 1; i < npoints; i ++) {

    d = dist2(px[i], py[i], x, y);
    if (d < mind) {
      mind = d;
      di = i;
    }
  }

  return di;
}

void test_position_map_nearest(void)
{
  #define NPOINTS 100
  #define SAMPLES 10000

  position_map2d_t *p;
  bbox2d_t bound;
  double px[NPOINTS];
  double py[NPOINTS];

  double x;
  double y;
  double nx;
  double ny;

  int li;
  int pi;

  int i;

  p = position_map2d_create(NPOINTS + 4, -1.0, 1.0, -1.0, 1.0);
  CU_ASSERT_FATAL(p != NULL);

  /*
   * Create points
   */
  /* printf("\n"); */
  for (i = 0; i < NPOINTS; i ++) {

    px[i] = rjmcmc_uniform() * 2.0 - 1.0;
    py[i] = rjmcmc_uniform() * 2.0 - 1.0;

    /* printf("%d %f %f\n", i, px[i], py[i]); */
    CU_ASSERT_FATAL(position_map2d_insert(p, px[i], py[i], &bound) >= 0);
  }
  
  CU_ASSERT_FATAL(position_map2d_validate(p) == 0);

  /*
   * 
   */
  for (i = 0; i < SAMPLES; i ++) {
    x = rjmcmc_uniform() * 2.0 - 1.0;
    y = rjmcmc_uniform() * 2.0 - 1.0;

    li = linear_nearest(px, py, NPOINTS, x, y);

    pi = position_map2d_nearest(p,
				x,
				y,
				&nx,
				&ny,
				0);
    if (li != (pi - 4)) {
      printf("%d %d (%f %f) (%f %f (%f)) (%f %f (%f))\n", 
	     li, pi, 
	     x, y, 
	     nx, ny, dist2(x, y, nx, ny),
	     px[li], py[li], dist2(x, y, px[li], py[li]));
    }

    CU_ASSERT_FATAL(pi >= 0);

    CU_ASSERT(li == (pi - 4));

    CU_ASSERT(nx == px[li]);
    CU_ASSERT(ny == py[li]);
  }
    


}
