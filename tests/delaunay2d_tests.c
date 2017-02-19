
#include <stdio.h>
#include <stdlib.h>

#include <CUnit/Basic.h>

#include <rjmcmc/delaunay2d.h>

int init_util_suite(void)
{
  return 0;
}

int clean_util_suite(void)
{
  return 0;
}

void test_circumcircle(void);
void test_point_in_triangle(void);

void test_delaunay2d_create(void);
void test_delaunay2d_basic_find(void);
void test_delaunay2d_basic_add(void);
void test_delaunay2d_edge_add_01(void);
void test_delaunay2d_edge_add_02(void);
void test_delaunay2d_random_add(void);
void test_delaunay2d_stress_add(void);

void test_delaunay2d_save_load(void);

void test_delaunay2d_delete_basic(void);
void test_delaunay2d_delete_internal(void);
void test_delaunay2d_delete_internal_concave(void);

void test_delaunay2d_stress_add_delete(void);
void test_delaunay2d_stress_add_delete_nve(void);

void test_delaunay2d_case_01(void);
void test_delaunay2d_case_02(void);
void test_delaunay2d_case_03(void);
void test_delaunay2d_case_04(void);
void test_delaunay2d_case_05(void);
void test_delaunay2d_case_06(void);

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

  pSuite = CU_add_suite("Delaunay2d Suite", init_util_suite, clean_util_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  ADD_TEST("Triangle Circumcircle", test_circumcircle);
  ADD_TEST("Triangle Point", test_point_in_triangle);

  ADD_TEST("Delaunay2d Create", test_delaunay2d_create);
  ADD_TEST("Delaunay2d Basic Find", test_delaunay2d_basic_find);
  ADD_TEST("Delaunay2d Basic Add", test_delaunay2d_basic_add);
  ADD_TEST("Delaunay2d Edge Add 1", test_delaunay2d_edge_add_01);
  ADD_TEST("Delaunay2d Edge Add 2", test_delaunay2d_edge_add_02);
  ADD_TEST("Delaunay2d Random Add", test_delaunay2d_random_add);
  ADD_TEST("Delaunay2d Stress Add", test_delaunay2d_stress_add);

  ADD_TEST("Delaunay2d Save/Load", test_delaunay2d_save_load);

  ADD_TEST("Delaunay2d Delete Basic", test_delaunay2d_delete_basic);
  ADD_TEST("Delaunay2d Delete Internal", test_delaunay2d_delete_internal);
  ADD_TEST("Delaunay2d Delete Internal Concave", test_delaunay2d_delete_internal_concave);

  ADD_TEST("Delaunay2d Stress Add/Delete", test_delaunay2d_stress_add_delete);
  ADD_TEST("Delaunay2d Stress Add/Delete Negative", test_delaunay2d_stress_add_delete_nve);

  ADD_TEST("Test Case 01", test_delaunay2d_case_01);
  ADD_TEST("Test Case 02", test_delaunay2d_case_02);
  ADD_TEST("Test Case 03", test_delaunay2d_case_03);
  ADD_TEST("Test Case 04", test_delaunay2d_case_04);
  ADD_TEST("Test Case 05", test_delaunay2d_case_05);
  ADD_TEST("Test Case 06", test_delaunay2d_case_06);

  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}

void test_circumcircle(void)
{
  double x1, y1;
  double x2, y2;
  double x3, y3;

  double cx;
  double cy;
  double r2;

  double dx;
  double dy;
  double pr2;

  x1 = 0.213938;
  y1 = -0.967399;
  x2 = 1.0;
  y2 = 1.0;
  x3 = 1.0;
  y3 = -1.0;
    
  CU_ASSERT_FATAL(triangle_circumcircle(x1, y1,
					x2, y2,
					x3, y3,
					&cx, &cy, &r2) == 0);

  dx = x1 - cx;
  dy = y1 - cy;
  pr2 = dx*dx + dy*dy;
  CU_ASSERT(fabs(pr2 - r2) < 1.0e-6);

  dx = x2 - cx;
  dy = y2 - cy;
  pr2 = dx*dx + dy*dy;
  CU_ASSERT(fabs(pr2 - r2) < 1.0e-6);

  dx = x3 - cx;
  dy = y3 - cy;
  pr2 = dx*dx + dy*dy;
  CU_ASSERT(fabs(pr2 - r2) < 1.0e-6);

  x1 = 1.0;
  y1 = 1.0;
  x2 = -1.0;
  y2 = -1.0;
  
}

void test_point_in_triangle(void)
{
  double x1, y1;
  double x2, y2;
  double x3, y3;

  x1 = -0.3;
  y1 = 0.5;

  x2 = 0.4;
  y2 = 0.6;

  x3 = 0.1;
  y3 = -0.7;

  CU_ASSERT(point_in_triangle(0.0, 0.0, x1, y1, x2, y2, x3, y3) < 0);
  CU_ASSERT(point_in_triangle(-0.5, -0.1, x1, y1, x2, y2, x3, y3) == 0);
  CU_ASSERT(point_in_triangle(0.5, 0.1, x1, y1, x2, y2, x3, y3) == 0);
  CU_ASSERT(point_in_triangle(0.1, 0.7, x1, y1, x2, y2, x3, y3) == 0);

}

void test_delaunay2d_create(void)
{
  delaunay2d_t *d;

  d = delaunay2d_create(10, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
  CU_ASSERT(delaunay2d_validate_delaunay(d) == 0);
  
  delaunay2d_destroy(d);

  d = delaunay2d_create(10, 
			-2.0,
			2.0,
			-1.0,
			1.0);
  CU_ASSERT(d != NULL);

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

}

void test_delaunay2d_basic_find(void)
{
  delaunay2d_t *d;

  int pa;
  int pb;
  int pc;

  double ba;
  double bb;
  double bc;

  d = delaunay2d_create(10, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_find_enclosing_triangle(d,
						     0,
						     -0.5,
						     0.5,
						     &pa,
						     &pb, 
						     &pc,
						     &ba,
						     &bb,
						     &bc) == 0);
  
  CU_ASSERT(pa == 0);
  CU_ASSERT(pb == 1);
  CU_ASSERT(pc == 2);

  CU_ASSERT_FATAL(delaunay2d_find_enclosing_triangle(d,
						     1,
						     -0.5,
						     0.5,
						     &pa,
						     &pb, 
						     &pc,
						     &ba,
						     &bb,
						     &bc) == 0);
  
  CU_ASSERT(pa == 0);
  CU_ASSERT(pb == 1);
  CU_ASSERT(pc == 2);

  CU_ASSERT_FATAL(delaunay2d_find_enclosing_triangle(d,
						     0,
						     0.5,
						     -0.5,
						     &pa,
						     &pb, 
						     &pc,
						     &ba,
						     &bb,
						     &bc) == 1);
  
  CU_ASSERT(pa == 0);
  CU_ASSERT(pb == 2);
  CU_ASSERT(pc == 3);

  CU_ASSERT_FATAL(delaunay2d_find_enclosing_triangle(d,
						     1,
						     0.5,
						     -0.5,
						     &pa,
						     &pb, 
						     &pc,
						     &ba,
						     &bb,
						     &bc) == 1);
  
  CU_ASSERT(pa == 0);
  CU_ASSERT(pb == 2);
  CU_ASSERT(pc == 3);

}

void test_delaunay2d_basic_add(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_create(10, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  /* delaunay2d_print_points(d); */
  /* delaunay2d_print_triangles(d); */
  
  CU_ASSERT_FATAL(delaunay2d_add(d, -0.5, 0.5, &bound) == 0);

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);


  CU_ASSERT(delaunay2d_validate_edges(d) == 0);

  printf("\n");
  delaunay2d_print_points(d);
  delaunay2d_print_triangles(d);
  delaunay2d_print_edges(d);
}

void test_delaunay2d_edge_add_01(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_create(10, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
  
  /* Edge add with no neighbour */
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.0, 1.0, &bound) == 0);
  
  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);
}

void test_delaunay2d_edge_add_02(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_create(10, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  /* Edge add with neighbour */
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.1, 0.1, &bound) == 0);

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);
}

void test_delaunay2d_random_add(void)
{
  delaunay2d_t *d;
  int i;
  double x;
  double y;
  bbox2d_t bound;

  d = delaunay2d_create(100, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  for (i = 0; i < 10; i ++) {

    x = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;
    y = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;

    CU_ASSERT(delaunay2d_add(d, x, y, &bound) == 0);

    CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

    CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);
  }
}

void test_delaunay2d_stress_add(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;
  int j;
  int i;
  double x;
  double y;

  for (j = 0; j < 100; j ++) {

    d = delaunay2d_create(100, 
			  -1.0,
			  1.0,
			  -1.0,
			  1.0);
    CU_ASSERT_FATAL(d != NULL);
    
    CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
    
    for (i = 0; i < 10; i ++) {

      x = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;
      y = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;
      
      CU_ASSERT(delaunay2d_add(d, x, y, &bound) == 0);
      
      CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
    }

    delaunay2d_destroy(d);
  }
}

void test_delaunay2d_save_load(void)
{
  delaunay2d_t *d;
  delaunay2d_t *dl;
  bbox2d_t bound;

  d = delaunay2d_create(15, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_add(d, -0.1, 0.3, &bound) == 0);   // 4
  CU_ASSERT_FATAL(delaunay2d_add(d, -0.6, 0.5, &bound) == 0);   // 5
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.5, 0.55, &bound) == 0);   // 6
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.3, 0.3, &bound) == 0);    // 7
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.1, -0.2, &bound) == 0);   // 8
  CU_ASSERT_FATAL(delaunay2d_add(d, -0.2, 0.0, &bound) == 0);   // 9
  CU_ASSERT_FATAL(delaunay2d_add(d, -0.4, -0.4, &bound) == 0);  // 10

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_save(d, "test_delaunay2d_save_load.d") == 0);
  
  dl = delaunay2d_load("test_delaunay2d_save_load.d");
  CU_ASSERT_FATAL(dl != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(dl) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_neighbours(dl) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(dl) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(dl) == 0);
}

void test_delaunay2d_delete_basic(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_create(10, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_add(d, -0.5, 0.5, &bound) == 0);

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);

  CU_ASSERT(delaunay2d_delete(d, 4, &bound) == 0);

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);

  CU_ASSERT(delaunay2d_validate_edges(d) == 0);
}

void test_delaunay2d_delete_internal(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_create(10, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_add(d, -0.1, 0.0, &bound) == 0); // Central point (to be deleted)

  CU_ASSERT_FATAL(delaunay2d_add(d, -0.5, 0.6, &bound) == 0);
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.6, 0.4, &bound) == 0);
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.5, -0.3, &bound) == 0);
  CU_ASSERT_FATAL(delaunay2d_add(d, -0.4, -0.5, &bound) == 0);

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);

  CU_ASSERT(delaunay2d_validate_edges(d) == 0);

  /* fprintf(stderr, "Pre delete\n"); */
  /* delaunay2d_print_points(d); */
  /* delaunay2d_print_triangles(d); */
  /* fprintf(stderr, "\n"); */

  CU_ASSERT(delaunay2d_delete(d, 4, &bound) == 0);

  /* printf("%f %f %f %f\n", bound.xmin, bound.xmax, bound.ymin, bound.ymax); */

  CU_ASSERT(bound.xmin == -0.5);
  CU_ASSERT(bound.xmax == 0.6);
  CU_ASSERT(bound.ymin == -0.5);
  CU_ASSERT(bound.ymax == 0.6);

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);

  CU_ASSERT(delaunay2d_validate_edges(d) == 0);
}

void test_delaunay2d_delete_internal_concave(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_create(15, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_add(d, -0.1, 0.3, &bound) == 0);   // 4
  CU_ASSERT_FATAL(delaunay2d_add(d, -0.6, 0.5, &bound) == 0);   // 5
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.5, 0.55, &bound) == 0);   // 6
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.3, 0.3, &bound) == 0);    // 7
  CU_ASSERT_FATAL(delaunay2d_add(d, 0.1, -0.2, &bound) == 0);   // 8
  CU_ASSERT_FATAL(delaunay2d_add(d, -0.2, 0.0, &bound) == 0);   // 9
  CU_ASSERT_FATAL(delaunay2d_add(d, -0.4, -0.4, &bound) == 0);  // 10

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);

  CU_ASSERT(delaunay2d_validate_edges(d) == 0);

  /* fprintf(stderr, "Pre delete\n"); */
  /* delaunay2d_print_points(d); */
  /* delaunay2d_print_triangles(d); */
  /* fprintf(stderr, "\n"); */

  CU_ASSERT(delaunay2d_delete(d, 4, &bound) == 0);

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);

  CU_ASSERT(delaunay2d_validate_edges(d) == 0);

}

void test_delaunay2d_stress_add_delete_nve(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;
  int i;
  int pi;
  double x;
  double y;
  double u;
  
  double xmin, xmax;
  double ymin, ymax;

  int np;

  int maxpoints;

  char filename[1024];

  maxpoints = 200;

  ymin = 10.0;
  ymax = 20.0;

  xmin = -30.0;
  xmax = -20.0;

  d = delaunay2d_create(maxpoints, 
			xmin,
			xmax,
			ymin,
			ymax);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  srand(0);

  for (i = 0; i < 10000; i ++) {

    u = (double)rand()/(double)RAND_MAX;
    np = delaunay2d_npoints(d);
    
    if ((np == 4 || u < 0.5) && (np < maxpoints)) {

      /* Add */

      x = (double)rand()/(double)RAND_MAX * (xmax - xmin) + xmin;
      y = (double)rand()/(double)RAND_MAX * (ymax - ymin) + ymin;
      /* fprintf(stderr, "Adding: %f %f\n", x, y); */
      
      CU_ASSERT_FATAL(delaunay2d_add(d, x, y, &bound) == 0);
  
      CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
      
      CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_edges(d) == 0);

    } else {

      /* Delete */

      pi = 4 + (int)((double)(np - 4) * (double)rand()/(double)RAND_MAX);

      /* fprintf(stderr, "Deleting: %d\n", pi); */
      /* fprintf(stderr, "***\n"); */
      /* delaunay2d_print_points(d); */
      /* delaunay2d_print_triangles(d); */
      /* fprintf(stderr, "***\n"); */

      CU_ASSERT_FATAL(delaunay2d_delete(d, pi, &bound) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
      
      CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_edges(d) == 0);

      if (delaunay2d_validate_nonintersecting(d) != 0) {
	delaunay2d_save_geo(d, "intersecting.geo");
      }
      CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);
    }

    /* sprintf(filename, "stress_add_delete_%d.geo", i); */
    /* if (delaunay2d_save_geo(d, filename) < 0) { */
    /*   fprintf(stderr, "error: failed to save goe\n"); */
    /*   exit(-1); */
    /* } */
       
  }
}

void test_delaunay2d_stress_add_delete(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;
  int i;
  int pi;
  double x;
  double y;
  double u;
  
  int np;

  int maxpoints;

  char filename[1024];

  maxpoints = 200;
  d = delaunay2d_create(maxpoints, 
			-1.0,
			1.0,
			-1.0,
			1.0);
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);

  srand(0);

  for (i = 0; i < 10000; i ++) {

    u = (double)rand()/(double)RAND_MAX;
    np = delaunay2d_npoints(d);
    
    if ((np == 4 || u < 0.5) && (np < maxpoints)) {

      /* Add */

      x = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;
      y = (double)rand()/(double)RAND_MAX * 2.0 - 1.0;
      /* fprintf(stderr, "Adding: %f %f\n", x, y); */
      
      CU_ASSERT_FATAL(delaunay2d_add(d, x, y, &bound) == 0);
  
      CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
      
      CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_edges(d) == 0);

    } else {

      /* Delete */

      pi = 4 + (int)((double)(np - 4) * (double)rand()/(double)RAND_MAX);

      /* fprintf(stderr, "Deleting: %d\n", pi); */
      /* fprintf(stderr, "***\n"); */
      /* delaunay2d_print_points(d); */
      /* delaunay2d_print_triangles(d); */
      /* fprintf(stderr, "***\n"); */

      CU_ASSERT_FATAL(delaunay2d_delete(d, pi, &bound) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
      
      CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);

      CU_ASSERT_FATAL(delaunay2d_validate_edges(d) == 0);

      if (delaunay2d_validate_nonintersecting(d) != 0) {
	delaunay2d_save_geo(d, "intersecting.geo");
      }
      CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);
    }

    /* sprintf(filename, "stress_add_delete_%d.geo", i); */
    /* if (delaunay2d_save_geo(d, filename) < 0) { */
    /*   fprintf(stderr, "error: failed to save goe\n"); */
    /*   exit(-1); */
    /* } */
       
  }
}

void test_delaunay2d_case_01(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;
  
  d = delaunay2d_load("test_failure_01_pre.d");
  CU_ASSERT_FATAL(d != NULL);

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_delete(d, 20, &bound) == 0);

  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

}

void test_delaunay2d_case_02(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_load("test_failure_02_pre.d");
  CU_ASSERT_FATAL(d != NULL);

  delaunay2d_save_geo(d, "case02_pre.geo");

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_add(d, 0.384927, -0.086972, &bound) == 0);

  delaunay2d_save_geo(d, "case02_post.geo");

  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

}

void test_delaunay2d_case_03(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_load("test_failure_03_pre.d");
  CU_ASSERT_FATAL(d != NULL);

  delaunay2d_save_geo(d, "case03_pre.geo");

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_delete(d, 132, &bound) == 0);

  delaunay2d_save_geo(d, "case03_post.geo");

  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

}

void test_delaunay2d_case_04(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_load("test_failure_04_pre.d");
  CU_ASSERT_FATAL(d != NULL);

  delaunay2d_save_geo(d, "case04_pre.geo");

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_delete(d, 105, &bound) == 0);

  delaunay2d_save_geo(d, "case04_post.geo");

  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

}

void test_delaunay2d_case_05(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_load("test_failure_05_pre.d");
  CU_ASSERT_FATAL(d != NULL);

  delaunay2d_save_geo(d, "case05_pre.geo");

  CU_ASSERT_FATAL(delaunay2d_validate_circumcircles(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_neighbours(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_nonintersecting(d) == 0);
  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_delete(d, 5, &bound) == 0);

  delaunay2d_save_geo(d, "case05_post.geo");

  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

}

void test_delaunay2d_case_06(void)
{
  delaunay2d_t *d;
  bbox2d_t bound;

  d = delaunay2d_load("test_failure_06_pre.d");
  CU_ASSERT_FATAL(d != NULL);

  delaunay2d_save_geo(d, "case06_pre.geo");
  delaunay2d_save_cc_geo(d, "cast06_pre_cc.geo");

  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);
  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);
  CU_ASSERT(delaunay2d_validate_nonintersecting(d) == 0);
  CU_ASSERT(delaunay2d_validate_delaunay(d) == 0);

  CU_ASSERT_FATAL(delaunay2d_delete(d, 5, &bound) == 0);

  delaunay2d_save_geo(d, "case06_post.geo");
  delaunay2d_save_cc_geo(d, "cast06_post_cc.geo");

  printf("\nCircum\n");
  CU_ASSERT(delaunay2d_validate_circumcircles(d) == 0);

  printf("\nNeighbours\n");
  CU_ASSERT(delaunay2d_validate_neighbours(d) == 0);

  printf("\nNon-inter\n");
  CU_ASSERT(delaunay2d_validate_nonintersecting(d) == 0);

  printf("\nDelaunay\n");
  CU_ASSERT_FATAL(delaunay2d_validate_delaunay(d) == 0);

}
