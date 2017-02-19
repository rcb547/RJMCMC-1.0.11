
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <rjmcmc/position_map2d.h>
#include <rjmcmc/position_map2d_linear.h>
#include <rjmcmc/position_map2d_delaunay.h>
#include <rjmcmc/position_map2d_quadtree.h>

#include "rjmcmc/rjmcmc_util.h"
#include "rjmcmc/rjmcmc_defines.h"

//struct _position_map2d {
//};

static struct position_map2d_map {
  
  position_map2d_t *(*create)(int, double, double, double, double);
  void (*destroy)(position_map2d_t *);
  void (*clone)(const position_map2d_t *, position_map2d_t *);
  int (*insert)(position_map2d_t *, double, double, bbox2d_t *);
  int (*delete)(position_map2d_t *, int, bbox2d_t *);
  int (*move)(position_map2d_t *, int, double, double, bbox2d_t *);
  int (*nearest)(position_map2d_t *, double, double, double*, double*, int);
  int (*enclosing_triangle)(position_map2d_t *, double, double, int*, int*, int*, double*, double*, double*);
  int (*position_of_index)(position_map2d_t *, int, double*, double*);
  int (*polygon_bound)(position_map2d_t*, int, bbox2d_t *);
  int (*validate)(position_map2d_t*);

} pmap[3] = {
  { position_map2d_linear_create,
    position_map2d_linear_destroy,
    position_map2d_linear_clone,
    position_map2d_linear_insert,
    position_map2d_linear_delete,
    position_map2d_linear_move,
    position_map2d_linear_nearest,
    position_map2d_linear_enclosing_triangle,
    position_map2d_linear_position_of_index,
    position_map2d_linear_polygon_bound,
    position_map2d_linear_validate },

  { position_map2d_delaunay_create,
    position_map2d_delaunay_destroy,
    position_map2d_delaunay_clone,
    position_map2d_delaunay_insert,
    position_map2d_delaunay_delete,
    position_map2d_delaunay_move,
    position_map2d_delaunay_nearest,
    position_map2d_delaunay_enclosing_triangle,
    position_map2d_delaunay_position_of_index,
    position_map2d_delaunay_polygon_bound,
    position_map2d_delaunay_validate },

  { position_map2d_quadtree_create,
    position_map2d_quadtree_destroy,
    position_map2d_quadtree_clone,
    position_map2d_quadtree_insert,
    position_map2d_quadtree_delete,
    position_map2d_quadtree_move,
    position_map2d_quadtree_nearest,
    position_map2d_quadtree_enclosing_triangle,
    position_map2d_quadtree_position_of_index,
    position_map2d_quadtree_polygon_bound,
    position_map2d_quadtree_validate }

};

static int map2d_method = 0;

void
position_map2d_set_type(position_map2d_type_t t)
{
  map2d_method = (int)t;
  if (map2d_method < 0 || 
      map2d_method > 2) {
    map2d_method = 0;
  }
}


position_map2d_t *
position_map2d_create(int max_partitions,
		      double xmin,
		      double xmax,
		      double ymin,
		      double ymax)
{
  return pmap[map2d_method].create(max_partitions, xmin, xmax, ymin, ymax);
}

void
position_map2d_destroy(position_map2d_t *p)
{
  pmap[map2d_method].destroy(p);
}

void
position_map2d_clone(const position_map2d_t *src,
		     position_map2d_t *dst)
{
  pmap[map2d_method].clone(src, dst);
}

int
position_map2d_insert(position_map2d_t *p,
		      double x, double y,
		      bbox2d_t *bound)
{
  return pmap[map2d_method].insert(p, x, y, bound);
}

int 
position_map2d_delete(position_map2d_t *p,
		      int iy,
		      bbox2d_t *bound)
{
  return pmap[map2d_method].delete(p, iy, bound);
}

int
position_map2d_move(position_map2d_t *p,
		    int iy,
		    double new_x, double new_y,
		    bbox2d_t *bound)
{
  return pmap[map2d_method].move(p, iy, new_x, new_y, bound);
}

int 
position_map2d_nearest(position_map2d_t *p,
		       double x,
		       double y,
		       double *nx,
		       double *ny,
		       int include_boundary_points)
{
  return pmap[map2d_method].nearest(p, x, y, nx, ny, include_boundary_points);
}

int
position_map2d_enclosing_triangle(position_map2d_t *p,
				  double x,
				  double y,
				  int *ta,
				  int *tb,
				  int *tc,
				  double *ba,
				  double *bb,
				  double *bc)
{
  return pmap[map2d_method].enclosing_triangle(p, x, y, ta, tb, tc, ba, bb, bc);
}

int 
position_map2d_position_of_index(position_map2d_t *p,
				 int iy,
				 double *x,
				 double *y)
{
  return pmap[map2d_method].position_of_index(p, iy, x, y);
}

int
position_map2d_polygon_bound(position_map2d_t *p,
			     int iy,
			     bbox2d_t *bound)
{
  return pmap[map2d_method].polygon_bound(p, iy, bound);
}

int
position_map2d_validate(position_map2d_t *p)
{
  return pmap[map2d_method].validate(p);
}

		     

