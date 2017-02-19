#ifndef position_map2d_quadtree_h
#define position_map2d_quadtree_h

#include <rjmcmc/bbox2d.h>

typedef struct _position_map2d_quadtree position_map2d_quadtree_t;

position_map2d_quadtree_t *
position_map2d_quadtree_create(int max_partitions,
			       double xmin,
			       double xmax,
			       double ymin,
			       double ymax);

void
position_map2d_quadtree_destroy(position_map2d_quadtree_t *p);

void
position_map2d_quadtree_clone(const position_map2d_quadtree_t *src,
			      position_map2d_quadtree_t *dst);

int
position_map2d_quadtree_insert(position_map2d_quadtree_t *p,
			       double x, double y,
			       bbox2d_t *bound);

int 
position_map2d_quadtree_delete(position_map2d_quadtree_t *p,
			       int iy,
			       bbox2d_t *bound);

int
position_map2d_quadtree_move(position_map2d_quadtree_t *p,
			     int iy,
			     double new_x, double new_y,
			     bbox2d_t *bound);

int 
position_map2d_quadtree_nearest(position_map2d_quadtree_t *p,
				double x,
				double y,
				double *nx,
				double *ny,
				int include_boundary_points);

int
position_map2d_quadtree_enclosing_triangle(position_map2d_quadtree_t *p,
					   double x,
					   double y,
					   int *ta,
					   int *tb,
					   int *tc,
					   double *ba,
					   double *bb,
					   double *bc);

int 
position_map2d_quadtree_position_of_index(position_map2d_quadtree_t *p,
					  int iy,
					  double *x,
					  double *y);

int
position_map2d_quadtree_polygon_bound(position_map2d_quadtree_t *p,
				      int iy,
				      bbox2d_t *bound);

int
position_map2d_quadtree_validate(position_map2d_quadtree_t *p);


#endif /* position_map2d_quadtree_h */
