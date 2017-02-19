#ifndef position_map2d_h
#define position_map2d_h

#include <rjmcmc/bbox2d.h>

typedef struct _position_map2d position_map2d_t;

typedef enum {
  POSITION_MAP2D_LINEAR = 0,
  POSITION_MAP2D_DELAUNAY,
  POSITION_MAP2D_QUADTREE
} position_map2d_type_t;

void
position_map2d_set_type(position_map2d_type_t t);
  
position_map2d_t *
position_map2d_create(int max_partitions,
		      double xmin,
		      double xmax,
		      double ymin,
		      double ymax);

void
position_map2d_destroy(position_map2d_t *p);

void
position_map2d_clone(const position_map2d_t *src,
		     position_map2d_t *dst);

int
position_map2d_insert(position_map2d_t *p,
		      double x, double y,
		      bbox2d_t *bound);

int 
position_map2d_delete(position_map2d_t *p,
		      int iy,
		      bbox2d_t *bound);

int
position_map2d_move(position_map2d_t *p,
		    int iy,
		    double new_x, double new_y,
		    bbox2d_t *bound);

int 
position_map2d_nearest(position_map2d_t *p,
		       double x,
		       double y,
		       double *nx,
		       double *ny,
		       int include_boundary_points);

int
position_map2d_enclosing_triangle(position_map2d_t *p,
				  double x,
				  double y,
				  int *ta,
				  int *tb,
				  int *tc,
				  double *ba,
				  double *bb,
				  double *bc);

int 
position_map2d_position_of_index(position_map2d_t *p,
				 int iy,
				 double *x,
				 double *y);

int
position_map2d_polygon_bound(position_map2d_t *p,
			     int iy,
			     bbox2d_t *bound);

int
position_map2d_validate(position_map2d_t *p);


#endif /* position_map2d_h */
