#ifndef bbox2d_h
#define bbox2d_h

/** \file bbox2d.h

\brief 2D Bounding Box routines

2D bounding boxes are used in the constrained delaunay code as well as
for 2D forward model codes.

*/

struct bbox2d {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};
typedef struct bbox2d bbox2d_t;

/**

\brief Set a bounding box

Sets a bounding box to given values. No validity checking is performed
(ie xmin < xmax).
*/
void
bbox2d_set(bbox2d_t *bound,
	   double xmin,
	   double xmax,
	   double ymin,
	   double ymax);

/**

\brief Initialise a bounding box to an initial point (creating an infinitesimal bound).

Sets the bound to an infinitesimal point, ie xmin = xmax = x and ymin
= ymax = y.
*/
void bbox2d_initialize(bbox2d_t *bound, 
		       double x, 
		       double y);
/**

\brief Expand a bounding box to include a point

Increases an existing bounding box to include the specified point. If the
point x, y is already in the bound, the bounding box will not be changed.
*/
void 
bbox2d_expand(bbox2d_t *bound,
	      double x,
	      double y);

/**

\brief Expand a bounding box to include another bounding box

Increases an existing bounding box to include the specified bounding
box. If the existing bound completely contains the specified bound, the
existing bound is unchanged.
*/
void
bbox2d_bound_expand(bbox2d_t *bound,
		    const bbox2d_t *exp_bound);

#endif /* bbox2d_h */
