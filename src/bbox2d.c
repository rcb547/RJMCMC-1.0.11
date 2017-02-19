
#include <rjmcmc/bbox2d.h>

void
bbox2d_set(bbox2d_t *bound,
	   double xmin,
	   double xmax,
	   double ymin,
	   double ymax)
{
  bound->xmin = xmin;
  bound->xmax = xmax;
  bound->ymin = ymin;
  bound->ymax = ymax;
}

void 
bbox2d_initialize(bbox2d_t *bound,
		  double x,
		  double y)
{
  bound->xmin = x;
  bound->xmax = x;

  bound->ymin = y;
  bound->ymax = y;
}

void 
bbox2d_expand(bbox2d_t *bound,
	      double x,
	      double y)
{
  if (x < bound->xmin) {
    bound->xmin = x;
  }

  if (x > bound->xmax) {
    bound->xmax = x;
  }

  if (y < bound->ymin) {
    bound->ymin = y;
  }

  if (y > bound->ymax) {
    bound->ymax = y;
  }
}

void
bbox2d_bound_expand(bbox2d_t *bound,
		    const bbox2d_t *exp_bound)
{
  if (exp_bound->xmin < bound->xmin) {
    bound->xmin = exp_bound->xmin;
  }

  if (exp_bound->xmax > bound->xmax) {
    bound->xmax = exp_bound->xmax;
  }

  if (exp_bound->ymin < bound->ymin) {
    bound->ymin = exp_bound->ymin;
  }

  if (exp_bound->ymax > bound->ymax) {
    bound->ymax = exp_bound->ymax;
  }
}
