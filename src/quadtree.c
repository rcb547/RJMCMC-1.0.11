
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <rjmcmc/quadtree.h>

#include <rjmcmc/rjmcmc_util.h>

struct _point {
  double x;
  double y;
};
typedef struct _point point_t;

struct _quadtree_leaf {
  int np;
  int *point;
};
typedef struct _quadtree_leaf quadtree_leaf_t;

typedef struct _quadtree_node quadtree_node_t;
struct _quadtree_node {
  bbox2d_t bound;
  double cx;
  double cy;

  int node;
  int grand_children;
  union {
    quadtree_node_t *child[4];
    quadtree_leaf_t leaf;
  } data;
};

struct _quadtree {
  bbox2d_t bound;
  
  int maxpoints;
  int np;
  point_t *point;

  quadtree_node_t *head;
};

static quadtree_node_t *
quadtree_node_create(int maxpoints,
		     int maxlevels,
		     double xmin,
		     double xmax,
		     double ymin,
		     double ymax);

static void
quadtree_node_destroy(quadtree_node_t *q);

static int quadtree_node_insert(quadtree_t *q,
				quadtree_node_t *qn,
				double x,
				double y,
				int pi);

static int quadtree_node_delete(quadtree_t *q, 
				quadtree_node_t *qn, 
				int pi, 
				bbox2d_t *bound);

static int quadtree_node_nearest(quadtree_t *q, 
				 quadtree_node_t *qn, 
				 double px, 
				 double py, 
				 int include_boundary_points,
				 int *minpi,
				 double *mindist2);

static int quadtree_node_clone(const quadtree_t *src, 
			       const quadtree_node_t *sn, 
			       quadtree_t *dst, 
			       quadtree_node_t *dn);

static int quadtree_node_shift_replace(quadtree_t *q,
				       quadtree_node_t *qn,
				       int pi);

quadtree_t *
quadtree_create(int maxpoints,
		int maxlevels,
		double xmin,
		double xmax,
		double ymin,
		double ymax)
{
  quadtree_t *q;

  q = malloc(sizeof(quadtree_t));
  if (q == NULL) {
    return NULL;
  }

  q->bound.xmin = xmin;
  q->bound.xmax = xmax;
  q->bound.ymin = ymin;
  q->bound.ymax = ymax;

  q->point = malloc(sizeof(point_t) * (maxpoints + 4));
  if (q->point == NULL) {
    return NULL;
  }

  q->maxpoints = maxpoints + 4;
  q->np = 4;

  q->point[0].x = xmin;
  q->point[0].y = ymin;

  q->point[1].x = xmin;
  q->point[1].y = ymax;

  q->point[2].x = xmax;
  q->point[2].y = ymax;

  q->point[3].x = xmax;
  q->point[3].y = ymin;

  q->head = quadtree_node_create(maxpoints, maxlevels, xmin, xmax, ymin, ymax);
  if (q->head == NULL) {
    return NULL;
  }

  if (quadtree_node_insert(q, q->head, xmin, ymin, 0) < 0) {
    return NULL;
  }

  if (quadtree_node_insert(q, q->head, xmin, ymax, 1) < 0) {
    return NULL;
  }

  if (quadtree_node_insert(q, q->head, xmax, ymax, 2) < 0) {
    return NULL;
  }

  if (quadtree_node_insert(q, q->head, xmax, ymin, 3) < 0) {
    return NULL;
  }

  return q;
}

void 
quadtree_destroy(quadtree_t *q)
{
  if (q != NULL) {
    
    free(q->point);
    quadtree_node_destroy(q->head);

    free(q);
  }
}

int
quadtree_clone(const quadtree_t *src,
	       quadtree_t *dst)
{
  int i;

  dst->np = src->np;
  for (i = 0; i < src->np; i ++) {
    dst->point[i].x = src->point[i].x;
    dst->point[i].y = src->point[i].y;
  }

  if (quadtree_node_clone(src, src->head, dst, dst->head) < 0) {
    return -1;
  }

  return 0;
}

int 
quadtree_add(quadtree_t *q, 
	     double x,
	     double y,
	     bbox2d_t *bound)
{
  int pi;
  int t;

  if (q == NULL ||
      q->np == q->maxpoints) {
    return -1;
  }

  pi = q->np;
  q->np ++;

  q->point[pi].x = x;
  q->point[pi].y = y;

  t = quadtree_node_insert(q, q->head, x, y, pi);
  if (t < 0) {
    return -1;
  }

  return pi;
}


int 
quadtree_delete(quadtree_t *q,
		int pi,
		bbox2d_t *bound)
{
  int i;

  if (q == NULL ||
      pi >= q->np ||
      pi < 4) {
    return -1;
  }

  for (i = pi + 1; i < q->np; i ++) {
    q->point[i - 1].x = q->point[i].x;
    q->point[i - 1].y = q->point[i].y;
  }
  q->np --;

  if (quadtree_node_delete(q, q->head, pi, bound) < 0) {
    return -1;
  }

  return 0;
}

int 
quadtree_shift_replace(quadtree_t *q,
		       int pi)
{
  int i;
  double px;
  double py;

  /*
   * Move the last point to the index pi (used for the move operation)
   */

  px = q->point[q->np - 1].x;
  py = q->point[q->np - 1].y;

  for (i = (q->np - 1); i > pi; i --) {
    q->point[i].x = q->point[i - 1].x;
    q->point[i].y = q->point[i - 1].y;
  }

  q->point[pi].x = px;
  q->point[pi].y = py;

  if (quadtree_node_shift_replace(q, q->head, pi) < 0) {
    return -1;
  }

  return 0;
}


int 
quadtree_nearest(quadtree_t *q,
		 int include_corners,
		 double px,
		 double py)
{
  int minpi;
  double min_dist2;

  minpi = -1;
  min_dist2 = DBL_MAX;
  
  if (quadtree_node_nearest(q, q->head, px, py, include_corners, &minpi, &min_dist2) < 0) {
    return -1;
  }

  return minpi;
}

int 
quadtree_point_of_index(const quadtree_t *q,
			int i,
			double *px,
			double *py)
{
  if (q == NULL ||
      i >= q->np) {
    return -1;
  }
  
  *px = q->point[i].x;
  *py = q->point[i].y;

  return 0;
}

int 
quadtree_index_of_point(const quadtree_t *q,
			double x,
			double y)
{
  int i;
  for (i = 0; i < q->np; i ++) {
    if (q->point[i].x == x &&
	q->point[i].y == y) {
      return i;
    }
  }

  return -1;
}

int
quadtree_npoints(const quadtree_t *q)
{
  if (q != NULL) {
    return q->np;
  }
}

int
quadtree_polygon_bound(const quadtree_t *q,
		       int i,
		       bbox2d_t *bound)
{
  bound->xmin = q->bound.xmin;
  bound->xmax = q->bound.xmax;
  bound->ymin = q->bound.ymin;
  bound->ymax = q->bound.ymax;

  return 0;
}

static quadtree_node_t *
quadtree_node_create(int maxpoints,
		     int maxlevels,
		     double xmin,
		     double xmax,
		     double ymin,
		     double ymax)
{
  quadtree_node_t *q;

  q = malloc(sizeof(quadtree_node_t));
  if (q == NULL) {
    return NULL;
  }

  q->bound.xmin = xmin;
  q->bound.xmax = xmax;
  q->bound.ymin = ymin;
  q->bound.ymax = ymax;

  q->cx = (xmin + xmax)/2.0;
  q->cy = (ymin + ymax)/2.0;

  if (maxlevels == 0) {
    q->node = 0;

    q->data.leaf.np = 0;
    q->data.leaf.point = rjmcmc_create_int_array_1d(maxpoints);

    if (q->data.leaf.point == NULL) {
      return NULL;
    }

  } else {

    q->node = -1;

    q->data.child[0] = quadtree_node_create(maxpoints, 
					    maxlevels - 1,
					    xmin,
					    q->cx,
					    q->cy,
					    ymax);

    q->data.child[1] = quadtree_node_create(maxpoints, 
					    maxlevels - 1,
					    q->cx,
					    xmax,
					    q->cy,
					    ymax);

    q->data.child[2] = quadtree_node_create(maxpoints, 
					    maxlevels - 1,
					    q->cx,
					    xmax,
					    ymin,
					    q->cy);
  
    q->data.child[3] = quadtree_node_create(maxpoints, 
					    maxlevels - 1,
					    xmin,
					    q->cx,
					    ymin,
					    q->cy);
    if ((q->data.child[0] == NULL) ||
	(q->data.child[1] == NULL) ||
	(q->data.child[2] == NULL) ||
	(q->data.child[3] == NULL)) {
      return NULL;
    }
  }

  q->grand_children = 0;

  return q;
}

static void
quadtree_node_destroy(quadtree_node_t *q) 
{
  if (q != NULL) {
    
    if (q->node) {
      quadtree_node_destroy(q->data.child[0]);
      quadtree_node_destroy(q->data.child[1]);
      quadtree_node_destroy(q->data.child[2]);
      quadtree_node_destroy(q->data.child[3]);
    } else {
      rjmcmc_destroy_int_array_1d(q->data.leaf.point);
    }

    free(q);
  }
}

static int quadtree_node_insert(quadtree_t *q,
				quadtree_node_t *qn,
				double x,
				double y,
				int pi)
{
  int qpi;

  if (qn->node) {

    qn->grand_children ++;

    if (x <= qn->cx) {
      if (y >= qn->cy) {
	return quadtree_node_insert(q, qn->data.child[0], x, y, pi);
      } else {
	return quadtree_node_insert(q, qn->data.child[3], x, y, pi);
      }
    } else {
      if (y >= qn->cy) {
	return quadtree_node_insert(q, qn->data.child[1], x, y, pi);
      } else {
	return quadtree_node_insert(q, qn->data.child[2], x, y, pi);
      }
    }

  } else {

    qpi = qn->data.leaf.np;
    qn->data.leaf.np ++;

    qn->data.leaf.point[qpi] = pi;

    return 0;
  }
}

static int quadtree_node_delete(quadtree_t *q, 
				quadtree_node_t *qn, 
				int pi, 
				bbox2d_t *bound)
{
  int i;
  int pii;
  int t;

  if (qn->node) {

    /*
     * We actually have to visit all nodes since we have to update the indicies
     */
    if (qn->grand_children > 0) {
      for (i = 0; i < 4; i ++) {
	t = quadtree_node_delete(q, qn->data.child[i], pi, bound);
	if (t < 0) {
	  return -1;
	} else if (t == 1) {
	  qn->grand_children --;
	  break;
	}
      }
    } 

  } else {

    pii = -1;
    for (i = 0; i < qn->data.leaf.np; i ++) {
      
      if (qn->data.leaf.point[i] > pi) {
	qn->data.leaf.point[i] --;
      } else if (qn->data.leaf.point[i] == pi) {
	pii = i;
      }
    }

    if (pii >= 0) {
      for (i = pii + 1; i < qn->data.leaf.np; i ++) {
	qn->data.leaf.point[i - 1] = qn->data.leaf.point[i];
      }
      qn->data.leaf.np --;
    }
  }

  return 0;
}

static int quadtree_node_nearest(quadtree_t *q, 
				 quadtree_node_t *qn, 
				 double px, 
				 double py, 
				 int include_corners,
				 int *minpi,
				 double *dist2)
{
  int i;
  int pi;
  double dx;
  double dy;
  double d2;

#define TEST_NEIGHBOURS   

  if (qn->node) {

    if (qn->grand_children == 0) {
      return 0;
    }

    if (px <= qn->cx) {
      if (py >= qn->cy) {

	if (quadtree_node_nearest(q, qn->data.child[0], px, py, include_corners, minpi, dist2) < 0) {
	  return -1;
	}

#if defined(TEST_NEIGHBOURS)
	dx = qn->cx - px;
	dx = dx*dx;

	if (dx < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[1], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}
	
	dy = qn->cy - py;
	dy = dy*dy;

	if (dy < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[3], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}

	d2 = dx+dy;
	if (d2 < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[2], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}
#endif 
	    
      } else {
	if (quadtree_node_nearest(q, qn->data.child[3], px, py, include_corners, minpi, dist2) < 0) {
	  return -1;
	}

#if defined(TEST_NEIGHBOURS)
	dx = qn->cx - px;
	dx = dx*dx;

	if (dx < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[2], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}
	
	dy = qn->cy - py;
	dy = dy*dy;

	if (dy < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[0], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}

	d2 = dx+dy;
	if (d2 < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[1], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}
#endif 

      }
    } else {
      if (py >= qn->cy) {
	if (quadtree_node_nearest(q, qn->data.child[1], px, py, include_corners, minpi, dist2) < 0) {
	  return -1;
	}

#if defined(TEST_NEIGHBOURS)
	dx = qn->cx - px;
	dx = dx*dx;

	if (dx < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[0], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}
	
	dy = qn->cy - py;
	dy = dy*dy;

	if (dy < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[2], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}

	d2 = dx+dy;
	if (d2 < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[3], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}
#endif
	
      } else {
	if (quadtree_node_nearest(q, qn->data.child[2], px, py, include_corners, minpi, dist2) < 0) {
	  return -1;
	}

#if defined(TEST_NEIGHBOURS)
	dx = qn->cx - px;
	dx = dx*dx;

	if (dx < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[3], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}
	
	dy = qn->cy - py;
	dy = dy*dy;

	if (dy < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[1], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}

	d2 = dx+dy;
	if (d2 < (*dist2)) {
	  if (quadtree_node_nearest(q, qn->data.child[0], px, py, include_corners, minpi, dist2) < 0) {
	    return -1;
	  }
	}
#endif
      }
    }

  } else {

    if (include_corners) {
      for (i = 0; i < qn->data.leaf.np; i ++) {
	pi = qn->data.leaf.point[i];
	
	dx = px - q->point[pi].x;
	dy = py - q->point[pi].y;
	
	d2 = dx*dx + dy*dy;
	
	if (d2 < (*dist2)) {
	  *minpi = pi;
	  *dist2 = d2;
	}
      }
    } else {
      for (i = 0; i < qn->data.leaf.np; i ++) {
	pi = qn->data.leaf.point[i];
	if (pi < 4) {
	  continue;
	}
	
	dx = px - q->point[pi].x;
	dy = py - q->point[pi].y;
	
	d2 = dx*dx + dy*dy;
	
	if (d2 < (*dist2)) {
	  *minpi = pi;
	  *dist2 = d2;
	}
      }
    }
  }	
}

static int quadtree_node_clone(const quadtree_t *src, 
			       const quadtree_node_t *sn, 
			       quadtree_t *dst, 
			       quadtree_node_t *dn)
{
  int i;

  if (sn->node) {
    if (!dn->node) {
      return -1;
    }

    dn->grand_children = sn->grand_children;
    for (i = 0; i < 4; i ++) {
      if (quadtree_node_clone(src, sn->data.child[i],
			      dst, dn->data.child[i]) < 0) {
	return -1;
      }
    }
  } else {
    if (dn->node) {
      return -1;
    }
    
    dn->data.leaf.np = sn->data.leaf.np;
    for (i = 0; i < sn->data.leaf.np; i ++) {
      dn->data.leaf.point[i] = sn->data.leaf.point[i];
    }
  }

  return 0;
}

static int quadtree_node_shift_replace(quadtree_t *q,
				       quadtree_node_t *qn,
				       int pi)
{
  int i;

  if (qn->node) {

    if (qn->grand_children > 0) {
      for (i = 0; i < 4; i ++) {
	if (quadtree_node_shift_replace(q, qn->data.child[i], pi) < 0) {
	  return -1;
	}
      }
    }

  } else {

    for (i = 0; i < qn->data.leaf.np; i ++) {

      if (qn->data.leaf.point[i] == (q->np - 1)) {
	qn->data.leaf.point[i] = pi;
      } else if (qn->data.leaf.point[i] >= pi) {
	qn->data.leaf.point[i] ++;
      }
    }
  }

  return 0;
}
