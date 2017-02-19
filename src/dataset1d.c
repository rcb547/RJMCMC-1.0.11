
#include <stdlib.h>
#include <stdio.h>

#include <float.h>

#include "rjmcmc/dataset1d.h"
#include "rjmcmc/rjmcmc_debug.h"

dataset1d_t *
dataset1d_create(int size)
{
  dataset1d_t *d;

  d = malloc(sizeof(dataset1d_t));
  if (d == NULL) {
    rjmcmc_fatal("dataset1d_create: unable to allocate memory\n");
    return NULL;
  }

  d->npoints = size;
  
  d->xmin = 0.0;
  d->xmax = 0.0;

  d->ymin = 0.0;
  d->ymax = 0.0;

  d->points = (point1d_t*)malloc(sizeof(point1d_t) * size);
  if (d->points == NULL) {
    rjmcmc_fatal("dataset1d_create: unable to allocate memory\n");
    return NULL;
  }

  d->lambdamin = 0.0;
  d->lambdamax = 0.0;
  d->lambdastd = 0.0;

  return d;
}

void
dataset1d_destroy(dataset1d_t *d)
{
  if (d != NULL) {
    free(d->points);
    free(d);
  }
}

dataset1d_t *
dataset1d_load_fixed(const char *filename, 
		     double n)
{
  FILE *fp;
  int i;
  int size;
  dataset1d_t *d;

  double x;
  double y;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    rjmcmc_error("dataset1d_load_fixed: unable to open file\n");
    return NULL;
  }

  i = 0;
  while (!feof(fp)) {
    if (fscanf(fp, "%lf %lf\n", &x, &y) != 2) {
      if (feof(fp)) {
	continue;
      } else {
	rjmcmc_error("dataset1d_load_fixed: unable to parse line\n");
	return NULL;
      }
    }

    i ++;
  }

  size = i;
  fseek(fp, 0, SEEK_SET);
  
  d = dataset1d_create(size);
  if (d == NULL) {
    return NULL;
  }
 
  d->xmin = DBL_MAX;
  d->xmax = -DBL_MAX;

  d->ymin = DBL_MAX;
  d->ymax = -DBL_MAX;

  for (i = 0; i < size; i ++) {
    if (fscanf(fp, "%lf %lf\n", &x, &y) != 2) {
      rjmcmc_error("dataset1d_load_fixed: unable to parse line\n");
      return NULL;
    }

    d->points[i].x = x;
    d->points[i].y = y;
    d->points[i].n = n;

    if (x < d->xmin) {
      d->xmin = x;
    }
    if (x > d->xmax) {
      d->xmax = x;
    }

    if (y < d->ymin) {
      d->ymin = y;
    } 
    if (y > d->ymax) {
      d->ymax = y;
    }
  }

  fclose(fp);

  return d;
}

dataset1d_t *
dataset1d_create_from_array(const double *x,
			    const double *y,
			    const double *n,
			    int size)
{
  dataset1d_t *d;
  int i;

  d = dataset1d_create(size);
  if (d == NULL) {
    return NULL;
  }
 
  d->xmin = DBL_MAX;
  d->xmax = -DBL_MAX;

  d->ymin = DBL_MAX;
  d->ymax = -DBL_MAX;

  for (i = 0; i < size; i ++) {
    d->points[i].x = x[i];
    d->points[i].y = y[i];
    d->points[i].n = n[i];

    if (x[i] < d->xmin) {
      d->xmin = x[i];
    }
    if (x[i] > d->xmax) {
      d->xmax = x[i];
    }

    if (y[i] < d->ymin) {
      d->ymin = y[i];
    } 
    if (y[i] > d->ymax) {
      d->ymax = y[i];
    }
  }

  return d;
}

dataset1d_t *
dataset1d_load_known(const char *filename)
{
  FILE *fp;
  int i;
  int size;
  dataset1d_t *d;

  double x;
  double y;
  double n;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    rjmcmc_error("dataset1d_load_known: failed to open file\n");
    return NULL;
  }

  i = 0;
  while (!feof(fp)) {
    if (fscanf(fp, "%lf %lf %lf\n", &x, &y, &n) != 3) {
      if (feof(fp)) {
	continue;
      } else {
	rjmcmc_error("dataset1d_load: failed to parse line\n");
	return NULL;
      }
    }

    i ++;
  }

  size = i;
  fseek(fp, 0, SEEK_SET);
  
  d = dataset1d_create(size);
  if (d == NULL) {
    return NULL;
  }
 
  d->xmin = DBL_MAX;
  d->xmax = -DBL_MAX;

  d->ymin = DBL_MAX;
  d->ymax = -DBL_MAX;

  for (i = 0; i < size; i ++) {
    if (fscanf(fp, "%lf %lf %lf\n", &x, &y, &n) != 3) {
      rjmcmc_error("dataset1d_load: failed to parse line\n");
      return NULL;
    }

    d->points[i].x = x;
    d->points[i].y = y;
    d->points[i].n = n;

    if (x < d->xmin) {
      d->xmin = x;
    }
    if (x > d->xmax) {
      d->xmax = x;
    }

    if (y < d->ymin) {
      d->ymin = y;
    } 
    if (y > d->ymax) {
      d->ymax = y;
    }
  }

  fclose(fp);

  return d;
}

static int cmp_point1d(const void *a, const void *b)
{
  point1d_t *pa = (point1d_t*)a;
  point1d_t *pb = (point1d_t*)b;

  return pa->x - pb->x;
}

void
dataset1d_sort(dataset1d_t *d)
{
  if (d == NULL) {
    rjmcmc_error("dataset1d_sort: NULL dataset\n");
    return;
  }

  qsort(d->points,
	sizeof(point1d_t),
	d->npoints,
	cmp_point1d);
}

int
dataset1d_range(const dataset1d_t *data,
		double xl,
		double xr,
		int *_xi,
		int *_xj)
{
  int xi;
  int xj;

  xi = 0;
  while (xi < data->npoints &&
	 data->points[xi].x < xl) {
    xi ++;
  }
  if (xi == data->npoints) {
    return -1;
  }
    
  xj = data->npoints - 1;
  while (xj > xi &&
	 data->points[xj].x > xr) {
    xj --;
  }
  
  *_xi = xi;
  *_xj = xj;

  return xj - xi + 1;
}

int 
dataset1d_mean_variance(const dataset1d_t *data,
			int xi,
			int xj,
			double *_mean,
			double *_variance)
{
  /* Algorithm adapted from Knuth Semi-numerical algorithms */
  int i;
  int n;

  double d;
  double mean;
  double m;

  double delta;

  n = 0;
  mean = 0.0;
  m = 0.0;

  if (data == NULL) {
    rjmcmc_error("dataset1d_mean_variance: NULL data\n");
    return -1;
  }

  if (xi < 0 ||
      xj >= data->npoints ||
      xi > xj) {
    rjmcmc_error("dataset1d_mean_variance: invalid start/end indices\n");
    return -1;
  }
      
  for (i = xi; i <= xj; i ++) {

    n ++;
    d = data->points[i].y;
    delta = d - mean;
    
    mean += delta/(double)n;
    m += delta * (d - mean);
  }

  if (n > 1) {
    *_mean = mean;
    *_variance = m/(double)(n - 1); /* Note the use of the sample variance */
  }

  return n;
}

