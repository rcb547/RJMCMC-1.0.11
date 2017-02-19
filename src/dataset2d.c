

#include <stdlib.h>
#include <stdio.h>

#include <float.h>

#include <rjmcmc/dataset2d.h>
#include <rjmcmc/rjmcmc_util.h>

dataset2d_t *
dataset2d_allocate(int npoints)
{
  dataset2d_t *d;
  int i;

  d = malloc(sizeof(dataset2d_t));
  if (d == NULL) {
    return NULL;
  }

  d->npoints = npoints;
  d->points = malloc(sizeof(point2d_t) * npoints);
  if (d->points == NULL) {
    return NULL;
  }
  for (i = 0; i < npoints; i ++) {
    d->points[i].x = 0.0;
    d->points[i].y = 0.0;
    d->points[i].z = 0.0;
    d->points[i].n = 0.0;
  }

  d->xmin = 0.0;
  d->xmax = 0.0;

  d->ymin = 0.0;
  d->ymax = 0.0;

  d->zmin = 0.0;
  d->zmax = 0.0;

  d->lambdamin = 0.0;
  d->lambdamax = 0.0;
  d->lambdastd = 0.0;

  return d;
}

void
dataset2d_destroy(dataset2d_t *d)
{
  if (d != NULL) {
    free(d->points);
    free(d);
  }
}

dataset2d_t *
dataset2d_load_fixed(const char *filename, double sigma)
{
  FILE *fp;
  int i;
  int size;
  dataset2d_t *d;

  double x;
  double y;
  double z;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return NULL;
  }

  i = 0;
  while (!feof(fp)) {
    if (fscanf(fp, "%lf %lf %lf\n", &x, &y, &z) != 3) {
      if (feof(fp)) {
	continue;
      } else {
	rjmcmc_error(
		"dataset2d_load_estimated: error: failed to parse line\n");
	return NULL;
      }
    }

    i ++;
  }

  size = i;
  fseek(fp, 0, SEEK_SET);
  
  d = dataset2d_allocate(size);
  if (d == NULL) {
    return NULL;
  }
 
  d->xmin = DBL_MAX;
  d->xmax = -DBL_MAX;

  d->ymin = DBL_MAX;
  d->ymax = -DBL_MAX;

  d->zmin = DBL_MAX;
  d->zmax = -DBL_MAX;

  for (i = 0; i < size; i ++) {
    if (fscanf(fp, "%lf %lf %lf\n", &x, &y, &z) != 3) {
      rjmcmc_error("dataset_load: error: failed to parse on 2nd scan\n");
      return NULL;
    }

    d->points[i].x = x;
    d->points[i].y = y;
    d->points[i].z = z;
    d->points[i].n = sigma;

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

    if (z < d->zmin) {
      d->zmin = z;
    } 
    if (z > d->zmax) {
      d->zmax = z;
    }
  }

  fclose(fp);

  return d;
}

dataset2d_t *
dataset2d_load_known(const char *filename)
{
  FILE *fp;
  int i;
  int size;
  dataset2d_t *d;

  double x;
  double y;
  double z;
  double n;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return NULL;
  }

  i = 0;
  while (!feof(fp)) {
    if (fscanf(fp, "%lf %lf %lf %lf\n", &x, &y, &z, &n) != 4) {
      if (feof(fp)) {
	continue;
      } else {
	rjmcmc_error(
		     "dataset2d_load_known: error: failed to parse line %d\n", i);
	return NULL;
      }
    }

    i ++;
  }

  size = i;
  fseek(fp, 0, SEEK_SET);
  
  d = dataset2d_allocate(size);
  if (d == NULL) {
    return NULL;
  }
 
  d->xmin = DBL_MAX;
  d->xmax = -DBL_MAX;

  d->ymin = DBL_MAX;
  d->ymax = -DBL_MAX;

  d->zmin = DBL_MAX;
  d->zmax = -DBL_MAX;

  for (i = 0; i < size; i ++) {
    if (fscanf(fp, "%lf %lf %lf %lf\n", &x, &y, &z, &n) != 4) {
      rjmcmc_error(
	      "dataset_load_known: error: failed to parse on 2nd scan\n");
      return NULL;
    }

    d->points[i].x = x;
    d->points[i].y = y;
    d->points[i].z = z;
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

    if (z < d->zmin) {
      d->zmin = z;
    } 
    if (z > d->zmax) {
      d->zmax = z;
    }
  }

  fclose(fp);

  return d;
}
