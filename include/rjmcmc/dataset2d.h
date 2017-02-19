#ifndef dataset2d_h
#define dataset2d_h

/** \file dataset2d.h

\brief 2D Dataset Storage

The ::dataset2d_t type is used for storing a 2D dataset with optional
noise parameters.
*/

struct _point2d {
  double x;
  double y;
  double z;
  double n;
};
typedef struct _point2d point2d_t;

struct _dataset2d {

  double xmin;
  double xmax;

  double ymin;
  double ymax;

  double zmin;
  double zmax;

  point2d_t *points;
  int npoints;

  /*
   * Estimated lambda for hierarchical datasets
   */
  double lambdamin;
  double lambdamax;
  double lambdastd;
};
typedef struct _dataset2d dataset2d_t;

/** \brief Create a new empty dataset

Allocate a new ::dataset2d_t object of a given size.

\param The maximum number of points.
*/
dataset2d_t *
dataset2d_allocate(int npoints);

/** \brief Destroy a dataset.

Deallocate an existing ::dataset2d_t object.

\param d The dataset to destroy.
*/
void
dataset2d_destroy(dataset2d_t *d);

/** \brief Load a dataset with a fixed noise parameter.

Loads a file containing x, y, z coordinates and applies a fixed
noise level to each data point. The expected format is a 
text file containing space separated x, y and z coordinates
on each line.

\param filename The file to load.
\param sigma The fixed noise level.
 */
dataset2d_t *
dataset2d_load_fixed(const char *filename, double sigma);

/** \brief Load a dataset with known noise.

Loads a file containing x, y, z coordinates and noise. 
The expected format is a text file containing space separated x, y and z
coordinates and noise on each line.

\param filename The file to load.
*/
dataset2d_t *
dataset2d_load_known(const char *filename);

#endif /* dataset2d_h */
