#ifndef dataset1d_h
#define dataset1d_h

/** \file dataset1d.h

\brief 1D Dataset Storage

The ::dataset1d_t type is used for storing a 1D dataset with optional
noise parameters.
*/

struct _point1d {
  double x;
  double y;
  double n;
};
typedef struct _point1d point1d_t;

struct _dataset1d {

  double xmin;
  double xmax;

  double ymin;
  double ymax;

  point1d_t *points;
  int npoints;

  /*
   * Estimated lambda for hierarchical datasets
   */
  double lambdamin;
  double lambdamax;
  double lambdastd;
};
typedef struct _dataset1d dataset1d_t;

/** \brief Create a new empty dataset

Allocate a new ::dataset1d_t object of a given size.

\param The maximum number of points.
*/
dataset1d_t *
dataset1d_create(int size);

/** \brief Destroy a dataset.

Deallocate an existing ::dataset1d_t object.

\param d The dataset to destroy.
*/
void
dataset1d_destroy(dataset1d_t *d);

/** \brief Load a dataset with a fixed noise parameter.

Loads a file containing x, y coordinates and applies a fixed
noise level to each data point. The expected format is a 
text file containing space separated x and y coordinates
on each line.

\param filename The file to load.
\param n The fixed noise level.
 */
dataset1d_t *
dataset1d_load_fixed(const char *filename, 
		     double n);

/** \brief Load a dataset with known noise.

Loads a file containing x, y coordinates and noise. 
The expected format is a text file containing space separated x and y 
coordinates and noise on each line.

\param filename The file to load.
*/
dataset1d_t *
dataset1d_load_known(const char *filename);

/** \brief Create a dataset from arrays.

Creates a ::dataset1d_t object from existing arrays of coordinates 
and noise. The arrays are copied into the ::dataset1d_t object.

\param x X coordinates
\param y Y coordinates
\param n Noise values per point
\int size The number of (X, Y, N) tuples.
*/
dataset1d_t *
dataset1d_create_from_array(const double *x,
			    const double *y,
			    const double *n,
			    int size);

/** \brief Reorder dataset by x coordinate

Sorts the points within the dataset to ascending in x-coordinate.

\param d The ::dataset1d_t object to sort.
*/
void
dataset1d_sort(dataset1d_t *d);

/** \brief Subset of a dataset

Outputs the indices of the subset of a dataset that is within a 
given range of x coordinates.

\param data The ::dataset1d_t object
\param xl The left x-coordinate of the range.
\param xr The right x-coordinate of the range.
\param xi The index of the first data point in the x-coordinate range.
\param xj The index of the last data point in the x-coordinate range.

\returns The number of data points in the x-coordinate range or -1 on error.
*/
int
dataset1d_range(const dataset1d_t *data,
		double xl,
		double xr,
		int *xi,
		int *xj);


/** \brief Compute the mean and variance

Computes the mean and variance of a set of data points in a given 
range of indices.

\param d The ::dataset1d_t object
\param xi The first point to include in the mean/variance
\param xj The last point to include in the mean/variance
\param mean The mean output
\param variance The variance output.

\returns The number of points in the range or -1 on error.

 *
 */
int 
dataset1d_mean_variance(const dataset1d_t *d,
			int xi,
			int xj,
			double *mean,
			double *variance);

		       
#endif /* dataset1d_h */
