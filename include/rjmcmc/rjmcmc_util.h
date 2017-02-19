#ifndef rjmcmc_util_h
#define rjmcmc_util_h

double *rjmcmc_create_array_1d(int i);
double **rjmcmc_create_array_2d(int i, int j);
double ***rjmcmc_create_array_3d(int i, int j, int k);

void rjmcmc_destroy_array_1d(double *a);
void rjmcmc_destroy_array_2d(int w, double **a);
void rjmcmc_destroy_array_3d(int w, int h, double ***a);

int *rjmcmc_create_int_array_1d(int i);
int **rjmcmc_create_int_array_2d(int i, int j);
int ***rjmcmc_create_int_array_3d(int i, int j, int k);
int ****rjmcmc_create_int_array_4d(int i, int j, int k, int l);

void rjmcmc_destroy_int_array_1d(int *a);
void rjmcmc_destroy_int_array_2d(int w, int **a);
void rjmcmc_destroy_int_array_3d(int w, int h, int ***a);
void rjmcmc_destroy_int_array_4d(int w, int h, int d, int ****a);

int rjmcmc_map_to_index(double v, double vmin, double vmax, int n);

double rjmcmc_mean(const double *v, int n);
double rjmcmc_mean_skip(const double *v, int skip, int n);

int rjmcmc_mean_variance(const double *v, int n, double *mean, double *variance);

int
rjmcmc_vector_to_histogram(int s, int n, const double *v,
			   int hn, double vmin, double vmax, int *hist);

double rjmcmc_median_from_histogram(int *hist, double vmin, double vmax, int n);
double rjmcmc_head_from_histogram(int *hist, double vmin, double vmax, int n, int drop);
double rjmcmc_tail_from_histogram(int *hist, double vmin, double vmax, int n, int drop);
double rjmcmc_mode_from_histogram(int *hist, double vmin, double vmax, int n);

/* 
 * Misc file output routines
 */
int rjmcmc_save_vector(const char *filename,
		       const double *v,
		       int n);

int rjmcmc_save_coords(const char *filename,
		       const double *x,
		       const double *y,
		       int n);

int rjmcmc_save_vector_as_histogram(const char *filename,
				    double minv,
				    double maxv,
				    int bins,
				    const double *v,
				    int n);

int rjmcmc_save_int_vector(const char *filename,
			   const int *i,
			   int n);

int rjmcmc_save_int_vector_as_histogram(const char *filename,
					int minv,
					int maxv,
					const int *i,
					int n);

int rjmcmc_save_int_coords(const char *filename,
			   const double *x,
			   const int *iy,
			   int n);

int rjmcmc_save_matrix(const char *filename,
		       const double **m,
		       int c,
		       int r);

int rjmcmc_save_int_matrix(const char *filename,
			   const int **m,
			   int c,
			   int r);

double rjmcmc_gaussian_probability(double phi,
				   double sigma);

double rjmcmc_log_gaussian_probability(double phi,
				       double sigma);

double rjmcmc_gaussian_interval_probability(double phi,
					    double sigma,
					    double deltasigma);

double rjmcmc_polynomial_value(const double *coeff,
			       int ncoeff,
			       double x);

int rjmcmc_fact(int n);
int rjmcmc_binomial(int n, int k);

int
rjmcmc_fill_coord_vector(double xmin,
			 double xmax,
			 int xsamples,
			 double *x);
			       


#endif /* rjmcmc_util_h */
