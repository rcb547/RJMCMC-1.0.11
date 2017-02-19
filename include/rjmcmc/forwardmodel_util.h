
#ifndef forwardmodel_util_h
#define forwardmodel_util_h


double forwardmodel_misfit_sigma_power(const double *phi,
				       int n,
				       double sigma,
				       double r);

double forwardmodel_log_determinant_sigma_power(int n,
						double sigma,
						double r);

#endif /* forwardmodel_util_h */
