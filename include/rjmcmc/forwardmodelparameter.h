
#ifndef forwardmodelparameter_h
#define forwardmodelparameter_h

struct _forwardmodelparameter {
  double fmin;
  double fmax;
  double fstd_value;
  double fstd_bd;
};
typedef struct _forwardmodelparameter forwardmodelparameter_t;

forwardmodelparameter_t *forwardmodelparameter_create(int nparameters);

void forwardmodelparameter_destroy(forwardmodelparameter_t *p);


#endif /* forwardmodelparameter_h */

