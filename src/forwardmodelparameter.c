
#include <stdlib.h>

#include <rjmcmc/forwardmodelparameter.h>

forwardmodelparameter_t *forwardmodelparameter_create(int nparameters)
{
  return (forwardmodelparameter_t *)malloc(sizeof(forwardmodelparameter_t) * nparameters);
}

void forwardmodelparameter_destroy(forwardmodelparameter_t *p)
{
  free(p);
}

