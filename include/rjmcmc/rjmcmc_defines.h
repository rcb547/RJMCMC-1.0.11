#ifndef RJMCMC_DEFINES_H
#define RJMCMC_DEFINES_H

#include <rjmcmc/rjmcmc_debug.h>

/* */
#define RJMCMC_INDEXCHECKVOID(i, len, msg) \
  if (i < 0 || i >= len) { \
  rjmcmc_error(msg); \
  return; \
  }

#define RJMCMC_INDEXCHECKPTR(i, len, msg) \
  if (i < 0 || i >= len) { \
  rjmcmc_error(msg); \
  return NULL; \
  }

#define RJMCMC_NULLCHECKVOID(p, msg) \
  if (p == NULL) { \
  rjmcmc_error(msg); \
  return; \
  }

#define RJMCMC_NULLCHECKINT(p, msg) \
  if (p == NULL) { \
  rjmcmc_error(msg); \
  return -1; \
  }

#define RJMCMC_NULLCHECKPTR(p, msg) \
  if (p == NULL) { \
  rjmcmc_error(msg); \
  return NULL; \
  }

#define RJMCMC_CONDITIONCHECKVOID(cond, msg) \
  if (cond) { \
  rjmcmc_error(msg); \
  return;\
  }

#define RJMCMC_CONDITIONCHECKINT(cond, msg) \
  if (cond) { \
  rjmcmc_error(msg); \
  return -1;\
  }

#define RJMCMC_CONDITIONCHECKPTR(cond, msg) \
  if (cond) { \
  rjmcmc_error(msg); \
  return NULL; \
  }

#define RJMCMC_INTCHECKPTR(i, msg) \
  if (i < 0) { \
  rjmcmc_error(msg); \
  return NULL; \
  }

#define RJMCMC_INTCHECKINT(i, msg) \
  if (i < 0) { \
  rjmcmc_error(msg); \
  return -1; \
  }


#endif /* RJMCMC_DEFINES_H */
