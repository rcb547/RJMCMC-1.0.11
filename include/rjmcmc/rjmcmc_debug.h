#ifndef rjmcmc_debug_h
#define rjmcmc_debug_h

#include <stdarg.h>
#include <rjmcmc/rjmcmc.h>

typedef enum {
  RJMCMC_FATAL = 0,
  RJMCMC_ERROR,
  RJMCMC_WARNING,
  RJMCMC_DEBUG
} rjmcmc_debug_level_t;

typedef void (*rjmcmc_debug_function_t)(rjmcmc_debug_level_t level,
					const char *fmt,
					va_list ap);

rjmcmc_debug_function_t 
rjmcmc_debug_set_output(rjmcmc_debug_function_t fn);

void 
rjmcmc_debug_set_default_output(void);

rjmcmc_debug_level_t
rjmcmc_debug_set_level(rjmcmc_debug_level_t new_level);

void rjmcmc_fatal(const char *fmt, ...);
void rjmcmc_error(const char *fmt, ...);
void rjmcmc_warning(const char *fmt, ...);
void rjmcmc_debug(const char *fmt, ...);

#endif /* rjmcmc_debug_h */
