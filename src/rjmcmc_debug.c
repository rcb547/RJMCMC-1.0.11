
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include <rjmcmc/rjmcmc_debug.h>

static void rjmcmc_out(rjmcmc_debug_level_t level,
		       const char *fmt,
		       va_list ap);

static rjmcmc_debug_function_t rd_default = rjmcmc_out;
static rjmcmc_debug_level_t rd_level = RJMCMC_ERROR;

rjmcmc_debug_function_t 
rjmcmc_debug_set_output(rjmcmc_debug_function_t fn)
{
  rjmcmc_debug_function_t old = rd_default;

  assert(fn != NULL);
  rd_default = fn;

  return old;
}

void rjmcmc_debug_set_default_output(void)
{
  rd_default = rjmcmc_out;
}

rjmcmc_debug_level_t
rjmcmc_debug_set_level(rjmcmc_debug_level_t new_level)
{
  rjmcmc_debug_level_t old = rd_level;
  
  rd_level = new_level;
  
  return old;
}

void rjmcmc_fatal(const char *fmt, ...)
{
  va_list ap;

  assert(rd_default != NULL);

  va_start(ap, fmt);
  rd_default(RJMCMC_FATAL, fmt, ap);
  va_end(ap);
}

void rjmcmc_error(const char *fmt, ...)
{
  va_list ap;

  assert(rd_default != NULL);

  if (rd_level >= RJMCMC_ERROR) {
    va_start(ap, fmt);
    rd_default(RJMCMC_ERROR, fmt, ap);
    va_end(ap);
  }
}

void rjmcmc_warning(const char *fmt, ...)
{
  va_list ap;

  assert(rd_default != NULL);

  if (rd_level >= RJMCMC_WARNING) {
    va_start(ap, fmt);
    rd_default(RJMCMC_WARNING, fmt, ap);
    va_end(ap);
  }
}

void rjmcmc_debug(const char *fmt, ...)
{
  va_list ap;

  assert(rd_default != NULL);

  if (rd_level >= RJMCMC_DEBUG) {
    va_start(ap, fmt);
    rd_default(RJMCMC_DEBUG, fmt, ap);
    va_end(ap);
  }
}

static void rjmcmc_out(rjmcmc_debug_level_t level,
		       const char *fmt,
		       va_list ap)
{
  switch(level) {
  case RJMCMC_DEBUG:
    fprintf(stderr, "debug:");
    break;
  case RJMCMC_WARNING:
    fprintf(stderr, "warning:");
    break;
  case RJMCMC_ERROR:
    fprintf(stderr, "error:");
    break;
  case RJMCMC_FATAL:
    fprintf(stderr, "fatal:");
    break;
  default:
    fprintf(stderr, "unknown:");
    break;
  }
      
  vfprintf(stderr, fmt, ap);
}
