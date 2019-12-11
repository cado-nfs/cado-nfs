#ifndef FACUL_FWD_HPP_
#define FACUL_FWD_HPP_

/* Some forward declarations to simplify dependencies */

typedef struct facul_strategy_s facul_strategy_t;
typedef struct facul_method_s {
  long method; /* Which method to use (P-1, P+1 or ECM) */
  void *plan;  /* Parameters for that method */
} facul_method_t;
class FaculModulusBase;

#endif
