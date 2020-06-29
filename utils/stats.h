#ifndef STATS_H_
#define STATS_H_

#include <stdint.h>
#include <stdio.h>

#define MIN_LOG_REPORT 10

struct stats_struct_s {
  uint64_t *followed_var;
  uint64_t last_report;
  double t0;
  uint8_t log_report;
  uint8_t max_log_report;
  FILE *out;
  const char *name;
  const char *abbrv;
  const char *verb;
  const char *outofname;
};

typedef struct stats_struct_s stats_data_t[1];
typedef struct stats_struct_s * stats_data_ptr;

#ifdef __cplusplus
extern "C" {
#endif
void stats_init (stats_data_t, FILE *, uint64_t *, uint8_t, const char *,
                 const char *, const char *, const char *);
int stats_test_progress (stats_data_t);
void stats_print_progress (stats_data_t, uint64_t, uint64_t, size_t, int);
#ifdef __cplusplus
}
#endif
#endif /* STATS_H_ */
