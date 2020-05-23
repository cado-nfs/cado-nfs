#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include "memusage.h"

/* Returns memory usage, in KB 
 * This is the VmSize field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure. When not supported, simply return 0.
 */
long
Memusage (void)
{
#if defined(__linux__) || defined(__linux)
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", (int) pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmSize: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
#else
  return 0;
#endif
}

/* same as above, for resident memory (column RES of top) */
long
Memusage2 (void)
{
#if defined(__linux__) || defined(__linux)
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", (int) pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmRSS: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
#else
  return 0;
#endif
}

/* Returns peak memory usage, in KB 
 * This is the VmPeak field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure.
 */
long
PeakMemusage (void)
{
#if defined(__linux__) || defined(__linux)
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", (int) pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmPeak: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
#else
  return 0;
#endif
}

