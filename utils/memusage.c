#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <unistd.h>
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#include "memusage.h"

/* Returns memory usage, in KiB (1024 bytes).
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
#elif defined(__APPLE__)
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  kern_return_t ret = task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count);
  if (ret != KERN_SUCCESS || count != MACH_TASK_BASIC_INFO_COUNT)
      return -1;
  /* seems that virtual_size is often inaccessible */
  return info.resident_size >> 10;
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
#elif defined(__APPLE__)
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  kern_return_t ret = task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count);
  if (ret != KERN_SUCCESS || count != MACH_TASK_BASIC_INFO_COUNT)
      return -1;
  return info.resident_size >> 10;
#else
  return 0;
#endif
}

/* Returns peak memory usage, in KiB (1024 bytes).
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
#elif defined(__APPLE__)
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  kern_return_t ret = task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count);
  if (ret != KERN_SUCCESS || count != MACH_TASK_BASIC_INFO_COUNT)
      return -1;
  return info.resident_size_max >> 10;
#else
  return 0;
#endif
}

