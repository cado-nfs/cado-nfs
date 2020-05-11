#include <hwloc.h>
#ifdef HWLOC_API_VERSION
#if HWLOC_API_VERSION < 0x00010400
#error "too old, never checked"
#endif
#endif
int main(void)
{
  return hwloc_get_api_version() == HWLOC_API_VERSION ? EXIT_SUCCESS : EXIT_FAILURE;
}
