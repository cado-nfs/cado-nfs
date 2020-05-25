#include "cado.h" // IWYU pragma: keep
#include <stdio.h>      // printf
#include <stdlib.h>        // for abort exit
#include <string.h>
#include <errno.h>
#include <inttypes.h>
#include "purgedfile.h"
#include "gzip.h"

/* Read all lines which begin with # in the input file, until we find one
 * which matches the desided format.
 */

void
purgedfile_read_firstline (const char *fname, uint64_t *nrows, uint64_t *ncols)
{
  FILE *f_tmp = fopen_maybe_compressed (fname, "rb");
  if (!f_tmp)
  {
    fprintf(stderr, "%s: %s\n", fname, strerror(errno));
    abort();
  }
  char buf[1024];
  while (fgets(buf, sizeof(buf), f_tmp)) {
      if (*buf != '#') {
          break;
      }
      int ret = sscanf(buf, "# %" SCNu64 " %" SCNu64 "", nrows, ncols);
      if (ret == 2) {
          fclose_maybe_compressed(f_tmp, fname);
          return;
      }
  }
  fprintf(stderr,
          "Parse error while reading %s: no header line with desired format\n", fname);
  fclose_maybe_compressed(f_tmp, fname);
  exit(EXIT_FAILURE);
}
