#!/bin/sh

set -ex

if ! [ -d flint2/.git ] ; then
    echo "Need flint2/.git" >&2
    exit 1
fi

(cd flint2/ ; git pull)

if git status -uno | grep -q flint-fft/ ; then
    echo "uncommitted changes in flint-fft">&2
    exit 1
fi

rm -rf flint-fft/
rsync -a flint2/fft/ flint-fft/
rm -rf flint-fft/profile
rm -rf flint-fft/test
rm -rf flint-fft/tune
cp flint2/fft_tuning64.in flint-fft/fft_tuning.h
cp flint2/fft.h flint-fft/
cp flint2/ulong_extras.h flint-fft/
cp flint2/ulong_extras/revbin.c flint-fft/
cp flint2/mpn_extras/mulmod_2expp1_basecase.c flint-fft/
cp flint2/flint.h flint-fft/
cp flint2/longlong.h flint-fft/
cp flint2/memory_manager.c flint-fft/
cp flint2/printf.c flint-fft/
cp flint2/scanf.c flint-fft/
cp flint2/exception.c flint-fft/
cp flint2/exception.h flint-fft/
cp flint2/clz_tab.c flint-fft/
mkdir flint-fft/doc
cp flint2/doc/source/fft.rst flint-fft/doc/

# This is a custom trimmed-down copy.

cat > flint-fft/mpn_extras.h <<EOF
#ifndef MPN_EXTRAS_H_
#define MPN_EXTRAS_H_

#include "gmp.h"
#include "flint.h"

#define BITS_TO_LIMBS(b) (((b) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

int flint_mpn_mulmod_2expp1_basecase(mp_ptr xp, mp_srcptr yp, mp_srcptr zp, int c,
    mp_bitcnt_t b, mp_ptr tp);

#endif	/* MPN_EXTRAS_H_ */
EOF

for f in $(find flint-fft -name '*.h' | xargs grep -l 'extern "C"') ; do
    ex -s $f >/dev/null <<EOF
/ifdef.*__cplusplus
s/^/\/\/ temp /

s/^/\/\/ temp /

s/^/\/\/ temp /
/ifdef.*__cplusplus
s/^/\/\/ temp /

s/^/\/\/ temp /

s/^/\/\/ temp /
wq
EOF
done

ex -s flint-fft/ulong_extras.h >/dev/null <<EOF
/n_primes_struct
?typedef struct
i
#if 0
.
/n_randlimb
i
#endif

.
/mp_limb_t.*p_lo/
s/p_lo/FLINT_SET_BUT_UNUSED(p_lo)/
wq
EOF

ex -s flint-fft/fft.h >/dev/null <<EOF
$
?ifdef.*__cplusplus
i
/* CADO-NFS addition */

/* mpn_mulmod_2expp1 is an internal function exposed by mpir, but the
 * real symbol is mpn_mulmod_2expp1_basecase anyway. The tarball we're
 * extracting here has the very same code (with a more permissive
 * license) as flint_mpn_mulmod_2expp1_basecase
 *
 * Bottom line: if we use gmp and not mpir, we may use the code we have
 * here.
 */
#define mpn_mulmod_2expp1 flint_mpn_mulmod_2expp1_basecase

/* end CADO-NFS addition */

.
/endif
+1
i

#include "transform_interface.h"
.
wq
EOF

sed -e '/thread_pool.h/ d' -i flint-fft/memory_manager.c
sed -e 's/^.*\(fmpz\|mpfr\)/\/\/ &/' -i flint-fft/memory_manager.c
ex -s flint-fft/memory_manager.c >/dev/null <<EOF
/global_thread_pool_initialized
i
/*
.
/}
+
i
*/
.
wq
EOF

sed -e 's/t = flint_malloc/t = (mp_ptr) flint_malloc/' -i flint-fft/fft.h
sed -e '/config.h/ d' -i flint-fft/flint.h
sed -e '/gmpcompat.h/ d' -i flint-fft/flint.h
sed -e '/<mpfr.h>/ d' -i flint-fft/flint.h
sed -e '/mpfr_struct/ d' -i flint-fft/flint.h
(echo '/MPFR_VERSION_MAJOR' ; echo .,+2d ; echo wq) | ex flint-fft/flint.h

find flint-fft -name '*.[ch]' | xargs -n 1 sed -e "s/FLINT_DLL *//g" -i
find flint-fft -name '*.[ch]' -o -name fft.rst | xargs -n 1 sed -e "s/fft_combine_bits/fft_addcombine_bits/g" -i
find flint-fft -name '*.[ch]' | xargs -n 1 sed -e "s/HAVE_OPENMP/defined(_OPENMP)/g" -i
find flint-fft -name '*.[ch]' | xargs -n 1 indent -kr -i4 -sc -fca -fc1 -lc78

find flint-fft -name '*.h' | xargs -n 1 sed -e "s/^ *\/\/ temp //g" -i
find flint-fft -type f  | xargs grep -l 'fmpz' | xargs -n 1 sed -e 's/^#include "fmpz/\/\/ #include "fmpz/g' -i

(cd flint-fft ; ctags -R . )

find flint-fft -name '*~' | xargs rm

perl ./parse-flint-fft-doc.pl

for f in \
        normmod_2expp1.c                \
        mulmod_2expp1.c                 \
        mul_fft_main.c                  \
        ifft_truncate_sqrt2.c           \
        ifft_mfa_truncate_sqrt2.c       \
        fft_mfa_truncate_sqrt2_inner.c  \
        fft_mfa_truncate_sqrt2.c        \
        split_bits.c                    \
        memory_manager.c                \
    ; do
ex -s flint-fft/$f > /dev/null <<EOF
1
i
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-parameter"
/* flint uses unprotected openmp pragmas every so often */
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
.
wq
EOF
done

#for f in `find flint-fft -name '*.c'` ; do
#ex -s $f > /dev/null <<EOF
#1
#i
##include "cado.h"
#.
#wq
#EOF
#done
#
git co HEAD flint-fft/transform_interface.c
git co HEAD flint-fft/transform_interface.h

(cd flint-fft ; ctags -R . )
