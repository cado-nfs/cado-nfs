#!/bin/sh

set -e

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

# This is a custom trimmed-down copy.

cat > flint-fft/mpn_extras.h <<EOF
#ifndef MPN_EXTRAS_H_
#define MPN_EXTRAS_H_

#include "mpir.h"
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

sed -e '/config.h/ d' -i flint-fft/flint.h
sed -e '/gmpcompat.h/ d' -i flint-fft/flint.h
sed -e '/<mpfr.h>/ d' -i flint-fft/flint.h
sed -e '/mpfr_struct/ d' -i flint-fft/flint.h
(echo '/MPFR_VERSION_MAJOR' ; echo .,+2d ; echo wq) | ex flint-fft/flint.h

find flint-fft -name '*.[ch]' | xargs -n 1 indent -kr -i4 -sc -fca -fc1 -lc78

find flint-fft -name '*.h' | xargs -n 1 sed -e "s/^ *\/\/ temp //g" -i
find flint-fft -type f  | xargs grep -l 'fmpz' | xargs -n 1 sed -e 's/^#include "fmpz/\/\/ #include "fmpz/g' -i
sed -e 's/^.*\(fmpz\|mpfr\)/\/\/ &/' -i flint-fft/memory_manager.c

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
    ; do
ex -s flint-fft/$f > /dev/null <<EOF
1
i
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
.
wq
EOF
done


(cd flint-fft ; ctags -R . )
