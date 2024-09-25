#!/bin/bash

# rsync -a --files-from=- /usr/local/ /usr/local/common/ <<EOF
# lib/cado-nfs-3.0.0/utils/libutils.so
# EOF

set -ex

rsync -a --files-from=- /usr/local/ /usr/local/factoring-client/ <<EOF
lib/cado-nfs-3.0.0/utils/libutils.so
lib/cado-nfs-3.0.0/sieve/las
lib/cado-nfs-3.0.0/polyselect/libpolyselect_common.so
lib/cado-nfs-3.0.0/polyselect/polyselect
lib/cado-nfs-3.0.0/polyselect/polyselect_ropt
bin/cado-nfs-client.py
EOF

rsync -a --files-from=- /usr/local/ /usr/local/discretelog-client/ <<EOF
lib/cado-nfs-3.0.0/utils/libutils.so
lib/cado-nfs-3.0.0/sieve/las
lib/cado-nfs-3.0.0/polyselect/libpolyselect_common.so
lib/cado-nfs-3.0.0/polyselect/dlpolyselect
bin/cado-nfs-client.py
EOF

rsync -a --files-from=- /usr/local/ /usr/local/common-client/ <<EOF
lib/cado-nfs-3.0.0/scripts/cadofactor/
EOF

rsync -a --files-from=- /usr/local/ /usr/local/common-server/ <<EOF
bin/cado-nfs.py
lib/cado-nfs-3.0.0/polyselect/polyselect3
lib/cado-nfs-3.0.0/filter/dup1
lib/cado-nfs-3.0.0/filter/dup2
lib/cado-nfs-3.0.0/filter/purge
lib/cado-nfs-3.0.0/scripts/cadofactor/
lib/cado-nfs-3.0.0/sieve/freerel
lib/cado-nfs-3.0.0/sieve/las
lib/cado-nfs-3.0.0/sieve/makefb
lib/cado-nfs-3.0.0/utils/antebuffer
lib/cado-nfs-3.0.0/utils/libutils.so
lib/cado-nfs-3.0.0/utils/libutils_with_io.so
share/cado-nfs-3.0.0/misc/cpubinding.conf
EOF

rsync -a --files-from=- /usr/local/ /usr/local/common-linalg/ <<EOF
lib/cado-nfs-3.0.0/utils/libutils.so
lib/cado-nfs-3.0.0/linalg/bwc/acollect
lib/cado-nfs-3.0.0/linalg/bwc/bwc.pl
lib/cado-nfs-3.0.0/linalg/bwc/cleanup
lib/cado-nfs-3.0.0/linalg/bwc/dispatch
lib/cado-nfs-3.0.0/linalg/bwc/gather
lib/cado-nfs-3.0.0/linalg/bwc/krylov
lib/cado-nfs-3.0.0/linalg/bblas/libbblas.so
lib/cado-nfs-3.0.0/linalg/bwc/libbwc_arith.so
lib/cado-nfs-3.0.0/linalg/bwc/libbwc_base.so
lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_common.so
lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_mf.so
lib/cado-nfs-3.0.0/linalg/bwc/mksol
lib/cado-nfs-3.0.0/linalg/bwc/prep
lib/cado-nfs-3.0.0/linalg/bwc/secure
EOF

rsync -a --files-from=- /usr/local/ /usr/local/factoring-linalg/ <<EOF
lib/cado-nfs-3.0.0/linalg/bwc/libarithmetic_b128.so
lib/cado-nfs-3.0.0/linalg/bwc/libarithmetic_b64.so
lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_b128_bucket.so
lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_b64_bucket.so
lib/cado-nfs-3.0.0/linalg/bwc/lingen_b64
EOF

# # This is wrong. We have multiple binary configurations per prime size, a
# # priori.
# rsync -a --files-from=- /usr/local/ /usr/local/discretelog-linalg/ <<EOF
# lib/cado-nfs-3.0.0/linalg/bwc/lingen_pz
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p1_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p2_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p3_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p4_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p5_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p6_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p7_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p8_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p9_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p10_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p11_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p12_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p13_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p14_zone.so
# lib/cado-nfs-3.0.0/linalg/bwc/libmatmul_p15_zone.so
# EOF
# 

rsync -a --files-from=- /usr/local/ /usr/local/factoring-server/ <<EOF
lib/cado-nfs-3.0.0/filter/merge
lib/cado-nfs-3.0.0/filter/replay
share/cado-nfs-3.0.0/factor/
lib/cado-nfs-3.0.0/linalg/characters
lib/cado-nfs-3.0.0/sqrt/sqrt
EOF

rsync -a --files-from=- /usr/local/ /usr/local/discretelog-server/ <<EOF
lib/cado-nfs-3.0.0/misc/descent_init_Fp
lib/cado-nfs-3.0.0/scripts/descent.py
lib/cado-nfs-3.0.0/sieve/las_descent
lib/cado-nfs-3.0.0/utils/numbertheory_tool
lib/cado-nfs-3.0.0/filter/merge-dl
lib/cado-nfs-3.0.0/filter/replay-dl
lib/cado-nfs-3.0.0/filter/reconstructlog-dl
lib/cado-nfs-3.0.0/filter/sm
lib/cado-nfs-3.0.0/filter/sm_simple
share/cado-nfs-3.0.0/dlp/
EOF

# rsync -a /usr/local/common-server/ /usr/local/discretelog-server/
rsync -a /usr/local/common-server/ /usr/local/factoring-server/
rsync -a /usr/local/common-linalg/ /usr/local/factoring-linalg/
rsync -a /usr/local/factoring-linalg/ /usr/local/factoring-server/
# rsync -a /usr/local/discretelog-linalg/ /usr/local/discretelog-server/
rsync -a /usr/local/common-client/ /usr/local/factoring-client/
# rsync -a /usr/local/common-client/ /usr/local/discretelog-client/

: unused <<EOF
lib/cado-nfs-3.0.0/misc/check_rels
lib/cado-nfs-3.0.0/linalg/bwc/blocklanczos
lib/cado-nfs-3.0.0/linalg/bwc/bwccheck
share/cado-nfs-3.0.0/polynomials
lib/cado-nfs-3.0.0/linalg/bwc/mf_bal
EOF
