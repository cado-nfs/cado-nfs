#!/bin/bash

if [ -z $1 ] ; then
  echo "Missing argument"
  echo "Usage: $0 <parameter_file>"
  exit 1
fi

source $1

check_error () {
  if [ $1 -ne 0 ] ; then
    echo "Error !"
    echo "Exited with $1"
    exit 1
  else
    echo "done"
  fi
}

POLY="${DATADIR}/${NAME}.poly"
BADIDEALS="${DATADIR}/${NAME}.magmanmbrthry.badideals"
BADINFO="${DATADIR}/${NAME}.magmanmbrthry.badidealinfo"
INPUTRELS="${DATADIR}/${NAME}.rels"
PURGED="${WDIR}/${NAME}.purge.purged"
DELETED="${WDIR}/${NAME}.purge.deleted"
RENUMBER="${WDIR}/${NAME}.renumber.gz"
FREERELS="${WDIR}/${NAME}.freerels.gz"
NODUP1="${WDIR}/${NAME}.nodup1"
NODUP2="${WDIR}/${NAME}.nodup2"
NODUP3="${WDIR}/${NAME}.nodup3"
HISFILE="${WDIR}/${NAME}.merge.his"
MERGEFORBIDDENCOLS="${WDIR}/${NAME}.merge.forbidden.cols"
MAT_TMP="${WDIR}/${NAME}.replay.matrix.txt"
IDEALSFILE_TMP="${WDIR}/${NAME}.replay.ideals"
INDEXFILE="${WDIR}/${NAME}.replay.relsets"
MAT="${WDIR}/${NAME}.rearrange.matrix.txt"
SIDEINFO="${WDIR}/${NAME}.rearrange.side.info.txt"
IDEALSFILE="${WDIR}/${NAME}.rearrange.ideals"
SMFILE_TMP="${WDIR}/${NAME}.sm"
SMFILE="${WDIR}/${NAME}.rearrange.sm"
KERFILE="${WDIR}/${NAME}.LA.ker"
FINALDLOG="${WDIR}/${NAME}.final.dlog"
BINMAT="${WDIR}/${NAME}.matrix.bin"

SM=$(printf ",%s" "${SM_ARRAY[@]}"); SM=${SM:1}
LPBS=$(printf ",%s" "${LPB_ARRAY[@]}"); LPBS=${LPBS:1}
NSM=`echo "$SM" | tr , + | bc`

declare -i i DUP1_NSLICES

DUP1_NSLICES_LOG="0"
DUP1_NSLICES=`echo "2^${DUP1_NSLICES_LOG}" | bc`

DUP1_SELECTOR="["
for ((i=0;i<DUP1_NSLICES;++i)); do DUP1_SELECTOR="${DUP1_SELECTOR}$i"; done
DUP1_SELECTOR="${DUP1_SELECTOR}]"

if [ "${#SM_ARRAY[@]}" -ne "${#LPB_ARRAY[@]}" ] ; then
  echo "SM_ARRAY and LPB_ARRAY must have the same size"
  exit 1
fi
if [ "${#SM_ARRAY[@]}" -ne "`grep -c ^poly[0-8] ${POLY}`" ] ; then
  echo "The size of SM_ARRAY and LPB_ARRAY must be the number of polynomials"
  echo "in ${POLY}"
  exit 1
fi

: ${KEEP:="$NSM"}
: ${coverNmax:="50.0"}

mkdir -p ${WDIR}

#### freerel
echo -n "freerel ... "
ARGS_FR="-poly ${POLY} -renumber ${RENUMBER} -out ${FREERELS} -lcideals\
         -lpbs ${LPBS} -badideals ${BADIDEALS}"
LOG_FR="${WDIR}/${NAME}.freerel.log"
${BUILDDIR}/sieve/freerel ${ARGS_FR} > ${LOG_FR} 2>&1
check_error "$?"

#### debug_renumber
echo -n "debug_renumber ... "
ARGS_DR="-poly ${POLY} -renumber ${RENUMBER} -check"
LOG1_DR="${WDIR}/${NAME}.debug_renumber.log"
LOG2_DR="${WDIR}/${NAME}.debug_renumber.check"
${BUILDDIR}/misc/debug_renumber ${ARGS_DR} > ${LOG1_DR} 2>${LOG2_DR}
check_error "$?"

#print some info
NCOLS=`grep -o "nprimes=[0-9]*$" ${LOG_FR} | cut -d= -f2`
echo "# Size of the renumbering table: ${NCOLS}"

#### dup1 (XXX useless for now but maybe we will need it later)
echo -n "dup1 ... "
mkdir -p ${NODUP1}/0 ${NODUP1}/1
ARGS_DUP1="-out ${NODUP1} -prefix ${NAME} -n ${DUP1_NSLICES_LOG} ${INPUTRELS}"
LOG_DUP1="${WDIR}/${NAME}.dup1.log"
${BUILDDIR}/filter/dup1 ${ARGS_DUP1} > ${LOG_DUP1} 2>&1
check_error "$?"

#### dup2
for ((i=0;i<DUP1_NSLICES;++i)); do
  echo -n "dup2 slice $i ... "
  NRELS_SLICE=`grep "slice $i rec" ${LOG_DUP1} | cut -d " " -f 5`
  mkdir -p ${NODUP2}/$i
  ARGS_DUP2="-renumber ${RENUMBER} -nrels ${NRELS_SLICE} -dl\
             -outdir ${NODUP2}/$i -badidealinfo ${BADINFO} ${NODUP1}/$i/*"
  LOG_DUP2="${WDIR}/${NAME}.dup2-$i.log"
  ${BUILDDIR}/filter/dup2 ${ARGS_DUP2} > ${LOG_DUP2} 2>&1
  check_error "$?"
done

#### split
echo -n "split ... "
mkdir -p ${WDIR}/${NAME}.nodup3/
ARGS_SP="-poly ${POLY} -renumber ${RENUMBER} -outprefix ${NAME}\
        -outdir ${NODUP3} ${NODUP2}/${DUP1_SELECTOR}/*"
LOG_SP="${WDIR}/${NAME}.split.log"
${BUILDDIR}/misc/split_renumbered_rels ${ARGS_SP} > ${LOG_SP} 2>&1
check_error "$?"

# print some info
grep "type 0" ${LOG_SP} | cut -d " " -f -6

# XXX (hack for p4dd5b example: too many 01 rels, even the number of 01 and 02
# rels)
TMP=`wc -l ${NODUP3}/${NAME}.02.rels | cut -d " " -f 1`
mv ${NODUP3}/${NAME}.01.rels ${NODUP3}/${NAME}.01.rels.bak
head -n ${TMP} ${NODUP3}/${NAME}.01.rels.bak > ${NODUP3}/${NAME}.01.rels

# print some info
wc -l ${NODUP3}/${NAME}.0[0-9].rels
NRELS_PURGE=`wc -l ${NODUP3}/${NAME}.0[0-9].rels | tail -n 1 | cut -d " " -f 3`

#### purge
echo -n "purge ... "
if [ -z ${NO_PURGE} ] ; then
  ARGS_PURGE="-out ${PURGED} -nrels ${NRELS_PURGE} -col-max-index ${NCOLS}\
              -col-min-index 0 -keep ${KEEP} -outdel ${DELETED}\
              ${NODUP3}/${NAME}.0[0-9].rels"
  LOG_PURGE="${WDIR}/${NAME}.purge.log"
  ${BUILDDIR}/filter/purge ${ARGS_PURGE} >${LOG_PURGE} 2>&1
  check_error "$?"
  
  # print some info
  tail -n 3 ${LOG_PURGE}
else
  echo " skipped"
  echo "# ${NRELS_PURGE} ${NCOLS} ${NCOLS}" > ${PURGED}
  cat ${NODUP3}/${NAME}.0[0-9].rels >> ${PURGED}
  touch ${DELETED}
fi

#### merge
echo -n "list_forbidden_cols_for_merge ... "
# Create list of columns to disable them during merge
ARGS_LIST="-purgedfile ${PURGED} -poly ${POLY} -renumber ${RENUMBER}\
           -outfile ${MERGEFORBIDDENCOLS}"
LOG_LIST="${WDIR}/${NAME}.forbidden.cols.log"
${BUILDDIR}/misc/list_forbidden_cols_for_merge ${ARGS_LIST} > ${LOG_LIST} 2>&1
check_error "$?"

#### merge
echo -n "merge ... "
if [ -z ${NO_MERGE} ] ; then

  ARGS_MERGE="-mat ${PURGED} -out ${HISFILE} -skip 0 -keep ${KEEP} -maxlevel 25\
              -forbw 3 -coverNmax ${coverNmax}\
              -forbidden-cols ${MERGEFORBIDDENCOLS}"
  LOG_MERGE="${WDIR}/${NAME}.merge.log"
  ${BUILDDIR}/filter/merge-dl ${ARGS_MERGE} > ${LOG_MERGE} 2>&1
  check_error "$?"

  # print some info
  grep "^Final matrix" ${LOG_MERGE}
else
  echo " skipped"
  touch ${HISFILE}
fi

#### replay
echo -n "replay ... "
ARGS_REPLAY="-purged ${PURGED} -his ${HISFILE} -out ${MAT_TMP}\
             -ideals ${IDEALSFILE_TMP} -index ${INDEXFILE}"
LOG_REPLAY="${WDIR}/${NAME}.replay.log"
${BUILDDIR}/filter/replay-dl ${ARGS_REPLAY} > ${LOG_REPLAY} 2>&1
check_error "$?"

# print some info
grep "^Sparse submatrix" ${LOG_REPLAY}

#### sm
echo -n "sm ... "
ARGS_SM="-poly ${POLY} -purged ${PURGED} -index ${INDEXFILE} -out ${SMFILE_TMP}\
         -ell ${ELL} -nsm ${SM} -t 4"
LOG_SM="${WDIR}/${NAME}.sm.log"
${BUILDDIR}/filter/sm ${ARGS_SM} > ${LOG_SM} 2>&1
check_error "$?"

#### rearrange the matrix by block
echo -n "rearrange ... "
ARGS_REA="-poly ${POLY} -renumber ${RENUMBER}\
          -ideals ${IDEALSFILE_TMP} -new-ideals ${IDEALSFILE}\
          -sm ${SMFILE_TMP} -new-sm ${SMFILE} -nsm ${SM}\
          -matrix ${MAT_TMP} -new-matrix ${MAT} -side-info ${SIDEINFO}"
LOG_REA="${WDIR}/${NAME}.rearrange.log"
${BUILDDIR}/misc/rearrange_MNFS_matrix ${ARGS_REA} > ${LOG_REA} 2>&1
check_error "$?"

#### LA (in magma)
echo -n "LA (in magma) ... "
ARGS_LA="ell:=${ELL} nmaps:=${NSM} sparsefile:=${MAT} smfile:=${SMFILE}\
         kerfile:=${KERFILE}"
LOG_LA="${WDIR}/${NAME}.LA.log"
CMD="magma ${ARGS_LA} ${BUILDDIR}/../../scripts/linalg.mag"
echo ${CMD} > ${LOG_LA}
${CMD} >> ${LOG_LA} 2>&1
check_error "$?"

#print some info
head ${KERFILE}

#### reconstructlog
echo -n "reconstructlog-dl ... "
if [ -z ${PARTIAL} ] ; then
  AP=""
else
  AP="-partial"
fi
ARGS_REC="-log ${KERFILE} -ell ${ELL} -out ${FINALDLOG} -poly ${POLY}\
          -renumber ${RENUMBER} -ideals ${IDEALSFILE} -nsm ${SM} -mt 4\
          -purged ${PURGED} -nrels ${NRELS_PURGE} -relsdel ${DELETED}"
LOG_REC="${WDIR}/${NAME}.reconstructlog.log"
${BUILDDIR}/filter/reconstructlog-dl ${ARGS_REC} ${AP} > ${LOG_REC} 2>&1
check_error "$?"

#### LA (with bwc)

echo -n "mf_scan (for bwc) ... "
# create binary file from ${MAT}
ARGS_MF_SCAN="--ascii-in mfile=${MAT} --binary-out ofile=${BINMAT} --freq\
              --withcoeffs"
LOG_MF_SCAN="${WDIR}/${NAME}.bwc.mf_scan.log"
CMD="${BUILDDIR}/linalg/bwc/mf_scan ${ARGS_MF_SCAN}" > ${LOG_MF_SCAN}
${CMD} >> ${LOG_MF_SCAN} 2>&1
check_error "$?"

echo -n "LA (with bwc) ... "
BWCWDIR="${WDIR}/${NAME}.bwc"
BWC_N=${NSM}
BWC_M=`expr 2 \* ${BWC_N}`
rm -rf ${BWCDIR}
mkdir -p ${BWCWDIR}
ARGS_BWC=":complete thr=4 m=${BWC_M} n=${BWC_N} nullspace=right matrix=${BINMAT}\
         rhs=${SMFILE} prime=${ELL} wdir=${BWCWDIR} mm_impl=basicp\
         allow_zero_on_rhs=1"
         #cpubinding=${BUILDDIR}/../../linalg/bwc/cpubinding.conf"
STDOUT_BWC="${WDIR}/${NAME}.bwc.stdout"
STDERR_BWC="${WDIR}/${NAME}.bwc.stderr"
${BUILDDIR}/linalg/bwc/bwc.pl ${ARGS_BWC} > ${STDOUT_BWC} 2> ${STDERR_BWC}
check_error "$?"

#### reconstructlog (with data from bwc)
echo -n "reconstructlog-dl ... "
if [ -z ${PARTIAL} ] ; then
  AP=""
else
  AP="-partial"
fi
BWCKERFILE="${BWCWDIR}/K.sols0-1.0.truncated.txt"
BWCFINALDLOG="${WDIR}/${NAME}.bwc.final.dlog"
ARGS_REC="-log ${BWCKERFILE} -ell ${ELL} -out ${BWCFINALDLOG} -poly ${POLY}\
          -renumber ${RENUMBER} -ideals ${IDEALSFILE} -nsm ${SM} -mt 4\
          -purged ${PURGED} -nrels ${NRELS_PURGE} -relsdel ${DELETED}"
LOG_REC="${WDIR}/${NAME}.bwc.reconstructlog.log"
${BUILDDIR}/filter/reconstructlog-dl ${ARGS_REC} ${AP} > ${LOG_REC} 2>&1
check_error "$?"

