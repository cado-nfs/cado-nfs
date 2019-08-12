#!/bin/sh

# Debug tools to track down ideals and relations in the filtering step.
# Examples are at the end of the file.

DIRNAME="/tmp/test.filter"
PREFIX="test."

PURGEDFILE="${DIRNAME}/${PREFIX}rels.purged.gz" # produced by purge
INDEXFILE="${DIRNAME}/${PREFIX}index" #produced by replay
MATRIXCOLS="${DIRNAME}/${PREFIX}ideals" #produced by replay
NODUPDIR="${DIRNAME}/${PREFIX}nodup"
RELSDIR="${DIRNAME}"

#should be create with print_table
RENUMBERTAB="${DIRNAME}/${PREFIX}table.pretty"


check_file_exists () {
  if [ ! -f "$1" ] ; then
    echo "$1 does not exist or is not a regular file."
    exit 1
  fi
}

check_dir_exists () {
  if [ ! -d "$1" ] ; then
    echo "$1 does not exist or is not a directory."
    exit 1
  fi
}

from_matrix_row_to_relation_set () {
  __row_number=$1

  check_file_exists ${INDEXFILE}
  check_file_exists ${PURGEDFILE}

  __two=2
  __rowp2=$(( __row_number + __two ))

  __line=`head -n ${__rowp2} ${INDEXFILE} | tail -n 1`

  echo "Row number ${__row_number} is the sum of ${__line%% *} relations"
  echo "These relations are:"
  for i in ${__line#* } ; do
    __i_dec=$(( 0x${i%%:*} ))
    __ip2=$(( __i_dec + __two ))

    __coeff=${i#*:}
    if [ "x${__coeff}x" = "x${i}x" ] ; then # there is no ':' in ${i}
      echo -n "  rel ${__i_dec}: " 
    else
      echo -n "  rel ${__i_dec} (coeff=${__coeff}): "
    fi
    zcat ${PURGEDFILE} | head -n ${__ip2} | tail -n 1
  done
}

from_matrix_col_to_ideal_index () {
  __col_number=$1

  check_file_exists ${MATRIXCOLS}

  __rep=`grep "^${__col_number} " ${MATRIXCOLS} | cut -d " " -f 2`

  echo ${__rep}
}

from_ideal_index_to_ideal_p_r () {
  __index=$1

  check_file_exists ${RENUMBERTAB}

  __rep=`grep "^i=${__index} " ${RENUMBERTAB} | cut -d " " -f 3-`

  echo ${__rep}
}

from_ideal_p_r_to_ideal_index () {
  __p=$1
  __r=$2

  check_file_exists ${RENUMBERTAB}

  if [ "x${__r}y" = "xraty" ] ; then
    __line=`grep " p=${__p} rat " ${RENUMBERTAB}`
  else
    __line=`grep " p=${__p} r=${__r} " ${RENUMBERTAB}`
  fi

  echo ${__line} | cut -d " " -f 1 | cut -d "=" -f 2
}

from_ideal_index_to_matrix_col () {
  __index=$1

  check_file_exists ${MATRIXCOLS}

  __rep=`grep " ${__index}$" ${MATRIXCOLS} | cut -d " " -f 1`

  echo ${__rep}
}

from_matrix_col_to_ideal_p_r () {
  __tmp=`from_matrix_col_to_ideal_index "$1"`
  from_ideal_index_to_ideal_p_r "${__tmp}"
}

from_ideal_p_r_to_matrix_col () {
  __tmp=`from_ideal_p_r_to_ideal_index "$1" "$2"`
  from_ideal_index_to_matrix_col "${__tmp}"
}


#### Examples:

# output the relation set corresponding to row 42 in the matrix
from_matrix_row_to_relation_set "42" #42 is in decimal
echo "####"

# output the index of the ideal corresponding to column 42 in the matrix
from_matrix_col_to_ideal_index "42" #42 is in decimal, output is in hexa
echo "####"

# output the column number corresponding to the ideal of index 31
from_ideal_index_to_matrix_col "31" #31 is in hexa, output is in decimal
echo "####"

# output (p,r) corresponding to the ideal of index 31
from_ideal_index_to_ideal_p_r "31" #31 is in hexa, output is in hexa
echo "####"

# output (p,r) corresponding to column 42 in the matrix
from_matrix_col_to_ideal_p_r "42" #42 is in decimal, output is in hexa
echo "####"

# output index of ideal corresponding to (p,r)
from_ideal_p_r_to_ideal_index "3b" "18" # args are in hexa, output is in hexa
from_ideal_p_r_to_ideal_index "13" "rat" # arg is in hexa, output is in hexa
echo "####"

# output the column number of ideal corresponding to (p,r)
from_ideal_p_r_to_matrix_col "3b" "18" # args are in hexa, output is in decimal
from_ideal_p_r_to_matrix_col "13" "rat" # arg is in hexa, output is in decimal
echo "####"
