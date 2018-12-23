#!/usr/bin/env bash

if [ ! -f cado.h ]
then
	echo "Script must be launched at the root path of CADO-NFS."
	echo "Something like bash nfs-hd/automatic_p6bd40.sh must work."
	exit
fi

PATH_BUILD=build/`hostname`
PATH_BUILD_NFS_HD=$PATH_BUILD/nfs-hd

echo "Modify CMakeLists.txt"
sed -i "s/# \(add_subdirectory (nfs-hd)\)/\1/g" CMakeLists.txt

if [ -f local.sh ]; then
	echo "Save previous local.sh in local.sh.old."
	cp local.sh local.sh.old
fi
echo 'CFLAGS="-O3 -DNDEBUG"' > local.sh
echo 'CXXFLAGS="-O3 -DNDEBUG"' >> local.sh

echo "Compile nfs-hd part of cado-nfs."
make cmake makefb-hd special-q_sieve

echo "Copy parameter file."
cp parameters/nfs-hd/p6bd40.poly $PATH_BUILD_NFS_HD/.
cd $PATH_BUILD_NFS_HD

echo "Launch makefb-hd (less than 10 minutes on a single 2GHz processor) for the p6bd40."
./makefb-hd -t 3 -fbb 524288,524288 -poly p6bd40.poly -lpb 23,23 -out sieve_base_p6bd40

echo "Launch special-q_sieve (less than 1 days on a single 2GHz processor) for the p6bd40."
./special-q_sieve -H 6,6,6 -fbb 524288,524288 -thresh 65,65 -poly p6bd40.poly -lpb 23,23 -q_side 1 -fb sieve_base_p6bd40 -q_range 524341,2354820 -out sieve_524341_2354820 -err /dev/null -gal autom6.1

cd ../../../nfs-hd/check_relations/
make build_tree=../../$PATH_BUILD_NFS_HD HEADERS_CADO_H=../.. HEADERS_CADO_CONFIG_H=../../$PATH_BUILD
cd ../../$PATH_BUILD_NFS_HD

echo "Check the relations."
grep -v "#" sieve_524341_2354820 > relations.raw
sort relations.raw | uniq > relations.uniq

echo "Needed relations: 1128604."
./check_relations relations.uniq relations.true p6bd40.poly 23,23 relations.err

cd ../../../
if [ -f local.sh.old ]; then
	echo "Restore previous local.sh."
	mv local.sh.old local.sh
fi
echo "Restore previous CMakeLists.txt"
sed -i "s/\(add_subdirectory (nfs-hd)\)/# \1/g" CMakeLists.txt
