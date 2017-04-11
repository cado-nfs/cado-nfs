#!/bin/bash
if [ -f local.sh ]; then
	length=$( wc -l local.sh | cut -d " " -f 1 )
	if [ $length -ne 2 ]; then
		echo "Save previous local.sh in local.sh.old."
		sleep 0.5
		cp local.sh local.sh.old
	fi
fi
echo 'CFLAGS="-O3 -DNDEBUG"' > local.sh
echo 'CXXFLAGS="-O3 -DNDEBUG"' >> local.sh
echo "Compile nfs-hd part of cado-nfs."
sleep 0.5
make makefb-hd special-q_sieve
cp parameters/nfs-hd/p6bd40.poly build/`hostname`/nfs-hd/.
cd build/`hostname`/nfs-hd
echo
echo "Launch makefb-hd (less than 10 minutes on a single 2GHz processor) for the p6bd40."
sleep 0.5
./makefb-hd -t 3 -fbb 524288,524288 -poly p6bd40.poly -lpb 23,23 -out sieve_base_p6bd40
echo "Launch special-q_sieve (less than 1 days on a single 2GHz processor) for the p6bd40."
sleep 0.5
./special-q_sieve -H 6,6,6 -fbb 524288,524288 -thresh 65,65 -poly p6bd40.poly -lpb 23,23 -q_side 1 -fb sieve_base_p6bd40 -q_range 524341,2354820 -out sieve_524341_2354820 -err sieve_524341_2354820.err -gal autom6.1
echo
echo "Check the relations"
sleep 0.5
cd ../../../nfs-hd/check_relations/
make build_tree=../../build/`hostname`/nfs-hd
cd ../../build/`hostname`/nfs-hd/
grep -v "#" sieve_524341_2354820 > relations.raw
sort relations.raw | uniq > relations.uniq
echo "Needed relations: 1128604"
./check_relations relations.uniq relations.true p6bd40.poly 23,23 relations.err
