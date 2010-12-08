#!/bin/sh

# run as
# ssh gcc64 < gcc-script.sh

remote_cmake_path=/home/zimmerma/bin
cado_tree=$HOME/cado-nfs-1.0-rc3

if [ -x $remote_cmake_path/cmake ] ; then
    echo "$remote_cmake_path/cmake found, good"
else
    echo "$remote_cmake_path/cmake missing" >&2
    exit 1
fi

if [ -d $cado_tree ] ; then
    echo "Clearing old $cado_tree"
    /bin/rm -rf $cado_tree
fi

echo "Checking out remote tree"
#svn co -q svn://scm.gforge.inria.fr/svn/cado-nfs/trunk $cado_tree
wget -q http://www.loria.fr/~zimmerma/cado-nfs-1.0-rc3.tar.gz -O- | tar xzf -

cd $cado_tree

export PATH=$remote_cmake_path:$PATH
echo "Compiling cado-nfs"
make -j 2
echo "Running factor.sh"
exec ./factor.sh 90377629292003121684002147101760858109247336549001090677693 -t 2
