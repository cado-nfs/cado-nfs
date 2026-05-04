move_to_cpp() {
    for f in "$@" ; do
        git mv $f ${f}pp
        grep "^$f\$" files.dist && sed -e "s,^$f\$,${f}pp,g" -i files.dist
        grep "^$f\$" files.nodist && sed -e "s,^$f\$,${f}pp,g" -i files.nodist
        b=$(basename "$f")
        find . -type f -name CMakeLists\*txt | \
            xargs -r grep -l "\\b$b\\b"  | \
            xargs -r -n 1 sed -e "s,\\b$b\\b,${b}pp,g" -i
        if [[ $b =~ ^(.*)\.h$ ]] ; then
            # for headers, we want to chase all possible uses
            git grep -l "\\b$b\\b"  | \
                grep -v files.dist | \
                xargs -r -n 1 sed -e "s,\\b$b\\b,${b}pp,g" -i
        fi
    done
    git commit -am "convert $* to c++"
}

rename_cpp_protectors() {
    git grep -l '_H_\?\b' '*.hpp'  | xargs -n1 sed -e 's/_H\b/_HPP/g' -i
}

rename_all_c_headers() {
git grep -l  '<\(std.*\|float\|limits\|math\|string\|time\|ctype\|inttypes\|errno\)\.h>' '*.[ch]pp' | xargs -rn1 sed -e 's/<\(std.*\|float\|limits\|math\|string\|time\|ctype\|inttypes\|errno\)\.h>/<c\1>/g' -i
}

offline_validate() {
    X=$PWD
    B=`git rev-parse --abbrev-ref HEAD`
    R=`git rev-parse --short HEAD`
    (cd /tmp ; git clone $X -b $B cado-nfs-$B-$R && cp -pf $X/local.sh  /tmp/cado-nfs-$B-$R)
    if [ "$1" ] ; then
        rsync -a "/tmp/cado-nfs-$B-$R/" "$1":"/tmp/cado-nfs-$B-$R/"
        ssh "$1" "cd /tmp/cado-nfs-$B-$R ; CLANG=1 make -j16 all all_test_dependencies  && CLANG=1 DOCKER_SAGEMATH_SILENT=1 DOCKER_SAGEMATH_NO_PULL=1 make check ARGS='-E builddep -j8'"
    else
        (cd /tmp/cado-nfs-$B-$R ; CLANG=1 make -j8 all all_test_dependencies && CLANG=1 DOCKER_SAGEMATH_SILENT=1 DOCKER_SAGEMATH_NO_PULL=1 make check ARGS='-E builddep -j4')
    fi
    if [ $? = 0 ] ; then
        rc="0 (ok)"
    else
        rc="$? (NOK NOK NOK)"
    fi
    echo "# Done checking commit $R on branch $B, return code is $rc"
}
