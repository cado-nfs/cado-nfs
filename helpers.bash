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
git grep -l  '<\(std.*\|float\|limits\|math\|string\|ctype\|inttypes\|errno\)\.h>' '*.[ch]pp' | xargs -rn1 sed -e 's/<\(std.*\|float\|limits\|math\|string\|ctype\|inttypes\|errno\)\.h>/<c\1>/g' -i
}
