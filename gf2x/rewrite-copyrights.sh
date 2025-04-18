#!/bin/bash
files=($(git ls-files '*.[ch]*' | xargs grep -l Copyright))
for f in "${files[@]}" ; do
    if [ -L "$f" ] ; then continue ; fi
    years=($(git log --format=format:%as "$f" | cut -d- -f1 | sort -nu))
    read -d  -r code <<EOF
    next unless /(?<=Copyright\\s)\\s*\\d{4}(?:-\\d{4})?(,\\s+\\d{4}(?:-\\d{4})?)*$/;
    my %years=();
    for my \$x (qw/${years[*]}/) {
        \$years{\$x}=1
    }
    while (s/(?<=Copyright\\s)\\s*(\\d{4})(?:-(\\d{4}))?(,\\s+)?//) {
        if (defined(\$2)) {
            for my \$x (\$1..\$2) {
                \$years{\$x}=1
            }
        } else {
            \$years{\$1}=1
        }
    }
    my \$allyears = join(", ", sort { \$a <=> \$b } keys %years);
    s/(?<=Copyright\\s)/\$allyears/g;
EOF
    perl -pe "$code" -i "$f"
done
