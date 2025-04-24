#!/bin/bash

re_grep='^#\(ifndef\|define\|endif\) *\(CADO_\)\?.*_H\(PP\)\?_\>\( \*\/\)\?$'
filter() { grep -v gf2x/ | grep -v embedded/ | grep -v flint-fft/ ; }
inputs() { git grep -l "$re_grep" ; }
perl_code='s/^(#(?:ifndef|define|endif))\s+(?:CADO_)?([A-Z0-9_]*_H(?:PP)?)_?\b(*pla:(?:\s+\/)?$)/$1 CADO_$2/ && s/([A-Z0-9_]{4,}_)\1/$1/;'

inputs | filter | xargs -n1 perl -pe "$perl_code" -i
