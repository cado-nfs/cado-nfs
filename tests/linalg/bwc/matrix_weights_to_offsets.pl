#!/usr/bin/env perl

use strict;
use warnings;

# Usage: matrix_weights_to_offsets.pl < "$base.cw.bin" > "$wdir/a.col_offsets"

open(CW, "<&STDIN");
open(C, ">&STDOUT");
binmode(CW);
binmode(C);
my $c = 0;
while (sysread(CW, my $x, 4)) {
    my $w = unpack("L", $x);
    $c += $w;
    syswrite(C, pack("L", $c), 4);
}
