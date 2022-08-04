#!/usr/bin/env perl

use strict;
use warnings;

# Usage:
#   matrix_bin2rows.pl < "$base.bin" > "$wdir/a.rows" 3> "$wdir/a.row_offsets" 4> "$wdir/a.bin"
#
# This is very slow, and intended for debugging only. Do not use in
# production, of course.

binmode(STDIN);
binmode(STDOUT);
my $r = 0;
open(R, ">&3"); binmode(R);
open(B, ">&4"); binmode(B);
while (sysread(STDIN, my $x, 4)) {
    my $w = unpack("L", $x);
    syswrite(B, $x, 4);
    $r += $w;
    syswrite(R, pack("L", $r), 4);
    my @s = ();
    for(my $i = 0 ; $i < $w ; $i++) {
        sysread(STDIN, my $x, 4);
        push @s, unpack("L", $x);
    }
    # Sort the .bin input, and write the sorted data to fd 4
    for my $j (sort { $a <=> $b } @s) {
        $x = pack("L", $j);
        syswrite(STDOUT, $x, 4);
        syswrite(B, $x, 4);
    }
}
