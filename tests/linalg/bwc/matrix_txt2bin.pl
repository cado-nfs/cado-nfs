#!/usr/bin/env perl

use strict;
use warnings;

# This very silly script takes a matrix in txt format, and writes the
# .rows .row_offsets .col_offsets and .bin file
#
# The matrix must be provided as a file name ending in .txt

sub usage {
    print STDERR <<EOF;
Usage: $0 [matrix file name]

matrix file name should be something like foo.txt

The files foo.bin foo.rows foo.row_offsets and foo.col_offsets are created
EOF
    print STDERR "\nError: @_\n" if @_;
    exit 1;
}

my $filename = shift @ARGV or usage;

$filename =~ /^(.*)\.txt$/ or usage "$filename does not end in .txt";

my $filename_rows = "$1.rows";
my $filename_row_offsets = "$1.row_offsets";
my $filename_col_offsets = "$1.col_offsets";
my $filename_bin = "$1.bin";

for my $f ($filename_rows, $filename_row_offsets, $filename_col_offsets, $filename_bin) {
    usage "$f already exists" if -e $f;
}

open F, "$filename" or die "$!";
open D, ">$filename_rows" or die "$!";
open R, ">$filename_row_offsets" or die "$!";
open C, ">$filename_col_offsets" or die "$!";
open B, ">$filename_bin" or die "$!";

binmode(D);
binmode(R);
binmode(C);
binmode(B);

my %c = ();
my $z = 0;
my $jmax = 0;
while (defined($_ = <F>)) {
    chomp($_);
    my @s = split(' ', $_);
    my $w = shift @s;
    die "$filename: wrong line $.: $_ [[$w $#s]]" unless $w == scalar @s;
    $z += $w;
    syswrite(R, pack("L", $z), 4);
    syswrite(B, pack("L", $w), 4);
    for my $j (@s) {
        syswrite(D, pack("L", $j), 4);
        syswrite(B, pack("L", $j), 4);
        $c{$j}+=1;
        $jmax = $j if $j >= $jmax;
    }
}

$z = 0;
for my $j (0..$jmax) {
    $z += $c{$j} || 0;
    syswrite(C, pack("L", $z), 4);
}
