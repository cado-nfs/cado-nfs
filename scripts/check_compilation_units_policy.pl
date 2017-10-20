#!/usr/bin/env perl

# This script makes sure that for all compilation units found in the
# source tree, we have '#include "cado.h"' near the top, as we should.
#
# Of course we have exceptions. See @path_exceptions below.
#       gf2x/
#       linalg/bwc/mpfq/
#       linalg/bwc/flint-fft/
# This script must be called from the top of the tree.
#
# 
# using vim, you may do 
#       :set makeprg=./scripts/check_compilation_units_policy.pl
# and cycle through the errors as if they were produced by gcc.

use strict;
use warnings;

die "Please call $0 from the top of the source tree" unless -f "cado.h";

my @path_exceptions=qw|
        gf2x/
        linalg/bwc/mpfq/
        linalg/bwc/flint-fft/
        config/
        misc/
        |;

my @all_files = `git ls-files`;

my $err=0;
FILE: for my $f (@all_files) {
    $f =~ s/^\s*//;
    $f =~ s/\s*$//;
    next unless $f =~ /\.[ch](?:pp)?$/;
    for my $p (@path_exceptions) {
        next FILE if $f =~ /^$p/;
    }
    my $is_header = ($f =~ /\.h(?:pp)?$/);

    open F, $f;
    
    my $contents = eval { undef $/; <F>; };
    $contents =~ s#/\*.*?\*/##sg;

    my @lines = split(/^/, $contents);
    my @numbered_lines = map { [1+$_, $lines[$_]] } (0..$#lines);

    my @includes = map {
            my ($lnum, $text) = @$_;
            chomp($text);
            my @x=();
            $text =~ /^\s*#\s*include\s*(\S+)/ && push @x,[$lnum,$1];
            @x;
        } @numbered_lines;

    if ($is_header) {
        my @include_cado = grep { $_->[1] =~ /cado\.h/ } @includes;
        if (@include_cado) {
            my $lnum = $include_cado[0]->[0];
            print STDERR "$f:$lnum: is a header file, it must not include cado.h\n";
            $err++;
        }
    } else {
        my $first = shift @includes;
        my @include_cado = grep { $_->[1] =~ /cado\.h/ } @includes;
        if ($first) {
            if ($first->[1] =~ /cado\.h/) {
                if ($first->[1] !~ /"cado\.h"/) {
                    print STDERR "$f:$first->[0]: cado.h should be included as \"cado.h\", not <cado.h>\n";
                    $err++;
                }
            } else {
                print STDERR "$f:$first->[0]: is a compilation unit, its first include file must be \"cado.h\"\n";
                $err++;
            }
            if (@include_cado) {
                my $lnum = $include_cado[0]->[0];
                print STDERR "$f:$lnum: there is no point in including cado.h twice\n";
                $err++;
            }
        } else {
            print STDERR "$f:1: is a compilation unit, it must include \"cado.h\"\n";
            $err++;
        }
    }
}
if ($err) {
    print STDERR "$err errors found. Please fix\n";
    exit 1;
}
