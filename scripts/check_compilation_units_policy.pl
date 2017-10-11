#!/usr/bin/env perl

# This script makes sure that for all compilation units found in the
# source tree, we have '#include "cado.h"' near the top, as we should.
#
# Of course we have exceptions. See @path_exceptions below.
#       gf2x/
#       linalg/bwc/mpfq/
#       linalg/bwc/flint-fft/
# This script must be called from the top of the tree.

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

    my @includes = map {
            chomp($_);
            my @x=();
            /^\s*#\s*include\s*(\S+)/ && push @x,$1;
            @x;
        } (split(/^/, $contents));
    if ($is_header) {
        if (grep /cado\.h/, @includes) {
            print STDERR "$f is a header file, it must not include cado.h\n";
            $err++;
        }
    } else {
        my $first = shift @includes;
        if ($first && $first !~ /cado\.h/) {
            print STDERR "$f is a compilation unit, its first include file must be cado.h\n";
            $err++;
        } elsif (grep /cado\.h/, @includes) {
            print STDERR "there is no point in including cado.h twice\n";
            $err++;
        }
    }
}
if ($err) {
    print STDERR "$err errors found. Please fix\n";
    exit 1;
}
