#!/usr/bin/env perl

use strict;
use warnings;

my $impl = shift @ARGV or die "Usage: $0 <impl name>";

(my $impl_filebase = $impl) =~ s/_/-/g;

for ("${impl_filebase}.h", "template.h") {
    die "No header file $_ found" unless -f $_;
}

open T, "template.h";

# First fetch the code blobs from the template.

my @blobs=();
my $cblob = undef;

while (defined($_=<T>)) {
    if (/BEGIN SECTION \d+: (.*) \*\//) {
        die if $cblob;
        push @blobs, [ $1, "" ];
        $cblob = \$blobs[$#blobs]->[1];
    } elsif (/END SECTION \d+/) {
        $cblob = undef;
    } elsif ($cblob) {
        $$cblob .= $_;
    }
}
close T;

print "Found ", scalar @blobs, " text blobs:\n";
print "\t$_->[0]\n" for (@blobs);

open F, "${impl_filebase}.h";
open G, ">${impl_filebase}.h.tmp";

my $blobnum = 0;
my $ttext;
my $ttitle;

while (defined($_=<F>)) {
    if (/The section below is automatically generated/) {
        my $x = shift @blobs;
        ($ttitle,$ttext) = @$x;
        $ttitle =~ s/XXX/$impl/g;
        print "Inserting at line $.: $ttitle\n";
        print G;
        defined($_=<F>) or die;
        my $params={ inline=>[], pod=>[] };
        while (/^\s*\/\* (inline|pod):((?:\s+\w+)+)\s\*\/\s*$/) {
            push @{$params->{$1}}, split(' ', $2);
            print G;
            defined($_=<F>) or die;
        }
        if ($blobnum == 0) {
            # print "Performing inline modifiations for ${impl}_info ; [" . join(", ", @{$params->{'inline'}}) . "\n";
            for my $func (@{$params->{'inline'}}) {
                $ttext =~ s/^((?:\w+) XXX_info_${func})\b/static inline $1/mg;
            }
            # sigh... https://gcc.gnu.org/ml/gcc-help/2014-06/msg00048.html
            $ttext =~ s/^((?!typedef)(?:\w+(?:\s*\*)*)) (XXX_info_\w+)\b([^;]*);/extern $1 $2$3 GF2X_FFT_EXPORTED;/mg;
        } elsif ($blobnum == 1) {
            # print "Performing inline modifiations for ${impl}_info ; [" . join(", ", @{$params->{'inline'}}) . "\n";
            for my $func (@{$params->{'inline'}}) {
                $ttext =~ s/^((?:\w+) XXX_${func})\b/static inline $1/mg;
            }
            $ttext =~ s/^((?!typedef)(?:\w+(?:\s*\*)*)) (XXX_\w+)\b([^;]*);/extern $1 $2$3 GF2X_FFT_EXPORTED;/mg;
        } elsif ($blobnum == 2 && @{$params->{'pod'}} && $params->{'pod'}->[0] eq 'yes') {
            # print "Performing pod modifiations for ${impl}_info\n";
            $ttext =~ s/^(\s*(?:inline\s*)?~XXX_info)\(\)[^\}]*\}/$1() = default;/sm;
            $ttext =~ s/^(\s*(?:inline\s*)?XXX_info)\((XXX_info const &) o\)[^\}]*\}/$1($2) = default;/sm;
            $ttext =~ s/^(\s*(?:inline\s*)?XXX_info)\(\)[^\}]*\}/$1() = default;/sm;
            $ttext =~ s/^(\s*(?:inline\s*)?XXX_info& operator=)\((XXX_info const &) o\)[^\}]*\}/$1($2) = default;/sm;
        }
        $ttext =~ s/XXX/$impl/g;

        print G $ttext;
        while (!/End of automatically generated section/) {
            defined($_=<F>) or die;
        }
        print G;
        $blobnum++;
        next;
    }
    print G;
}
close F or die;
close G or die;
rename "${impl_filebase}.h.tmp", "${impl_filebase}.h";
