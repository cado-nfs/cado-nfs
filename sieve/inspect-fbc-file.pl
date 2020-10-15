#!/usr/bin/env perl

use strict;
use warnings;

use Fcntl;
use Fcntl 'SEEK_SET';

# This script lists the contents of an fbc file. 
# XXX this should be kept in sync with ./sieve/fb.cpp

sub usage {
    print <<EOF;
Usage: ./inspect-fbc-file.pl [options]

Options understood:
    -fbc <fbc file> (mandatory) fbc file to parse. Multiple files can be
    parsed in a single script call.

EOF
    die @_ if @_;
}

my @fbcs=();

while (defined($_=shift @ARGV)) {
    /^-fbc$/ && scalar @ARGV && do {
        push @fbcs, shift @ARGV;
        next;
    };
    (/^--help$/ || /^-h$/) && do { usage; exit 0; };
    usage "ERROR: Unexpected argument: $_";
}

if (!scalar @fbcs) {
    usage "ERROR: No -fbc argument provided";
}

for my $fbc (@fbcs) {
    sysopen F, "$fbc", O_RDONLY or die "$fbc: $!";
    my $size = (stat $fbc)[7];
    printf "$fbc has size $size (0x%x)\n", $size;
    for(my $pos = 0 ; $pos < $size ; ) {
        sysseek F, $pos, SEEK_SET;
        sysread F, my $header, 4096;

        my ($version, $length, $deg, $poly, $lim, $powlim);
        my ($v_offset, $w_offset, $entrybytes, $nentries);
        my @tokens = split(/\s+/m, $header);
        ($version, @tokens) = @tokens or die;
        print '-' x 70, "\n";
        printf "FBC HEADER BLOCK version $version at file position $pos (0x%x)\n", $pos;
        die unless $version == 1;
        ($length, @tokens) = @tokens or die;
        printf "\tfbc block length $length (0x%x)\n", $length;
        ($deg, $poly, @tokens) = @tokens or die;
        print "\tdegree-$deg polynomial: $poly\n";
        ($lim, $powlim, @tokens) = @tokens or die;
        printf "\tFB bound $lim (0x%x, 2^%.1f)\n", $lim, log($lim)/log(2);
        printf "\tpower limit $powlim (0x%x, 2^%.1f)\n", $powlim, log($powlim)/log(2);
        ($v_offset, $w_offset, $nentries, $entrybytes, @tokens) = @tokens;
        my $ppv = $pos + $v_offset;
        my $ppw = $pos + $w_offset;
        printf "\tGeneral entries ($nentries, $entrybytes bytes each) start at file position $ppv (0x%x) ; corresponding weights start at file position $ppw (0x%x)\n", $ppv, $ppw;
        for(my $i = 0 ; $i <= $deg ; $i++) {
            ($v_offset, $w_offset, $nentries, $entrybytes, @tokens) = @tokens or die;
            $ppv = $pos + $v_offset;
            $ppw = $pos + $w_offset;
            printf "\tEntries with $i roots ($nentries, $entrybytes bytes each) start at file position $ppv (0x%x) ; corresponding weights start at file position $ppw (0x%x)\n", $ppv, $ppw;
        }
        print "\n";

        $pos += $length;
    }

    close F;
}
