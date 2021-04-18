#!/usr/bin/env perl

use strict;
use warnings;
use JSON;
use Data::Dumper;

# This script can be used to split a given relation file, or set of
# concatenated relation files, into a hierarchy of directories for easier
# random access. This is particularly useful fr test filtering runs based
# on an existing computation, notably to experiment with the impact of
# the shrink factor.
#
# Resulting files are stored uncompressed.
#
# A metadata file dirinfo.json is written, and reflects the configuration.
# 
# Unless otherwise stated by the command-line option -a, the output
# directory must be either an existent empty directory, or not exist at
# all. With -a, an existing directory also works, provided that it has a
# compatible dirinfo.json file.
#
# Usage example:
#
# localhost ~ $ scripts/estimate_matsize/build_relation_cache.pl -o /tmp/split -s 1000,10000 /tmp/cado*/*.upload/*gz
# localhost ~ $ ls  /tmp/split/
# 08  09  10  11  12  13  14
# localhost ~ $ ls  /tmp/split/11
# 00000-5000   10000-15000  20000-25000  30000-35000  40000-45000  50000-55000  60000-65000  70000-75000  80000-85000  90000-95000
# 05000-10000  15000-20000  25000-30000  35000-40000  45000-50000  55000-60000  65000-70000  75000-80000  85000-90000  95000-100000
#
#
# Note about append mode (-a): all files passed on the command line are
# processed no matter what. So if you use append, for performance your
# best bet is to make sure that you pass distinct sets of files to the
# different passes.

my $json = JSON->new->pretty(1);

my $outputdir = "./split";
my @splits = (10000);
my $append_ok;

sub deep_compare
{
    my ($u,$v)=@_;
    return 0 if (!defined($u) && defined($v));
    return 0 if (defined($u) && !defined($v));
    return 1 if (!defined($u) && !defined($v));
    return 0 if ref $u ne ref $v;
    if (ref $u eq 'HASH') {
        my @ku = sort { $a cmp $b } keys %$u;
        my @kv = sort { $a cmp $b } keys %$v;
        return 0 if scalar @ku != scalar @kv;
        do { return 0 if $ku[$_] ne $kv[$_] } for (0..$#ku);
        do { return 0 unless deep_compare($u->{$_}, $v->{$_}); } for @ku;
    } elsif (ref $u eq 'ARRAY') {
        return 0 if scalar @$u != scalar @$v;
        do { return 0 if $u->[$_] ne $v->[$_] } for (0..$#$u);
    } else {
        return $u eq $v;
    }
    return 1;
}

while (defined($_ = $ARGV[0])) {
    /^(-o|--output-dir(?:ectory)?)$/ && do {
        shift @ARGV;
        $outputdir = shift @ARGV;
        next;
    };
    /^(-s|--splits)$/ && do {
        shift @ARGV;
        @splits = split(',',shift @ARGV);
        for my $x (@splits) {
            die "Split \"1\" doesn't really make sense" if $x == 1;
            if (length($x) == length($x-1)) {
                die "Only splits at powers of 10 are supported\n";
            }
        }
        next;
    };
    /^(-h|--help)$/ && do {
        my $s = join(',',@splits);
        print STDERR <<EOF;
Usage: $0 [options] [files...]
Splits the input files in the split directory.
according to the number of split digits.
Options:
    -o, --output-dir, --output-directory (default: $outputdir)
        place files under this directory (the directory must exist)
    -s, --splits (default: $s)
        comma separated numbers of items (see below).

The splitting rule is as follows. Within the files that are read (which
possibly contain several concatenated las runs), every line that matches
the (perl) regexp:
    ^#.*-q0 (\\d+).*-q1 (\\d+)
starts a new output file, whose location in the directory hierarchy below
the output directory depends on the value q0.  The length of the list
passed to -s determines the depth of the directory tree in which output
files are stored. An empty -s string means that all output files are
created right under the output directory. If the list passed to -s has n
values, then the depth below the output directory is n.
If the list of splits is specified as n2,n1, then a part of the output
that has q0=Q and q1=Q' goes into the file
    N2/N1/N0-N'0,
with
    (N2*n2+N1)*n1+N0 = Q 
and (N2*n2+N1)*n1+N'0 = Q'.
Therefore the qrange must be a divisor of n1. Note that the program
refuses to create more than 1000 directories at a given level, which
means that all numbers in the -s list except the last one must be at most
1000.
EOF
        exit 0;
    };
    /^-a$/ && do { $append_ok=1; shift @ARGV; next; };
    last;
}

sub is_empty_directory {
    opendir my $dh, $outputdir or die "$outputdir: $!";
    my @x = grep { !/^\.{1,2}$/ } readdir $dh;
    closedir $dh;
    return scalar @x == 0;
}

@splits = map { 0 + $_ } @splits;

my $dirinfo = {   
    "splits" => \@splits,
    "compression" => "none",
};


if (!-d $outputdir || is_empty_directory) {
    if (!-d $outputdir) {
        mkdir $outputdir or die "$outputdir: $!";
    }
    open F, ">$outputdir/dirinfo.json" or die "$outputdir/dirinfo.json: $!";
    print F $json->encode($dirinfo);
} else {
    die "$outputdir is not empty" unless $append_ok;
    open F, "$outputdir/dirinfo.json" or die "$outputdir/dirinfo.json: $!";
    my $j = eval { local $/=undef; <F>; };
    unless (deep_compare($dirinfo, $json->decode($j))) {
        die "cache configuration in $outputdir/dirinfo.json is not consistent with command line, cannot append\n" .
        Dumper({ commandline => $dirinfo, ondisk => $json->decode($j) }) . "\n";
    }
}

for my $filename (@ARGV) {
    my $openpath = $filename;

    if ($filename =~ /\.gz$/) {
        $openpath = "zcat $filename |";
    } elsif ($filename =~ /\.bz2$/) {
        $openpath = "bzcat $filename |";
    } elsif ($filename =~ /\.xz$/) {
        $openpath = "xzcat $filename |";
    };

    open IN, $openpath or die "$filename: $!";

    my $fh;

    while (defined($_=<IN>)) {
        /^# .*-q0 (\d+).*-q1 (\d+)/ && do {
            close $fh if defined $fh;
            my $q0 = $1;
            my $q1 = $2;
            my $oq0 = $1;
            my $oq1 = $2;
            my @qq0;
            my @formats = ();
            for my $i (0..$#splits) {
                my $nn = $splits[$#splits-$i];
                unshift @qq0, $q0 % $nn; $q0 = int($q0/$nn);
                unshift @formats, "%0" . length($nn-1) . "d";
            }
            die "$oq0 is too large for current settings" if $q0;
            my $qq1 = sprintf $formats[$#qq0], $qq0[$#qq0]+($oq1-$oq0);
            for my $i (0..$#splits) {
                $qq0[$i] = sprintf $formats[$i], $qq0[$i];
            }
            my $d = "";
            for my $i (0..$#qq0-1) {
                die "Refusing to create more than 1000 subdirs below a given level" if $qq0[$i] >= 1000;
                $d .= "/$qq0[$i]";
                mkdir "$outputdir$d" unless -d "$outputdir/$d";
            }
            my $f = "$outputdir$d/$qq0[$#qq0]-$qq1";
            open $fh, ">$f" or die "$f: $!";
        };
        die "Found data with no active file to write to: $_\n" unless $fh;
        print $fh $_;
    }
    close $fh if defined $fh;

    close IN;
}
