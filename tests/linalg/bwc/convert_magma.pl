#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

# Converts one bwc file as given on stdin to magma-parseable text on
# stdout

my $mode = shift @ARGV;

if ($mode eq 'matrix') {
    die unless defined($_=<>);
    die unless /^(\d+) (\d+)$/;
    my $nr=$1;
    my $nc=$2;
    my $i = 0;
    my $ne = 0;
    my $needcomma = 0;
    print "var:=[\n";
    while (<>) {
        s/^(\d+)// or die;
        my $l = $1;
        if ($l == 0) { $i++; next; }
        print ",\n" if $needcomma;
        $needcomma=0;
        for(my $j = 0;$j < $l;$j++) {
            s/^\s+(\d+)// or die;
            print ", " if $needcomma;
            $ne++;
            my $i1 = $i + 1;
            my $j1 = $1 + 1;
            print "<$i1,$j1,1>";
            $needcomma=1;
        }
        die unless /^\s*$/;
        $i++;
    }
    print "];\n";
    print "var:=#var eq 0 select SparseMatrix($nr,$nc) else SparseMatrix($nr,$nc,var);\n";
    if ($i != $nr) {
        print STDERR "Matrix nrows was wrong (read $i, expected $nr)\n";
    }
    exit;
}

if ($mode eq 'bmatrix') {
    my $nr=0;
    my $nc=0;
    my $curr=0;
    my @coeffs=();
    while(sysread(STDIN, my $x, 4)) {
        my $v=unpack("L",$x);
        if ($curr-- == 0) {
            $nr++;
            $curr=$v;
        } else {
            $nc = $v+1 if $v >= $nc;
            push @coeffs, [$nr, $v+1];
        }
    }
    print "var:=[\n";
    print join(", ", map { "<$_->[0], $_->[1],1>" } @coeffs);
    print "];\n";
    print "var:=#var eq 0 select SparseMatrix($nr,$nc) else SparseMatrix($nr,$nc,var);\n";
    die unless $curr==0;
    exit;
}

if ($mode =~ /^bpmatrix(?:_(\d+)_(\d+))?$/) {
    my $nr0=$1;
    my $nc0=$2;
    my $nr=0;
    my $nc=0;
    my $curr=0;
    my @coeffs=();
    while(sysread(STDIN, my $x, 4)) {
        my $v = unpack("L", $x);
        if ($curr-- == 0) {
            $nr++;
            $curr=$v;
        } else {
            $nc = $v + 1 if $v >= $nc;
            sysread(STDIN, my $x, 4) or die;
            my $w=unpack("l", $x);
            push @coeffs, [$nr, $v+1, $w];
        }
    }
    if (defined($nr0) && defined($nc0)) {
        if ($nr > $nr0 || $nc > $nc0) {
            die "Unexpected: $mode is incompatible with seeing $nr rows and $nc cols";
        }
        $nr = $nr0;
        $nc = $nc0;
    }
    print "var:=[\n";
    print join(", ", map { "<$_->[0], $_->[1], $_->[2]>" } @coeffs);
    print "];\n";
    print "var:=#var eq 0 select SparseMatrix($nr,$nc) else SparseMatrix($nr,$nc,var);\n";
    die unless $curr==0;
    exit;
}

if ($mode eq 'balancing') {
    sysread(STDIN, my $x, 8);
    my ($zero, $magic) = unpack("L2", $x);
    die unless $zero == 0;
    die unless $magic == 0xba1a0000;
    sysread(STDIN, $x, 40);
    my ($nh,$nv,$nr,$nc,$nzr,$nzc,$ncoeffs,$checksum,$flags) = unpack("L6QLL", $x);
    sysread(STDIN, $x, 8);
    my ($pa, $pb) = unpack("L2", $x);
    sysread(STDIN, $x, 8);
    my ($pai, $pbi) = unpack("L2", $x);
    $checksum = sprintf("%04x", $checksum);
    my $txflags="";
    my $colperm = $flags & 1;
    my $rowperm = $flags & 2;
    $txflags .= ", colperm" if $colperm;
    $txflags .= ", rowperm" if $rowperm;
    $txflags .= ", replicate" if $flags & 8;
    # pad.
    my $pad = sub {
        my ($n,$K,$b)=@_;
        my $x = int(($n + $K-1)/$K);
        while ($x % $b) { $x++; }
        return $x*$K;
    };
    my $tr = &$pad($nr,$nh*$nv,8);
    my $tc = &$pad($nc,$nh*$nv,8);
    if ($flags & 8) {
        $tr = $tc = $tr > $tc ? $tr : $tc;
    }
    print "nr:=$tr; // originally $nr\n";
    print "nc:=$tc; // originally $nc\n";
    print "nzr:=$nzr;\n";
    print "nzc:=$nzc;\n";
    print "nr_orig:=$nr;\n";
    print "nc_orig:=$nc;\n";
    my $s = $nh * $nv;
    while ($tr % $s) { $tr++; }
    while ($tc % $s) { $tc++; }
    print "// $nr rows $nc cols, split ${nh}x${nv}, checksum $checksum$txflags\n";
    print "tr:=$tr; // originally $nr\n";
    print "tc:=$tc; // originally $nc\n";
    print "nh:=$nh;\n";
    print "nv:=$nv;\n";
    print "preshuf:=func<x|x le $nc select 1+(($pa*(x-1)+$pb) mod $nc) else x>;\n";
    print "preshuf_inv:=func<x|x le $nc select 1+(($pai*(x-1)+$pbi) mod $nc) else x>;\n";
    if ($rowperm) {
        my @p=();
        for(my $i = 0 ; $i < $tr ; $i++) {
            die unless sysread(STDIN, my $x, 4);
            push @p, 1+unpack("L",$x);
        }
        print "rowperm:=[", join(", ", @p), "];\n";
    }
    if ($colperm) {
        my @p=();
        for(my $i = 0 ; $i < $tc ; $i++) {
            die unless sysread(STDIN, my $x, 4);
            push @p, 1+unpack("L",$x);
        }
        print "colperm:=[", join(", ", @p), "];\n";
    }
    exit;
}


if ($mode =~ /^(permutation|weights|pvector32)$/) {
    # Dump a list of unsigned ints
    my $add1 = $mode eq 'permutation';
    my @p=();
    while(sysread(STDIN, my $x, 4)) {
        my $v = unpack("L",$x);
        push @p, $v+$add1;
    }
    print "var:=[",join(',',@p),"];\n";
    exit;
}

if ($mode =~ /^(spvector32)$/) {
    # Dump a list of SIGNED ints
    my $add1 = $mode eq 'permutation';
    my $delim="var:=[";
    while(sysread(STDIN, my $x, 4)) {
        my $v = unpack("l",$x);
        print $delim, $v+$add1;
        $delim=", ";
    }
    print "];\n";
    exit;
}

if ($mode =~ /^(spvector64)$/) {
    # Dump a list of SIGNED ints
    my $add1 = $mode eq 'permutation';
    my @p=();
    while(sysread(STDIN, my $x, 8)) {
        my $v = unpack("q",$x);
        push @p, $v+$add1;
    }
    print "var:=[",join(',',@p),"];\n";
    exit;
}

if ($mode eq 'x') {
    die unless defined($_=<>);
    die unless /^(\d+)$/;
    my $nx = 1;
    my @p=();
    while(<>) {
        my @xs = split ' ', $_;
        @xs = map { $_+1; } @xs;
        my $xstring=join(",",@xs);
        push @p,"[$xstring]";
    }
    print "var:=[",join(',',@p),"];\n";
    exit;
}

if ($mode eq 'vector') {
#    # Dump a list of uint64_t's ; this version is robust on 32-bits.
#    my @p=();
#    while(sysread(STDIN, my $x, 8)) {
#        my @v = unpack("L2",$x);
#        my $sv = join(',',@v);
#        my $v = "Seqint([$sv],2^32)";
#        push @p, $v;
#    }
#    print "var:=[",join(',',@p),"];\n";
#    exit;
#   get rid of 32-bit clutter.
    my @p=();
    while(sysread(STDIN, my $x, 8)) {
        push @p, unpack("Q",$x); 
    }
    print "var:=[",join(',',@p),"];\n";
    exit;
}

