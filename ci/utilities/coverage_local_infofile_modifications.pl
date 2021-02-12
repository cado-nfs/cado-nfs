#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd qw/getcwd realpath/;

# this is absolute.
my $pwd = getcwd();

my $discard = 0;
my $kept = 0;

my $bdir;

if ($#ARGV > 1 && $ARGV[0] eq '-d') {
    shift @ARGV;
    $bdir = shift @ARGV;
}

my $f = shift @ARGV;
open F, $f or die "$f: $!";
open G, ">$f.tmp" or die "$f.tmp: $!";

sub print_record {
    do { $discard++; return; } if $_[0] =~ m{^SF:/} ||  $_[1] =~ m{^SF:/};
    do { $discard++; return; } if $_[1] =~ m{^SF:.*<built-in>$};
    do { $discard++; return; } if $_[1] =~ m{^SF:.*CompilerId.c(?:pp)?$};
    $kept++;
    print G $_ for @_;
}

my @current;

while (defined($_=<F>)) {
    # make paths relative
    chomp($_);
    s,^SF:$pwd/,SF:,g;
    (my $file = $_) =~ s/^SF:(.*)$/$1/;
    if (defined($bdir) && ! -f $file && -l "$bdir/gf2x/$file" && -e "$bdir/gf2x/$file") {
        my $rpath = realpath("$bdir/gf2x/$file");
        $rpath =~ s,^$pwd/?,,;
        print STDERR "Rewrite path $file --> $bdir/gf2x/$file --> $rpath\n";
        $_="SF:$rpath"; # relative
    }
    # dereference symlinks found in build directory
    if (/^SF:(build\/.*)$/) {
        my $f = (-e $1) ? realpath($1) : $1;   # absolute
        $f =~ s,^$pwd/?,,;
        $_="SF:$f"; # relative
    }
    push @current, $_ . "\n";
    if (/end_of_record/) {
        print_record @current;
        @current=();
    }
}

close F;
close G;
rename "$f.tmp", $f or die "rename: $!";

if ($discard) {
    print STDERR "Discarded $discard entries, kept $kept\n";
}
