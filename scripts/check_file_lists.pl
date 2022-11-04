#!/usr/bin/env perl

# This script makes sure that the files.dist and files.nodist cover the 
# set of files in the source tree (as defined by the SCM system of course).
# Note that nothing is there to make sure that the files which are _required_
# for compiling the tree are actually there. That's admittedly unfortunate,
# but doing it the right way would require some cooperation from the build
# system. At the moment, CMake has no functionality to build source package
# (worse, it's got a very lame excuse for it called CPack).

# This script must be called from the top of the tree.

use strict;
use warnings;

die "Please call $0 from the top of the source tree" unless -f "cado.h";

my @files_scm;

open SCM, "git ls-files |";
@files_scm = <SCM>;
close SCM;
die "where am I ?" unless @files_scm;

chomp for @files_scm;

die "No SCM found ???" unless @files_scm;

my $known_lists=[];

my $error=0;

my $allow_rewrite=0;

while (defined($_ = shift @ARGV)) {
    /^-w$/ && do { $allow_rewrite = 1; next; };
    die"Unexpected command line argument: $_";
}

sub add_known_list {
    my $n = shift;
    my $f = {};
    my $d = {};
    my @comments;
    my @data;
    open my $fh, $n or return;
    while (defined($_=<$fh>)) {
        do { push @comments, $_; next; } unless /^[^#]/;
        push @data, $_;
        chomp($_);
        m{/$} && do {
            if (exists $d->{$_}) {
                print STDERR "Pattern $_ appears more than once in $n\n";
                $error++;
            }
            $d->{$_}=0;
            next;
        };
        if (exists $f->{$_}) {
            print STDERR "Pattern $_ appears more than once in $n\n";
            $error++;
        }
        $f->{$_}=0;
    }
    close $fh;
    my @sorted_data = sort { $a cmp $b } @data;
    my $i = 0;
    my $sort_ok = 1;
    while ($sort_ok && $i <= $#sorted_data) {
        if ($sorted_data[$i] ne $data[$i]) {
            $sort_ok=0;
        } else {
            $i++;
        }
    }
    $sort_ok = $sort_ok && $i > $#data;
    if (!$sort_ok) {
        if ($allow_rewrite) {
            print STDERR "$n is not sorted. Sorting in place\n";
            open my $fh, ">$n" or die;
            print $fh $_ for @comments;
            print $fh $_ for @sorted_data;
            close $fh;
        } else {
            print STDERR "$n is not sorted. Run scripts/check_file_lists.pl -w to sort in-place\n";
            $error++;
        }
    }

    if (scalar keys %$d == 0 && scalar keys %$f == 0) {
        return;
    }
    return [ $n, $f, $d ];
}

push @$known_lists, &add_known_list('files.dist');
push @$known_lists, &add_known_list('files.nodist');
if (defined(my $x = &add_known_list('files.unknown'))) {
    warn "Warning: files.unknown is not empty. Please arrange so that every file is covered by a pattern in either files.dist or files.nodist";
    push @$known_lists, $x;
}

for my $file (@files_scm) {
    # Check whether this matches in any way.
    my @matches=();
    for my $k (@$known_lists) {
        my ($n,$f,$d) = @$k;
        if (exists($f->{$file})) {
            $f->{$file}++;
            push @matches, "$n:$file";
        }
        for my $p (keys %$d) {
            if ($file =~ /^$p/) {
                $d->{$p}++;
                push @matches, "$n:$p";
            }
        }
    }
    my $nm = scalar @matches;
    next if ($nm == 1);
    if ($nm == 0) {
        print STDERR "File $file is not covered by any pattern\n";
        $error++;
    } else {
        print STDERR "File $file is covered by $nm patterns: @matches\n";
        $error++;
    }
}

for my $k (@$known_lists) {
    my ($n,$f,$d) = @$k;
    for (keys %$f) {
        next if $f->{$_};
        print STDERR "Pattern $n:$_ matches nothing\n";
        $error++;
    }
    for (keys %$d) {
        next if $d->{$_};
        print STDERR "Pattern $n:$_ matches nothing\n";
        $error++;
    }
}
if ($error) {
    die "$error errors found. Please fix files.dist and files.nodist";
}

