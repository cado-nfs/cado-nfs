#!/usr/bin/env perl
#  This file is part of the gf2x library.
#
#  Copyright 2007, 2008, 2009, 2010, 2012, 2013, 2015, 2016
#  Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of either:
#   - If the archive contains a file named toom-gpl.c (not a trivial
#     placeholder), the GNU General Public License as published by the
#     Free Software Foundation; either version 3 of the License, or (at
#     your option) any later version.
#   - If the archive contains a file named toom-gpl.c which is a trivial
#     placeholder, the GNU Lesser General Public License as published by
#     the Free Software Foundation; either version 2.1 of the License, or
#     (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
#  
#  You should have received a copy of the GNU General Public License as
#  well as the GNU Lesser General Public License along with this program;
#  see the files COPYING and COPYING.LIB.  If not, write to the Free
#  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
#  02110-1301, USA.


use warnings;
use strict;

sub usage {
    die "Usage: ./tune-lowlevel.pl <list of test programs>\n";
}

my %sizes=();

for my $x (@ARGV) {
    $x =~ /^tune_/ or die "Bad test program $x";
    $x =~ /mul(\d)/ or die "Bad test program $x";
    my $s=$1;
    if (!exists $sizes{$s}) { $sizes{$s}=[]; }
    push @{$sizes{$s}}, $x;
}

sub mysys {
    print STDERR "@_\n";
    system @_;
}



my $make;
if (!defined($make=$ENV{'MAKE'})) {
    $make="$make";
}

# make sure the upper build is complete. It might be bad to do this
# check, but since this script is bound to call make anyway...
mysys "cd .. ; $make";

my @summary = ();

for my $s (sort { $a <=> $b } keys %sizes) {
    print STDERR "Tuning for $s word(s)\n";
    # Now we check everything, always.
#    if (scalar @{$sizes{$s}} == 1 && !defined($ENV{'BENCH'})) {
#        (my $x = $sizes{$s}->[0]) =~ s/^tune_//;
#        print STDERR "Only one possibility ($x) -- check skipped\n";
#        push @summary, "mul$s -> $x.c (only choice)\n";
#        next;
#    }
    my @results;
    print STDERR "Building the library afresh with current code selection.\n";
    mysys "cd .. ; $make clean";
    mysys "cd .. ; $make";
    for my $p (@{$sizes{$s}}) {
        mysys "$make $p";
        my $r = `./$p`;
        chomp($r);
        print STDERR "$r\n";
        $r =~ /^(.*)\s:\s([\d\.]+)\sns$/;
        die "$p returned no results" unless defined($1) && defined($2);
        push @results, [$1, $2];
    }
    @results = sort { $a->[1] <=> $b->[1] } @results;
    my $best = $results[0];
    $best->[0] =~ /tune_(.*)$/ or die;
    my $selected="$1.c";
    print STDERR "Selected $selected\n";
    my $keep_already_selected;
    my $ltarget="already_tuned/tuned/gf2x_mul$s.h";
    my $slot="gf2x/gf2x_mul$s.h";
    my $link_target_in_already_tuned_subdir;
    for my $r (@results) {
        $r->[0] =~ /tune_(.*)$/ or die;
        my $cfile = $1 . ".c";
        my $msg = "mul$s -> $cfile [ $r->[1] ns ]";
        if ($r == $best) {
            $msg .= " **BEST**";
        }
        my $xxlink_target_in_already_tuned_subdir = "../../lowlevel/$cfile";
        # Arrange so that $cfile is something reachable for us, and also
        # so that $xxlink_target_in_already_tuned_subdir is reachable from the
        # already_tuned/tuned subdirectory.
        if (! -f $cfile) {
            my $e;
            if (defined($e=$ENV{'abs_srcdir'}) && -f "$e/../lowlevel/$cfile") {
                # We are building out of source, so we resort to putting an
                # absolute path in the link target.
                $xxlink_target_in_already_tuned_subdir="$e/../lowlevel/$cfile";
                $cfile = "$e/../lowlevel/$cfile";
            } else {
                die "Cannot find $cfile anywhere !"
            }
        }
        # my $prepared="ready_gf2x_mul$s.c";
        # mysys "sed -e s///g $selected > $prepared";
        # print "diff ../$slot $cfile\n";
        my $rc=system "diff ../$slot $cfile > /dev/null";
        if ($rc == 0) {
            $msg .= " (previous)";
            if ($r == $best) {
                $msg .= " -> no change";
                $keep_already_selected = 1;
            }
        }
        if ($r == $best) {
            $selected = $cfile;
            $link_target_in_already_tuned_subdir = $xxlink_target_in_already_tuned_subdir;
        }
        push @summary, $msg . "\n";
    }


    if ($keep_already_selected) {
        next;
    }
    mkdir "../already_tuned" unless -d "../already_tuned";
    mkdir "../already_tuned/tuned" unless -d "../already_tuned/tuned";
    # Show the commands at the same time as we execute them.
    mysys "rm -f ../$slot";
    mysys "rm -f ../$ltarget";
    if ($selected =~ /gen_/) {
        # generated file: do a copy, not a link.
        mysys "cp -f $selected ../$ltarget";
    } else {
        mysys "ln -sf $link_target_in_already_tuned_subdir ../$ltarget";
    }
    mysys "ln -sf ../$ltarget ../$slot";
}

print STDERR "Summary of tune-lowlevel results\n";
for my $m (@summary) {
    print STDERR $m;
}
