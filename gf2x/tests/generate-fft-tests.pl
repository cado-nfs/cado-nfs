#!/usr/bin/env perl
#  This file is part of the gf2x library.
#
#  Copyright 2007, 2008, 2009, 2010, 2011, 2012, 2015
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

my @lines;
my @pattern;
my @arguments;
my @sed;
my $target;

my $in_pat=0;
my $in_arg=0;
my $in_sed=0;
my $in_gen=0;

while (<>) {
    if (!$in_pat && /^# -- begin pattern --$/) {
        $in_pat=1;
        print;
        next;
    }
    if (!$in_arg && /^# -- begin arguments --$/) {
        $in_arg=1;
        print;
        next;
    }
    if (!$in_arg && /^# -- begin sed --$/) {
        $in_sed=1;
        print;
        next;
    }
    if (!$in_gen && /^# -- begin generated code --$/) {
        $in_gen=1;
        print;
        next;
    }
    print unless $in_gen;
    if ($in_pat && /^# -- end pattern --$/) {
        $in_pat = 0;
        next;
    }
    if ($in_arg && /^# -- end arguments --$/) {
        $in_arg = 0;
        next;
    }
    if ($in_arg && /^# -- end sed --$/) {
        $in_sed = 0;
        next;
    }
    if ($in_gen && /^# -- end generated code --$/) {
        # generate everything now.
        my $footer=$_;
        for my $arg (map { s/^#\s//; $_ } @arguments) {
            if ($arg =~ /^(if|endif)/) {
                print "$arg";
                next;
            }
            my @args=split(' ', $arg);
            my @locpat = @pattern;
            for (@locpat) {
                s/^#\s//;
                my $i = 0;
                for my $a (@args) {
                    s/ARG$i/$a/g;
                    $i++;
                }
                for my $s (@sed) {
                    $s =~ s/^#\s//;
                    eval $s;
                }
                print;
            }
        }
        print $footer;
        $in_gen=0;
    }
    next if $in_gen;
    push @pattern, $_ if $in_pat;
    push @arguments, $_ if $in_arg;
    push @sed, $_ if $in_sed;
}
print join('',@lines);
