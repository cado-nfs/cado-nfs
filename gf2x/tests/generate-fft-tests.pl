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
my $target;

while (<>) {
    push @lines, $_;
    if (/^# -- begin generated code --$/) {
        while (<>) {
            if (/^# -- end generated code --$/) {
                push @lines, '';
                $target=\$lines[$#lines];
                last;
            }
        }
        push @lines, $_;
        next;
    }
    if (/^# -- begin pattern --$/) {
        while (<>) {
            push @lines, $_;
            last if /^# -- end pattern --$/;
            s/^#\s//;
            push @pattern, $_;
        }
        next;
    }
    if (/^# -- begin arguments --$/) {
        while (<>) {
            push @lines, $_;
            last if /^# -- end arguments --$/;
            chomp($_);
            s/^#\s//;
            push @arguments, $_;
        }
        next;
    }
}

my $pat=join('',@pattern);
for my $arg (@arguments) {
    if ($arg =~ /^(if|endif)/) {
        $$target .= "$arg\n";
        next;
    }
    my @args=split(' ', $arg);
    my $i = 0;
    my $y = $pat;
    while (my $a = shift @args) {
        $y =~ s/ARG$i/$a/g;
        $i++;
    }
    $$target .= $y;
}

select(STDOUT);

print join('',@lines);
