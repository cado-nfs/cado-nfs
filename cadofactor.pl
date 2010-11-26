#!/usr/bin/perl -w
#
# General script for factoring integers with Cado-NFS.
#
# Copyright 2008, 2009, 2010 Pierrick Gaudry, Emmanuel Thome, Paul Zimmermann,
#                            Jeremie Detrey, Lionel Muller
#
# This file is part of CADO-NFS.
#
# CADO-NFS is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CADO-NFS; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# 
# NB:
# ===
#
# For a small factorization, it is possible to use the factor.sh script
# (easier to use) which runs the factorization only on the local machine.
# Using this strict is recommended for larger factorizations involving
# several machines.
#
# Usage:
# ======
#
# cadofactor.pl param=<paramfile> wdir=<...> n=<...> ...
# 
# where:
#   <paramfile> is a file describing the choice of parameters for the
#               algorithm. It depends essentially only on the size of 
#               the input number, but it also contain non-algorithmic
#               information, like the nice level you want, and the name
#               of the file where to find the list of machines to use.
#               The distribution comes with a list of default paramfiles
#               that you can use.
#   <wdir>      a directory where cadofactor.pl will put all the
#               intermediate information (relation files, log files, etc)
#   <n>         the number you want to factor. 
#
# At the end of the command line, you can add parameters, in the form
# x=my_choice_for_x. If they were defined in <paramfile>, they will
# be overriden. It is useful if you choose for <paramfile> one of the 
# default files in the params/ directory. You can then override some
# parameter without editing the file.
#
# The possible parameters in a <paramfile> are not listed here. The 
# file params/params.c91 is an example where a lot of comments have been
# added, describing the role of each parameter.
#
# An important parameter, though, is machines=/path/to/mach_desc. It
# tells cadofactor.pl the list of available computing ressources for this
# factorization. If you just want to run it on your local computer, the
# minimal mach_desc file looks like:
#   tmpdir=/tmp/
#   bindir=/tmp/cado-nfs/build/localhost
#   localhost cores=4
# (this file is created automatically by the script factor.sh).
# The format of mach_desc is described in the example file: 
#   params/mach_desc
# It allows in particular to tell exactly how many cores there are to be
# used on each computer.
#
# Communications between nodes:
# =============================
#
# The communications between cadofactor.pl and the machines are done with
# SSH. You must have an automatic way to authenticate yourself:  having
# an ssh-agent running is recommended. See README for a few more things
# about SSH configuration.
#
# The communications between machines during the linear algebra step is
# done using MPI (if cado-nfs was compiled with MPI suppport). 
#
# What happens when it runs?
# ==========================
#
# After parsing (and checking) the parameters and the command line, the
# script starts running jobs of available ressources. It periodically
# checks if the tasks are finished and if it is the case, it import the 
# resulting data and start new jobs.
#
# The filtering phase is currently fully sequential, and is run on the 
# host machine that run cadofactor.pl. The same is true for the square
# root task.
#
# The linear algebra can be run in parallel if MPI is activated at
# compilation and if the mach_desc file contains machines that are
# explicitly marked as supporting mpi. Otherwise it is also run on the
# host machine.
# 
# Therefore, for non-trivial factorization, the host machine that runs
# cadofactor.pl must have enough memory.
#
# The cadofactor.pl script write in $wdir all its intermediate data and
# some diagnostics. Notice in particular the two following files:
#   - $wdir/$name.cmd  
#     It contains all the shell commands that are run by the script
#     with full arguments, so that it is easy to reproduce part of
#     the computation.
#   - $wdir/*.log
#     Contains log files of various steps.
#   - $wdir/$bwc/
#     This directory is specific to the linear algebra step.
#
# cadofactor.pl is supposed to resist to a crash of the host or remote
# machines: if it is run again with the same parameters, it will (try to)
# fix the problems and continue the computation.
#
# The script uses special (empty) files of the form
#     $wdir/xxxxx_done
# to remember that the step xxxxx is finished. It is therefore possible
# to force cadofactor.pl to redo some step by removing the corresponding
# file. Conversely, creating such a file can be used to skip some steps,
# if the corresponding data is produced by a third-party software or
# by hand. It is used for instance for an SNFS computation (see below).
#  
#
# Example: factorization of an integer of 155 digits (on several machines):
# =========================================================================
# 
# Before starting the factorization, you must configure ssh (cf README)
# for all the machines.
#
# $ mkdir $HOME/c155
# $ cd $HOME/c155
# $ cp $CADO_DIR/params/mach_desc .
#   Edit the file mach_desc and configure it.
# $ cp $CADO_DIR/params/params.c155 .
#   It is not necessary to edit the paramfile but you can choose to modify
#   some parameters in this file or on the cadofactor.pl command line.
#   If the paramfile does not exist then you must create it.
# $ $CADO_DIR/cadofactor.pl params.c155 wdir=$HOME/c155 name=rsa155 \
#    machines=mach_desc \
#    n=10941738641570527421809707322040357612003732945449205990913842131476349984288934784717997257891267332497625752899781833797076537244027146743531593354333897
#
#
# To factorize an integer with SNFS (on several machines):
# ==============================================================================
# 
# Before starting the factorization, you must configure ssh (cf README)
# for all the machines.
#
# $ mkdir $HOME/snfs
# $ cd $HOME/snfs
# $ cp $CADO_DIR/params/mach_desc .
#   Edit the file mach_desc and configure it.
# $ cp $CADO_DIR/params/params.c<size> .  
#   [substitute <size> by the integer size]
#   It is not necessary to edit the paramfile but you can choose to modify
#   some parameters in this file or the cadofactor.pl command line.
#   If the paramfile does not exist then you must create it.
# $ cp <polynomial choosen in the cado format> snfs<size>.poly
#   The polynomial file must also contain the sieve parameters.
# $ $CADO_DIR/cadofactor.pl params.c<size> wdir=$HOME/snfs name=snfs<size> \
#    machines=mach_desc n=<n>
# $ touch snfs<size>.polysel_done
# $ $CADO_DIR/cadofactor.pl params.c<size> wdir=$HOME/snfs name=snfs<size> \
#    machines=mach_desc n=<n>

use Cwd qw(abs_path);
use File::Basename;
use lib abs_path(dirname($0));
use cadofct;
use strict;
use warnings;

select(STDERR); $| = 1; # always flush stderr
select(STDOUT); $| = 1; # always flush stdout

print "$0 @ARGV\n";

do_init();
do_task("sqrt");

banner("All done!");

open FILE, "$param{'prefix'}.allfactors"
or die "Cannot open `$param{'prefix'}.allfactors' for reading: $!.\n";
my @l = <FILE>;
close FILE;
chomp for @l;
@l = sort {$a <=> $b} @l;
print "@l\n";

if (defined(my $e = $param{'expected_factorization'})) {
	my @exp=split(',',$e);
	@exp = sort {$a <=> $b} @exp;
	my $ok = "@exp" eq "@l";
	if ($ok) {
		print "Factorization matches expected result\n";
	} else {
		die "Factorization does not match expected result\n";
	}
}

