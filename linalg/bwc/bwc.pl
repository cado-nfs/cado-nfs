#!/usr/bin/env perl
use warnings;
use strict;

use POSIX qw/getcwd/;
use File::Basename;
use File::Temp qw/tempdir tempfile mktemp/;
use File::Spec;
use List::Util qw/max/;
# do not include Data::Dumper except for debugging.
# use Data::Dumper;
use Fcntl;
use Carp;

# This companion program serves as a helper for running bwc programs.
#
# It should be rewritten someday, perhaps in python. At least this has to
# be more modular (to be honest, it used to be much worse).

# MPI_BINDIR=/opt/mpich2-1.1a2-1.fc10.x86_64/usr/bin ./bwc.pl :complete matrix=c72 mn=64 ys=0..64 mpi=4x4

# {{{ usage
sub usage {
    print STDERR <<EOF;
Usage: ./bwc.pl <action> <parameters>

# 20140901: some revamping in progress. The instructions below may be
# partly inaccurate.

The action to be performed is either:
- the path to a program, in which case the command line eventually called
  is quite similar to ``<action> <parameters>'', with some substitutions
  performed. Depending on the program name, some options are discarded
  because they are not meaningful.
- one of the special actions defined by the script. These actions are
  prepended by a colon.
    :complete -- do a complete linear system solving. The solution ends
                 up in \$wdir/W.
    :bench    -- does prep, secure, and krylov with
                 an exceedingly large finish bound in order to perform
                 some timings.

Parameters are specified in the form <key>=<value>.
All parameters having a meaning for the bwc programs are accepted. Some
parameters control the wrapping script itself.

mpi thr              same meaning as for bwc programs.
mn m n prime         same meaning as for bwc programs.
nullspace            same meaning as for bwc programs.
ys                   same meaning as for bwc programs (krylov)
solutions            same meaning as for bwc programs (mksol, gather)
simd                 SIMD width to be used for krylov (ys=) and mksol,
                     gather (solutions=) programs. Defaults to 64 if
                     prime=2 or 1 otherwise
lingen_mpi           like mpi, but for lingen only

matrix               input matrix file. Must be a complete path. Binary, 32-b LE
rhs                  rhs file for inhomogeneous systems (ascii with header)
                     Note that inhomogeneous systems over GF(2) are not
                     supported yet (but we're not that far), and the file
                     format for the RHS file is not completely decided
                     yet.
wdir                 working directory.
                     If unset, derived from 'matrix', 'mpi' 'thr'
mpiexec              Command to invoke for mpiexec
hosts                Comma-separated list of hosts. 'hosts' may be supplied
                     several times, so that lists concatenate. Thus one
                     may use the shell to supply hosts=node{1,3,5,7}.
hostfile             path to a host list file to be understood by mpiexec

Amongst parameters, some options with leading dashes may be set:

-d  Do a dry run ; show the commands that would be executed.
    (-d is relevant to the script only)
-v  Increase verbosity level ; understood by bwc, but meaningless at the moment
-h  Show this help.

Some environment variables have an impact on this script.

MPI_BINDIR   The path to the directory holding the mpiexec program (the
             mpiexec= parameter allows to specify the same thing. Note that
             MPI alone serves as an alias for MPI_BINDIR, but this may
             become legacy. Normally mpiexec is obtained by cmake
             substitution, at least this is the default.).
BWC_BINDIR   The directory holding the bwc binaries. Defaults to the
             directory holding the current script.

EOF

    if (scalar @_) {
        print STDERR "Error messages:\n", join("\n", @_), "\n";
    }

    exit 1;
}
# }}}

# $program denotes either a simple command, or something prepended by ':'
# which indicates something having a special meaning for the script.
my $main = shift @ARGV or usage;

my $my_cmdline="$0 $main @ARGV";

# ----- cmake substituted variables -----
## mpiexec is substituted by cmake in case mpi has been used for the
## compilation. NOTE that this means that a priori, mpiexec _must_ be
## used for running all programs.
my $mpiexec='@MPIEXEC@';
$mpiexec='' if $mpiexec =~ m{^\@.*\@$};
my $mpiexec_extra_stanzas='@MPIEXEC_EXTRA_STANZAS@';
$mpiexec_extra_stanzas='' if $mpiexec_extra_stanzas =~ m{^\@.*\@$};


my $bindir;
if (!defined($bindir=$ENV{'BWC_BINDIR'})) {
    $bindir=$0;
    $bindir =~ s{/[^/]*$}{};
    if ($bindir eq $0) {
        $bindir='.';
    }
}

# It's only used so that we can keep a non-zero reference count on
# temporary filenames which we would like to keep until program exit.
my @tempfiles;


##################################################
# Parse the command line, and fill the following:
#  - %$param hash, with special rules for the hosts argument.
#  - $show_only, @extra_args
#  - obey -h
my @extra_args=();
my $show_only=0;
my $param={};
my $param_defaults={
    # some defaults.
    prime => 2,         # for factoring.
    thr => '1x1',
    mpi => '1x1',
    lingen_mpi => '1x1',
};
# {{{
while (defined($_ = shift @ARGV)) {
    if ($_ eq '--') {
        push @extra_args, splice @ARGV;
        last;
    }
    # -d will never be found as a bw argument, it is only relevant here.
    if (/^(-d|--show|--dry-run)$/) { $show_only=1; next; }
    if (/^(-h|--help)$/) { usage; }
    my ($k,$v);
    if (/^--([^=]+)$/ && scalar @ARGV) { $k=$1; $v=shift(@ARGV); }
    elsif (/^(-.*)$/) { $k=$1; $v=undef; }
    elsif (/^([^=]+)=(.*)$/) { $k=$1; $v=$2; }
    if (!defined($k)) {
        usage "Garbage not understood on command line: $_";
    }
    if ($k eq 'bwc_bindir') { $bindir=$v; next; }
    if (!defined($param->{$k})) {
        $param->{$k}=$v;
        next;
    }
    if ($k eq 'hosts') {
        if (ref $param->{$k} eq '') {
            $param->{$k} = [$param->{$k}, $v];
        } elsif (ref $param->{$k} eq 'ARRAY') {
            $param->{$k} = [@{$param->{$k}}, $v];
        } else {
            die "\$param->{$k} has gone nuts.\n";
        }
    } else {
        usage "parameter $k may not be specified more than once";
    }
}
for my $k (keys %$param_defaults) {
    next if exists $param->{$k};
    $param->{$k} = $param_defaults->{$k};
}


# }}}

# {{{ global variables. Much resides in $param, but we do have globals
# too.
my @main_args;

my @mpi_split=(1,1);
my @lingen_mpi_split=(1,1);
my @thr_split=(1,1);
my $nh;
my $nv;

my $wdir;
my $matrix;
my $random_matrix;
my $rhs;
my $nrhs=0;

## The mpi_extra_args argument is used to pass information to the mpiexec 
## command. The idea is that mpiexec, the mpi driver program, may need
## additional info to properly setup communications between jobs. As an 
## example, if libvirt-bin is installed on a linux system, openmpi will
## need the following additions.
##    
## ./build/x86_64/linalg/bwc/bwc.pl :complete matrix=/tmp/mats/c90b wdir=/tmp/c59 mn=64 mpi=2x1 thr=2x4 hosts=node0,node1 mpi_extra_args='--mca btl_tcp_if_exclude lo,virbr0'
##    
## This tells mpi not to try routing traffic through either the lo or the 
## virbr0 interface. For the former, it's already openmpi's default
## behaviour (unless /etc/hosts is screwed up, which this file does tend to
## be sometimes). For virbr0, it's again some local mess, but each machine 
## having this class C network defined, openmpi can't really tell whether 
## they're connected or not -- and of course they're not.
##
## Note that we also obey the $MPI_EXTRA_ARGS environment variable.
my $mpi_extra_args;
 
my ($m, $n);
my $prime;
my $splitwidth;
my $simd;

# my $force_complete;
my $stop_at_step;

my $hostfile;
my @hosts=();
my @mpi_path=();
my $mpi_ver;
my $needs_mpd;
# }}}
#
# Some more command line interpretation. Here we build @main_args, taking
# out the arguments which are only of interest to us.
# Here we read (from $param) the
# following things.
#  - $hostfile $mpi_extra_args
#  - ## broken ## $force_complete
#  - @hosts $m $n
#  - @mpi_split @thr_split
#  - $nh $nv
#    --> dimensions of the split of the matrix (horizontal chunks,
#    vertical chunks).
#  - $prime
#    --> Linear system is over GF(p) for that prime. For factoring, which
#    is the default, we have prime==2
#  - $splitwidth
#    --> Number of vectors which are put together in a file. This is only
#    inferred from $prime, and set to 64 for prime==2, 1 otherwise.
#
# {{{
sub set_mpithr_param { # {{{ utility
    my $v = shift @_;
    my @s = (1,$v);
    if ($v=~/^(\d+)$/) {
        my $nthreads = $1;
        my $s = int(sqrt($nthreads));
        for (; $s >= 1; $s--) {
            next if $nthreads % $s > 0;
            @s = ($s, $nthreads / $s);
            last;
        }
        die unless $s[0]*$s[1] == $v;
    } elsif ($v=~/(\d+)x(\d+)$/) {
        @s = ($1, $2);
    } else {
        usage "bad splitting value '$v'";
    }
    return @s;
} # }}}

my $my_verbose_flags = {
    cmdline => 1,
    checks => 1,
    sections => 1,
};

while (my ($k,$v) = each %$param) {
    # Some parameters are relevant to us just like they're relevant to
    # the bwc programs, so we'll set a variable based on their value, and
    # also copy them to the @main_args array.
    if ($k eq 'matrix') { $matrix=$v; }
    if ($k eq 'random_matrix') { $random_matrix=$v; }
    # Some parameters which are simply _not_ relevant to the bwc programs.
    if ($k eq 'hostfile') { $hostfile=$v; next; }
    if ($k eq 'simd') { $simd=$v; next; }
    if ($k eq 'mpi_extra_args') { $mpi_extra_args=$v; next; }
    ## if ($k eq 'force_complete') { $force_complete=$v; next; }
    if ($k eq 'stop_at_step') {
        die "stop_at_step requires :complete" unless $main =~ /^(?::srun)?:complete/;
        $stop_at_step=$v;
        next;
    }
    if ($k eq 'hosts') {
        $v=[$v] if (ref $v eq '');
        for (@$v) { push @hosts, split(',',$_); }
        next;
    }
    # Some of the command-line arguments are modified before being put
    # into the main argument list.
    if ($k eq 'mpi') { @mpi_split = set_mpithr_param $v; $param->{$k} = $v = join "x", @mpi_split;}
    if ($k eq 'thr') { @thr_split = set_mpithr_param $v; $param->{$k} = $v = join "x", @thr_split; }
    if ($k eq 'lingen_mpi') { @lingen_mpi_split = set_mpithr_param $v; $param->{$k} = $v = join "x", @lingen_mpi_split;}
    if ($k eq 'verbose_flags') {
        my @heritage;
        for my $f (split(',', $v)) {
            my $w=1;
            my $f0=$f;
            $w = 0 if $f =~ s/^[\^!]//;
            if ($f =~ s/^perl-//) {
                $my_verbose_flags->{$f}=$w;
            } else {
                push @heritage, $f0;
            }
        }
        $v=join(',',@heritage);
    }

    # The rest is passed to subprograms, unless explicitly discarded on a
    # per-program basis.
    if (defined($v)) {
        push @main_args, "$k=$v";
    } else {
        # This is for switches (like -v) which don't have a "value"
        # associated with them
        push @main_args, "$k";
    }

    # Yet it does make sense to read some parameters, as we are
    # interested by their value.

    if ($k eq 'wdir') { $wdir=$v; next; }
    if ($k eq 'mn') { $m=$n=$v; next; }
    if ($k eq 'n') { $n=$v; next; }
    if ($k eq 'm') { $m=$v; next; }
    if ($k eq 'prime') { $prime=$v; next; }
    if ($k eq 'rhs') { $rhs=$v; next; }

    # Ok, probably there's a fancy argument we don't care about.
}

$nh = $mpi_split[0] * $thr_split[0];
$nv = $mpi_split[1] * $thr_split[1];
$splitwidth = ($prime == 2) ? 64 : 1;

if (!defined($simd)) {
    $simd = $splitwidth;
}

# }}}

print "$my_cmdline\n" if $my_verbose_flags->{'cmdline'};

if ($main =~ /^:(?:mpi|s)run(?:_single)?$/) {
    # ok, this is really an ugly ugly hack. We have some mpi detection
    # magic in this script, which we would like to use. So the :mpirun
    # meta-command is just for that. Of course the argument requirements
    # are mostly waived in this case.
    $matrix=$param->{'matrix'}=$0;
    $param->{'prime'}=2;
    $param->{'simd'}=$simd=64;
    $m=$n=64;
    $wdir=$param->{'wdir'}="/";
    if ($main =~ /^:srun/) {
        # Setting this will trigger mpi_needed, mpi detection check,
        # and ultimately passing -n
        $mpiexec='srun';
    }
}

# {{{ Some important argument checks
{
    my @miss = grep { !defined $param->{$_}; } (qw/prime/);
    die "Missing argument(s): @miss" if @miss;
    if (!defined($param->{'matrix'}) && !defined($param->{'random_matrix'})) {
        die "Missing parameter: matrix (or random_matrix)";
    }
    if (!defined($m) || !defined($n)) {
        die "Missing parameters: m and/or n";
    }
}

if (defined($param->{'matrix'})) {
    $param->{'matrix'} =~ s/~/$ENV{HOME}/g;
    die "$param->{'matrix'}: $!" unless -f $param->{'matrix'};
}

if ($prime == 2 && (($m % 64 != 0) || ($n % 64 != 0))) {
    die "Currently for p=2 bwc supports only block sizes which are multiples of 64";
}
if (($m % $splitwidth != 0) || ($n % $splitwidth != 0)) {
    die "Currently bwc supports only block sizes which are multiples of $splitwidth";
}
if ($prime == 2 && defined($rhs)) {
    die "inhomogeneous systems currently supported only for p>2";
}
if (defined($rhs)) {
    # Try to read the first line.
    my $header = eval { open F, $rhs or die "$rhs: $!"; <F>; };
    if ($header !~ /^(\d+) (\d+) (\d+)$/) {
        die "$rhs does not seem to be ascii with header, we can't proceed.\n";
    }
    $param->{'nrhs'} = $nrhs = $2;
    print "$rhs is ascii file with header. Getting nrhs=$nrhs from there\n";
}
if ($nrhs > $n || $nrhs > $m) {
    die "nrhs > n or m is not supported";
}

if ($prime == 2 && ($lingen_mpi_split[0] != 1 || $lingen_mpi_split[1] != 1)) {
    die "Current binary lingen code does not support mpi";
}

if ($lingen_mpi_split[0] != $lingen_mpi_split[1]) {
    die "lingen_mpi ($lingen_mpi_split[0]x$lingen_mpi_split[1]) must be a square split";
}

if ($m % $lingen_mpi_split[0] != 0 || $n % $lingen_mpi_split[0] != 0) {
    die "lingen_mpi must divide gcd(m,n)"
}

# }}}

# {{{ wdir preliminary (pre-mpi) checks.
#  - Check that wdir is defined (or define it to some default);
#  - Check that it exists if it has to exist (:complete is entitled to
#  create $wdir

# Selection of the working directory -- there are three relevant options.
# The 'matrix' option gives the file name of the matrix to be used.
# The working directory used is either specified as a parameter on its
# own named 'wdir', or inferred from the basename of the matrix file,
# appended with some information on the splitting, e.g.  c160-4x4.

# Note that when doing random matrices, we don't even need a wdir.

if (!defined($random_matrix)) {
    if (!defined($param->{'wdir'}) && defined($param->{'matrix'})) {
        # ok, we're doing this a bit ahead of time
        $wdir=$param->{'wdir'}="$param->{'matrix'}-${nh}x${nv}";
        push @main_args, "wdir=$wdir";
    }

    if ($main !~ /^(?::srun)?:(?:complete|bench)$/ && !-d $param->{'wdir'}) {
        die "$param->{'wdir'}: no such directory";
    }
}
# }}}

# Some quick shortcuts in case the user failed to specify some arguments
#  - ys=
#    --> this is valid only for krylov, for *one* sequence run.
#    In general, this must correspond to an interval defining a sequence
#    (i.o.w. two consecutive multiples of simd width).
#    For interleaving, this is a bit trickier, since two
#    contiguous intervals are to be treated.
# The code below is just setting sensible defaults.
# {{{
##################################################
# Some more default values to set. Setting separately ys= and solutions= 
# is quite awkward in the case of a single-site test.
if ((!defined($m) || !defined($n))) {
    usage "The parameters m and n must be set";
}

if (!defined($param->{'solutions'})) {
    $param->{'solutions'}=["0-$simd"];
    if ($param->{'interleaving'}) {
        $param->{'solutions'} = [ "0-" . (2*$simd) ];
    }
} else {
    # make that a list.
    $param->{'solutions'} = [ split(',',$param->{'solutions'}) ];
    if (scalar(@{$param->{'solutions'}}) == 1 and $splitwidth == 1) {
        my ($s0, $s1) = split('-',$param->{'solutions'}->[0]);
        $param->{'solutions'} = [];
        for(my $i = 0 ; $s0 + $i < $s1 ; $i += $splitwidth) {
            push @{$param->{'solutions'}}, join('-', ($s0 + $i, $s0 + $i + $splitwidth ));
        }
        print STDERR "## Note: rewriting solutions=$s0-$s1 into solutions=" . join(",", @{$param->{'solutions'}} ) . "\n";
    }
}

for (@{$param->{'solutions'}}) {
    /^\d+-\d+$/ or die "the 'solutions' parameter must match \\d+-\\d+";
}

# Default settings for ys= --> see krylov / mksol / gather.
# }}}

# Done playing around with what we've been given.

##################################################
# System interaction.
# - set environment variables with my_setenv (which does a bit more)
# - run commands with dosystem
# - ssh to nodes with ssh_program ; this is MPI-backend-dependent.
# {{{
my $env_strings = "";
sub my_setenv
{
    return if (exists($ENV{$_[0]}) && $ENV{$_[0]} eq $_[1]);
    $env_strings .= "$_[0]=$_[1] ";
    $ENV{$_[0]}=$_[1];
}

sub dosystem
{
    my $nofatal=0;
    my $prg = shift @_;
    if ($prg eq '-nofatal') {
        $nofatal=1;
        $prg=shift @_;
    }
    my @args = @_;
    print STDERR '#' x 77, "\n" if $my_verbose_flags->{'sections'};
    my $msg = "$env_strings$prg " . (join(' ', @args)) . "\n";

    if ($show_only) {
        print $msg;
        return 0;
    }

    print STDERR $msg if $my_verbose_flags->{'cmdline'};
    my $rc = system $prg, @args;
    return if $rc == 0;
    if ($rc == -1) {
        print STDERR "Cannot execute $prg\n";
    } elsif ($rc & 127) {
        my $sig = $rc & 127;
        my $coreinfo = ($rc & 128) ? 'with' : 'without';
        print STDERR "$prg: died with signal $sig, $coreinfo coredump\n";
    } else {
        my $ret = $rc >> 8;
        print STDERR "$prg: exited with status $ret\n";
    }
    die "aborted on subprogram error";
}

sub ssh_program
{
    my $ssh;
    if ($ssh = $ENV{'SSH'}) {
        return $ssh;
    }
    $ssh='ssh';
    if (exists($ENV{'OAR_JOBID'})) {
        $ssh="/usr/bin/oarsh";
    }
    return $ssh
}
# }}}

if (!defined($random_matrix) && !-d $wdir) {        # create $wdir on script node.
    if ($show_only) {
        print "mkdir $wdir\n";
    } else {
        mkdir $wdir;
    }
}

##################################################
# MPI-specific detection code. Beyond probing the environment variables
# to see which MPI middleware we're running, this code sets the following
# important argument lists:
#  - @mpi_precmd @mpi_precmd_single
#    --> Simply put, the former is prepended before all important
#    mpi-level programs (secure krylov mksol gather), while the
#    other is of course for leader-node-only programs.
#  - @mpi_precmd_lingen
#    --> use lingen_mpi (no lingen_thr exists at this point)
#  - $mpi_needed
#    --> to be used to check whether we want mpi. This is important
#    before calling dosystem.
# {{{

# If we've been built with mpi, then we _need_ mpi for running. Otherwise
# we run into shared libraries mess.

my $mpi_needed = $mpiexec ne '';

# @mpi_precmd_single is something silly; we want provision for the case
# where mpi is used for running non-mpi jobs. It's something which does
# not officially work, yet it always does. And some programs do turn out
# to be compiled with mpi, so we need the mpi libraries at runtime... So
# short of a more accurate solution, this is a hack.
#
# Note also that @mpi_precmd_single is absolutely needed when we run a
# non-mpi binary, yet that does load all of the mpi infrastructure, in a
# context where that mpi infrastructure smartly reads the environment
# variables and detects that we need X jobs (use case: bwccheck does
# bw_common_init_mpi out of convenience, but really isn't an mpi program).
my @mpi_precmd;
my @mpi_precmd_single;
my @mpi_precmd_lingen;


# This now has $mpi a local variable, and the array @mpi_path
#
sub detect_mpi {
    my $mpi;
    if (defined($_=$ENV{'MPI_BINDIR'}) && -x "$_/mpiexec") {
        $mpi=$_;
    } elsif (defined($_=$ENV{'MPI'}) && -x "$_/mpiexec") { 
        $mpi=$_;
    } elsif (defined($_=$ENV{'MPI'}) && -x "$_/bin/mpiexec") { 
        $mpi="$_/bin";
    } elsif ($mpiexec && -x $mpiexec) {
        my($basename, $dirname) = fileparse($mpiexec);
        $mpi=$dirname;
    } else {
        my @path = split(':',$ENV{'PATH'});
        for my $d (@path) {
            if (-x "$d/mpiexec") {
                $mpi=$d;
                print STDERR "Auto-detected MPI_BINDIR=$mpi\n";
                last;
            }
        }
    }

    my $maybe_mvapich2=1;
    # mpich versions implementing the mpi-2 standard were named mpich2.
    # From standard mpi-3 on, the name has returned to mpich.
    my $maybe_mpich=1;
    my $maybe_openmpi=1;
    my $maybe_intelmpi=1;

    if (defined($mpi)) {
        SEVERAL_CHECKS: {
            # first check the alternatives system, which is fairly
            # commonplace.
            my $mpiexec = "$mpi/mpiexec";
            while (-l $mpiexec) {
                my $target=readlink($mpiexec);
                print STDERR "readlink($mpiexec)->$target\n";
                if ($target =~ m{^/}) {
                    $mpiexec = $target;
                } else {
                    $mpiexec =~ s{[^/]+$}{$target};
                }
                if ($mpiexec =~ /openmpi/) {
                    print STDERR "Auto-detecting openmpi based on alternatives\n";
                    $maybe_mvapich2=0;
                    $maybe_mpich=0;
                    $maybe_intelmpi=0;
                    last;
                } elsif ($mpiexec =~ /mpich2/) {
                    print STDERR "Auto-detecting mpich2(old) based on alternatives\n";
                    $maybe_mpich='mpich2';
                    $maybe_mvapich2=0;
                    $maybe_openmpi=0;
                    $maybe_intelmpi=0;
                    last;
                } elsif ($mpiexec =~ /mpich/) {
                    print STDERR "Auto-detecting mpich based on alternatives\n";
                    $maybe_mvapich2=0;
                    $maybe_openmpi=0;
                    $maybe_intelmpi=0;
                    last;
                } elsif ($mpiexec =~ /hydra/) {
                    # Newer mvapich2 uses hydra as well...
                    print STDERR "Auto-detecting mpich or mvapich2 (hydra) or intel mpi (hydra) based on alternatives\n";
                    $maybe_mvapich2='hydra';
                    $maybe_mpich='hydra';
                    $maybe_intelmpi='hydra';
                    $maybe_openmpi=0;
                    last;
                } elsif ($mpiexec =~ /mvapich2/) {
                    print STDERR "Auto-detecting mvapich2 based on alternatives\n";
                    $maybe_mpich=0;
                    $maybe_openmpi=0;
                    $maybe_intelmpi=0;
                    last;
                }
            }
            CHECK_MVAPICH2: {
                if ($maybe_mvapich2 && -x "$mpi/mpiname") {
                    my $v = `$mpi/mpiname -n -v`;
                    chomp($v);
                    if ($v =~ /MVAPICH2\s+([\d\.]+)((?:\D\w*)?)/) {
                        # Presently all versions of mvapich2 up
                        # until 1.6rc3 included need mpd daemons.
                        # Released version 1.6 uses hydra.
                        $mpi_ver="mvapich2-$1$2";
                        if (($1 < 1.6) || ($1 == 1.6 && $2 =~ /^rc\d/)) {
                            $needs_mpd = 1;
                        } else {
                            $mpi_ver .= "+hydra" unless $mpi_ver =~ /hydra/;
                        }
                        last SEVERAL_CHECKS;
                    }
                }
            }
            CHECK_MPICH_VERSION: {
                if ($maybe_mpich =~ /^(hydra|mpich2)/ && -x "$mpi/mpich2version") {
                    my $v = `$mpi/mpich2version -v`;
                    chomp($v);
                    if ($v =~ /MPICH2 Version:\s*(\d.*)$/) {
                        $mpi_ver="mpich2-$1";
                        # Versions above 1.3 use hydra.
                        $needs_mpd=($mpi_ver =~ /^mpich2-(0|1\.[012])/);
                    } else {
                        $mpi_ver="mpich2-UNKNOWN";
                        $needs_mpd=1;
                    }
                    $v = `$mpi/mpich2version -c`;
                    chomp($v);
                    if ($v =~ /--with-pm=hydra/) {
                        $mpi_ver .= "+hydra";
                        $needs_mpd=0;
                    }
                    if ($maybe_mpich eq 'hydra') {
                        $mpi_ver .= "+hydra" unless $mpi_ver =~ /hydra/;
                        $needs_mpd=0;
                    }
                    last SEVERAL_CHECKS;
                } elsif ($maybe_mpich && -x "$mpi/mpichversion") {
                    my $v = `$mpi/mpichversion -v`;
                    chomp($v);
                    if ($v =~ /MPICH Version:\s*(\d.*)$/) {
                        $mpi_ver="mpich-$1";
                        # Only antique mpich-1.x versions don't use
                        # hydra.
                        $needs_mpd=($mpi_ver =~ /^mpich-[01]\./);
                    } else {
                        $mpi_ver="mpich-UNKNOWN";
                        $needs_mpd=1;
                    }
                    last SEVERAL_CHECKS;
                } elsif ($maybe_mpich && -x "$mpi/mpivars.sh") {
                    my $v = `unset I_MPI_ROOT ; . $mpi/mpivars.sh  ; echo \$I_MPI_ROOT`;
                    if ($v) {
                        $mpi_ver="Intel MPI";
#       comment out, as this requires mpicc on the compute nodes (at least for debugging), and that isn't always available...
#			my $v = `$mpi/mpicc -v`;
#			chomp($v);
#			if ($v =~ /Intel \(R\) MPI Library (.*?) for/) {
#			    $mpi_ver="$1";
#			    $mpi_ver=~/^(\d+)/;
#                            my $mpi_year=$1;
#                            $mpi_ver=~s/\s+/_/g;
#			    $mpi_ver="intel-$mpi_ver";
#                            $needs_mpd=0;
#			} else {
#                            warn "Could not correctly recognize mpi version (perhaps intel ?)";
#                        }
                        $needs_mpd=0;
                    }
                    last SEVERAL_CHECKS;
                }
            }
            CHECK_OMPI_VERSION: {
                if ($maybe_openmpi && -x "$mpi/ompi_info") {
                    my @v = `$mpi/ompi_info`;
                    my @vv = grep { /Open MPI:/; } @v;
                    last CHECK_OMPI_VERSION unless scalar @vv == 1;
                    $needs_mpd=0;
                    if ($vv[0] =~ /Open MPI:\s*(\d\S*)$/) {
                        $mpi_ver="openmpi-$1";
                        last SEVERAL_CHECKS;
                    } else {
                        $mpi_ver="openmpi-UNKNOWN";
                        last SEVERAL_CHECKS;
                    }
                }
            }
            CHECK_INTELMPI_VERSION: {
                # there's an impi_info binary, but it seems to be mostly
                # useless
                if ($maybe_intelmpi && -x "$mpi/mpicc") {
                    my @v = `$mpi/mpicc -v`;
                    my @vv = grep { /Intel.*MPI.*Library/i } @v;
                    last CHECK_INTELMPI_VERSION unless scalar @vv == 1;
                    if ($vv[0] =~ /Intel.*MPI.*Library\s*(\S+)/i) {
                        $mpi_ver="intelmpi-$1";
                        last SEVERAL_CHECKS;
                    } else {
                        $mpi_ver="intelmpi-UNKNOWN";
                        last SEVERAL_CHECKS;
                    }
                }
            }
        }

        if (defined($mpi_ver)) {
            print STDERR "Using $mpi_ver, MPI_BINDIR=$mpi\n";
        } else {
            print STDERR "Using UNKNOWN mpi, MPI_BINDIR=$mpi\n";
            if (defined($needs_mpd=$ENV{NEEDS_MPD})) {
                warn "Assuming needs_mpd=$needs_mpd as per env variable.\n";
            } else {
                $needs_mpd=1;
                warn "Assuming needs_mpd=$needs_mpd ; " .
                    "modify env variable \$NEEDS_MPD to " .
                    "change fallback behaviour\n";
            }
        }
    }

    if (!defined($mpi)) {
        print STDERR <<EOMSG;
***ERROR*** No mpi library was detected. Arrange for mpiexec to be in
***ERROR*** your PATH, or set the MPI_BINDIR environment variable.
EOMSG
        exit 1;
    }

    # Some final tweaks.
    if ($mpi_ver =~ /^mvapich2/) {
        print "## setting MV2_ENABLE_AFFINITY=0 (for mvapich2)\n";
        $ENV{'MV2_ENABLE_AFFINITY'}=0;
    }

    @mpi_path=($mpi);

    if ($mpi_ver =~ /openmpi/ && $mpi =~ m{^/*bin/*}) {
        print STDERR "Rewriting MPI prefix to just nothing\n";
        @mpi_path=();
    }
}

# Starting daemons for mpich2 1.[012].x and mvapich2 ; we're assuming
# this works the same for older mpich2's, although this has never been
# checked.
sub check_mpd_daemons
{
    return unless $needs_mpd;

    my $ssh = ssh_program();

    my $mpdtrace = File::Spec->catfile(@mpi_path, "mpdtrace");
    my $mpdboot = File::Spec->catfile(@mpi_path, "mpdboot");
    my $rc = system "$mpdtrace > /dev/null 2>&1";
    if ($rc == 0) {
        print "mpi daemons seem to be ok\n";
        return;
    }

    if ($rc == -1) {
        die "Cannot execute $mpdtrace";
    } elsif ($rc & 127) {
        die "$mpdtrace died with signal ", ($rc&127);
    } else {
        print "No mpi daemons found by $mpdtrace, restarting\n";
    }

    open F, $hostfile;
    my %hosts;
    while (<F>) {
        $hosts{$_}=1;
    }
    close F;
    my $n = scalar keys %hosts;
    print "Running $n mpi daemons\n";
    dosystem "$mpdboot -n $n -r $ssh -f $hostfile -v";
}

sub get_mpi_hosts_oar {
    my @x = split /^/, eval {
		local $/=undef;
		open F, "$ENV{OAR_NODEFILE}";
		<F> };
    @hosts=();
    my %h=();
    for (@x) {
        push @hosts, $_ unless $h{$_}++;
    }
}

sub get_mpi_hosts_slurm {
    @hosts=();
    my %h=();
    for (`scontrol show hostname`) {
        chomp($_);
        push @hosts, $_ unless $h{$_}++;
    }
    if (!scalar @hosts) {
        die "no hosts found with \"scontrol show hostname\"; is scontrol in \$PATH ?";
    }
}

sub get_mpi_hosts_torque {
    my @x = split /^/, eval {
		local $/=undef;
		open F, "$ENV{PBS_NODEFILE}";
		<F> };
    my $nth = $thr_split[0] * $thr_split[1];
    @hosts=();
    while (scalar @x) {
        my @z = splice @x, 0, $nth;
        my %h=();
        $h{$_}=1 for @z;
        die "\$PBS_NODEFILE not consistent mod $nth\n" unless scalar
			keys %h == 1;
        my $c = $z[0];
        chomp($c);
        push @hosts, $c;
    }
}

sub get_mpi_hosts_sge {
    print STDERR "Building hosts file from $ENV{PE_HOSTFILE}\n";
    my @x = split /^/, eval { local $/=undef; open F, "$ENV{PE_HOSTFILE}"; <F> };
    my $cores_on_node = {};
    for my $line (@x) {
        my ($node, $ncores, $toto, $tata) = split ' ', $line;
        print STDERR "$node: +$ncores cores\n";
        $cores_on_node->{$node}+=$ncores;
    }
    my $values_for_cores_on_node = {};
    local $_;
    $values_for_cores_on_node->{$_}=1 for (values %$cores_on_node);

    die "Not always the same number of cores obtained on the different nodes, as per \$PE_HOSTFILE" if keys %$values_for_cores_on_node != 1;

    my $ncores_obtained = (keys %$values_for_cores_on_node)[0];
    my $nnodes = scalar keys %$cores_on_node;
    print STDERR "Obtained $ncores_obtained cores on $nnodes nodes\n";

    my $nthr = $thr_split[0] * $thr_split[1];
    my $nmpi = $mpi_split[0] * $mpi_split[1];

    die "Not enough cores ($ncores_obtained) obtained: want $nthr\n" if $nthr > $ncores_obtained;

    @hosts=();
    push @hosts, $_ for keys %$cores_on_node;

    die "Not enough mpi nodes ($nnodes): want $nmpi\n" if $nmpi > $nnodes;
}

# {{{ utilities
sub version_ge {
    my ($a, $b) = @_;
    my @a = split(/[\.-]/, $a);
    my @b = split(/[\.-]/, $b);
    while (@a && @b) {
        $a = shift @a;
        $b = shift @b;
        $a =~ /^(\d*)(.*)$/; my ($an, $at) = ($1, $2);
        $b =~ /^(\d*)(.*)$/; my ($bn, $bt) = ($1, $2);
        $an=0 if length($an)==0;
        $bn=0 if length($bn)==0;
        return 1 if $an > $bn;
        return 0 if $an < $bn;
        return 1 if $at gt $bt;
        return 0 if $at lt $bt;
    }
    return 0 if @b;
    return 1;
}

# }}}

open STDIN,  '<',  '/dev/null' or die "can't redirect STDIN: $!\n";

if ($mpi_needed) {
    # This is useful for debugging in case we see new MPI environments.
    # print STDERR "Inherited environment:\n";
    # print STDERR "$_=$ENV{$_}\n" for keys %ENV;
    detect_mpi;

    # quirk. If called with srun, then we want to avoid mpiexec anyway.
    if ($main =~ /^:srun/) {
        push @mpi_precmd, 'srun', '--cpu-bind=verbose,none', '-u';
	if ($mpi_ver =~ /^openmpi/) {
            # srun might like that we give it the hint that we're running
            # openmpi.
            push @mpi_precmd, "--mpi=openmpi";
	}
        # quirk
        $main =~ s/^:srun://;
        $main =~ s/^complete/:complete/;
    } else {
        # Otherwise we'll start via mpiexec, and we need to be informed
        # on the list of nodes.
        push @mpi_precmd, File::Spec->catfile(@mpi_path, "mpiexec");

        my $auto_hostfile_pattern="/tmp/cado-nfs.hosts_XXXXXXXX";

        # Need hosts. Put that to the list @hosts first.
        if ($main =~ /^:srun/) {
            print STDERR "srun environment detected, not detecting hostfile.\n";
        } elsif (exists($ENV{'OAR_JOBID'}) && !defined($hostfile) && !scalar @hosts) {
            print STDERR "OAR environment detected, setting hostfile.\n";
            get_mpi_hosts_oar;
            $auto_hostfile_pattern="/tmp/cado-nfs.hosts.$ENV{'OAR_JOBID'}.XXXXXXX";
        } elsif (exists($ENV{'PBS_JOBID'}) && !defined($hostfile) && !scalar @hosts ) {
            print STDERR "Torque/OpenPBS environment detected, setting hostfile.\n";
            get_mpi_hosts_torque;
            $auto_hostfile_pattern="/tmp/cado-nfs.hosts.$ENV{'PBS_JOBID'}.XXXXXXX";
        } elsif (exists($ENV{'PE_HOSTFILE'}) && exists($ENV{'NSLOTS'}) && !defined($hostfile) && !scalar @hosts) {
            print STDERR "Oracle/SGE environment detected, setting hostfile.\n";
            get_mpi_hosts_sge;
            # } elsif (exists($ENV{'SLURM_STEP_NODELIST'})) {
            # not clear that I want this detection to be done _only_ within a
            # job step. I do see cases where mpiexec from the batch itself seems
            # to make sense
        } elsif (exists($ENV{'SLURM_JOBID'}) && !defined($hostfile) && !scalar @hosts) {
            get_mpi_hosts_slurm;
        }

        if (scalar @hosts) {
            # I think that when doing so, the file will get deleted at
            # program exit only.
            my $fh = File::Temp->new(TEMPLATE=>$auto_hostfile_pattern);
            $hostfile = $fh->filename;
            push @tempfiles, $fh;
            for my $h (@hosts) { print $fh "$h\n"; }
            close $fh;
        }

        if (defined($hostfile)) {
            if ($needs_mpd) {
                # Assume daemons do the job.
            } elsif ($mpi_ver =~ /^\+hydra/) {
                my_setenv 'HYDRA_HOST_FILE', $hostfile;
                # I used to have various setups using --hostfile for openmpi,
                # --file in some other cases and so on. I think that
                # -machinefile is documented in the published standard, so
                # it's better to stick to it.
                #        } elsif ($mpi_ver =~ /^openmpi/) {
                #            push @mpi_precmd, "--hostfile", $hostfile;
                #        } else {
                #            push @mpi_precmd, "-file", $hostfile;
            } else {
                push @mpi_precmd, "-machinefile", $hostfile;
            }

        }
        if (!defined($hostfile)) {
            # At this point we're going to run processes on localhost.
            if ($mpi_ver =~ /^\+hydra/) {
                my_setenv 'HYDRA_USE_LOCALHOST', 1;
            }
            # Otherwise we'll assume that the simple setup will work fine.
        }
        check_mpd_daemons();
        if (!$needs_mpd) {
            # Then we must configure ssh
            # Note that openmpi changes the name of the option pretty
            # frequently.
            if ($mpi_ver =~ /^openmpi-1\.2/) {
                push @mpi_precmd, qw/--mca pls_rsh_agent/, ssh_program();
            } elsif ($mpi_ver =~ /^openmpi-1\.[34]/) {
                push @mpi_precmd, qw/--mca plm_rsh_agent/, ssh_program();
            } elsif ($mpi_ver =~ /^openmpi-1\.[56]/) {
                push @mpi_precmd, qw/--mca orte_rsh_agent/, ssh_program();
            } elsif ($mpi_ver =~ /^openmpi/) {
                push @mpi_precmd, qw/--mca plm_rsh_agent/, ssh_program();
                if (version_ge($mpi_ver, "openmpi-1.8")) {
                    # This is VERY important for bwc ! The default policy for
                    # openmpi 1.8 seems to be by slot, which obviously
                    # schedules most of the desired jobs on only one node...
                    # which doesn't work too well.
                    if (!$param->{'only_mpi'}) {
                        # with only_mpi=1, the default policy works fine.
                        push @mpi_precmd, qw/--map-by node/;

                        # I don't exactly when the --bind-to argument
                        # appeared.
                        push @mpi_precmd, qw/--bind-to none/;
                    }
                }
            } elsif ($mpi_ver =~ /^mpich2/ || $mpi_ver =~ /^mvapich2/) {
                # Not older mpich2's, which need a daemon.
                push @mpi_precmd, qw/-launcher ssh -launcher-exec/, ssh_program();
            }
        }
        # End of section where we try to set up the proper options so
        # that mpiexec reaches the desired nodes correctly.
    }
    if (defined($mpi_extra_args)) {
        my @a = split(' ', $mpi_extra_args);
        my @b;
        while (defined($_=shift(@a))) {
            if (/^--(map-by|bind-to)$/) {
                shift @a;
                next;
            }
            if (/^plm_rsh_agent$/) {
                pop @b;
                shift @a;
                next;
            }
            push @b, $_;
        }
        # filter out stuff that we already provide (there are cases where
        # we supply this information redundantly).
        push @mpi_precmd, @b;
    }
    push @mpi_precmd, split(' ', $mpiexec_extra_stanzas);
    push @mpi_precmd, split(' ', $ENV{'MPI_EXTRA_ARGS'}) if $ENV{'MPI_EXTRA_ARGS'};

    @mpi_precmd_single = @mpi_precmd;
    @mpi_precmd_lingen = @mpi_precmd;
    if (!$param->{'only_mpi'}) {
        push @mpi_precmd, '-n', $mpi_split[0] * $mpi_split[1];
    } else {
        push @mpi_precmd, '-n', $nh * $nv;
    }
    push @mpi_precmd_lingen, '-n', $lingen_mpi_split[0] * $lingen_mpi_split[1];
    push @mpi_precmd_single, '-n', 1;
}

# }}}

if ($mpi_needed) {
    if ($ENV{'DISPLAY'}) {
        print "## removing the DISPLAY environment variable, as it interacts badly with MPI startup\n";
        delete $ENV{'DISPLAY'};
    }
    if ($ENV{'SSH_AUTH_SOCK'}) {
        if ($ENV{'OAR_JOBID'}) {
            print "## removing the SSH_AUTH_SOCK environment variable, as it interacts badly with MPI startup, and OAR does not need it.\n";
            delete $ENV{'SSH_AUTH_SOCK'};
        } else {
            print "## WARNING: the environment variable SSH_AUTH_SOCK is set. If it so happens that you do *not* need it, you should unset it for faster job startup.\n";
        }
    }

    # openmpi seems to properly propagate the path to mpirun as the path
    # to be forwarded on the remote nodes, so no manual propagation of
    # LD_LIBRARY_PATH or --prefix is needed.

    # mkdir must not be marked fatal, because if the command terminates
    # without having ever tried to join in an mpi collective like
    # mpi_init(), there's potential for the mpirun command to complain.
    dosystem('-nofatal', @mpi_precmd, split(' ', "mkdir -p $wdir"))
        unless defined($random_matrix) || $wdir eq '/';
}

if ($main eq ':mpirun') {
    # we don't even put @main_args in, because we're tinkering with it
    # somewhat.
    dosystem(@mpi_precmd, @extra_args);
    exit 0;
}

if ($main eq ':mpirun_single') {
    # we don't even put @main_args in, because we're tinkering with it
    # somewhat.
    dosystem(@mpi_precmd_single, @extra_args);
    exit 0;
}

if ($main eq ':srun') {
    # we don't even put @main_args in, because we're tinkering with it
    # somewhat.
    s/mpiexec/srun/ for (@mpi_precmd);
    dosystem(@mpi_precmd, @extra_args);
    exit 0;
}

if ($main eq ':srun_single') {
    # we don't even put @main_args in, because we're tinkering with it
    # somewhat.
    s/mpiexec/srun/ for (@mpi_precmd);
    dosystem(@mpi_precmd_single, @extra_args);
    exit 0;
}

##################################################
### ok -- now @main_args is something relatively useful.
# print "main_args:\n", join("\n", @main_args), "\n";

push @main_args, splice @extra_args;

# Now we have one function for each of the following macroscopic steps of
# the program
#
# prep
# krylov
# lingen
# mksol
# gather

# {{{ Some pretty-printing
my $current_task;

my $terminal_colors = {
    BLACK	=> "\e[01;30m",
    RED 	=> "\e[01;31m",
    GREEN	=> "\e[01;32m",
    YELLOW	=> "\e[01;33m",
    BLUE	=> "\e[01;34m",
    VIOLET	=> "\e[01;35m",
    black	=> "\e[00;30m",
    red         => "\e[00;31m",
    green	=> "\e[00;32m",
    yellow	=> "\e[00;33m",
    blue	=> "\e[00;34m",
    violet	=> "\e[00;35m",
    normal      => "\e[0m",
};
$terminal_colors = {} if !defined($ENV{'TERM'}) || $ENV{'TERM'} !~ /^(xterm|screen|linux)/;

sub task_begin_message {
    my $blue = $terminal_colors->{'BLUE'} || '';
    my $normal = $terminal_colors->{'normal'} || '';
    print "## Entering task: ${blue}$current_task${normal}\n";
}

sub task_check_message {
    my $status = shift;
    return unless $my_verbose_flags->{'checks'} || $status eq 'error';
    my $normal = $terminal_colors->{'normal'} || '';
    my $color = {
        'ok' => $terminal_colors->{'green'} || '',
        'missing' => $terminal_colors->{'YELLOW'} || '',
        'error' => $terminal_colors->{'RED'} || '',
    };
    my @lines;
    push @lines, split(/^/m, $_) for @_;
    my $head = shift @lines;
    chomp($head);
    for (@lines) {
        chomp($_);
        $_ = "\t" . $_;
    }
    unshift @lines, $head;
    my @printme;
    push @printme, "## $color->{$status}Check for $current_task$normal: $_\n" for @lines;
    die @printme if $status eq 'error';
    print @printme;
}
# }}}

# Inspection code for what is currently in $wdir

# Some of the inspection code correspond to queries which are done
# several times, so it's better to cache the results as long as they are
# expected to remain valid.
# {{{ cache hash, and corresponding queries
my $cache = {
    # This corresponds to data which isn't supposed to expire with
    # immediate effect.
    # Here we'll find things which are inferred from reading the data
    # files in $wdir, but would be cumbersome to find otherwise.
};
sub expire_cache_entry {
    my $k = shift;
    delete $cache->{$k};
}
sub store_cached_entry {
    my ($k, $v) = @_;
    if (defined($cache->{$k}) && $cache->{$k} != $v) {
        die "Fatal error: conflict for cache entry $k: previously had $cache->{$k}, now detected $v";
    }
    $cache->{$k} = $v;
}


sub get_cached_leadernode_filelist {
    my $key = 'leadernode_filelist';
    my $opt = shift;
    my @x;
    if (defined(my $z = $cache->{$key})) {
        @x = @$z;
    } else {
        if (!$hostfile) {
            print STDERR "Listing files in $wdir via readdir\n";
            # We're running locally, thus there's no need to go to a remote
            # place just to run find -printf.
            my $dh;
            opendir $dh, $wdir;
            for (readdir $dh) {
                push @x, [$_, (stat "$wdir/$_")[7]];
            }
            closedir $dh;
        } else {
            my $foo = join(' ', @mpi_precmd_single, "find $wdir -maxdepth 1 -follow -type f -a -printf '%s %p\\n'");
            print STDERR "Listing files in $wdir via $foo\n";
            for my $line (`$foo`) {
                $line =~ s/^\s*//;
                chomp($line);
                $line =~ s/^(\d+)\s+//;
                my $s = $1;
                push @x, [basename($line), $s];
            }
        }
        $cache->{$key}=\@x;
    }
    if ($opt && $opt eq 'HASH') { # we could also use wantarray
        my $h = {};
        $h->{$_->[0]} = $_->[1] for @x;
        return $h;
    } else {
        return @x;
    }
}

sub rename_file_on_leader {
    my ($old, $new) = @_;
    if (!$hostfile) {
        # We're running locally, thus there's no need to go to a remote
        # place just to run find -printf.
        rename $old, $new;
    } else {
        system(join(' ', @mpi_precmd_single, "mv $old $new"));
    }
}


# {{{ get_cached_bfile -> check for balancing file.
sub get_cached_bfile {
    my $key = 'balancing';
    return undef if defined($random_matrix);
    if ($param->{$key}) {
        $cache->{$key}=$param->{$key};
    }
    if (defined(my $z = $cache->{$key})) {
        return $z;
    }
    # We're checking on the leader node only. Because the other nodes
    # don't really mind if they don't see the balancing file.
    # sense.
    my $pat;
    my $x = $matrix;
    $x =~ s{^(?:.*/)?([^/]+)$}{$1};
    $x =~ s/\.(?:bin|txt)$//;
    my $bfile = "$wdir/$x.${nh}x${nv}/$x.${nh}x${nv}.bin";
    if (!-f $bfile) {
        print STDERR "$bfile: not found\n";
        return undef;
    }
    $cache->{$key} = $bfile;
    return $bfile;
}
# }}}

sub get_cached_balancing_header {
    return undef if defined $random_matrix;
    my $key = 'balancing_header';
    if (defined(my $z = $cache->{$key})) { return @$z; }
    defined(my $balancing = get_cached_bfile) or confess "\$balancing undefined";
    sysopen(my $fh, $balancing, O_RDONLY) or die "$balancing: $!";
    sysread($fh, my $bhdr, 24);
    my @x = unpack("LLLLLL", $bhdr);
    my $zero = shift @x;
    die "$balancing: no leading 32-bit zero" unless $zero == 0;
    my $magic = shift @x;
    die "$balancing: bad file magic" unless $magic == 0xba1a0000;
    $cache->{$key} = \@x;
    close($fh);
    return @x;
}
sub get_nrows_ncols {
    my $key = 'nrows_ncols';
    if (defined(my $z = $cache->{$key})) {
        return @$z;
    }
    if (defined($random_matrix)) {
        my $nrows;
        my $ncols;
        my @tokens = split(',', $random_matrix);
        if ($tokens[0] =~ /^(\d+)/) {
            $nrows = $1;
            $ncols = $1;
        }
        if ($tokens[1] =~ /^(\d+)/) {
            $ncols = $1;
        }
        die "parameter random_matrix does not give nrows and ncols ??" unless defined($nrows) && defined($ncols);
        $cache->{$key} = [ $nrows, $ncols ];
        return @{$cache->{$key}};
    }


    if (defined(my $z = $cache->{'balancing_header'})) {
        my @x = @$z;
        shift @x;
        shift @x;
        $cache->{$key} = \@x;
        return @x;
    }
    if (defined(my $balancing = get_cached_bfile)) {
        my ($bnh, $bnv, $bnrows, $bncols) = get_cached_balancing_header;
        my @x = ($bnrows, $bncols);
        $cache->{$key} = \@x;
        return @x;
    }
    (my $mrw = $matrix) =~ s/(\.(?:bin|txt))$/.rw$1/;
    (my $mcw = $matrix) =~ s/(\.(?:bin|txt))$/.cw$1/;
    if ($mrw ne $matrix && $mcw ne $matrix && -f $mrw && -f $mcw) {
        my $nrows = ((stat $mrw)[7] / 4);
        my $ncols = ((stat $mcw)[7] / 4);
        my @x = ($nrows, $ncols);
        $cache->{$key} = \@x;
        return @x;
    }
    confess "Cannot find nrows and ncols ???";
}
# }}}

# {{{ List the V files in $wdir -- output is given as a hash y-range =>
# iterations found. We also return the file size, which has to be
# constant.
sub list_files_generic {
    my $files = {};
    my $filesize;
    my ($n, $pattern) = @_;
    for my $fileinfo (get_cached_leadernode_filelist) {
        my ($file, $size) = @$fileinfo;
        $file =~ /^$pattern$/ or next;
        my $n0 = scalar @-;
        die "Found $n0 < $n matches. Bad pattern $pattern ??\n" if $n0 < $n;
        # See man perlvar
        my @matches = map { 0+substr $file, $-[$_], $+[$_] - $-[$_] } (1..$n0-1);
        my @kmatches = splice(@matches, 0, $n);
        push @{$files->{join("..", @kmatches)}}, \@matches;
        $filesize = $size if !defined($filesize);
        if ($filesize != $size) {
            task_check_message 'error', "Inconsistency detected for the sizes of the ${pattern} files. We have seen at least $filesize and $size (last seen: $file, $size). Please fix.\n";
        }
    }
    return $files, $filesize;
}
sub list_vfiles {
    my ($f, $filesize) = list_files_generic(2, qr/V(\d+)-(\d+)\.(\d+)/);
    if ($filesize) {    # take the occasion to store it.
        for my $k (keys %$f) {
            $k =~ /^(\d+)\.\.(\d+)$/ or die "Bad key $k returned by list_files_generic";
            if ($2-$1 != $splitwidth) {
                die "Problem with the width of V files: $1-$2 is not good";
            }
        }
        # Note that it might happen that we haven't computed the
        # balancing file yet.
        my ($bnrows, $bncols) = get_nrows_ncols;
        my $N = $bncols > $bnrows ? $bncols : $bnrows;
        eval {store_cached_entry('nbytes_per_splitwidth', $filesize / $N);};
        die "Problem with the size of V files ($filesize bytes, $N rows):\n$@" if $@;
    }
    my $flatten = sub { local $_; map { shift @$_ } @_; };
    $f->{$_} = [&$flatten(@{$f->{$_}})] for keys %$f;
    @{$f->{$_}} = sort { $a <=> $b } @{$f->{$_}} for keys %$f;
    return $f;
}
sub lexcmp {
    my $a = shift;
    my $b = shift;
    for my $k (0..$#$a) {
        my $z = $a->[$k] <=> $b->[$k];
        return $z if $z;
    }
    return 0;
}

sub list_afiles {
    my ($f, $filesize) = list_files_generic(2, qr/A(\d+)-(\d+)\.(\d+)-(\d+)/);
    if ($filesize) {    # take the occasion to store it.
        for my $k (keys %$f) {
            $k =~ /^(\d+)\.\.(\d+)$/ or die "Bad key $k returned by list_files_generic";
            # We can tolerate $2-$1 == $n, because we still have acollect
            # around.
            if ($2-$1 != $splitwidth && $2-$1 != $n) {
                die "Problem with the width of A files: $1-$2 is not good";
            }
            my $length = $m * ($2-$1) * ($f->{$k}->[0]->[1] - $f->{$k}->[0]->[0]);
            eval { store_cached_entry('nbytes_per_splitwidth', $filesize / ($length / $splitwidth));};
            die "Problem with the size of A files:\n$@" if $@;
        }
    }
    @{$f->{$_}} = sort { lexcmp($a, $b) } @{$f->{$_}} for keys %$f;
    return $f;
}

sub list_sfiles {
    my ($f, $filesize) = list_files_generic(2, qr/S\.sols(\d+)-(\d+)\.(\d+)-(\d+)/);
    if ($filesize) {    # take the occasion to store it.
        for my $k (keys %$f) {
            $k =~ /^(\d+)\.\.(\d+)$/ or die "Bad key $k returned by list_files_generic";
            if ($2-$1 != $splitwidth) {
                die "Problem with the width of S files: $1-$2 is not good";
            }
        }
        my ($bnrows, $bncols) = get_nrows_ncols;
        my $N = $bncols > $bnrows ? $bncols : $bnrows;
        my $length = $N;
        eval { store_cached_entry('nbytes_per_splitwidth', $filesize / $length); };
        die "Problem with the size of S files:\n$@" if $@;
    }
    # Because our S pattern now includes the iteration range, we pick
    # $_->[1] for the identifier.
    # my $flatten = sub { local $_; map { $_->[1] } @_; };
    # $f->{$_} = [&$flatten(@{$f->{$_}})] for keys %$f;
    #
    my $flatten = sub { local $_; return map { $_->[0] => $_->[1] } @_; };
    for (keys %$f) {
        my %x = &$flatten(@{$f->{$_}});
        $f->{$_} = \%x;
    }
    # @{$f->{$_}} = sort { $a <=> $b } @{$f->{$_}} for keys %$f;
    return $f;
}
# }}}

# {{{ how many bytes per ($splitwidth) element ?
sub get_cached_nbytes_per_splitwidth {
    my $key = 'nbytes_per_splitwidth';
    if (defined(my $z = $cache->{$key})) { return $z; }
    # This function should practically never go beyond this point, as
    # presumably the list_* functions below, which all fill the cache as
    # a side effect, should have been run beforehand.
    list_vfiles;
    list_afiles;
    list_sfiles;
    if (defined(my $z = $cache->{$key})) { return $z; }
    die "Cannot find the number of bytes needed per finite field element";
}
# }}}


# {{{ inferring max iteration indices.
sub max_krylov_iteration {
    # read matrix dimension from the balancing header.
    my ($bnrows, $bncols) = get_nrows_ncols;
    my $length = $bncols > $bnrows ? $bncols : $bnrows;
    $length = int(($length+$m-1)/$m) + int(($length+$n-1)/$n);
    $length += 2 * int(($m+$n-1)/$m);
    $length += 2 * int(($m+$n-1)/$n);
    $length += 10;
    return $length;
}
# This one is of course much harder to guess.
sub max_mksol_iteration {
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    my $nbytes_per_splitwidth = get_cached_nbytes_per_splitwidth;
    my @x;
    for my $file (keys %$leader_files) {
        $file =~ /^F\.sols(\d+)-(\d+)\.(\d+)-(\d+)$/ or next;
        if ($2-$1 != $splitwidth) {
            die "Problem with the width of F files: $1-$2 is not good";
        }
        if ($4-$3 != $splitwidth) {
            die "Problem with the width of F files: $1-$2 is not good";
        }
        my $size = $leader_files->{$file};
        return $size / (($2-$1)*($4-$3)/$splitwidth*$nbytes_per_splitwidth);
    }
    die "can't find any F file, cannot infer the mksol max iteration";
} # }}}

sub interval_default {
    my $krylov_length = max_krylov_iteration;
    my $interval = 1;
    while($interval * $interval < $krylov_length) {
        $interval *= 2;
    }
    if ($interval < 64) { $interval = 64; }
    return $interval;
}

# {{{ task_common_run is just a handy proxy.
sub task_common_run {
    my $program = shift @_;
    expire_cache_entry 'leadernode_filelist';
    # Some arguments are relevant only to some contexts.
    #
    # We start with lingen, because it's slightly specific
    # take out the ones we don't need (and acollect shares some
    # peculiarities).
    @_ = grep !/^(skip_bw_early_rank_check|rebuild_cache|cpubinding|balancing.*|interleaving|matrix|mm_impl|mpi|thr)?=/, @_ if $program =~ /(lingen|acollect$)/;
    if ($program =~ /lingen/) {
        @_ = map { s/^lingen_mpi\b/mpi/; $_; } @_;
    } else {
        @_ = grep !/^lingen_mpi?=/, @_;
    }
    # @_ = grep !/lingen_threshold/, @_ unless $program =~ /lingen/;
    # @_ = grep !/lingen_mpi_threshold/, @_ unless $program =~ /lingen/;
    @_ = grep !/allow_zero_on_rhs/, @_ unless $program =~ /^lingen/;
    @_ = grep !/^save_submatrices?=/, @_ unless $program =~ /^(prep|krylov|mksol|gather)$/;
    # are we absolutely sure that lingen needs no matrix ?
    @_ = grep !/^ys=/, @_ unless $program =~ /(krylov|dispatch)$/;
    @_ = grep !/^solutions=/, @_ unless $program =~ /(?:mksol|gather)$/;
    @_ = grep !/^rhs=/, @_ unless $program =~ /(?:prep|gather|lingen.*|mksol)$/;
    @_ = grep !/(?:precmd|tolerate_failure)/, @_;

    if ($program =~ /bwccheck$/) {
        my @x_;
        while (defined($_=shift @_)) {
            if ($_ eq '--') {
                push @x_, $_, splice @_;
                last;
            }
            next unless /^(?:[mn]|prime|wdir)=/;
            push @x_, $_;
        }
        @_ = splice @x_;
    }

    $program="$bindir/$program";
    unshift @_, $program;

    if ($param->{'precmd'}) {
        unshift @_, split(' ', $param->{'precmd'});
    }

    if ($mpi_needed) {
        if ($program =~ /\/lingen[^\/]*$/) {
            unshift @_, @mpi_precmd_lingen;
        } elsif ($program =~ /\/(?:split|acollect|lingen|cleanup|bwccheck)$/) {
            unshift @_, @mpi_precmd_single;
        } elsif ($program =~ /\/(?:prep|secure|krylov|mksol|gather|dispatch)$/) {
            unshift @_, @mpi_precmd;
        } else {
            die "Don't know the parallel status of program $program ... ?";
        }
    } else {
        @_ = grep { !/^(mpi)=/ } @_;
    }

    eval { dosystem @_; };

    if ($@) {
        if (defined(my $tol = $param->{'tolerate_failure'})) {
            my $re = qr/$tol/;
            if ($program =~ /$re/) {
                print STDERR "Not aborting because $program matches tolerate_failure regexp $tol\n";
                return;
            }
        } else {
            die;
        }
    }
}
# }}}

# {{{ SUBTASKS are used by one or several tasks.

# {{{ subtask_krylov_todo - the hard checkpoint recovery work
sub subtask_krylov_todo {
    # For each sequence, this finds the "most advanced" V file. Being
    # advanced takes two things. First this means an existing file, and
    # the code here checks for this. But we also request that some extra
    # function returns true, and this check is provided by the caller.
    my $length = shift;
    my $morecheck = shift || sub{1;};

    my @ys;
    if (defined($random_matrix)) {
        my @ys=();
        my $dx = $simd;
        $dx +=  $simd if $param->{'interleaving'};
        for(my $x = 0; $x < $n ; $x += $dx) {
            my $y = $x + $dx;
            push @ys, "ys=$x..$y";
        }
        return @ys;
    }

    # We unconditionally do a check of the latest V checkpoint found in
    # $wdir, for all vectors.
    my @all_ys = map { [ $_*$splitwidth, ($_+1)*$splitwidth ] } (0..$n/$splitwidth - 1);

    my $vfiles = list_vfiles;
    my $vstarts = {};

    print "## complete range list for krylov: ".join(" ", map {"$_->[0]-$_->[1]";} @all_ys)."\n";

    for my $y (@all_ys) {
        my $yrange = $y->[0] . ".." . $y->[1];
        print "## V files for $yrange:";
        my @for_this_y =
            grep { $_ == 0 || &$morecheck($yrange, $_); } @{$vfiles->{$yrange}}
            or do { print " none\n"; next; };
        $vstarts->{$yrange} = $for_this_y[$#for_this_y];
        $vstarts->{$yrange} .= " (DONE!)" if $vstarts->{$yrange} >= $length;
        print " last $current_task checkpoint is $vstarts->{$yrange}\n";
    }

    my @todo;

    my $ys_comment='';
    if (@ys = grep { /^ys/ } @main_args) {
        @ys = map { /^(?:ys=)?(\d+)\.\.(\d+)$/; $_=[$1, $2]; } @ys;
        $ys_comment = ' (from command line)';
        # Do we also have a specification for the starting checkpoint ?
        if (my @ss = grep { /^start/ } @main_args) {
            my ($start) = grep { s/^start=(\d+)/$1/; } @ss;
            my $yrange = $ys[0]->[0] . ".." . $ys[0]->[1];
            my %ok = map { $_ => 1 } grep { $_ == 0 || &$morecheck($yrange, $_); } @{$vfiles->{$yrange}};
            if ($ok{$start}) {
                print "## Forcing ys=$yrange start=$start as requested on command line\n";
                push @todo, "ys=$yrange start=$start";
                return @todo;
            } else {
                my @avail = sort { $a <=> $b } keys %ok;
                task_check_message 'error', "Cannot start at ys=$yrange start=$start: missing checkpoints (available: @avail)";
            }
        }
    } else {
        my $y;
        for(my $x = 0; $x < $n ; $x = $y) {
            $y = $x + $simd;
            $y += $simd if $param->{'interleaving'};
            $y = $n if $y > $n;
            push @ys, [ $x, $y];
        }
        $ys_comment = "simd=$simd" if $simd > $splitwidth;
        $ys_comment .= ", interleaving" if $param->{'interleaving'};
        $ys_comment = " ($ys_comment)" if $ys_comment;
    }

    print "## range list for krylov: ".join(" ", map {"$_->[0]..$_->[1]";} @ys)."$ys_comment\n";

    my @impossible;
    for my $ab (@ys) {
        # most ys are double ranges, except maybe the last one.
        my ($a, $b) = @{$ab};
        my $yrange = $a . ".." . $b;
        my $start;
        # now all saved vectors are for width = $splitwidth ; given that
        # we have one single start= argument, we need to make sure that
        # all the sub-vectors have reached the same checkpoint.
        my $discard;
        for(my $x = $a; $x < $b ; $x += $splitwidth) {
            my $r = $x."..".($x + $splitwidth);
            my $s = $vstarts->{$r};
            if (!defined($s)) {
                task_check_message 'warning',
                "No starting vector for range $r"
                . " (sub-range from $yrange)\n";
                $discard = 1;
            } elsif (!defined($start)) {
                $start = $s;
            } elsif ($s != $start) {
                task_check_message 'warning',
                "Inconsistent start vectors for the "
                . " sub-ranges from $yrange"
                . " (found both $start and $s)\n";
                $discard = 1;
            }
        } 
        undef $start if $discard;
        if (!defined($start)) {
            print STDERR "## Can't schedule any work for $yrange, as *NO* checkpoint is here\n";
            push @impossible, $yrange;
            next;
        }
        next if $start =~ /DONE/;
        push @todo, "ys=$yrange start=$start";
    }
    if (!@todo && @impossible) {
        task_check_message 'error', "Cannot schedule remaining work. No checkpoints for ranges:\n", @impossible;
    }
    return @todo;
}
# }}}
# }}}

# {{{ prep

sub task_prep_missing_output_files {
    my @starts;
    for my $i (0..$n/$splitwidth-1) {
        my $x0 = $i*$splitwidth;
        my $x1 = ($i+1)*$splitwidth;
        push @starts, "V${x0}-${x1}.0";
    }
    my @missing;
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    unshift @starts, "X";
    for my $file (@starts) {
        push @missing, $file unless exists $leader_files->{$file};
    }
    if (@missing == @starts) { return "all", @missing; }
    elsif (@missing) { return "some", @missing; }
    else { return "none" };
}

sub task_prep {
    task_begin_message;
    # This prepares the starting vectors for block Wiedemann.
 
    # input files:
    #   - output files from previous tasks: dispatch
    #
    # output files, created if needed.
    #   - X
    #   - V${x0}-${x1}.0 ; the different starting vectors.

    my ($status, @missing) = task_prep_missing_output_files;

    if ($status eq 'none') {
        task_check_message 'ok', "All output files for $current_task have been found, good";
        return;
    }

    if ($status eq 'some') {
        task_check_message 'error', "Missing output files for $current_task", @missing, "We don't fix this automatically. Please investigate.";
    }

    task_check_message 'missing', "none of the output files for $current_task have been found, need to run $current_task now. We want to create files:", @missing;

    if ($prime == 2) {
        task_common_run('prep', @main_args);
    } else {
        # The prime case is somewhat different. We'll generate random
        # data, that's it. We might prepend the RHS if it exists.
        task_common_run('prep', @main_args, "ys=0..$splitwidth");
    }

} # }}}


# {{{ secure -- this is now just a subtask of krylov or mksol
sub subtask_secure {
    return if $param->{'skip_online_checks'};
    my $wanted_stops={};
    if (defined(my $x = $param->{'check_stops'})) {
        my @x = split(',', $x);
        $wanted_stops->{$_}=1 for @x;
        if (defined(my $x = $param->{'interval'})) {
            $wanted_stops->{$x} = 1;
        }
    } else {
        # We're taking a pretty aggressive approach of always storing
        # extra checkpoints that correspond to the interval value + a
        # small offset See #30025
        my $x = $param->{'interval'};
        $x = interval_default unless defined $x;
        $wanted_stops->{$x}=1;
        my $a = int($x/2);
        $a = 16 if $a > 16;
        if ($a) {
            $wanted_stops->{$x+$a}=1;
            $wanted_stops->{$a}=1;
        }
    }
    my $leader_files = get_cached_leadernode_filelist 'HASH';

    my $mustrun = 0;
    if (scalar keys %$wanted_stops) {
        $wanted_stops->{0}=1;
        for (keys %$wanted_stops) {
            my @wnames=("Cv0-$splitwidth.$_");
            push @wnames, "Cd0-$splitwidth.$_" if $_ > 0;
            for my $name (@wnames) {
                next if $leader_files->{$name};
                task_check_message 'missing', "missing check vector $name\n";
                $mustrun = 1;
            }
        }
    } else {
        # we can only check for the existence of _some_ check vector.
        unless (grep { /^Cv0-(\d+)\.\d+$/ && ($1 == $splitwidth) && !/^Cv0-\d+\.0$/ } keys %$leader_files) {
            task_check_message 'missing', "no check vector (Cv) found\n";
            $mustrun = 1;
        }
        unless (grep { /^Cd0-(\d+)\.\d+$/ && ($1 == $splitwidth) } keys %$leader_files) {
            task_check_message 'missing', "no check vector (Cd) found\n";
            $mustrun = 1;
        }
    }
    for my $name ("Cr0-$splitwidth.0-$splitwidth", "Ct0-$splitwidth.0-$m") {
        next if $leader_files->{$name};
        task_check_message 'missing', "missing check vector $name\n";
        $mustrun = 1;
    }
    if ($mustrun) {
        my $cs = "check_stops=" . join(",", sort { $a <=> $b } keys %$wanted_stops);
        task_common_run('secure', @main_args, $cs);
    } else {
        my $x = join(", ", sort { $a <=> $b } keys %$wanted_stops);
        task_check_message 'ok', "All auxiliary files for checkpointing are here, good (wanted: $x).\n";
    }
}
# }}}

sub task_dispatch {
    task_begin_message;
    return if defined $random_matrix;
    if (defined(get_cached_bfile)) {
        task_check_message 'ok', "All balancing files are present, good.\n";
        return;
    }
    # we do need to keep the ys, because those control the width that get
    # passed to the cache building procedures.
    task_common_run('dispatch', @main_args);
}

sub task_secure {
    task_begin_message;
    if (defined($random_matrix)) {
        task_check_message 'ok', "No checkpoint verification with random_matrix\n";
    } else {
        subtask_secure;
    }
}

sub task_safety_check_krylov {
    task_begin_message;
    my @args = @_;

    my $leader_files = get_cached_leadernode_filelist 'HASH';

    my @Cvs = sort { $a <=> $b } map { /\.(\d+)$/ && $1 } grep { /^Cv0-(\d+)\.(\d+)$/ && ($1 == $splitwidth) && $2 } keys %$leader_files;

    # refuse to check against C*.0, this gives false negatives (see
    # #30025).
    shift @Cvs if (@Cvs && $Cvs[0] == 0);

    if (!@Cvs) {
        task_check_message 'warning', "No check files present, cannot verify that the starting point is consistent";
    }
    
    my $ys;    for (@args) { $ys=$1 if    /^(?:ys=)?(\d+)/;    }
    my $start; for (@args) { $start=$1 if /^(?:start=)?(\d+)/; }

    my $doable_checks = {};

    my $yrange = $ys . ".." . ($ys + $splitwidth);
    my $vfiles = list_vfiles;
    my $v_found_checkpoints = {};
    $v_found_checkpoints->{$_}=1 for @{$vfiles->{$yrange}};
    my $V = sprintf("V%d-%d.%d", $ys, $ys+$splitwidth, $start);
    for(my $i = 0 ; $i < scalar @Cvs ; $i++) {
        my $Ci = sprintf("Cv%d-%d.%d", 0, $splitwidth, $Cvs[$i]);
        for(my $j = $i + 1 ; $j < scalar @Cvs ; $j++) {
            my $Cj = sprintf("Cv%d-%d.%d", 0, $splitwidth, $Cvs[$j]);
            my $b = $start + $Cvs[$i] - $Cvs[$j];
            next unless $v_found_checkpoints->{$b};
            my $Vx = sprintf("V%d-%d.%d", $ys, $ys+$splitwidth, $b);
            print "## doable check against iteration $b, using check stops $Cvs[$i] and $Cvs[$j]\n";
            $doable_checks->{$V} = 1;
            $doable_checks->{$Vx} = 1;
            $doable_checks->{$Ci} = 1;
            $doable_checks->{$Cj} = 1;
        }
    }

    if (!keys %$doable_checks) {
        task_check_message 'warning', "Found no way to verify that the starting point is consistent";
        return;
    }

    @args = grep { !/^ys/ && !/^start/ } @args;
    task_common_run('bwccheck', @args, '--', keys %$doable_checks);
    task_check_message 'ok', "Consistency check for starting vector $V was successful";
}


# {{{ krylov
sub task_krylov {
    task_begin_message;

    # input files:
    #   - output files from previous tasks: prep
    #     More accurately, we need V${x0}-${x1}.$i for some iteration i.
    #
    # if ys is found in the command line, then we focus on that sequence
    # specifically. Otherwise, we process all sequences one after
    # another.
    #
    # side-effect files, which have no impact on whether this step gets
    # re-run or not:
    #   - A${x0}-${x1}.* ; dot product files which are fed to lingen.
    #   These are write-only as far as this task is concerned.

    subtask_secure unless defined $random_matrix;

    my $length = max_krylov_iteration;
    print "## krylov max iteration is $length\n";

    my @todo = subtask_krylov_todo $length;

    if (!@todo) {
        # Note that we haven't checked for the A files yet !
        task_check_message 'ok',
                    "All krylov tasks are done, good.\n";
        return;
    }

    task_check_message 'missing',
                "Pending tasks for krylov:\n" . join("\n", @todo);
    for my $t (@todo) {
        # take out ys and start from main_args, put the right ones in place if
        # needed. (note that if both were specified, this is essentially
        # putting the very same parameters back in place, but with the
        # benefit of having performed a check inbetween).
        my @args = grep { !/^ys/ && !/^start/ && !/^max_ram/ } @main_args;
        push @args, split(' ', $t);
        task_safety_check_krylov @args, split(' ', $t);
        task_common_run 'krylov', @args;
    }
} # }}}

# {{{ lingen

sub task_lingen_input_errors {
    # krylov_length is not documented. It's just here to fool lingen and
    # let it believe that it does have the full data set.
    my $h = shift;
    my $length = $param->{'krylov_length'} || max_krylov_iteration;
    my $afiles = list_afiles;
    my $astrings = {};
    my $ok = 1;
    my $nb_afiles = 0;
    my $maxz1;
    my $minz1;
    for my $k (keys %$afiles) {
        my @strings;
        my ($z0, $z1);
        for my $seg (@{$afiles->{$k}}) {
            $nb_afiles++;
            if (defined($z1)) {
                if ($seg->[0] == $z1) {
                    $z1 = $seg->[1];
                } else {
                    push @strings, [$z0, $z1];
                    $z0 = undef;
                    $z1 = undef;
                }
                next;
            }
            ($z0, $z1) = @$seg;
        }
        $ok = 0 unless scalar @strings == 0 && $z0 == 0 && $z1 >= $length;
        push @strings, [$z0, $z1] if defined $z1;
        $maxz1 = $z1 unless defined($maxz1) && $z1 < $maxz1;
        $minz1 = $z1 unless defined($minz1) && $z1 > $minz1;
        $astrings->{$k} = \@strings;
    }

    $ok = 0 if $minz1 != $maxz1;
    my $rlength = $minz1;

    if (!$ok) {
        my @errors;
        push @errors, "Incomplete set of A files found";

        for my $k (sort {$a cmp $b} keys %$afiles) {
            push @errors, "$k, computed range(s): " . join(" ", map { $_->[0]."..".$_->[1]; } @{$astrings->{$k}}) . "\n";
        }
        return @errors;
    }
    $h->{'concatenated_A'} = "A0-$n.0-$rlength";
    $h->{'need_acollect'} = $nb_afiles > 1;
    return ();
}

sub task_lingen {
    task_begin_message;
    
    # input files:
    #   - A.*-*.*-* ; all A files created by krylov.
    #     -> this task accumulates them in a single file.
    # output_file
    #   - F.sols*-*.*-* ; all polynomial entries of the generating
    #     matrix.

    # The function task_lingen_input_errors has a few side-effects,
    # returned back in the hash $h.
    my $h = {};
    if (my @errors = task_lingen_input_errors $h) {
        task_check_message 'error', @errors;
    }

    my $concatenated_A = $h->{'concatenated_A'};

    if ($h->{'need_acollect'}) {
        task_check_message 'missing', "Running acollect to create $concatenated_A";
        task_common_run("acollect", @main_args, "--remove-old");
    }

    # We expect lingen to split its output into pieces of predictable
    # size.
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    my @missing;
    my @expected;
    for my $i (0..$n/$splitwidth-1) {
        my $sol = sprintf("sols%d-%d", $i*$splitwidth, ($i+1)*$splitwidth);
        for my $j (0..$n/$splitwidth-1) {
            my $j0 = $splitwidth * $j;
            my $j1 = $splitwidth + $j0;
            my $f = "F.$sol.${j0}-${j1}";
            push @expected, $f;
            push @missing, $f unless exists $leader_files->{$f};
            if ($j1 <= $nrhs) {
                $f .= ".rhs";
                push @expected, $f;
                push @missing, $f unless exists $leader_files->{$f};
            }
        }
    }
    if (@missing == 0) {
        task_check_message 'ok', "lingen result found, good.";
        return;
    } elsif (@missing < @expected) {
        task_check_message 'error', "Incomplete lingen output found. Missing files:\n" , @missing;
    }

    task_check_message 'missing', "lingen has not run yet. Running now.";
    # Now run lingen itself. Which binary we'll run is not totally
    # obvious though.
    # NOTE: It may be worthwhile to run specifically this step, but
    # with adapted mpi and thr parameters.
    my @args;
    push @args, "split-output-file=1";
    push @args, "afile=$concatenated_A";
    push @args, "ffile=F";
    push @args, grep { /^(?:mn|m|n|wdir|prime|rhs)=/ || /allow_zero_on_rhs/ || /max_ram/ } @main_args;
    if (!$mpi_needed && ($lingen_mpi_split[0]*$lingen_mpi_split[1] != 1)) {
        print "## non-MPI build, avoiding multi-node lingen\n";
        # We keep thr=
        @args = grep { !/^(mpi)=/ } @args;
    }
    push @args, grep { /^verbose_flags=/ } @main_args;
    if (! -f "$wdir/$concatenated_A.gen") {
        if ($prime == 2) {
            push @args, "tuning_thresholds=recursive:128,ternary:6400,cantor:6400,notiming:0";
            task_common_run("lingen_b64", @args);
        } else {
            push @args, "tuning_thresholds=recursive:10,flint:10,notiming:0";
            task_common_run("lingen_pz", @args);
        }
    } else {
        task_check_message 'ok', "lingen already has .gen file, good.";
    }
}
# }}}

# {{{ mksol
sub task_mksol {
    task_begin_message;

    # input files:
    #   - V*-*.<some iteration>, e.g. maybe the starting one generated by
    #     prep, or any other.
    #   - F.sols*-*.*-* ; all polynomial entries of the generating
    #     matrix.
    #
    # Note that the computation will start at the first iteration number
    # for which a V file is here, and no S file yet (for any of the
    # solutions considered, which is a number limited by $nrhs if
    # specified).
    #
    # side-effect files, which have no impact on whether this step gets
    # re-run or not:
    #   - S.sols*-*.*-*.* ; partial sum which are fed to gather.
    #   These are write-only as far as this task is concerned.

    subtask_secure unless defined $random_matrix;

    # We choose to require the S files for considering checkpoints as
    # valid. This is because in all likelihood, the V files from krylov
    # will still be there, so even though this is a hint as to what can
    # be started right now, we can't really use it...
    my $sfiles = list_sfiles;
    # $sfiles->{$_} = eval { my %x=map {$_=>1;} @{$sfiles->{$_}};\%x;} for keys %$sfiles;
    #
    my $length = eval { max_mksol_iteration; };
    if ($@) {
        task_check_message 'error', "Lingen output files missing", $@, "Please run lingen first.";
        die "abort";
    }

    print "## mksol max iteration is $length\n";

    my @solutions=@{$param->{'solutions'}};
    print "## main solution ranges for mksol ".join(" ", @solutions)."\n";

    my @all_solutions = map
                        { $_=sprintf("%d-%d",
                                $_*$splitwidth,
                                ($_+1)*$splitwidth);
                        } (0..$n/$splitwidth-1);

    my %solutions_importance;
    $solutions_importance{$_}=0 for (@all_solutions);
    for (@solutions) {
        /^(\d+)-(\d+)$/ or die;
        my $x=$1;
        my $y;
        my $z=$2;
        for( ; $x < $z ; $x = $y) {
            $y = $x + $splitwidth;
            my $s = "$x-$y";
            die "No solution file defined for $s: we have only @all_solutions" unless defined $solutions_importance{$s};
            $solutions_importance{$s}=1;
        }
    }
    my @todo;
    my @only_start = grep(/^start=\d+$/, @main_args);

    if (@only_start) {
        die "start says: @only_start" unless scalar @only_start == 1;
        die "solution says: @solutions" unless scalar @solutions == 1;
        @todo = "@only_start solutions=@solutions";
        task_check_message 'ok', "Command line imposes one specific subtask @todo";
    } else {
        for my $s (@solutions) {
            $s =~ /^(\d+)-(\d+)$/ or die;
            my $x=$1;
            my $y;
            my $z=$2;
            my $all_optional=1;
            my $n_common;
            # scan all subranges.
            for( ; $x < $z ; $x = $y) {
                $y = $x + $splitwidth;
                my $s = "$x-$y";
                my $optional = ($solutions_importance{$s} == 0);
                $all_optional = 0 unless $optional;
                # find latest checkpoint for that subrange
                my $n = 0;
                $s =~ /^(\d+)-(\d+)$/ or die;
                my $graph = $sfiles->{"$1..$2"};
                while (defined(my $e = $graph->{$n})) {
                    $n = $e;
                }
                my $msg = "## S files for $s --> last mksol checkpoint is $n";
                $msg .= " (DONE!)" if $n >= $length;
                $msg .= " (optional)" if $optional;
                print "$msg\n";
                if (!defined($n_common)) {
                    $n_common = $n;
                } elsif ($n_common != $n) {
                    task_check_message 'warning',
                    "Inconsistent latest vectors for the "
                    . " sub-ranges from $s"
                    . " (found both $n_common and $n)\n";
                }
            }
            if ($n_common < $length) {
                push @todo, "solutions=$s start=$n_common" unless $all_optional;
            }
        }
        if ($@) {
            task_check_message 'error', "Failure message while checking $current_task files";
        }
        if (!@todo) {
            # Note that we haven't checked for the A files yet !
            task_check_message 'ok', "All $current_task tasks are done, good.\n";
            return;
        }

        task_check_message 'missing',
                    "Pending tasks for $current_task:\n" . join("\n", @todo);
    }

    for my $t (@todo) {
        # take out ys from main_args, put the right one in place if
        # needed.
        # print "main_args: @main_args\n";
        my @args = grep { !/^(ys|n?rhs|start|max_ram)/ } @main_args;
        push @args, split(' ', $t);

        task_common_run 'mksol', @args;
    }
}
# }}}

# {{{ gather
sub task_gather {
    task_begin_message;

    # input files:
    #   - S.sols*-*.*-*.* ; partial sums computed by mksol

    my @missing;
    my @todo;
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    my @solutions=@{$param->{'solutions'}};
    my @all_solutions = map
                        { $_=sprintf("%d-%d",
                                $_*$splitwidth,
                                ($_+1)*$splitwidth);
                        } (0..$n/$splitwidth-1);
    my %solutions_importance;
    $solutions_importance{$_}=0 for (@all_solutions);
    for my $s (@solutions) {
        $s =~ /^(\d+)-(\d+)$/ or die;
        my @loc_found;
        my @loc_notfound;
        my $x=$1;
        my $y;
        my $z=$2;
        for( ; $x < $z ; $x = $y) {
            $y = $x + $splitwidth;
            my $ls = "$x-$y";
            die "No solution file defined for $s: we have only @all_solutions" unless defined $solutions_importance{$ls};
            $solutions_importance{$ls}=1;

            my $kfile = "K.sols$ls.0";
            if (exists $leader_files->{$kfile}) {
                push @loc_found, $kfile;
            } else {
                push @loc_notfound, $kfile;
            }
        }
        if (@loc_found && @loc_notfound) {
            task_check_message 'error', "not all files are consistent for solutions=$s : found @loc_found, missing @loc_notfound";
        }
        next if @loc_found;
        push @missing, @loc_notfound;
        if ($param->{'interleaving'}) {
            # Gather doesn't support interleaving at all, so let's do two
            # distinct executions.
            $s =~ /^(\d+)-(\d+)$/ or die;
            my $x=$1;
            my $y=int(($1+$2)/2);
            my $z=$2;
            push @todo, "solutions=$x-$y";
            push @todo, "solutions=$y-$z";
        } else {
            push @todo, "solutions=$s";
        }
    }
    if (@missing == 0) {
        task_check_message 'ok', "All solution files produced by gather seem to be present, good.";
        return;
    }
    task_check_message 'missing', "Need to run gather to create the files: @missing\n";
    @missing=();
    # {{{ Print the number of files found in a matrix.
    my $sfiles = list_sfiles;
    my $maxmksol = eval { max_mksol_iteration; };
    if ($@) {
        task_check_message 'error', "Lingen output files missing", $@, "Please run lingen first.";
    }
    print "## mksol max iteration is $maxmksol\n";

    my $cmat = {};
    for my $k (keys %$sfiles) {
        $k =~ /^(\d+)\.\.(\d+)$/ or die;
        my $graph = $sfiles->{$k};
        my $key = "sols$1-$2";
        $cmat->{$key}= scalar keys %$graph;
        # find max advanced point.
        my $n = 0;
        while (defined(my $e = $graph->{$n})) {
            $n = $e;
        }
        $cmat->{$key}.= "*" if $n < $maxmksol;
    }
    my $c0 = 4;
    my $c = 0;

    for my $s (@solutions) {
        $s =~ /^(\d+)-(\d+)$/ or die;
        my $x=$1;
        my $y;
        my $z=$2;
        for( ; $x < $z ; $x = $y) {
            $y = $x + $splitwidth;
            my $ls = "$x-$y";
            my $key = "sols$ls";
            my $optional = ($solutions_importance{$ls} == 0);
            my $l = 4+length($ls); $c0 = $l if $l > $c0;
            my $n = $cmat->{$key} || 'NONE'; $cmat->{$key} = $n;
            $l = length($n); $c = $l if $l > $c;
            next if $optional;
            push @missing, "S.sols$s" if $n eq 'NONE' || $n =~ /\*$/; 
        }
    }
    print "## Number of S files found:\n";
    for my $s (@solutions) { 
        $s =~ /^(\d+)-(\d+)$/ or die;
        my $x=$1;
        my $y;
        my $z=$2;
        for( ; $x < $z ; $x = $y) {
            $y = $x + $splitwidth;
            my $ls = "$x-$y";
            my $optional = ($solutions_importance{$ls} == 0);
            print "##    " . sprintf("%${c0}s","$ls") . " | "
                . sprintf("%${c}s", $cmat->{"sols$ls"})
                . ($optional ? " (optional)" : "")
                . "\n";
        }
    }
    # }}}

    my $rhs_companion;
    if ($param->{'rhs'}) {
        for my $s (@solutions) {
            for my $j (0..$n/$splitwidth - 1) {
                my $y0 = $j * $splitwidth;
                my $y1 = $y0 + $splitwidth;
                next if $y0 >= $nrhs;
                my $f = "F.sols$s.$y0-$y1.rhs";
                push @missing, $f unless exists($leader_files->{$f});
            }
        }
    }

    if (@missing) {
        task_check_message 'error', "Missing files for $current_task:", @missing;
    }

    task_check_message 'ok', "All required files for gather seem to be present, good.";

    my @args = grep { !/^ys/ } @main_args;
    for my $t (@todo) {
        task_common_run 'gather', @args, $t;
    }
}
# }}}

# {{{ cleanup -- For p=2, this extra step produces a RREF solution.
sub task_cleanup {
    # This is only for p=2 for the moment.
    return if $prime ne 2;
    task_begin_message;

    my @missing;
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    my $per_solution_files = {};
    my $err;
    my @todo;
    my @solutions=@{$param->{'solutions'}};
    for my $sol (@solutions) {
        $sol =~ /^(\d+)-(\d+)$/ or die;
        my $x=$1;
        my $y;
        my $z=$2;
        for( ; $x < $z ; $x = $y) {
            $y = $x + $splitwidth;
            my $s = "$x-$y";

            my $wfile = "W.sols$s";
            if (exists($leader_files->{$wfile})) {
                task_check_message 'ok', "Solution $s is already computed in file $wfile, good.\n";
                next;
            }

            my @ks =
                sort {
                    basename($a)=~/\.(\d+)$/; my $xa=$1;
                    basename($b)=~/\.(\d+)$/; my $xb=$1;
                    $xa <=> $xb;
                }
                grep {
                    /^(.*)\.\d+$/ && $1 eq "K.sols$s";
                }
                (keys %$leader_files);

            if (!@ks) {
                task_check_message 'error', "No files created by gather for solution $sol";
                $err=1;
            } else {
                task_check_message 'missing', "Will run cleanup to compute $wfile from the files:", @ks;
            }
            push @todo, [$s, $wfile, @ks];
        }
    }
    for my $klist (@todo) {
        my $sol = shift @$klist;
        my @x = map { "$wdir/$_"; } @$klist;
        my $wfile = shift @x;
        $sol =~ /^(\d+)-(\d+)$/ or die "$sol: could not parse";
        my $nsols = $2-$1;
        task_common_run('cleanup', "--ncols", $nsols, "--out", $wfile, @x);
    }
    if (scalar @todo == 1 && (!-f"$wdir/W" || ((stat "$wdir/W")[7] lt (stat "$wdir/W.sols$solutions[0]")[7]))) {
        print STDERR "## Providing $wdir/W as an alias to $wdir/W.sols$solutions[0]\n";
        symlink "W.sols$solutions[0]", "$wdir/W";
    }
}
# }}}

my $tasks = {
    prep	=> \&task_prep,
    krylov	=> \&task_krylov,
    # also keep "secure" and "dispatch", but just as handy proxies.
    secure	=> \&task_secure,
    dispatch	=> \&task_dispatch,
    lingen	=> \&task_lingen,
    mksol	=> \&task_mksol,
    gather	=> \&task_gather,
    cleanup	=> \&task_cleanup,
};

my @complete = qw(
        prep
        krylov
        lingen
        mksol
        gather
        cleanup
    );

my @tasks_todo=();

if ($main eq ':complete') {
    @tasks_todo=@complete;
} else {
    @tasks_todo=($main);
}

for (@tasks_todo) {
    $current_task = $_;
    if ($stop_at_step && $current_task eq $stop_at_step) {
        print "Exiting early, because of stop_at_step=$stop_at_step\n";
        last;
    }
    if (defined($tasks->{$current_task})) {
        &{$tasks->{$current_task}}(@main_args);
    } else {
        die "No task by that name: $current_task";
    }
}
