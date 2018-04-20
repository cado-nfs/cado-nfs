#!/usr/bin/perl -w
#
use strict;
use warnings;

# example code:
# machine=grisou ; for c in `git rev-list b074914913bff3732222ca8b13b5c80afa8e43da^..HEAD | tac` ; do c=`git rev-parse --short $c`; ./output_to_tikz_timechart.pl `ls -rt  rsa220/cado.$c.rsa220.gantt*${machine}*out | tail -1`  > chart-out.tex && pdflatex chart && \mv chart.pdf chart-$machine-$c.pdf ; done ; pdfjoin $(for c in `git rev-list b074914913bff3732222ca8b13b5c80afa8e43da^..HEAD | tac` ; do c=`git rev-parse --short $c`; echo chart-$machine-$c.pdf ; done) -o charts-$machine-by-thread.pdf

my $slicings={};

my $by_bucket = 0;
my $cropfraction=400;
my $merge_all_fib_levels=0;

sub display_slicing {
    my ($k,$S)=@_;
    my ($ns, $nsp, $ni, $w) = @{$S};
    print "Slicing for $k: $ns slices:\n";
    print "\tpart 0: $ni->[0] ideals, weight $w->[0]\n";
    print "\tpart 1: $ni->[1] ideals, weight $w->[1], $nsp->[1] slices\n";
    print "\tpart 2: $ni->[2] ideals, weight $w->[2], $nsp->[2] slices\n";
}
sub parse_new_slicing {
    my ($side, $K) = @_;
    my $ideals = [(0)x3];
    my $weight = [(0)x3];
    my $nslices = [(0)x3];
    my $nslices_total = 0;
    my $cpart;
    while (<>) {
        next if /^# Number of primes/;
        /^# Number of prime ideals.*part (\d+) = (\d+)/ && do {
            $ideals->[$1] = $2;
            next;
        };
        /^# Weight of primes.*part (\d+) = ([\d\.]+)/ && do {
            $weight->[$1] = $2;
            next;
        };
        /^# slices for side-\d part (\d+)/ && do { $cpart=$1; next; };
        /^# \[side-\d (\d+)\]/ && do {
            die unless $nslices_total == $1;
            $nslices->[$cpart]++;
            $nslices_total++;
            next;
        };
        next if /slice.*overflows.*Trying/;
        next if /to make sure/;
        my $S = [$nslices_total, $nslices, $ideals, $weight];
        $K = "side $side, $K";
        $slicings->{$K}=$S;
        # display_slicing($K, $S);
        return;
    }
}

sub ship_chart {
    my ($desc, $tab) = @_;
    my ($sq,$kind,$level,$side,$narr) = @$desc;
    my ($tmin, $tmax);
    return unless @$tab;
    my $NS = scalar @$tab;
    my $dt=0;
    my $thrmax=0;
    for my $S (@$tab) {
        my ($side, $aidx, $slice, $thr, $t0, $t1, $w, $hits, $time) = @$S;
        $tmin = $t0 unless defined($tmin) && $t0 > $tmin;
        $tmax = $t1 unless defined($tmax) && $t1 < $tmax;
        $thrmax = $thr unless $thr < $thrmax;
        $dt += $t1 - $t0;
    }
    my $nthr = $thrmax + 1;
    $dt = $dt / $NS;
    # We want to plot for a fraction of the time that is no less than
    # $cropfraction times the average per-slice time.
    my $tmax_cap = $tmin + $cropfraction * $dt;
    if ($tmax_cap > $tmax) { $tmax_cap = $tmax; }
    print "\\newpage\n";
    my $tdesc = "level-$level buckets";
    my $arr_desc = "$narr arrays";
    if ($level == 2) {
        $tdesc .= " on side $side";
    } else {
        $tdesc .= " on both sides";
        $arr_desc = "2*$arr_desc";
    }
    print "\n\n\\noindent\\textbf{$tdesc ($sq; $nthr threads, $NS slices on $arr_desc)}\n\n";
    print "\\begin{center}\n";
    print "\\begin{adjustbox}{max totalsize={\\textwidth}{.9\\textheight},center}\n";
    my $ncol;
    if ($by_bucket) {
        $ncol = $nthr;
    } else {
        $ncol = $narr;
        if ($level == 1) {
            $ncol += $narr;
        }
    }
    print "\\resetcolorseries[$ncol]{rainbow}\n";
    print "\\begin{tikzpicture}\n";
    for my $S (@$tab) {
        my ($cside, $aidx, $slice, $thr, $t0, $t1, $w, $hits, $time) = @$S;
        next if $t0 >= $tmax_cap;
        my $x0 = 20 * ($t0 - $tmin) / ($tmax_cap - $tmin);
        my $x1 = 20 * ($t1 - $tmin) / ($tmax_cap - $tmin);
        # our current value for microseconds() is useless, it's
        # counted globally...
        # # my $pcpu = 0;
        # # if ($t1 > $t0) { $pcpu = 1.0e-6*$time / ($t1-$t0); }
        # # $pcpu = int(10*$pcpu+0.5)*0.1;
        $x0 = sprintf("%.3f", $x0);
        $x1 = sprintf("%.3f", $x1);
        if ($by_bucket) {
            if ($level == 1) {
                $aidx += $cside * $narr;
            }
            my $y0 = -$aidx;
            my $y1 = -$aidx - 1;
            $y0 = $y0/2;
            $y1 = $y1/2;
            print "\\fill[fill={rainbow!![$thr]}] ($x0,$y0) rectangle ($x1,$y1);\n";
        } else {
            if ($level == 1) {
                $aidx += $cside * $narr;
            }
            my $y0 = -$thr;
            my $y1 = -$thr - 1;
            $y0 = $y0/2;
            $y1 = $y1/2;
            print "\\fill[fill={rainbow!![$aidx]}] ($x0,$y0) rectangle ($x1,$y1);\n";
        }
    }
    print "\\end{tikzpicture}\n\\end{adjustbox}\n\\end{center}\n\n";
    printf "Total time taken: %.2f real (\$t_{\\text{max}}-t_{\\text{min}}\$). Time displayed here: %.2f\n\\bigskip\n\n", $tmax-$tmin, $tmax_cap-$tmin;
}

sub parse_a_special_q {
    my $sq = shift @_;
    $sq=~s/;//;
    my $stack;
    my @ship=();
    while (<>) {
        /diagnosis for (\d)(\w) buckets on side (\d) \((\d+) arrays/ && do {
            # depending on the choice we made, $narr is either the same
            # as the number of threads, or one more, or even more.
            my $kind = "$1$2";
            my $level = $1;
            my $side = $3;
            my $narr = $4;
            my $tab;
            if (defined($stack)) {
                if ($merge_all_fib_levels) {
                    $tab = $stack->[1];
                } elsif ($stack->[0]->[2] != $level) {
                    push @ship, $stack;
                    $stack = undef;
                } else {
                    $tab = $stack->[1];
                }
            }
            if (!defined($stack)) {
                $tab=[];
                $stack=[[$sq,$kind,$level,$side,$narr],$tab];
            }
            my ($tmin, $tmax);
            for my $aidx (0..$narr-1) {
                defined($_=<>) or die;
                /side-$side array #$aidx processed (\d+).*(\d+) updates/ or die;
                my $nslices_arr = $1;
                my $hits_arr = $1;
                for my $sidx (0..$nslices_arr-1) {
                    # sidx is essentially meaningless.
                    defined($_=<>) or die;
                    /#\s*$kind side $side B $aidx/ or die;
                    my ($slice, $thr, $t0, $t1, $w, $hits, $time);
                    /slice (\d+)/ && do { $slice=$1; };
                    /thr(?:ead)? (\d+)/ && do { $thr=$1; };
                    /t0 ([\d\.]+)/ && do { $t0=$1; };
                    /t1 ([\d\.]+)/ && do { $t1=$1; };
                    /ecost ([\d\.]+)/ && do { $w=$1; };
                    /hits (\d+)/ && do { $hits=$1; };
                    /time ([\d\.]+)/ && do { $time=$1; };
                    push @$tab, [$side, $aidx, $slice, $thr, $t0, $t1, $w, $hits, $time];
                }
            }
            # get something to eat.
            next;
        };
        next if /^# small sieve/;
        next if /^# i=/;
        next if /^# DUP/;
        next if /^-?\d+/;       # relations.
        if (defined($stack)) {
            push @ship, $stack;
            $stack = undef;
        }
        next if /^# I=/;
        next if /^# f_\d'/;
        next if /^# (Reserving|Allocating|Reallocating)/;
        next if /^#\s+\d+\.\d+/;
        if (/redoing/) {
            print STDERR "Canceling $sq (overfull buckets)\n";
            @ship=();
            while (<>) {
                last if /^# Sieving side/;
            }
            next;
        }
        print STDERR "Broke loop at $_";
        last;
    }
    ship_chart(@$_) for @ship;
}

while (<>) {
        /cado-source\.([0-9a-f]+)/ && do {
                print "\\newpage\n";
                print "\\begin{verbatim}\n";
                system "git show --stat $1 | cat";
                print "\\end{verbatim}\n";
                next;
            };
    while (defined($_) && /^# Creating new slicing on side (\d+) for (scale=.*)$/) {
        parse_new_slicing($1, $2);
    }
    if (/^# Sieving (side-\d q=\d+; rho=\d+).*/) {
        parse_a_special_q($1);
    }
}

# for my $k (keys %$slicings) { display_slicing($k, $slicings->{$k}); }

