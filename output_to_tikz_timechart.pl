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
    my ($sq,$nthreads,$tab) = @_;
    my ($tmin, $tmax);
    return unless @$tab;
    my $NS = scalar @$tab;
    my $dt=0;
    my $idxmax={};
    @$tab = sort { $a->[3] <=> $b->[3] } @$tab;
    for my $S (@$tab) {
        my ($kind, $thr, $idx, $t0, $t1, $time) = @$S;
        $tmin = $t0 unless defined($tmin) && $t0 > $tmin;
        $tmax = $t1 unless defined($tmax) && $t1 < $tmax;
        (my $mkind = $kind) =~ s/([A-Z]*).*$/$1/g;
        $idxmax->{$mkind} = $idx unless defined($idxmax->{$mkind}) && $idxmax->{$mkind} > $idx;
        $dt += $t1 - $t0;
    }
    my $nidx = {};
    for my $k (keys %$idxmax) { $nidx->{$k}=1+$idxmax->{$k}; }
    $dt = $dt / $NS;
    # We want to plot for a fraction of the time that is no less than
    # $cropfraction times the average per-slice time.
    my $tmax_cap = $tmin + $cropfraction * $dt;
    if ($tmax_cap > $tmax) { $tmax_cap = $tmax; }
    print "\\newpage\n";
    print "\n\n\\noindent\\textbf{$sq; $nthreads threads}\n\n";
    # we want to print everything, but based on moving windows of width
    # $tmax_cap-$tmin.
    
    my $ngraphs = 1 + int(($tmax-$tmin)/($tmax_cap-$tmin));
    my @lists=map { [] } (0..$ngraphs-1);
    for my $S (@$tab) {
        my ($kind, $thr, $idx, $t0, $t1, $time) = @$S;
        my $i0 = int(($t0-$tmin)/($tmax_cap-$tmin));
        my $i1 = int(($t1-$tmin)/($tmax_cap-$tmin));
        push @{$lists[$i0]}, $S;
        push @{$lists[$i1]}, $S if $i1 != $i0;
    }

    print "\\begin{center}\n";
    if ($ngraphs > 16) { $ngraphs=16; }
    my $ckind='';
    for(my $j = 0 ; $j < $ngraphs ; $j++) {
        print "\\begin{adjustbox}{max totalsize={\\textwidth}{.9\\textheight},center}\n";
        print "\\begin{tikzpicture}\n";
        my $tscale = 20;
        my $Y = -$nthreads/2;
        print "\\draw[thin] (0,0) rectangle ($tscale,$Y);\n";
        my $YY = $Y-1;
        print "\\path[use as bounding box] (0,0) rectangle ($tscale,$YY);\n";
        $YY--;
        print "\\clip (0,1) rectangle ($tscale,$YY);\n";
        my $ccol=0;
        for my $S (@{$lists[$j]}) {
            my ($kind, $thr, $idx, $t0, $t1, $time) = @$S;
            (my $mkind = $kind) =~ s/([A-Z]*).*$/$1/g;
            my $ncol=$nidx->{$mkind};
            if ($ncol != $ccol) {
                print "\\resetcolorseries[$ncol]{rainbow}\n";
                $ccol=$ncol;
            }
            my $x0 = ($t0 - $tmin) / ($tmax_cap - $tmin) - $j;
            my $x1 = ($t1 - $tmin) / ($tmax_cap - $tmin) - $j;
            $x0 = $tscale * $x0;
            $x1 = $tscale * $x1;
            if ($kind ne $ckind) {
                $ckind=$kind;
                print "\\draw[thick] ($x0,0.25) -- ($x0,$Y) node[right,rotate=-90,scale=.7] {\\tiny $kind};\n";
            }
            # # my $pcpu = 0;
            # # if ($t1 > $t0) { $pcpu = 1.0e-6*$time / ($t1-$t0); }
            # # $pcpu = int(10*$pcpu+0.5)*0.1;
            $x0 = sprintf("%.3f", $x0);
            $x1 = sprintf("%.3f", $x1);
            my $y0 = -$thr;
            my $y1 = -$thr - 1;
            $y0 = $y0/2;
            $y1 = $y1/2;
            print "\\draw[fill={rainbow!![$idx]}] ($x0,$y0) rectangle ($x1,$y1); % $kind\n";
        }
        print "\\end{tikzpicture}\n\n";
        print "\\end{adjustbox}\n";
        my $cumul = (1+$j)*($tmax_cap-$tmin);
        if ($tmax - $tmin < $cumul) { $cumul = $tmax - $tmin; }
        printf "Total time %.2f (\$t_{\\text{max}}-t_{\\text{min}}\$). Time above: %.2f (cumulated %.2f)\n\\bigskip\n\n",
        $tmax-$tmin, $tmax_cap-$tmin, $cumul;
    }
    print "\\end{center}\n\n";
}

sub parse_a_special_q {
    my $sq = shift @_;
    $sq=~s/;//;
    my $stack;
    while (<>) {
        /time chart for (\d+) threads/ && do {
            my $nthreads = $1;
            my $tab = [];
            for(my $i = 0 ; $i < $nthreads ; $i++) {
                defined($_=<>) or die;
                /time chart has (\d+) entries/ or die;
                my $nentries = $1;
                for(my $i = 0 ; $i < $nentries ; $i++) {
                    defined($_=<>) or die;
                    my ($kind, $thr, $idx, $t0, $t1, $time);
                    /thread (\d+)/ && do { $thr=$1; };
                    /FIB side (\d+) level (\d+) B (\d+) slice (\d+)/ && do {
                        $kind = "FIB$2s.$1";
                        $idx = $3;
                        # slice unused.
                    };
                    /DS side (\d+) level (\d+) B (\d+) slice (\d+)/ && do {
                        $kind = "FIB$2l.$1";
                        $idx = 0;
                        # slice unused.
                    };
                    /PBR M (\d+) B (\d+)/ && do {
                        $kind = "PBR$1";
                        $idx = $thr;
                        # slice unused.
                    };
                    /PCLAT side (\d+) level (\d+)/ && do {
                        $kind = "PCLAT$2.$1";
                        $idx = $thr;
                        # slice unused.
                    };
                    /t0 ([\d\.]+)/ && do { $t0=$1; };
                    /t1 ([\d\.]+)/ && do { $t1=$1; };
                    /time ([\d\.]+)/ && do { $time=$1; };
                    push @$tab, [$kind, $thr, $idx, $t0, $t1, $time];
                }
            }
            # get something to eat.
            ship_chart($sq, $nthreads, $tab);
            return;
        };
    }
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

