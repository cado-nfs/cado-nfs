#!/usr/bin/perl -w
#
use strict;
use warnings;
use POSIX qw/getcwd/;
use File::Temp;
use File::Copy;

# example code:
# machine=grisou ; for c in `git rev-list b074914913bff3732222ca8b13b5c80afa8e43da^..HEAD | tac` ; do c=`git rev-parse --short $c`; ./output_to_tikz_timechart.pl `ls -rt  rsa220/cado.$c.rsa220.gantt*${machine}*out | tail -1`  > chart-out.tex && pdflatex chart && \mv chart.pdf chart-$machine-$c.pdf ; done ; pdfjoin $(for c in `git rev-list b074914913bff3732222ca8b13b5c80afa8e43da^..HEAD | tac` ; do c=`git rev-parse --short $c`; echo chart-$machine-$c.pdf ; done) -o charts-$machine-by-thread.pdf

my $slicings={};

my $by_bucket = 0;
my $cropfraction=300;
my $merge_all_fib_levels=0;
my $color_by_step=1;

sub display_slicing {
    my ($fh, $k,$S)=@_;
    my ($ns, $nsp, $ni, $w) = @{$S};
    print $fh "Slicing for $k: $ns slices:\n";
    print $fh "\tpart 0: $ni->[0] ideals, weight $w->[0]\n";
    print $fh "\tpart 1: $ni->[1] ideals, weight $w->[1], $nsp->[1] slices\n";
    print $fh "\tpart 2: $ni->[2] ideals, weight $w->[2], $nsp->[2] slices\n";
}
sub parse_new_slicing {
    my ($in, $fh, $side, $K) = @_;
    my $ideals = [(0)x3];
    my $weight = [(0)x3];
    my $nslices = [(0)x3];
    my $nslices_total = 0;
    my $cpart;
    while (defined($_=<$in>)) {
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
        # display_slicing($fh, $K, $S);
        return;
    }
}

sub ship_chart {
    my ($fh,$sq,$nthreads,$tab) = @_;
    print STDERR "Shipping chart with ", scalar @$tab, " entries\n";
    my ($tmin, $tmax);
    return unless @$tab;
    my $NS = scalar @$tab;
    my $dt=0;
    my $idxmax={};
    my %kinds=();
    @$tab = sort { $a->[3] <=> $b->[3] } @$tab;
    my @chrono_kinds;
    for my $S (@$tab) {
        my ($kind, $thr, $idx, $t0, $t1, $time) = @$S;
        $tmin = $t0 unless defined($tmin) && $t0 > $tmin;
        $tmax = $t1 unless defined($tmax) && $t1 < $tmax;
        (my $mkind = $kind) =~ s/([A-Z]*).*$/$1/g;
        $idxmax->{$mkind} = $idx unless defined($idxmax->{$mkind}) && $idxmax->{$mkind} > $idx;
        $dt += $t1 - $t0;
        $kind =~ s/PBR\d+/PBR/g;
        push @chrono_kinds, $kind unless $kinds{$kind}++;
    }
    my $nc_chrono = scalar @chrono_kinds;
    $kinds{$chrono_kinds[$_]}=$_ for (0..$nc_chrono-1);

    my $nidx = {};
    for my $k (keys %$idxmax) { $nidx->{$k}=1+$idxmax->{$k}; }
    $dt = $dt / $NS;
    # We want to plot for a fraction of the time that is no less than
    # $cropfraction times the average per-slice time.
    my $tmax_cap = $tmin + $cropfraction * $dt;
    if ($tmax_cap > $tmax) { $tmax_cap = $tmax; }
    print $fh "\\newpage\n";
    print $fh "\n\n\\noindent\\textbf{$sq; $nthreads threads}\n\n";
    # we want to print everything, but based on moving windows of width
    # $tmax_cap-$tmin.
    
    if ($color_by_step) {
        print $fh "\\resetcolorseries[$nc_chrono]{rainbow}\n";
        for my $i (0..$nc_chrono-1) {
            my $style="fill={rainbow!![$i]}";
            print $fh "\\par\\noindent\\tikz\\draw[$style] (0,0) rectangle ++ (2em,2ex);\\quad {\\small\\sffamily $chrono_kinds[$i]}\n";
        }
    }

    my $ngraphs = 1 + int(($tmax-$tmin)/($tmax_cap-$tmin));
    my @lists=map { [] } (0..$ngraphs-1);
    for my $S (@$tab) {
        my ($kind, $thr, $idx, $t0, $t1, $time) = @$S;
        my $i0 = int(($t0-$tmin)/($tmax_cap-$tmin));
        my $i1 = int(($t1-$tmin)/($tmax_cap-$tmin));
        push @{$lists[$i0]}, $S;
        push @{$lists[$i1]}, $S if $i1 != $i0;
    }
    my $last_step=-1;
    print $fh "\\begin{center}\n";
    if ($ngraphs > 200) { $ngraphs=200; }
    my $ckind='';
    for(my $j = 0 ; $j < $ngraphs ; $j++) {
        print $fh "\\begin{adjustbox}{max totalsize={\\textwidth}{.9\\textheight},center}\n";
        print $fh "\\begin{tikzpicture}\n";
        my $tscale = 20;
        my $Y = -$nthreads/2;
        print $fh "\\draw[thin] (0,0) rectangle ($tscale,$Y);\n";
        my $YY = $Y-1;
        print $fh "\\path[use as bounding box] (0,0) rectangle ($tscale,$YY);\n";
        $YY--;
        print $fh "\\clip (0,1) rectangle ($tscale,$YY);\n";
        my $ccol=0;
        for my $S (@{$lists[$j]}) {
            my ($kind, $thr, $idx, $t0, $t1, $time) = @$S;
            (my $mkind = $kind) =~ s/([A-Z]*).*$/$1/g;
            my $ncol=$nidx->{$mkind};
            if (!$color_by_step && $ncol != $ccol) {
                print $fh "\\resetcolorseries[$ncol]{rainbow}\n";
                $ccol=$ncol;
            }
            my $x0 = ($t0 - $tmin) / ($tmax_cap - $tmin) - $j;
            my $x1 = ($t1 - $tmin) / ($tmax_cap - $tmin) - $j;
            $x0 = $tscale * $x0;
            $x1 = $tscale * $x1;
            my $style;
            if ($color_by_step) {
                (my $akind = $kind) =~ s/PBR\d+/PBR/g;
                my $i = $kinds{$akind};
                if ($i > $last_step) {
                    $last_step = $i;
                    print $fh "\\draw[thick] ($x0,0.25) -- ($x0,$Y) node[right,rotate=-90,scale=.7] {\\tiny $kind};\n";
                }
                $style="fill={rainbow!![$i]}";
            } else {
                if ($kind ne $ckind) {
                    $ckind=$kind;
                    print $fh "\\draw[thick] ($x0,0.25) -- ($x0,$Y) node[right,rotate=-90,scale=.7] {\\tiny $kind};\n";
                }
                $style="fill={rainbow!![$idx]}";
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
            print $fh "\\draw[$style] ($x0,$y0) rectangle ($x1,$y1); % $kind\n";
        }
        print $fh "\\end{tikzpicture}\n\n";
        print $fh "\\end{adjustbox}\n";
        my $cumul = (1+$j)*($tmax_cap-$tmin);
        if ($tmax - $tmin < $cumul) { $cumul = $tmax - $tmin; }
        printf $fh "Total time %.2f (\$t_{\\text{max}}-t_{\\text{min}}\$). Time above: %.2f (cumulated %.2f)\n\\bigskip\n\n",
        $tmax-$tmin, $tmax_cap-$tmin, $cumul;
    }
    print $fh "\\end{center}\n\n";
}

sub parse_a_special_q {
    my ($in,$fh,$sq) = @_;
    $sq=~s/;//;
    my $stack;
    while (defined($_=<$in>)) {
        /time chart for (\d+) threads/ && do {
            my $nthreads = $1;
            my $tab = [];
            { # for(my $i = 0 ; $i < $nthreads ; $i++) {
                defined($_=<$in>) or die;
                /time chart has (\d+) entries/ or die;
                my $nentries = $1;
                for(my $i = 0 ; $i < $nentries ; $i++) {
                    defined($_=<$in>) or die;
                    my ($kind, $thr, $idx, $t0, $t1, $time);
                    /thread (\d+)/ && do { $thr=$1; };
                    if (/FIB side (\d+) level (\d+) B (\d+) slice (\d+)/) {
                        $kind = "FIB$2s.$1";
                        $idx = $3;
                        # slice unused.
                    } elsif (/DS side (\d+) level (\d+) B (\d+)/) {
                        $kind = "DS$2l.$1";
                        $idx = $3;
                        # slice unused.
                    } elsif (/PBR M (\d+) B (\d+)/) {
                        $kind = "PBR$1";
                        $idx = $thr;
                        # slice unused.
                    } elsif (/PCLAT side (\d+) level (\d+)/) {
                        $kind = "PCLAT$2.$1";
                        $idx = $thr;
                        # slice unused.
                    } elsif (/ECM/) { $kind = "ECM"; $idx = $thr;
                    } elsif (/INIT/) { $kind = "INIT"; $idx = $thr;
                    } elsif (/BOTCHED/) { $kind = "BOTCHED"; $idx = $thr;
                    } elsif (/SKEWGAUSS/) { $kind = "SKEWGAUSS"; $idx = $thr;
                    } elsif (/SLICING/) { $kind = "SLICING"; $idx = $thr;
                    } elsif (/ADJUST/) { $kind = "ADJUST"; $idx = $thr;
                    } else {
                        die "Unexpected chart line: $_";
                    }
                    /t0 ([\d\.]+)/ && do { $t0=$1; };
                    /t1 ([\d\.]+)/ && do { $t1=$1; };
                    /time ([\d\.]+)/ && do { $time=$1; };
                    push @$tab, [$kind, $thr, $idx, $t0, $t1, $time];
                }
            }
            # get something to eat.
            ship_chart($fh, $sq, $nthreads, $tab);
            return;
        };
    }
}

my @files=();
my $unlink=1;

while (defined($_=shift @ARGV)) {
    if (/--by-thread/) {
        $by_bucket = 0;
        next;
    }
    if (/--by-bucket/) {
        $by_bucket = 1;
        next;
    }
    if (/--color-by-step/) {
        $color_by_step=1;
        next;
    }
    if (/--cropfraction/) {
        $cropfraction = shift @ARGV;
        next;
    }
    if (/--keep/) {
        $unlink=0;
        next;
    }
    chomp($_);
    /^(.*?)\.(txt|out)$/ or die "input files must be named *.out or *.txt [got $_]";
    my $base = $1;
    die "$_: $!" unless -f $_;
    push @files, [$_, "$base.pdf"];
}

scalar @files or do { print "nothing to do\n"; exit 0; };

my $dir = File::Temp->newdir(UNLINK=>$unlink, CLEANUP=>$unlink);

print STDERR "Created temporary directory $dir\n";
print STDERR "Will keep temporary files\n" unless $unlink;

open(my $fh, ">", "$dir/chart.tex") or die $!;
print $fh <<'EOF';
\documentclass[12pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage{a4wide}
\usepackage[utf8]{inputenc}
\usepackage{amsmath}
% \usepackage[landscape]{geometry}
\usepackage[rgb]{xcolor}
\usepackage{tikz}
\usepackage{adjustbox}
\usetikzlibrary{calc,positioning}
\pagestyle{empty}
\definecolorseries{rainbow}{hsb}{grad}[hsb]{.575,1,1}{.987,-.234,0}
\begin{document}
\input{chart-out.tex}
\end{document}
EOF
close $fh;

for my $pair (@files) {
    my ($file, $dst) = @$pair;
    print "\n\nCreating $dst from $file\n\n";
    open($fh, ">", "$dir/chart-out.tex");
    $file=~/cado(?:-source)?\.([0-9a-f]+)/ && do {
            my $commit=$1;
            print $fh "\\newpage\n";
            print $fh "\\begin{verbatim}\n";
            open(my $oldout, ">&STDOUT");
            open(STDOUT, ">&", $fh);
            system "git show --stat $commit | cat";
            open(STDOUT, ">&", $oldout);
            print $fh "\\end{verbatim}\n";
        };
    open my $in, "<$file";
#    while (defined($_=<$in>)) {
#        while (defined($_) && /^# Creating new slicing on side (\d+) for (scale=.*)$/) {
#            parse_new_slicing($in, $fh, $1, $2);
#        }
#        if (/^# Sieving (side-\d q=\d+; rho=\d+).*/) {
#            parse_a_special_q($in, $fh, $1);
#        }
#    }
    parse_a_special_q($in, $fh, "merge special-q's");
    close $in;
    close $fh;
    my $opwd = getcwd;
    chdir $dir;
    system "pdflatex chart";
    chdir $opwd;
    print "$dir/chart.pdf -> $dst\n";
    move("$dir/chart.pdf", $dst) or die $!;
}

# for my $k (keys %$slicings) { display_slicing($k, $slicings->{$k}); }

