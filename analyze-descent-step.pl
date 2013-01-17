#!/usr/bin/perl -w

use strict;
use warnings;
# use Time::HiRes qw/time/;
use Data::Dumper;
use POSIX qw/uname/;

(my $node = (uname)[1]) =~ s/\..*$//;
my $bins="./build/$node";
my $config='p6w3';

my $lpbr_base = 35;
my $lpba_base = 35;

sub blen {
    my $x = shift;
    my $l = length($x);
    my $head = {
        0=>0, 1=>1, 2=>2, 3=>2,
        4=>3, 5=>3, 6=>3, 7=>3,
        8=>4, 9=>4, a=>4, b=>4,
        c=>4, d=>4, e=>4, f=>4,
    };
    return 4*($l-1)+$head->{substr($x,0,1)};
}

my $table = {
};


sub dump_table {
    for my $k (sort { $a cmp $b } keys %$table) {
        print "$k " . join(" ", @{$table->{$k}}) . "\n";
    }
}

sub finish {
    if ($@) {
        print STDERR "ERROR: $@\n";
    }
    print "\nTable at program exit:\n\n";
    dump_table;
    exit 1;
}

# $SIG{__DIE__}=\&finish;
$SIG{'INT'}=\&finish;

sub feed_to_table {
    my $feed = shift;

    for my $line (split /^/m, $feed) {
        my ($prefix, @x) = split ' ', $line;
        $table->{$prefix} = \@x;
    }
}

sub tryit {
    my $qbits = shift;
    $qbits =~ /^(\d+)([ra])$/ or die;
    $qbits = $1;
    my $rat = $2;
    my $prefix = "$qbits$rat";

    my @wo_lpb=();
    my @pass2las=();
    my ($lpbr_max,$lpba_max);
    my @fbparams=();
    {
        my @temp=@_;
        while (defined(my $v = shift @temp)) {
            if ($v =~ /^(\d+),(\d+),(\d+),([\d\.]+)$/) {
                push @fbparams, { lim=>$1, lpb=>$2, mfb=>$3, lambda=>$4 };
            } else {
                push @wo_lpb, $v;
                push @pass2las, $v;
            }
        }
    }

    die "fbparams has wrong size" unless scalar @fbparams == 2;

    push @pass2las, "-ratq" if $rat eq 'r';

    {
        my @sides=(qw/r a/);
        for my $h (@fbparams) {
            my $side = shift @sides;
            my $mfb = $h->{'mfb'};
            my $prod = int($h->{'lpb'}*$h->{'lambda'});
            if ($mfb == 0) {
                print "$prefix Setting mfb$side=$prod\n";
                $h->{'mfb'} = $prod;
            } else {
                my $c = $mfb / $prod;
                if ($c < 0.9 || $c > 1.1) {
                    printf STDERR "$prefix Warning, mfb$side (%d) should be closer to %d\n",
                        $h->{'mfb'}, int($h->{'lpb'}*$h->{'lambda'});
                }
            }
        }
        $lpbr_max = $fbparams[0]->{'lpb'};
        $lpba_max = $fbparams[1]->{'lpb'};
        push @pass2las, "lpbr=$fbparams[0]->{'lpb'}";
        push @pass2las, "lpba=$fbparams[1]->{'lpb'}";
        push @pass2las, "mfbr=$fbparams[0]->{'mfb'}";
        push @pass2las, "mfba=$fbparams[1]->{'mfb'}";
        push @pass2las, "rlim=$fbparams[0]->{'lim'}";
        push @pass2las, "alim=$fbparams[1]->{'lim'}";
        push @pass2las, "rlambda=$fbparams[0]->{'lambda'}";
        push @pass2las, "alambda=$fbparams[1]->{'lambda'}";

    }

    # Create a random $qbits bits number, together with that number +
    # 2^20
    my ($q0, $q1);
    {
        my $n = $qbits-32;
        my $lo = 2**($n-2);
        my $hi = 2**($n-1);
        my $head = sprintf("%x", $hi + $lo + int(rand($lo)));
        $q0 = "0x" . $head . "00000000";
        $q1 = "0x" . $head . "00010000";
    }

    my $splittings;
    my $n=0;
    my $nok=0;

    my $args_common="-poly $config.poly -fb $config.roots $config.params -v -v -v -timed";
    my $cmd = "$bins/sieve/las $args_common -q0 $q0 -q1 $q1";
    $cmd .= " $_" for @pass2las;

    my @results=();

    print "$prefix: $cmd\n";
    open F, "$cmd |";
    my ($mnorm_alg, $mnorm_rat);
    my $sbr;
    my $nrel;
    ALLQ: while (<F>) {
        /^# Alg.*log2\(maxnorm\)=([\d\.]+)/ && do { $mnorm_alg=$1; next; };
        /^# Rat.*log2\(maxnorm\)=([\d\.]+)/ && do {
            $mnorm_rat=$1;
            printf "$prefix Norms: %d %d\n", int($mnorm_rat), int($mnorm_alg);
            next;
        };
        /^# Sieving q=(\d+); rho=(\d+)/ && do {
            $splittings=[];
            $sbr = 0;
            $nrel = 0;
            last ALLQ if $n >= 40;
            $n++;
            next;
        };
        if (/^#\s(\d+)\sin sbr table/) { $sbr += $1; }
        if (/^#\s(\d+)\srelation/) { $nrel += $1; }
        next unless defined $splittings;
        if (/^# Time for this special-q: ([\d\.]+)s \[norm [\d\.]+\+[\d\.]+, sieving ([\d\.]+) \([\d\.]+ \+ [\d\.]+\), factor ([\d\.]+)\]$/) {
            my $tfac = $3;
            my $tsieve = $2;
            my $ratio = $2/$1;
            # Time to flush the data for the current prime
            # What is the best splitting obtained for this special q ?
            my @s = sort { $a->[0] <=> $b->[0] } @$splittings;
            if (defined(my $best = shift @s)) {
                $nok++;
                push @results, $best;
                my ($x, $self, $p, $f) = @$best;
                my $sub = $x - $self;
                my $tstring = sprintf '%.3f', $x;
                if ($f) {
                    $tstring = sprintf "%.3f=%.3f+%.3f [$f]", $x, $self, $sub;
                }
                printf "$prefix #$nok/$n (%.2f%%) [$nrel, sieve=$tsieve, %.1f%%]: $tstring\n",
                    100.0*$nok/$n, 100.0*$ratio;
            }
            next;
        }
        next if /^#/;
        my $rel = $_;
        /^\(([\d\.e\+-]+)\)\s(-?\d+,-?\d+):([a-f\d,]+):([a-f\d,]+)/ or die "$_";
        my $x0 = $1;
        my $ab = $2;
        my $sub = {};
        my @fr = grep { $_ > $lpbr_base; } map { blen($_); } split(',',$3);
        my @fa = grep { $_ > $lpba_base; } map { blen($_); } split(',',$4);
        $sub->{$_ . 'r'}++ for @fr;
        $sub->{$_ . 'a'}++ for @fa;
        $sub->{$prefix}--;
        my $y = 0;
        my $pr = 1;
        for (keys %$sub) {
            next unless $sub->{$_};
            my $l = $table->{$_} or die "Undefined cost for $_ ??? ; offending rel:\n$rel\n";
            my ($cost, $prob) = @$l;
            $y += $cost * $sub->{$_};
            $pr *= $prob;
        }
        my $factype = join(' ', grep { $sub->{$_} > 0 } keys %$sub);

        if ($pr < 0.6) {
            print STDERR "Discarding split $factype, pr=$pr\n";
            next;
        }
        push @$splittings, [$x0 + $y, $x0, $pr, $factype];
    }

    # Now we need to decide which are the proper lpb bounds for getting most
    # of the results.
    my $lpbr = $lpbr_base;
    my $lpba = $lpba_base;
    my $accepted = {};
    my $accept_bound=999.99;
    my @restricted_results;
    my $restricted_par;
    while ($lpbr <= $lpbr_max || $lpba <= $lpba_max) {
        my $par="lpbr=$lpbr lpba=$lpba";
        $fbparams[0]->{'lpb'} = $lpbr;
        $fbparams[1]->{'lpb'} = $lpba;
        # How many of the results satisfy these bounds ?
        my @nres = ();
        my $t1=0;
        my $p=0;
        SCAN_RESULTS: for my $r (@results) {
            for my $next (split(' ', $r->[3])) {
                next SCAN_RESULTS unless exists $accepted->{$next};
            }
            $t1 += $r->[0];
            $p += $r->[2];
            push @nres, $r;
        }
        if (!scalar @nres) {
            print "$prefix No results above $par\n";
            last;
        }
        my $m = $t1 / scalar @nres;
        if ($m > $accept_bound) {
            printf "$prefix Not taking results with $par\n";
            last;
        }
        $p = $p/$n;
        if ($p > 0.5) {
            if (scalar @restricted_results < scalar @nres) {
                # There are two ways to count some sort of a weighted
                # time.
                # 1) for a time t which is reached with some probability
                # p, consider that the event where this occurs can be
                # reached with probability p, whence an expected time
                # could be t/p. This is not really proper, beacuase it is
                # not possible of course to randomize all descent steps
                # independently of eachother. Therefore, such a measure
                # is overoptimistic.
                # 2) another reasoning might be to consider that if we
                # fail, then we have to resort to other, possibly more
                # tolerant, factorizations. Therefore, with probability
                # 1-p, we are likely to spend *at least* 2t, maybe more.
                # Say k*t for some k. k=2 is perhaps optimistic, we don't
                # know. Based on this, the weighted time can be something
                # like t*(k-p).
                # my $weighted = $m / $p;
                my $k = 2.5; # Here we fix k=2.5
                my $weighted = $p*$m + $k*(1-$p)*$m;
                printf "$prefix Can do t=%.3f prob=%.2f%%, weighted=%.3f ($par)\n", $m, 100.0*$p, $weighted;
                $accept_bound = $m * 1.5;
                @restricted_results = @nres;
                $restricted_par = $par;
            }
        }
        if ($lpba == $lpbr) {
            $lpba++;
            $accepted->{$lpba . 'a'}=1
        } else {
            $lpbr++;
            $accepted->{$lpbr . 'r'}=1
        }
    }
    if (scalar @restricted_results) {
        @results = @restricted_results;
    } else {
        print "$prefix Cannot reach a success probability above 50%\n";
    }
    my $t1=0;
    my $t2=0;
    for (@results) {
        my ($x, $self, $p, $f) = @$_;
        $t1 += $x;
        $t2 += $x * $x;
    }
    close F;
    my $line;
    if ($nok) {
        my $m = $t1 / $nok;
        my $s = ($t2 / $nok - $m * $m);
        my $prob = 1.0 * $nok/$n;
        my $cost = $m;
        # $extra .= " $restricted_par";
        $line = sprintf "$prefix %.3f %.2f", $cost, $prob;
        $line .= " $_" for @wo_lpb;
        $line .= " " . join(',', map { $fbparams[0]->{$_} } qw/lim lpb mfb lambda/);
        $line .= " " . join(',', map { $fbparams[1]->{$_} } qw/lim lpb mfb lambda/);
    } else {
        $line = "$prefix 999.999 0.00 FAILED";
    }
    feed_to_table($line);
    print "\n$line\n\n\n";
}

# This is on magret.loria.fr (Core i5-2500)
sub readfile {
    my $f = shift;
    local $/=undef;
    open my $x, "<$f";
    my $r=<$x>;
    close $x;
    return $r;
};
my $feed_magret = readfile("analyze-descent-step-data.magret.txt");
my $feed_truffe = readfile("analyze-descent-step-data.truffe.txt");

# feed_to_table($feed_truffe);

# The lpbr/lpba bounds indicated below are *indicative*. Better bounds, based
# on the results obtained, are returned by tryit().

# tryit(split(' ', "42a I=10 400000,41,0,1.8 400000,42,0,2.0"));

do { my $i = 36; my $j = $i - 1;
tryit(split(' ', "${i}a I=10  400000,$j,0,2.0  400000,$j,0,2.0"));
tryit(split(' ', "${i}r I=11 2000000,$j,0,1.8 1000000,$i,0,1.8"));
};

for my $i (37..40) {
my $j = $i - 1;
tryit(split(' ', "${i}a I=11  400000,$j,0,2.0  400000,$j,0,2.0"));
tryit(split(' ', "${i}r I=11 2000000,$j,0,1.8 1000000,$i,0,1.8"));
}
for my $i (41..46) {
my $j = $i - 1;
tryit(split(' ', "${i}a I=10  400000,$j,0,1.8  400000,$j,0,2.0"));
tryit(split(' ', "${i}r I=11 2000000,$j,0,1.8 1000000,$i,0,1.8"));
}
for my $i (47..51) {
my $j = $i - 1;
tryit(split(' ', "${i}a I=10  400000,$j,0,1.4  400000,$j,0,2.0"));
tryit(split(' ', "${i}r I=10 1800000,$j,0,1.8 1000000,$i,0,1.8"));
}
for my $i (52..53) {
my $j = $i - 1;
tryit(split(' ', "${i}a I=10  400000,$j,0,1.4  400000,$j,0,2.0"));
tryit(split(' ', "${i}r I=10 1800000,$j,0,1.8 1000000,$i,0,1.8"));
}
for my $i (54..64) {
my $j = $i - 1;
tryit(split(' ', "${i}a I=10  600000,$j,0,1.4  600000,$j,0,1.4"));
tryit(split(' ', "${i}r I=10 1800000,$j,0,1.8 1000000,$i,0,1.8"));
}
# do { my $i = 58; my $j = $i - 1;
# tryit(split(' ', "${i}r I=8 1800000,$j,0,1.8 1000000,$i,0,1.8"));
# exit 0;};
dump_table;
