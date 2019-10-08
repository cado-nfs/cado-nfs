#!/usr/bin/env perl

use strict;
use warnings;
use File::Temp qw/ tempdir /;
use IO::Handle;

$|=1;

my $verbose = 1;

# both flags below default to $verbose.
# 0 means nothing is printed.
# 1 means that failed (or pending) outputs are printed.
#   Passed outputs are not printed.
# 2 means that all outputs are printed.
my $print_outputs_on_normal_exit;
my $print_outputs_on_interrupt;

my $use_color = 1;

my $tmpdir = tempdir( CLEANUP => 1 );

my %archive=();
my @live=();

my $tests_pending = 0;
my $tests_passed = 0;
my $tests_failed = 0;

my $csi = {
    alert => "\033[01;31m",     # bold red
    failure => "\033[01;31m",   # bold red
    success => "\033[01;32m",   # bold green
    starting => "\033[01;33m",  # bold yellow
    normal => "\033[00;30m",    # black
};

while (defined($_ = shift @ARGV)) {
    if (/^-nc$/) { $csi->{$_}='' for keys %$csi; }
    elsif (/^-v$/) { $verbose++; }
    elsif (/^-q$/) { $verbose--; }
    else {
        die "unexpected arg\n";
    }
}

if (!defined($print_outputs_on_interrupt)) {
    $print_outputs_on_interrupt = $verbose;
}
if (!defined($print_outputs_on_normal_exit)) {
    $print_outputs_on_normal_exit = $verbose;
}


sub sformat {
    my $n = shift @_;
    if (!defined($n)) {
        # return sprintf("%3s  %12s", "", "");
        return "   ";
    }
    my $title = $archive{$n}->{'title'};
    # $title = substr($title, -12);
    return sprintf("%3d: %s", $n, $title);
}

sub short_global_status {
    my @tcounts = map {
                    $_->[1] >= 1 ?
                    $csi->{$_->[0]} .
                    sprintf("%3d", $_->[1]) .
                    $csi->{'normal'}
                    :
                    sprintf("%3d", $_->[1])
                    } (
                        ['starting', $tests_pending],
                        ['success', $tests_passed],
                        ['failure', $tests_failed]
                    );
    return join(" ", @tcounts);
}

sub print_live {
    my $n = shift @_;
    my $color = shift @_;
    my $ex='';
    my @txs=();
    for (@live) {
        my $tx = defined($_) ? sprintf("%3d", $_) : "   ";
        if (defined($_) && $_ == $n && defined($color)) {
            $ex = "$csi->{$color}$tx: $archive{$n}->{'title'}$csi->{'normal'}";
            $tx = "$csi->{$color}$tx$csi->{'normal'}";
        }
        push @txs, $tx;
    }
    my $tests = "[" . short_global_status . "]";
    $tests .= " " . join(" ", @txs);
    $tests .= "   $ex" if $ex;
    print "$tests\n" unless ($verbose < 0 && defined($color) && $color eq 'starting');
}

sub put_to_live_slot {
    my $n = shift @_;
    my $color = shift @_;
    my $i;
    for($i = 0 ; $i < scalar @live ; $i++) {
        if (!defined($live[$i])) {
            $live[$i]=$n;
            last;
        }
    }
    if ($i == scalar @live) {
        push @live, $n;
    }
    print_live($n, $color);
}

sub remove_from_live_slot {
    my $n = shift @_;
    my $color = shift @_;
    my $i;
    for($i = 0 ; $i < scalar @live ; $i++) {
        if (defined($live[$i]) && $live[$i] == $n) {
            last;
        }
    }
    die "unexpected test number $i, not currently running:\n$_" if $i == scalar @live;
    print_live($n, $color);
    $live[$i]=undef;
}

sub print_all_outputs {
    my $interrupt = scalar @_;
    return if $interrupt && $print_outputs_on_interrupt == 0;
    return if !$interrupt && $print_outputs_on_normal_exit == 0;
    for my $n (sort { $a <=> $b } keys %archive) {
        my $outcome=$archive{$n}->{'outcome'};
        my $passed = defined($outcome) && $outcome =~ /Passed/;
        if ($passed) {
            next if $interrupt && $print_outputs_on_interrupt <= 1;
            next if !$interrupt && $print_outputs_on_normal_exit <= 1;
        }
        my $pending = ($interrupt && !defined($outcome));
        next if ($verbose >= 3 && !$pending);
        open F, "$tmpdir/$n.out";
        while (defined($_=<F>)) {
            print;
        }
        close F;
        if ($pending) {
            print "Note: test $n was interrupted\n";
        }
        print "\n";
    }
}

sub sigint {
    my $sig = shift @_;
    print "\n$csi->{'alert'}Testing was interrupted$csi->{'normal'}\n";
    print_all_outputs($sig);
    print "\n$csi->{'starting'}Note: testing was interrupted$csi->{'normal'}\n";
    print "[" . short_global_status . "]\n";
    my $all_pending = '';
    my $all_failed = '';
    for my $n (keys %archive) {
        my $outcome=$archive{$n}->{'outcome'};
        my $title=$archive{$n}->{'title'};
        my $passed = defined($outcome) && $outcome =~ /Passed/;
        next if $passed;
        if ($outcome) {
            $all_failed .= sprintf("%3d", $n) . ": $title\n";
        } else {
            $all_pending .= sprintf("%3d", $n) . ": $title\n";
        }
    }
    if ($all_pending ne '') {
        print "$csi->{'starting'}=== $tests_pending interrupted tests ===\n$all_pending$csi->{'normal'}";
    }
    if ($all_failed ne '') {
        print "$csi->{'failure'}=== $tests_failed failed tests ===\n$all_failed$csi->{'normal'}";
    }
    exit 1;
}

$SIG{'INT'} = \&sigint;

while (<>) {
    next if /^test \d+$/;
    next if /^$/;
    if (/^\s+Start\s+(\d+): (.*)$/) {
        my $n = $1;
        my $title = $2;
        die "unexpected test number $n, already running:\n$_" if $archive{$n};
        open(my $fh, ">", "$tmpdir/$n.out");
        $archive{$n} = {
            title=>$title,
            data=>"",
            fd=>$fh
        };
        print $fh $_;
        print $fh "\n";
        $fh->flush();
        $tests_pending++;
        put_to_live_slot($n, 'starting');
        next;
    } elsif (/^(\d+):\s?/) {
        my $n=$1;
        die "unexpected data from stdin:\n$_" unless $archive{$1};
        my $fh = $archive{$n}->{'fd'};
        print $fh $_;
        $fh->flush();
    } elsif (/^\s*\d+\/\d+ Test\s+\#(\d+):\s+(\S+) \.*(.*)/) {
        my $n = $1;
        my $title = $2;
        my $outcome = $3;
        die "unexpected test number $n, not currently running:\n$_" unless $archive{$n};
        die "unexpected test title, does not match: \"$archive{$n}->{'title'}\"\n$_" unless $archive{$n}->{'title'} eq $title;
        $archive{$n}->{'outcome'}=$outcome;
        my $fh = $archive{$n}->{'fd'};
        print $fh "\n";
        print $fh $_;
        close  $fh;
        delete $archive{$n}->{'fd'};
        if ($verbose >= 3) {
            # print the full test contents.
            open F, "$tmpdir/$n.out";
            while (defined($_=<F>)) {
                print;
            }
            close F;
        }
        $tests_pending--;
        if ($outcome =~ /Passed/) {
            $tests_passed++;
        } else {
            $tests_failed++;
        }
        remove_from_live_slot($n, ($outcome =~ /Passed/ ? 'success' : 'failure'));
    } elsif (/^\s*\d+% tests/ || /The following/) {
        die "unexpected final text, some tests are still running: @live\n$_" if scalar grep { defined($_) } @live;
        my $tailmsg = $_;
        print_all_outputs;
        if ($verbose >= 0) {
            print "\n";
            print $tailmsg;
            while (<>) {
                print;
            }
        } else {
            # consume output, so that we avoid the dreaded broken pipe
            while (<>) {}
        }
        last;
    } elsif (scalar @live == 0) {
        next;
    } else {
        die "unexpected data from stdin:\n$_";
    }
}
