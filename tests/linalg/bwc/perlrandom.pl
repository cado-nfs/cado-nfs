#!/usr/bin/env perl

use Digest::MD5 qw(md5);
# Generates $1 bytes of random crap, based on a seed.
#
# if 0 is passed as a length, produce an infinite amount od crap.
# 
# the optional second argument serves as a random seed.

sub usage {
    die "Usage: $0 <nbytes> [<integer seed>]\n";
}

defined(my $nbytes = shift @ARGV) or usage;
my $seed = shift @ARGV || 0;
my $infinite = $nbytes == 0;
$|=1;

$ctx = Digest::MD5->new;
$ctx->add($seed);

while ($infinite || $nbytes > 0) {
    my $crap = $ctx->clone->digest;
    my $take = ($infinite || $nbytes >= 16) ? 16 : $nbytes;
    syswrite(STDOUT, $crap, $take) or die;
    $nbytes -= $take;
    $ctx->add($crap);
}
