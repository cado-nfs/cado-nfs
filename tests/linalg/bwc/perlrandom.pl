#!/usr/bin/env perl

use Digest::MD5 qw(md5);
# Generates $1 bytes of random crap, based on a seed.

sub usage {
    die "Usage: $0 <nbytes> [<integer seed>]\n";
}

my $nbytes = shift @ARGV or usage;
my $seed = shift @ARGV || 0;

$ctx = Digest::MD5->new;
$ctx->add($seed);

while ($nbytes > 0) {
    my $crap = $ctx->clone->digest;
    my $take = ($nbytes >= 16) ? 16 : $nbytes;
    syswrite(STDOUT, $crap, $take);
    $nbytes -= $take;
    $ctx->add($crap);
}
