#!/usr/bin/perl

# do not include Data::Dumper except for debugging.
# use Data::Dumper;

my $current = undef;
my $all = {};


open F, "<flint-fft/doc/fft.rst";

while (defined($_=<F>)) {
    chomp($_);
    # print STDERR "$_\n";
    if (/^\.\.\s+function::\s+((\S+)\s+(\w+)\s*\(.*\))$/) {
        $current = undef;
        $key = $3;
        $all->{$key}="";
        $current = \$all->{$key};
        die unless defined($_=<F>);
        die unless /^\s*$/;
    } elsif (/^\S/) {
        $current = undef;
        next;
    } else {
        $$current .= "$_\n";
    }
}
close F;

my $tags={};
open F, "<flint-fft/tags";
while (defined($_=<F>)) {
    chomp($_);
    next if /^!/;
    /^(\S+)\011(\S+)\011(.*);".*$/ or die;
    push @{$tags->{$1}}, [$2, $3];
}
close F;

# print Dumper($tags);
# print Dumper($all);

for my $func (keys %$all) {
    my $t = $tags->{$func} or die "No tag for $func";
    if (@$t != 1) {
        die "Multiple tags for $func";
    }
    my ($file, $position) = @{$t->[0]};

    # create a comment block.
    my $doc = $all->{$func};
    my @doclines = split(/^/m, $doc);
    my $prefix = 999;
    for (@doclines) {
        next if /^\s*$/;
        /^(\s*)/;
        my $l = length($1);
        if ($l < $prefix) { $prefix = $l; }
    }

    for (@doclines) {
        s/^\s{$prefix}//;
        s/^/ * /gm;
    }
    $doc = "/*\n" . join("",@doclines) . " */\n";

    print "Adding documentation to flint-fft/$file\n";
    open G, "| ex -s flint-fft/$file > /dev/null";
    print G "set nomagic\n";
    # print G "set noai\n";
    print G "$position\n";
    print G "i\n";
    print G $doc;
    print G ".\n";
    print G "wq\n";
}

