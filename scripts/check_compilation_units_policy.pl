#!/usr/bin/env perl

# This script makes sure that for all compilation units found in the
# source tree, we have '#include "cado.h"' near the top, as we should.
#
# Of course we have exceptions. See @path_exceptions below.
#       gf2x/
#       linalg/bwc/mpfq/
#       linalg/bwc/flint-fft/
# This script must be called from the top of the tree.
#
# 
# using vim, you may do 
#       :set makeprg=./scripts/check_compilation_units_policy.pl
# and cycle through the errors as if they were produced by gcc.

use strict;
use warnings;

die "Please call $0 from the top of the source tree" unless -f "cado.h";

my @all_files;

if (scalar @ARGV) {
    @all_files = @ARGV;
} else {
    @all_files = `git ls-files`;
}


my $ENABLE_CHECK_FOR_CADO_H_IN_COMPILATION_UNITS = 1;
my $ENABLE_CHECK_FOR_NO_CADO_H_IN_HEADERS = 1;
my $ENABLE_CHECK_FOR_PORTABILITY_H_WHEN_NEEDED = 1;
my $ENABLE_CHECK_FOR_HASH_IFNDEF_GUARDS_IN_HEADERS = 1;
my $ENABLE_CHECK_FOR_CXX_GUARDS_IN_C_HEADERS = 1;
my $ENABLE_CHECK_FOR_OWN_HEADER_IN_COMPILATION_UNITS = 1;

my @path_exceptions=qw|
        gf2x/
        linalg/bwc/mpfq/
        linalg/bwc/flint-fft/
        config/
        misc/
        utils/embedded/fmt/
        ci/coverity_model.c
        |;

my @header_guard_check_skip_patterns = (
    qr{^cado\.h$},
);

my @include_portability_skip_patterns = (
    qr{^cado\.h$},
    qr{^portability\.h$},
);

my $err=0;
FILE: for my $f (@all_files) {
    $f =~ s/^\s*//;
    $f =~ s/\s*$//;
    next unless $f =~ /\.[ch](?:pp)?$/;
    for my $p (@path_exceptions) {
        next FILE if $f =~ /^$p/;
    }
    my $is_header = ($f =~ /\.h(?:pp)?$/);
    my $is_cxx = ($f =~ /pp?$/);

    (my $self_include = $f) =~ s/c(pp|)$/h$1/;
    if (-f "$self_include") {
        $self_include =~ s,^.*/,,g;
    } else {
        undef $self_include;
    }

    open F, $f;

    my @includes;
    my @include_cado;
    my @include_self;
    my $guard_score=0;
    my $cxxguard_score=0;
    my $portability_score=0;
    my @needs_portability_h;

    while (defined($_ = <F>)) {
        my $lnum = $.;
        my $aa = [ $lnum, $_];
    
        chomp($_);
        my @x=();

        if ($is_header) {
            if ($guard_score == 0 && /^#ifndef/) { $guard_score++; }
            if ($guard_score == 1 && /^#define/) { $guard_score++; }
            if ($guard_score == 2 && /^#endif/) { $guard_score++; }
            if (/pragma multi include/) { $guard_score=3; }
            if ($cxxguard_score == 0 && /__cplusplus/) { $cxxguard_score=1; }
            if (/pragma no prototypes/) { $cxxguard_score=2; }
            if (/generated automatically/i) { $cxxguard_score=2; }
            if (/automatically generated/i) { $cxxguard_score=2; }
            if ($ENABLE_CHECK_FOR_NO_CADO_H_IN_HEADERS && /^\s*#\s*include\s*.*cado\.h/) {
                print STDERR "$f:$lnum: is a header file, it must not include cado.h\n";
                $err++;
            }
        } else {
            if (/^\s*#\s*include/) {
                push @includes, [ $lnum, $_ ];
                if (defined($self_include) && /\"$self_include\"/) {
                    push @include_self, $aa;
                }
                if (/.*cado\.h/) { push @include_cado, $aa; }
                if (/.*portability\.h/) { $portability_score++; }
            }
            for my $func (qw/strlcat strlcpy asprintf realpath pagesize sleep strdup strndup/) {
                last if @needs_portability_h;
                if (/\b$func\b/) { push @needs_portability_h, [@$aa, $func]; }
            }
        }
    }

    if ($is_header) {
        my $skip = 0;
        for (@header_guard_check_skip_patterns) {
            if ($f =~ $_) {
                $skip = 1;
                last;
            }
        }
        if ($ENABLE_CHECK_FOR_HASH_IFNDEF_GUARDS_IN_HEADERS && !$skip && $guard_score != 3) {
            my $guard_name = $f;
            $guard_name =~ s,^.*/,,g;
            $guard_name =~ tr/a-z/A-Z/;
            $guard_name =~ s,[^A-Z],_,g;
            $guard_name = "CADO_${guard_name}_";
            my $guard_construct = <<EOF;
#ifndef $guard_name
#define $guard_name

/* stuff */

#endif /* $guard_name */
EOF
            print STDERR "$f:1: is a header, it must be protected by a guard construct such as\n" . $guard_construct;
            $err++;
        }

        if ($ENABLE_CHECK_FOR_CXX_GUARDS_IN_C_HEADERS && !$skip && !$is_cxx && !$cxxguard_score) {
            my $guard_construct = <<EOF;
#ifdef __cplusplus
extern "C" {
#endif

/* prototypes, etc. typedefs and structs are ok. #include are not */

#ifdef __cplusplus
}
#endif
EOF
            print STDERR "$f:1: is a C header. If it exposes any prototypes, these must be protected by the following construct:\n" .  $guard_construct;
            $err++;
        }
    } else {
        my $first = shift @include_cado;
        if ($ENABLE_CHECK_FOR_CADO_H_IN_COMPILATION_UNITS) {
            if ($first) {
                if ($first->[1] =~ /cado\.h/) {
                    if ($first->[1] !~ /"cado\.h"/) {
                        print STDERR "$f:$first->[0]: cado.h should be included as \"cado.h\", not <cado.h>\n";
                        $err++;
                    }
                } else {
                    print STDERR "$f:$first->[0]: is a compilation unit, its first include file must be \"cado.h\"\n";
                    $err++;
                }
                if (@include_cado) {
                    my $lnum = $include_cado[0]->[0];
                    print STDERR "$f:$lnum: there is no point in including cado.h twice\n";
                    $err++;
                }
            } else {
                print STDERR "$f:1: is a compilation unit, it must include \"cado.h\"\n";
                $err++;
            }
        }

        if ($ENABLE_CHECK_FOR_OWN_HEADER_IN_COMPILATION_UNITS && defined($self_include)) {
            unless (@include_self) {
                print STDERR "$f:1: is a compilation unit, and $self_include exists. It seems that #include \"$self_include\" is missing\n";
                $err++;
            }
        }
    }
    if ($ENABLE_CHECK_FOR_PORTABILITY_H_WHEN_NEEDED && @needs_portability_h && $portability_score == 0) {
        my $first = shift @needs_portability_h;
        print STDERR "$f:$first->[0]: uses function $first->[2], which is not always present. You must #include \"portability.h\" // $first->[2]\n";
        $err++;
    }

    close F;
}
if ($err) {
    print STDERR "$err errors found. Please fix\n";
    exit 1;
}
