#!/usr/bin/perl -w

# Usage: cadosieve.pl $polyfile $paramfile1 $paramfile2 ...
#        read polyfile and paramfiles and estimate time for sieve
#        the paramfiles must contain: different name!, machines, params sieve...
# note: - it's not necessary than the polyfile contain the parameters of sieve.
#         For each paramfiles, an other polyfile is created and factors
#         base is calculated.
#       - you can to choose the interval of sieve for one job. (by
#         default range=4000 for Imin)
#
# TODO: generate input

use Cwd qw(abs_path);
use File::Basename;
use lib abs_path(dirname(dirname($0)));
use cadofct;
use strict;
use warnings;

my $range = 4000; # for Imin ($range = $range / 2 for Imin+1 ...)
my $int_num = 2; # approximate interval numbers
my $old_range = $range;

# input
print "$0 @ARGV\n";
my $poly = shift @ARGV;
my @files;
while (@ARGV) {
	push @files, shift @ARGV;
}
@files = sort @files;

# recover
my @done = <*.sieve_done>;
my @recover = <*.sieve_jobs>;

# sieve
my %link_name_params;
sub plus_range {
  my $x = shift;
  return ( $x + $range );
}
my %max_rels = ( # lpbr => [Imin, approximate total interval for Imin,
                 #          nrels max with ~20% duplicates]
                 25 => [ 11,  4000000,  5000000],
                 26 => [ 11,  4000000,  7500000],
                 27 => [ 11,  8000000, 13000000], 
                 28 => [ 12, 80000000, 26000000],
                 29 => [ 12,160000000, 52000000],
                 30 => [ 12,160000000,105000000],
                 31 => [ 13,400000000,210000000], 
                 32 => [ 13,800000000,420000000]);
my $file;
while (@files) {
  $file = shift @files;
  read_param(\%param, { strict => 1 }, "$file");
  $range = $old_range;
  $link_name_params{"$param{'name'}"} = $file;
  next if (@done);
  if (@recover) {
    next if ($recover[0] ne basename("$param{'prefix'}.sieve_jobs"));
    shift @recover;
  }
  if (basename($poly) ne basename("$param{'prefix'}.poly")) {
    open FILE, "< $poly"
      or die "Cannot open `$poly' for reading: $!.\n";
    open FILE2, "> $param{'prefix'}.poly"
      or die "Cannot open `$param{'prefix'}.poly' for writing: $!.\n";
    while (<FILE>) {
      print FILE2 "$_" 
        unless /^(rlim|alim|lpbr|lpba|mfbr|mfba|rlambda|alambda)/;
    }
    close FILE;
    print FILE2 "$_: $param{$_}\n"
      for qw(rlim alim lpbr lpba mfbr mfba rlambda alambda);
    close FILE2;
  }
  read_machines();
  do_factbase() unless (-f "$param{'prefix'}.roots");
	
  #init @t
  my @t = @{$max_rels{"$param{'lpbr'}"}};
  my $i = $param{'I'};
  while ( $i > $t[0] ) {
    $i--;
    $t[1] = int ( $t[1] / 2 );
    $range = int ( $range / 2 );
  }
  $t[1] = int ( $t[1] / $int_num );
		
  my @r;
  $r[0] = plus_range($param{'qmin'});
  $param{'qrange'} = $range;
  $r[1] = $param{'qmin'} + $t[1];
  while ( $r[1] < 30 * $t[2] ) {
    open FILE, "> $param{'prefix'}.rels.$r[0]-$r[1]"
      or die "Cannot open `$param{'prefix'}.rels.$r[0]-$r[1]' for writing: $!.\n";
    close FILE;
    $r[0] += $t[1];
    $r[1] += $t[1];
  }
  if (@files) {
    do_sieve_bench(\@t);
  } else {
    do_sieve_bench(\@t, 1);
    open FILE, "> $param{'prefix'}.sieve_done"
      or die "Cannot open `$param{'prefix'}.sieve_done' for writing: $!.\n";
   close FILE;
  }
}

# info rels
banner "Info rels";
my @rels_files;

foreach my $name (sort keys %link_name_params) {
  my $nrels = 0;
  my $reports;
  my $time_by_rels;
  my $estimate_time = 0;
  
  opendir DIR, "."
    or die "Cannot open directory `.': $!\n";
  @rels_files = grep /^$name\.rels\.[\de.]+-[\de.]+\.gz$/,
                     readdir DIR;
  closedir DIR;
  my @int_files = map { /^$name\.rels\.([\de.]+)-/; $1 } @rels_files;
  @int_files = sort ( {$a <=> $b } @int_files );

  $file = $link_name_params{"$name"};
  read_param(\%param, { strict => 1 }, "$file");
  $range = $old_range;

  #init @t
  my @t = @{$max_rels{"$param{'lpbr'}"}};
  my $i = $param{'I'};
  while ( $i > $t[0] ) {
    $i--;
    $t[1] = int ( $t[1] / 2 );
    $range = int ( $range / 2 );
  }
  $t[1] = int ( $t[1] / $int_num );

  my @rels_files = map { /^(\d+)$/; "$name.rels.$1-".plus_range($1).".gz" } @int_files;

  info "Params: rlim $param{'rlim'} alim $param{'alim'} lpbr $param{'lpbr'} ".
       "lpba $param{'lpba'} mfbr $param{'mfbr'} mfba $param{'mfba'}\n".
       "        rlambda $param{'rlambda'} alambda $param{'alambda'} ".
       "I $param{'I'} qmin $param{'qmin'} sieve_max_threads $param{'sieve_max_threads'}\n";
  foreach my $f (@rels_files) {
    if ( last_line($f) =~ /# Total (\d+) reports \[(\S+)s/ ) {
      info "$f: ".last_line($f)."\n";
      $reports = $1;
      $time_by_rels = $2; 
      my $rels = $t[1] * $reports / $range;
      $nrels += $rels;
      if ($nrels < $t[2]) {
        $estimate_time += $rels * $time_by_rels;
      } else {
        $estimate_time += ( $t[2] - ($nrels - $rels) ) * $time_by_rels;
        last;
      }
    }
  }
	
  $tab_level++;
  info "Estimate time for configuration $name : \033[01;31m$estimate_time"."s, ".
       format_dhms($estimate_time)."\033[01;00m\n";
  print "\n";
  $tab_level--;
}
	
