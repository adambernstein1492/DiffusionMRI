#!/usr/bin/perl -w
# Generate gradient vectors for dpfg simulation, as cartesian product of
# entries from a tensor.dat file of directions.
# Also generates DTI sequences, along with corresponding QUEST files for
# analysis.
# Author: David Rideout <drideout@math.ucsd.edu>
# Date: March 2017
use strict;

# Parse options
use Getopt::Std;
my %opts;
# our ($opt_b, $opt_s);
getopts('bsm:d:n:', \%opts);
# print "b=$opts{b} s=$opt_s\n";
# print "[$opts{m}]\n";

# Parse command line
die "usage: generate-dpfg-tensor-file.pl [OPTIONS] "
    ."<DifSim/DifSimGui/tensor.dat> > <dpfg gradient file>\n"
    ."  -b generates output in 'Benni format'\n"
    ."  -d '<list of ndirs>'  number of directions for each shell\n"
    ."  -m '<list of shells>' list b values\n"
    ."  -n <name> (used in filenames)\n"
    ."  -s generates output for SPFG\n" unless @ARGV>0;
my $benni = $opts{b};
my $spfg = $opts{s};
my $name = $opts{n};
my $tensorfile = shift;
my @ndirs = split /\s/, $opts{d};
# my $ndirs = shift;
my @bvals = split /\s/, $opts{m};
my $nshells = @bvals;
my $readem;
my @directions;
my $maxb = 0;

# Check parameters
die "Number of directions (-d) must match number of shells (-m)\n" unless @ndirs==@bvals;
print STDERR "Sorry I probably broke the non-benni format!\n" unless $benni;
open DIRS, "<$tensorfile" or die "error opening $tensorfile: $!\n";

my $mag;
while (<DIRS>) {

  # Ignore comments
  next if /^\#/;

  # Is this the beginning of a new block?
  if (/^(\d+)$/) {
    $readem = 0;
    for (my $shell=0; $shell<$nshells; ++$shell) {
      if ($1==$ndirs[$shell]) {
	print STDERR "reading shell $shell\n";
	$readem = 1;
	$mag = sqrt($bvals[$shell]);
	$maxb = $mag if $mag>$maxb;
      }
    }
    next;
  }
  next unless $readem;

  # Strip trailing zeros
  s/[~0]0+(\s)/$1/g;
  s/\.(\s)/$1/g;

  # Grab directions
  die "[$_]" unless /([-\d.]+) ([-\d.]+) ([-\d.]+)/;
#   chomp;
#   s/ /,/g if $benni;
  push @directions, ($1*$mag, $2*$mag, $3*$mag);
}
print STDERR "Generating DPFG sequence from $nshells shells of ", @directions/3, " directions from $tensorfile\n";


# Output gradient file
if ($spfg) {
  # bval.txt
  open QVAL, ">$name-bvals.txt" or die "error opening QUEST bvals file: $!\n";
  print QVAL "0\n0\n"; # two leading zeros?
#   foreach (@bvals) {
#     print QUEST "$_\n"; # join('\n',@bvals);
#   }
#   close QUEST;

  # bvec.txt
  open QVEC, ">$name-bvecs.txt" or die "error opening QUEST bvecs file: $!\n";
  print QVEC "0 0 0\n0 0 0\n"; # two leading zeros?
  
  if ($benni) {
    printf("# Set of %03d directions\n", @directions/3);
    print "# additional header info\n";
    print "[directions=", @directions/3, "]\n";
    print "CoordinateSystem = xyz\n";
    print "Normalization = none\n";
  }
  my $bvali=0;
  my $bdiri=0;
  for (my $i=0; $i<@directions; $i+=3) {
    if (++$bdiri > $ndirs[$bvali]) {
      ++$bvali;
      $bdiri=1;
    }
    print QVAL "$bvals[$bvali]\n"; # $bdiri $bvali\n";
    print "Vector[", $i/3, "] = ( " if $benni;
    print $directions[$i]/$maxb, ",", $directions[$i+1]/$maxb, ",", $directions[$i+2]/$maxb;
    print QVEC "$directions[$i] $directions[$i+1] $directions[$i+2]\n";
    print " )" if $benni;
    print "\n";
  }

} else {
# if ($benni) {
#   printf("# Set of %03d directions\n", $ndirs*$ndirs*$nshells*$nshells);
#   print "# additional header info\n";
#   print "[directions=", $ndirs*$ndirs*$nshells*$nshells, "]\n";
#   print "CoordinateSystem = xyz\n";
#   print "Normalization = none\n";
# }
# for (my $i = 0; $i < $ndirs; ++$i) {
#   for (my $j = 0; $j < $ndirs; ++$j) {
#     for (my $shelli = 0; $shelli < $nshells; ++$shelli) {
#       my $scalei = ($shelli+1)/$nshells;
#       for (my $shellj = 0; $shellj < $nshells; ++$shellj) {
# 	my $scalej = ($shellj+1)/$nshells;
# 
# 	print "Vector[", $i*$ndirs*$nshells*$nshells + $j*$nshells*$nshells + $shelli*$nshells + $shellj, "] = ( " if $benni;
#     #print "$directions[$i], $directions[$j]";
# 	print $scalei*$directions[$i*3], ",", $scalei*$directions[$i*3+1], ",", $scalei*$directions[$i*3+2], ", ", $scalej*$directions[$j*3],",", $scalej*$directions[$j*3+1], ",", $scalej*$directions[$j*3+2];
# 	print " )" if $benni;
# 	print "\n";
#       }
#     }
#   }
# }
}
