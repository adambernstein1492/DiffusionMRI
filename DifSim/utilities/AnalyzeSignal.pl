#!/usr/bin/perl -w
# Analyze output signal
# author: David Rideout
# version 0.1
use strict;

print "Generic analysis script for DifSim\n";

# Parse command line
# print scalar(@ARGV), " $ARGV[0]\n";
die "usage: AnalyzeSignal.pl <signal file> <tensor file> q=<q>\n" unless
    @ARGV==3 && $ARGV[2] =~ /q=([\d.e]+)/;
my $q = $1;
my $signalfile = shift;
my $tensorfile = shift;

# Read signal file
open SIG, "<$signalfile" or die "error opening signal file $signalfile: $!\n";
open TEN, "<$tensorfile" or die "error opening tensor file $tensorfile: $!\n";
print "Warning: ignoring all but first column!\n";
# my @sigs;
my ($b0, $same, $nsame, $diff, $ndiff);
while (<SIG>) {
  next if /^#/;
  die "unfamiliar line [$_]\n" unless /^([\d.]+) /;
  my $signal = $1;
#     push @sigs, $_;
  unless ($b0) {
    $b0 = $signal;
    next;
  }

  my $vecs = <TEN>;
  die "[sig = $_  vecs = $vecs]\n" unless $vecs =~ /^([\d.-]+ [\d.-]+ [\d.-]+) ([\d.-]+ [\d.-]+ [\d.-]+)/;
  if ($1 eq $2) {
#     print "same\n";
    $same += $signal;
    ++$nsame;
  } else {
#     print "diff\n";
    $diff += $signal;
    ++$ndiff;
  }
}

# Output epsilon
print "epsilon = ", (log($same/$nsame) - log($diff/$ndiff))/$q**4, "\n";
