#!/usr/bin/perl -w
# Run a series of difsim jobs, for different parameter settings
# For now it just varies q, and sits in the cwd etc
# author: David Rideout
# version 0.1
use strict;

# Tunable parameters
my $exe = "~/projects/neuroviz/DifSim/difsim";
my $minq = .05;
my $maxq = 1.0;
my $nruns = 9;
my $delta = .025; # 25 microsec pulse width -- only used for display now

# Parse command line
die "usage: RunSeries.pl <parameter file>\n" unless @ARGV==1;

# Read parameter file
my $parfile;
while (<>) {
  next if /^\/\//;
  $parfile .= $_;
}
# Replace signal file in parameter file
die "signal file parse error\n" unless $parfile =~ s/signal file = [\w.]+.*/signal file = signaljunk.asc\n/;
# Grab tensor file from parameter file
die "failed extracting tensor file\n" unless $parfile =~ /tensor file = (.*)/;
my $tensorfile = $1;
print "using tensor file $tensorfile\n";

# Set up output file
# die "[$parfile]" unless $parfile =~ s/par/asc/;
open OF, ">eps_vs_q.asc" or die "error opening file for output: $!\n";
print OF "# q\tepsilon\tsigma\terr\n";

# Main DifSim loop
my $seed = 0;
for (my $q = $minq; $q <= $maxq+.001; $q += .05) { # not sure why I have to add epsilon to the upper bound

  # Compute run parameters
  printf "Running DifSim $nruns times for q = %.2f...", $q;
  my $G = $q/267.513e6/$delta*1e11;
  printf " G =%5.0f\n", $G;

  # Replace parameters in parfile
  die "gradient strength parse error at [$G]\n" unless $parfile =~ s/gradient strength = [\d.]+\n/gradient strength = $G\n/;

  my @vals;
  for (my $run=0; $run<$nruns; ++$run) {
    # Place new seed into parameter file
    ++$seed;
    if ($parfile =~ /random seed/) {
      die "random seed parse error at [$G]\n"
	  unless $parfile =~ s/random seed = \d+\n/random seed = $seed\n/;
    } else {
      $parfile .= "random seed = $seed\n";
    }

    # Write parfile
#   print $parfile;
    open FP, ">junk.par" or die "error writing parameter file";
    print FP $parfile;
    close FP;

    # Run DifSim
    #system "../DifSim/difsim junk.par";
    my $difsim_output = `$exe junk.par`;
    # need some error checking here!!
    #print STDERR "Is this the return value of difsim?? : ${^CHILD_ERROR_NATIVE}\n"; # http://beerpla.net/2008/04/29/how-do-i-get-both-the-return-value-and-text-in-perl-backticks-vs-system-perl-510/
    die "Error ${^CHILD_ERROR_NATIVE} in difsim\n" if ${^CHILD_ERROR_NATIVE};
    open DIFSIM_OUT, ">difsim_output";
    print DIFSIM_OUT $difsim_output;
    close DIFSIM_OUT;

    # Analyze output
#     my $analysis_output = `~/projects/neuroviz/dde_72_signal_analysis.py signaljunk.asc $G $delta 0`;
    my $analysis_output = `~/projects/neuroviz/DifSim/utilities/AnalyzeSignal.pl signaljunk.asc $tensorfile q=$q`;
#   print $analysis_output;
    die "error parsing analysis output [$analysis_output]\n" unless $analysis_output =~ /epsilon = ([\-e\d.]+)/;
    #my $epsilon = $1;
    push @vals, $1; #epsilon;

  } # nrun runs

  # Compute statistics
  print join(' ', @vals), "\n";
  my $mean = 0;
  foreach (@vals) {
    $mean += $_;
  }
  $mean /= $nruns;
  die "too few runs\n" unless @vals>1;
  my $sigma = 0;
  foreach (@vals) {
    $sigma += ($_-$mean)*($_-$mean);
  }
  $sigma = sqrt($sigma/($nruns-1));
  my $err = $sigma/sqrt($nruns);
  
  #my $epsilon = @vals;
#  print OF "$q\t$epsilon\n";
  print OF "$q\t$mean\t$sigma\t$err\n";
}

close OF;
