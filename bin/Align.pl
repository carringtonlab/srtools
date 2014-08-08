#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Config::Tiny;
use Cwd;
use FindBin qw($Bin);

############## Begin variables ##############

my (%opt, $seqFile, $confFile, $species, $verbose, $method);
our (%stats, $conf, $outfile);
my $defaultConf = "$Bin/../include/default.conf";
getopts('f:o:c:s:a:vh', \%opt);
var_check();

#Initialize stats hash, stores summary statistics for the run
# hits = genome hits
# hit_reads = reads that have hits
# nohits = sequences without hits
# nohit_reads = reads without hits
$stats{'hits'} = 0;
$stats{'hit_reads'} = 0;
$stats{'nohits'} = 0;
$stats{'nohit_reads'} = 0;

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
$conf = $Conf->{$species};

############## End variables ##############

############## Begin main program ##############

# Parse method input to determine the subroutine to use to align the input file
if ($method eq 'bowtie') {
  bowtie($seqFile);
} else {
  print STDERR " Alignment method type $method is not currently supported.\n\n";
  var_error();
}

# Print out the run statistics
while (my ($k, $v) = each(%stats)) {
  print $k." => ".$v."\n";
}

exit;

############## End main program ##############

############## Begin subroutines ##############

sub bowtie {
  my $file = shift;
  my $db = $conf->{'bowtieDB'};
  if (!$db) {
    print STDERR " Bowtie database $db is not defined in $confFile.\n\n";
    var_error();
  }
  open (OUT, ">$outfile") or die " Cannot open file $outfile: $!\n\n";
  # Run bowtie, parse output
  #open BOW, "$Conf->{'PIPELINE'}->{'bowtie'} -f -n 0 -a -S --best $db $file |";
  # Options
  # -f    Input file is FASTA
  # -v 0  Report ungapped alignments with zero mismatches
  # -a    Report all alignments, limited by -v 0
  # -S    Output alignments in SAM format
  # -p 4  Use four CPUs
  #open BOW, "$Conf->{'PIPELINE'}->{'bowtie'} -f -v 0 -a -S -p 4 $db $file |";
  open BOW, "$Conf->{'PIPELINE'}->{'bowtie'} -f -v 0 -a -S -p 1 $db $file |";
  while (my $row = <BOW>) {
    if (substr($row,0,1) eq '@') {
      print OUT $row;
    } else {
      my @fields = split /\t/, $row;
      my ($id, $reads) = split /:/, $fields[0];
      if ($fields[2] eq '*') {
        $stats{'nohits'}++;
        $stats{'nohit_reads'} += $reads;
      } elsif ($fields[5] !~ /^\d+M$/) {
        print $row;
        exit;
      } else {
        $stats{'hits'}++;
        $stats{'hit_reads'} += $reads;
        print OUT $row;
      }
    }
  }
  close BOW;
  close OUT;
}

# var_check parses command-line options
# Activates var_error if required options missing, sets defaults
sub var_check {
  if ($opt{'h'}) {
    var_error();
  }
  if ($opt{'f'}) {
    $seqFile = $opt{'f'};
  } else {
    var_error();
  }
  if ($opt{'v'}) {
    $verbose = 1;
  } else {
    # Verbose output is off by default
    $verbose = 0;
  }
  if ($opt{'o'}) {
    $outfile = $opt{'o'};
  } else {
    var_error();
  }
  if ($opt{'c'}) {
    $confFile = $opt{'c'};
  } else {
    # Default conf file
    $confFile = $defaultConf;
  }
  if ($opt{'s'}) {
    $species = $opt{'s'};
  } else {
    var_error();
  }
  if ($opt{'a'}) {
    $method = $opt{'a'};
  } else {
    $method = 'bowtie';
  }
}

# var_error prints out command-line options, defaults, optional settings
sub var_error {
  print STDERR "\n\n";
  print STDERR " This script will align reads to a reference genome and process the hits.\n";
  print STDERR " Usage: Align.pl -f <sequence file> -o\n\n";
  print STDERR "   -f     Input sequence file\n\n";
  print STDERR "   -o     Output file name.\n\n";
  print STDERR "   -s     Species code. Example: A_THALIANA\n\n";
  print STDERR "   -c     Configuration file. Default = $defaultConf\n\n";
  print STDERR "   -a     Alignment method. Default = bowtie.\n";
  print STDERR "               Current methods: bowtie\n\n";
  print STDERR "   -v     Verbose output. Prints program status to terminal. Default is quiet.\n\n";
  print STDERR "   -h     Print this menu\n\n";
  exit 1;
}
