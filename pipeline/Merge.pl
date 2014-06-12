#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Config::Tiny;
use Cwd;
use DBI;
use Env qw(HOME);
use lib "$HOME/lib/perl";
use CommonFunctions qw(parseListToArray);

############## Begin variables ##############

my (%opt, $file, $verbose, $confFile, $species, $samFile);
my $defaultConf = '/home/nfahlgren/hts_data/sRNAmp.conf';
my $readsTable = 'reads.csv';
my $sequencesTable = 'sequences.csv';
getopts('f:r:o:l:c:s:vh', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};
my $cwd = getcwd();

############## End variables ##############

############## Begin main program ##############

# Load reads and sequences tables into database
print STDERR " Loading data into sequences... " if ($verbose);
`sqlite3 -separator ',' $conf->{'db'} '.import $cwd/$sequencesTable sequences'`;
print STDERR "done\n" if ($verbose);
print STDERR " Loading data into reads... " if ($verbose);
`sqlite3 -separator ',' $conf->{'db'} '.import $cwd/$readsTable reads'`;
print STDERR "done\n" if ($verbose);

# If there is already a BAM database, convert it to SAM
# Then concatenate the input SAM file onto the existing database
if (-e $conf->{'bam'}) {
  $samFile = 'tmp.sam';
  if (-e $samFile) {
    unlink($samFile);
  }
  `$sam view -h $conf->{'bam'} > $samFile`;
  `cat $file | grep -v '\^\@' >> $samFile`;
} else {
  $samFile = $file;
}

# Create BAM alignment database
if (-e $conf->{'bam'}) {
  unlink($conf->{'bam'});  
}
if (-e 'tmp.bam') {
  unlink('tmp.bam');
}
#`$sam view -bt $conf->{'fai'} $cwd/$samFile > tmp.bam`;
`$sam view -bS $samFile -o tmp.bam`;

# Sort the temporary BAM database
my $prefix = $conf->{'bam'};
$prefix =~ s/\.bam$//;
`$sam sort -m 53687091200 tmp.bam $prefix`;

# Index the BAM database
unlink('tmp.bam');
unlink('tmp.sam');
if (-e "$conf->{'bam'}.bai") {
  unlink("$conf->{'bam'}.bai");
}
`$sam index $conf->{'bam'}`;

## Create BAM alignment database
## Convert the input SAM file to a temporary BAM file
#`$sam view -bt $conf->{'fai'} $file > tmp.bam`;

#if (! -e $conf->{'bam'}) {
#  # Since the BAM database does not exist, sort the temporary
#  # BAM database and save it as the official database
#  my $outBAM = $conf->{'bam'};
#  $outBAM =~ s/\.bam$//;
#  `$sam sort tmp.bam $outBAM`;
#} else {
#  # Sort the temporary BAM database
#  `$sam sort tmp.bam tmp.sorted`;
#  # Merge the sorted BAM files
#  `$sam merge -f $conf->{'bam'} $conf->{'bam'} tmp.sorted.bam`;
#  unlink('tmp.sorted.bam');
#}
# Remove the temporary BAM file
#unlink('tmp.bam');
## Index the new BAM database
#`$sam index $conf->{'bam'}`;

exit;

############## End main program ##############

############## Begin subroutines ##############

# var_check parses command-line options
# Activates var_error if required options missing, sets defaults
sub var_check {
  if ($opt{'h'}) {
    var_error();
  }
  if ($opt{'f'}) {
    $file = $opt{'f'};
  } else {
    var_error();
  }
  if ($opt{'v'}) {
    $verbose = 1;
  } else {
    # Verbose output is off by default
    $verbose = 0;
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
}

# var_error prints out command-line options, defaults, optional settings
sub var_error {
  print STDERR "\n\n";
  print STDERR " This script will load new sequences and reads into the SQLite database and create a merged BAM hits database.\n";
  print STDERR " Usage: Merge.pl -f <SAM file> -s <species>\n\n";
  print STDERR "   -f     Input sequence file\n\n";
  print STDERR "   -s     Species code. Example: A_THALIANA\n\n";
  print STDERR "   -c     Configuration file. Default = $defaultConf\n\n";
  print STDERR "   -v     Verbose output. Prints program status to terminal. Default is quiet.\n\n";
  print STDERR "   -h     Print this menu\n\n";
  exit 1;
}
