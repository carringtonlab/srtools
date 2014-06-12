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

my (%opt, %stats, $file, $verbose, $confFile, $species, $libID, $outfile, %checked, %new);
my $defaultConf = '/home/nfahlgren/hts_data/sRNAmp.conf';
my $readsTable = 'reads.csv';
my $sequencesTable = 'sequences.csv';
my $newBAM = 0;
getopts('f:r:o:l:c:s:vh', \%opt);
var_check();

# Initialize stats hash, stores summary statistics for the run
# hits = total genome hits
# seqs = total unique sequences
# add_seq = sequences added to sequences table
# add_hit = hits to add
$stats{'hits'} = 0;
$stats{'seqs'} = 0;
$stats{'add_seq'} = 0;
$stats{'add_hit'} = 0;

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

# Get max sid for incrementing purposes
my $sth = $dbh->prepare('SELECT MAX(`sid`) FROM `sequences`');
$sth->execute();
my $maxSid = $sth->fetchrow_array;
if (!$maxSid) {
  $maxSid = 0;
}

############## End variables ##############

############## Begin main program ##############

print STDERR " Processed hits... 0\r" if ($verbose);

open (RD, ">$readsTable") or die " Cannot open file $readsTable: $!\n\n";
open (SEQ, ">$sequencesTable") or die " Cannot open file $sequencesTable: $!\n\n";

$sth = $dbh->prepare('SELECT * FROM `sequences` WHERE `seq` = ?');
open (OUT, ">$outfile") or die " Cannot open file $outfile: $!\n\n";
open (IN, $file) or die " Cannot open file $file: $!\n\n";
while (my $line = <IN>) {
  # We don't need the SQ lines of the SAM file
  if (substr($line,0,1) eq '@') {
		print OUT $line;
		next;
	}
  # These are the alignment lines
  chomp $line;
  my $sid;
  my ($name, $flag, $accession, $start, $mapq, $cigar, $mrnm, $mStart, $insertLen, $seq, $qual, $opt) = split /\t/, $line;
	next if ($accession eq '*');
  my ($id, $reads) = split /:/, $name;
  $stats{'hits'}++;
  print STDERR " Processed hits... $stats{'hits'}\r" if ($verbose && $stats{'hits'} % 100000 == 0);
  if (exists($checked{$id})) {
    # We have already dealt with this sequence, assume that the hits already exist
    $sid = $checked{$id};
  } else {
    # We haven't done anything with this sequence yet, add sequences and hits
    $stats{'seqs'}++;
    my $querySeq = $seq;
    # If hit is on reverse strand, reverse complement
    if ($flag == 16) {
      $querySeq = reverse($querySeq);
      $querySeq =~ tr/ATGC/TACG/;
    }
    # Have we seen this sequence in previous runs?
    $sth->execute($querySeq);
    my $match = $sth->fetchrow_hashref;
    if ($match) {
      # Yes
      $sid = $match->{'sid'};
    } else {
      # No
      $maxSid++;
      $sid = $maxSid;
      # Add to sequences table
      print SEQ "$sid,$querySeq\n";
      $stats{'add_seq'}++;
      $new{$id} = 1;
    }
    # Add to reads table
    print RD "$sid,$libID,$reads\n";
    $checked{$id} = $sid;
  }
  if (exists($new{$id})) {
    # Add hits to BAM database
    print OUT join("\t", ($sid,$flag,$accession,$start,$mapq,$cigar,$mrnm,$mStart,$insertLen,$seq,$qual,$opt))."\n";
    $stats{'add_hit'}++;
  }
}
close IN;
close OUT;
close RD;
close SEQ;
print STDERR " Processed hits... $stats{'hits'}\n" if ($verbose);

# Print out the run statistics
while (my ($k, $v) = each(%stats)) {
  print $k." => ".$v."\n";
}

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
  if ($opt{'o'}) {
    $outfile = $opt{'o'};
  } else {
    var_error();
  }
  if ($opt{'l'}) {
    $libID = $opt{'l'};
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
}

# var_error prints out command-line options, defaults, optional settings
sub var_error {
  print STDERR "\n\n";
  print STDERR " This script will populate the SQLite database with library and sequence information and reformat the SAM alignment file.\n";
  print STDERR " Usage: PopulateDB.pl -f <sequence file> -l <library ID> -o <output SAM file> -s <species>\n\n";
  print STDERR "   -f     Input sequence file\n\n";
  print STDERR "   -o     Output SAM file.\n\n";
  print STDERR "   -l     Library ID.\n\n";
  print STDERR "   -s     Species code. Example: A_THALIANA\n\n";
  print STDERR "   -c     Configuration file. Default = $defaultConf\n\n";
  print STDERR "   -v     Verbose output. Prints program status to terminal. Default is quiet.\n\n";
  print STDERR "   -h     Print this menu\n\n";
  exit 1;
}
