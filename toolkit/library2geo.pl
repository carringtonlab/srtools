#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Config::Tiny;
use DBI;

my (%opt, $library, $gff, $conf, $species, %reads);
getopts('g:l:c:s:h', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($conf);
my $sam = $Conf->{'PIPELINE'}->{'sam'};
my $config = $Conf->{$species};
my $refdb = $config->{'refseq'};
my $source = $config->{'source'};
my $type = $config->{'type'};
my $mol = $type;
$mol =~ s/_/ /g;

if (!$refdb || !$source || !$type) {
  print STDERR " No reference database, source or type information found in $conf.\n";
  print STDERR " Please add refseq, source and type fields before running this program.\n\n";
  exit 1;
}

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$config->{'db'}","","");

open (GFF, ">$gff") or die " Cannot open $gff: $!\n\n";
print GFF "##gff-version 3\n\n";
print GFF "##Reference sequence database=$refdb\n";
print GFF "##COLUMN1=<Name=seqid,Type=String,Description=\"An identifier for the reference sequence landmark.\">\n";
print GFF "##COLUMN2=<Name=source,Type=String,Description=\"Source sample/tissue.\">\n";
print GFF "##COLUMN3=<Name=type,Type=String,Description=\"Sample type/fraction.\">\n";
print GFF "##COLUMN4=<Name=start,Type=Integer,Description=\"The start of the feature, in 1-based integer coordinates, relative to the landmark given in COLUMN1. Start is always less than end.\">\n";
print GFF "##COLUMN5=<Name=start,Type=Integer,Description=\"The end of the feature, in 1-based integer coordinates, relative to the landmark given in COLUMN1.\">\n";
print GFF "##COLUMN6=<Name=score,Type=Integer,Description=\"Read frequency, the number of times the feature was observed in the given library.\">\n";
print GFF "##COLUMN7=<Name=strand,Type=String,Description=\"Strand, relative to the landmark given in COLUMN1. Positive '+' or minus '-'.\">\n";
print GFF "##COLUMN8=<Name=phase,Type=String,Description=\"No phase given, indicated with a '.'.\">\n";
print GFF "##COLUMN9=<Name=attributes,Type=String,Description=\"List of tag-value pairs separated by ';'.\">\n";
print GFF "##ATTRIBUTE: Name=<Type=Integer,Description=\"An integer assigned to each unique sequence within the given sample.\">\n";
print GFF "##ATTRIBUTE: Seq=<Type=String,Description=\"The nucleotide sequence (5'->3') of the $mol.\">\n";

my $sth = $dbh->prepare("SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `library_id` = $library");
$sth->execute();
while (my $row = $sth->fetchrow_hashref) {
  $reads{$row->{'sid'}}->{'reads'} = $row->{'reads'};
  $reads{$row->{'sid'}}->{'seq'} = $row->{'seq'};
}

open SAM, "$sam view $config->{'bam'} |";
while (my $hit = <SAM>) {
  my @fields = split /\t/, $hit;
	my $sid = $fields[0];
	my $strand = $fields[1];
	my $chrom = $fields[2];
	my $start = $fields[3];
	my $seq = $fields[9];
	next if (!exists($reads{$sid}));
	# Convert strand bitwise operator to + or - strand
	$strand = ($strand == 0) ? '+' : '-';
	# Use length of sequence to determine end position
	my $end = $start + length($seq) - 1;
	my $score = $reads{$sid}->{'reads'};
  my $phase = '.';
	my $attributes = 'Name='.$sid.';Seq='.$reads{$sid}->{'seq'};
  
  print GFF $chrom."\t";
  print GFF $source."\t";
  print GFF $type."\t";
  print GFF $start."\t";
  print GFF $end."\t";
  print GFF $score."\t";
  print GFF $strand."\t";
  print GFF $phase."\t";
  print GFF $attributes."\n";
}
close SAM;
close GFF;
exit;

sub var_check {
  if ($opt{'h'}) {
    var_error();
  }
  if ($opt{'l'}) {
    $library = $opt{'l'};
  } else {
    var_error();
  }
  if ($opt{'g'}) {
    $gff = $opt{'g'};
  } else {
    var_error();
  }
  if ($opt{'c'}) {
    $conf = $opt{'c'};
  } else {
    var_error();
  }
  if ($opt{'s'}) {
    $species = $opt{'s'};
  } else {
    var_error();
  }
  
}

sub var_error {
  print STDERR " This script will output a library to GFF3 for GEO submissions.\n\n";
  print STDERR " Usage: library2geo.pl -l <Library ID> -g <GFF output file> -c <conf file> -s <species>\n\n";
  print STDERR "     -l     Library ID.\n\n";
  print STDERR "     -g     GFF3 output file.\n\n";
  print STDERR "     -c     Configuration file.\n\n";
  print STDERR "     -s     Species.\n\n";
  print STDERR "     -h     Print this menu.\n\n";
  exit 1;
}
