#!/usr/bin/perl
use strict;
use warnings;
use CommonFunctions qw(parseFileList);
use Getopt::Std;

my (%opt, @files, $bamFile, $species);
my $sam = 'tmp.sam';
my $bam = 'tmp.bam';

getopts('s:f:o:h',\%opt);
var_check();

my %species;
$species{'Athaliana'} = '/shares/jcarrington_share/gbrowse/Arabidopsis_thaliana/gff_files/A_THALIANA.gff3';
$species{'TuMV'} = '/shares/jcarrington_share/gbrowse/Turnip_mosaic_virus/gff_files/TUMV.gff3';

if (!exists($species{$species})) {
  print STDERR " $species is not a valid species identifier. See makeSAMfeatureDB.pl -h for options.\n\n";
  exit 1;
}

open (SAM, ">$sam") or die " Cannot open $sam: $!\n\n";
print SAM '@HD'."\t".'VN:1.0'."\t".'SO:unsorted'."\n";

open (GFF, $species{$species}) or die " Cannot open $species{$species}: $!\n\n";
while (my $row = <GFF>) {
  chomp $row;
  my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $row;
  if ($type eq 'chromosome') {
    print SAM '@SQ'."\t".'SN:'.$chr."\t".'LN:'.$end."\n";
  }
}
close GFF;

foreach my $file (@files) {
  open (IN, $file) or die " Cannot open $file: $!\n\n";
  while (my $row = <IN>) {
    next if (substr($row,0,1) eq '#' || $row =~ /^\s*$/);
    chomp $row;
    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $row;
    my %tags;
    my @attributes = split /;/, $attributes;
    foreach my $attribute (@attributes) {
      my ($tag, $value) = split /=/, $attribute;
      $tags{$tag} = $value;
    }
    if ($strand eq '+') {
      $strand = 0;
    } elsif ($strand eq '-') {
      $strand = 16;
    } elsif ($strand eq '.') {
      $strand = 0;
    } else {
      $strand = 0;
    }
    my $length = $end - $start + 1;
    my $id;
    if ($type eq 'gene' || $type eq 'mRNA' || $type eq 'pseudogene' || $type eq 'pseudogenic_transcript' || $type eq 'transposable_element' || $type eq 'stem' || $type eq 'TAS' || $type eq 'transposable_element_gene' || $type eq 'transposable_element_transcript') {
      if (exists($tags{'ID'})) {
        $id = $tags{'ID'}.':'.$type;
      } else {
        gff_error($file, $row, 9);
      }
    } elsif ($type eq 'five_prime_UTR' || $type eq 'CDS' || $type eq 'exon' || $type eq 'three_prime_UTR' || $type eq 'pseudogenic_exon' || $type eq 'transposable_element_exon') {
      if (exists($tags{'Parent'})) {
        $id = $tags{'Parent'}.':'.$type;
      } else {
        gff_error($file, $row, 9);
      }
    } elsif ($type eq 'miRNA' || $type eq 'star') {
      if (exists($tags{'ID'})) {
        $id = $tags{'ID'}.':'.$type;
        $id =~ s/\*//;
      } else {
        gff_error($file, $row, 9);
      }
    } elsif ($type eq 'peptide') {
      if (exists($tags{'Name'})) {
        $id = $tags{'Name'};
      } else {
        gff_error($file, $row, 9);
      }
    } else {
      gff_error($file, $row, 3);
    }
    print SAM $id."\t".$strand."\t".$chr."\t".$start."\t255\t".$length."M\t*\t0\t0\t*\t*\tXA:i:0\n";
  }
  close IN;
}

close SAM;

`samtools view -bS $sam -o $bam`;
`samtools sort -m 53687091200 $bam $bamFile`;
`samtools index $bamFile.bam`;

exit;

sub gff_error {
  my $file = shift;
  my $row = shift;
  my $column = shift;
  print STDERR " Malformed GFF3 file: $file\n";
  print STDERR " At this line, column $column => $row\n\n";
  exit 1;
}

sub var_check {
  if ($opt{'h'}) {
    var_error();
  }
  if ($opt{'f'}) {
    @files = parseFileList($opt{'f'});
  } else {
    var_error();
  }
  if ($opt{'o'}) {
    $bamFile = $opt{'o'};
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
  print STDERR " This script will generate a feature annotation database indexed in SAM/BAM format.\n";
  print STDERR " Usage: makeSAMfeatureDB.pl -f <GFF file(s)> -o <Output BAM file prefix>\n\n";
  print STDERR "      -f     List of single or multiple GFF3 files.\n\n";
  print STDERR "      -o     Output BAM file prefix\n\n";
  print STDERR "      -s     Species. Available options (Athaliana, TuMV)\n\n";
  print STDERR "      -h     Print this menu.\n\n";
  exit 1;
}
