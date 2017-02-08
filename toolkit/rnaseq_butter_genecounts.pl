#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my (%opt, $bamFile, $outfile, @genes, %genes);
getopts('b:o:h', \%opt);
arg_check();

open HEAD, "samtools view -H $bamFile |";
while (my $line = <HEAD>) {
  chomp $line;
  my @cols = split /\t/, $line;
  if ($cols[0] eq '@SQ') {
    my $gene = $cols[1];
    $gene =~ s/SN://;
    $gene =~ s/\.\d+//;
    next if (exists($genes{$gene}));
    push @genes, $gene;
    $genes{$gene} = 0;
  }
}
close HEAD;

open SAM, "samtools view -F 4 $bamFile |";
while (my $line = <SAM>) {
  my @fields = split /\t/, $line;
  if ($fields[1] == 16) {
    my $gene = $fields[2];
    $gene =~ s/\.\d+//;
    $genes{$gene}++;
  }
}
close SAM;

open(OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
print OUT "Gene\t$bamFile\n";
foreach my $gene (@genes) {
  print OUT $gene."\t".$genes{$gene}."\n";
}
close OUT;
exit;

sub arg_check {
  if ($opt{'h'}) {
    arg_error();
  }
  if ($opt{'b'}) {
    $bamFile = $opt{'b'};
  } else {
    arg_error('No BAM file was provided!');
  }
  if ($opt{'o'}) {
    $outfile = $opt{'o'};
  } else {
    arg_error('No output file was provided!');
  }
}

sub arg_error {
  my $error = shift;
  if ($error) {
    print STDERR $error."\n";
  }
  my $usage = "
usage: rnaseq_butter_genecounts.pl -b BAMFILE -o OUTFILE [-h]

Extracts gene-specific read counts from a BUTTER output BAM file aligned against the transcriptome.

arguments:
  -d BAMFILE            BAM file from BUTTER.
  -p OUTFILE            Output counts file.
  -h                    Show this help message and exit.

  ";
  print STDERR $usage;
  exit 1;
}