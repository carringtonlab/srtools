#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use threads;
use Thread::Queue;
use CommonFunctions qw(parseFileList);

my (%opt, $outfile);
our @files : shared;
getopts('f:o:h', \%opt);
var_check();

my $threads = scalar(@files);
my $workq = Thread::Queue->new();
my $retq = Thread::Queue->new();

for (my $i = 1; $i <= $threads; $i++) {
  $workq->enqueue($i);
}
for (my $i = 1; $i <= $threads; $i++) {
  $workq->enqueue('EXIT');
}
for (my $i = 1; $i <= $threads; $i++) {
  threads->create("stats");
}

open (OUT, ">$outfile") or die " Cannot open $outfile: $!\n\n";

while (threads->list(threads::running)) {
  while ($retq->pending) {
    my $result = $retq->dequeue;
    print OUT $result."\n";
  }
  sleep 1;
}
while ($retq->pending) {
  my $result = $retq->dequeue;
  print OUT $result."\n";
}

close OUT;
exit;

sub stats {
  while (my $todo = $workq->dequeue()) {
    last if ($todo eq 'EXIT');
    my $file = $files[$todo - 1];
    my %list;
    my $total = 0;
    open (SAM, $file) or die " Cannot open file $file: $!\n\n";
    while (my $line = <SAM>) {
      next if (substr($line,0,1) eq '@');
      my @tmp = split /\t/, $line;
      next if ($tmp[2] eq '*');
      $list{$tmp[0]} = 1;
    }
    close SAM;
    foreach my $name (keys(%list)) {
      my ($sid, $reads) = split /:/, $name;
      $total += $reads;
    }
    $retq->enqueue("$file $total");
  }
  threads->detach;
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
    $outfile = $opt{'o'};
  } else {
    var_error();
  }
}

sub var_error {
  print STDERR " This script will calculate how many reads were mapped to the genome.\n";
  print STDERR " Usage: samstats.pl -f <SAM file(s)> -o <output file>\n\n";
  print STDERR "      -f     SAM file(s)\n\n";
  print STDERR "      -o     Output log file\n\n";
  print STDERR "      -h     Print this menu\n\n";
  exit 1;
}