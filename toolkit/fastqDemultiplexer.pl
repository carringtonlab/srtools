#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Cwd;
use FileHandle;
use Env qw(HOME);
use lib "$HOME/lib/perl";
use CommonFunctions qw(parseFileList);

############## Begin variables ##############

my (%opt, $prefix, $file, $verbose, @indexList, %fh);
getopts('f:i:o:vh', \%opt);
var_check();

############## End variables ##############

############## Begin main program ##############

my $indexLength = 0;
foreach my $tag (@indexList) {
	$indexLength = length($tag);
	$fh{$tag} = FileHandle->new;
	open ($fh{$tag}, ">".$prefix."_".$tag.".fastq") or die " Cannot open ".$prefix."_".$tag.".fastq\n";
}

if ($file =~ /\.gz$/) {
	open IN, "zcat $file |";	
} else {
	open (IN, $file) or die " Cannot open file $file: $!\n\n";	
}

while (1) {
	my $head = <IN>;
	my $seq = <IN>;
	my $head2 = <IN>;
	my $qual = <IN>;
	chomp $seq;
	
	my $parsed = 0;
  for (my $i = length($seq) - $indexLength; $i >= 17; $i--) { # Modban is 17, so the index has to be at the 18th position or higher
    foreach my $tag (@indexList) {
      if (substr($seq,$i,$indexLength) eq $tag) {
				print {$fh{$tag}} $head;
				print {$fh{$tag}} $seq."\n";
				print {$fh{$tag}} $head2;
				print {$fh{$tag}} $qual;
				$parsed = 1;
				last;
      }
    }
    last if ($parsed == 1);
  }	
	last if (eof(IN));
}
close IN;

foreach my $tag (@indexList) {
	close $fh{$tag};
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
	if ($opt{'o'}) {
    $prefix = $opt{'o'};
  } else {
    var_error();
  }
	if ($opt{'i'}) {
		@indexList = parseFileList($opt{'i'});
	} else {
		var_error();
	}
}

# var_error prints out command-line options, defaults, optional settings
sub var_error {
  print STDERR "\n\n";
  print STDERR " This script will separate reads into library sets based on adaptor sequences.\n";
  print STDERR " Usage: fastqDemultiplexer.pl -f <sequence file> -o <prefix> -i <index sequences>\n\n";
	print STDERR " REQUIRED:\n";
  print STDERR "   -f     Input sequence file.\n\n";
  print STDERR "   -o     Output file prefix.\n\n";
	print STDERR "   -i     Index sequences. List the index sequences in a comma-separated list.\n\n";
	print STDERR " OPTIONAL:\n";
  print STDERR "   -v     Verbose output. Prints program status to terminal. Default is quiet.\n\n";
  print STDERR "   -h     Print this menu\n\n";
  exit 1;
}
