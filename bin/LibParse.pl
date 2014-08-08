#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use CommonFunctions qw(parseFileList);

############## Begin variables ##############

my (%opt, $type, $prefix);
our (%seqs, %stats, @files, $verbose, @startList, @endList, @indexList, $minInsertSize, $maxInsertSize, $readLength, $adaptLen, $scanLen, $trim, $failed, $log);
getopts('f:t:e:E:m:M:r:o:l:c:s:a:vh', \%opt);
var_check();

# Initialize stats hash, stores summary statistics for the run
# qc = failed quality control, sequence contains 1+ N character
# adapt = failed adaptor parsing. Either no adaptor or insert length is shorter/longer than allowed
# total = total reads in file
# adatpor tags are stored to count library-specific reads
$stats{'qc'} = 0;
$stats{'adapt'} = 0;
$stats{'total'} = 0;
if (scalar(@startList) > 0) {
	foreach my $tag (@startList) {
		$stats{$tag} = 0;
	}
	$adaptLen = length($startList[0]);
	$scanLen = $readLength - $adaptLen;
	$trim = '5prime';
} elsif (scalar(@indexList) > 0) {
	foreach my $tag (@indexList) {
		$stats{$tag} = 0;
	}
	$adaptLen = length($indexList[0]);
	$scanLen = $readLength - $adaptLen;
	$trim = 'trim2';
} elsif (scalar(@endList) > 0) {
	foreach my $tag (@endList) {
		$stats{$tag} = 0;
	}
	$adaptLen = length($endList[0]);
	$scanLen = $readLength - $adaptLen;
	if ($maxInsertSize && $maxInsertSize < $scanLen) {
		$scanLen = $maxInsertSize;
	}
	$trim = '3prime';
} else {
	$trim = 'none';
}

if ($failed) {
	open (FAIL, ">$failed") or die " Cannot open file $failed: $!\n\n";
}
if ($log) {
	open (LOG, ">$log") or die " Cannot open file $log: $!\n\n";
}

############## End variables ##############

############## Begin main program ##############

# Parse type input to determine the subroutine to use to parse the input file
if ($type eq 'fasta') {
  ReadFasta(@files);
} elsif ($type eq 'fastq') {
  ReadFastq(@files);
} elsif ($type eq 'eland') {
  ReadEland(@files);
} elsif ($type eq 'raw') {
	ReadScarf(@files);
} elsif ($type eq 'seqonly') {
	ReadSeqOnly(@files);
} elsif ($type eq 'fastq_wustl') {
	ReadFastq_wustl(@files);
} else {
  print STDERR " File type $type is not a recognized format.\n\n";
  var_error();
}

# Print out the run statistics
if ($log) {
	while (my ($k, $v) = each(%stats)) {
	  print LOG $k." => ".$v."\n";
	}
}

# This section will get sequence IDs (sid) for each sequence
# sids, reads and sequences will be output to files per library
# Output is FASTA-format, for input to bowtie or other aligner
if ($trim eq '5prime') {
	foreach my $tag (@startList) {
		my $count = 1;
		open (OUT, ">$prefix\_$tag.parsed.fasta") or die " Cannot open file $prefix\_$tag.parsed.fasta: $!\n\n";
		while (my ($seq, $reads) = each(%{$seqs{$tag}})) {
			print OUT ">$count:$reads\n";
			print OUT "$seq\n";
			$count++;
		}
		close OUT;
	}
} elsif ($trim eq '3prime') {
	foreach my $tag (@endList) {
		my $count = 1;
		open (OUT, ">$prefix\_$tag.parsed.fasta") or die " Cannot open file $prefix\_$tag.parsed.fasta: $!\n\n";
		while (my ($seq, $reads) = each(%{$seqs{$tag}})) {
			print OUT ">$count:$reads\n";
			print OUT "$seq\n";
			$count++;
		}
		close OUT;
	}
} elsif ($trim eq 'trim2') {
	foreach my $tag (@indexList) {
		my $count = 1;
		open (OUT, ">$prefix\_$tag.parsed.fasta") or die " Cannot open file $prefix\_$tag.parsed.fasta: $!\n\n";
		while (my ($seq, $reads) = each(%{$seqs{$tag}})) {
			print OUT ">$count:$reads\n";
			print OUT "$seq\n";
			$count++;
		}
		close OUT;
	}
} else {
	my $count = 1;
	open (OUT, ">$prefix.parsed.fasta") or die " Cannot open file $prefix.parsed.fasta: $!\n\n";
	while (my ($seq, $reads) = each(%{$seqs{'none'}})) {
		print OUT ">$count:$reads\n";
		print OUT "$seq\n";
		$count++;
	}
	close OUT;
}

if ($failed) {
	close FAIL;
}
if ($log) {
	close LOG;
}

exit;

############## End main program ##############

############## Begin subroutines ##############

# ReadFasta parses FASTA-formatted files
# Currently assumes reads are on a single line
sub ReadFasta {
  my @file_list = @_;
	foreach my $file (@file_list) {
		print STDERR " Processing FASTA-formatted $file... \n" if ($verbose);
		print STDERR "   Reads processed... 0\r" if ($verbose);
		if ($file =~ /\.gz$/) {
			open IN, "zcat $file |";	
		} else {
			open (IN, $file) or die " Cannot open file $file: $!\n\n";	
		}
		while (my $line = <IN>) {
			next if (substr($line,0,1) eq '>');
			$stats{'total'}++;
			print STDERR "   Reads processed... $stats{'total'}\r" if ($verbose && $stats{'total'} % 100000 == 0);
			chomp $line;
			trim($line);
		}
		close IN;
		print STDERR "   Reads processed... $stats{'total'}\n" if ($verbose);
	}
}

# ReadFastq parses Illumina FASTQ-formatted files
# Currently assumes reads are on a single line
sub ReadFastq {
	my @file_list = @_;
	foreach my $file (@file_list) {
		print STDERR " Processing FASTQ-formatted $file... \n" if ($verbose);
		print STDERR "   Reads processed... 0\r" if ($verbose);
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
			$stats{'total'}++;
			print STDERR "   Reads processed... $stats{'total'}\r" if ($verbose && $stats{'total'} % 100000 == 0);
			chomp $seq;
			trim($seq);
			last if (eof(IN));
		}
		close IN;
		print STDERR "   Reads processed... $stats{'total'}\n" if ($verbose);
	}
}

sub ReadFastq_wustl {
	my @file_list = @_;
	foreach my $file (@file_list) {
		my $ifile = $file;
		print STDERR " Processing FASTQ(wustl)-formatted $file... \n" if ($verbose);
		print STDERR "   Reads processed... 0\r" if ($verbose);
		# Wustl Fastq files should have reads in file R1
		# Index sequences are in file R2
		if ($file =~ /R1/) {
			$ifile =~ s/R1/R2/;
		} else {
			print STDERR " The file $file does not appear to be in Wustl Fastq format.\n";
			print STDERR " Reads should be in a file called lane#_NoIndex_L00#_R1_00#.fastq\n";
			exit 1;
		}
		unless (-e $ifile) {
			print STDERR " The file $ifile does not exist. Files in Wustl Fastq format should have a paired index file.\n";
			print STDERR " Reads should be in a file called lane#_NoIndex_L00#_R1_00#.fastq\n";
			print STDERR " Index sequences should be in a file called lane#_NoIndex_L00#_R2_00#.fastq\n";
			exit 1;
		}
		if ($file =~ /\.gz$/) {
			open SEQ, "zcat $file |";
			open INDEX, "zcat $ifile |";
		} else {
			open (SEQ, $file) or die " Cannot open file $file: $!\n\n";
			open (INDEX, $ifile) or die " Cannot open file $ifile: $!\n\n";
		}
		while (1) {
			# Get first sequence from the reads file (R1)
			my $headR1 = <SEQ>;
			my $seqR1 = <SEQ>;
			my $head2R1 = <SEQ>;
			my $qualR1 = <SEQ>;
			
			# Get the first sequence from the index file (R2)
			my $headR2 = <INDEX>;
			my $seqR2 = <INDEX>;
			my $head2R2 = <INDEX>;
			my $qualR2 = <INDEX>;
			
			# Make sure the two files line up
			my ($idR1, $pairR1) = split /\s/, $headR1;
			my ($idR2, $pairR2) = split /\s/, $headR2;
			unless ($idR1 eq $idR2) {
				print STDERR " Reads in $file and $ifile do not match up.\n";
				print STDERR " Failed on $idR1 and $idR2\n";
				exit 1;
			}
			
			$stats{'total'}++;
			print STDERR "   Reads processed... $stats{'total'}\r" if ($verbose && $stats{'total'} % 100000 == 0);
			chomp $seqR1;
			chomp $seqR2;
			trim($seqR1.$seqR2);
			last if (eof(SEQ));
		}
		close SEQ;
		close INDEX;
	}
}

# ReadEland parses ELAND-formatted files
sub ReadEland {
  my @file_list = @_;
	foreach my $file (@file_list) {
		print STDERR " Processing ELAND-formatted $file... \n" if ($verbose);
		print STDERR "   Reads processed... 0\r" if ($verbose);
		if ($file =~ /\.gz$/) {
			open IN, "zcat $file |";	
		} else {
			open (IN, $file) or die " Cannot open file $file: $!\n\n";	
		}
		while (my $line = <IN>) {
			$stats{'total'}++;
			print STDERR "   Reads processed... $stats{'total'}\r" if ($verbose && $stats{'total'} % 100000 == 0);
			my @tmp = split /\t/, $line;
			trim($tmp[1]);
		}
		close IN;
		print STDERR "   Reads processed... $stats{'total'}\n" if ($verbose);
	}
}

# SeqOnly parses files with sequences separated by newlines
sub ReadSeqOnly {
	my @file_list = @_;
	foreach my $file (@file_list) {
		print STDERR " Processing sequence file $file... \n" if ($verbose);
		print STDERR "   Reads processed... 0\r" if ($verbose);
		if ($file =~ /\.gz$/) {
			open IN, "zcat $file |";	
		} else {
			open (IN, $file) or die " Cannot open file $file: $!\n\n";	
		}
		while (my $line = <IN>) {
			$stats{'total'}++;
			print STDERR "   Reads processed... $stats{'total'}\r" if ($verbose && $stats{'total'} % 100000 == 0);
			chomp $line;
			trim($line);
		}
		close IN;
		print STDERR "   Reads processed... $stats{'total'}\n" if ($verbose);
	}
}

# Scarf parses Illumina sequence.txt files
sub ReadScarf {
	my @file_list = @_;
	foreach my $file (@file_list) {
		print STDERR " Processing sequence file $file... \n" if ($verbose);
		print STDERR "   Reads processed... 0\r" if ($verbose);
		if ($file =~ /\.gz$/) {
			open IN, "zcat $file |";	
		} else {
			open (IN, $file) or die " Cannot open file $file: $!\n\n";	
		}
		while (my $line = <IN>) {
			$stats{'total'}++;
			print STDERR "   Reads processed... $stats{'total'}\r" if ($verbose && $stats{'total'} % 100000 == 0);
			my @tmp = split /:/, $line;
			trim($tmp[5]);
		}
		close IN;
		print STDERR "   Reads processed... $stats{'total'}\n" if ($verbose);
	}
}

# Trimming controller
sub trim {
	my $seq = shift;
	if ($trim eq '5prime') {
		startParse($seq);
	} elsif ($trim eq '3prime') {
		endParse($seq);
	} elsif ($trim eq 'trim2') {
		endParse2($seq);
	} else {
		if ($seq =~ /N/) {
			$stats{'qc'}++;
		} elsif (length($seq) < $minInsertSize || length($seq) > $maxInsertSize) {
			$stats{'adapt'}++;
		} else {
			indexSeq('none',$seq);
		}
	}
	return;
}

# startParse locates the 5' adaptor/barcode sequence in the input read
sub startParse {
	my $seq = shift;
	my $insert;
	my $parsed = 0;
	foreach my $tag (@startList) {
		if (substr($seq,0,$adaptLen) eq $tag) {
			$parsed = 1;
			$insert = substr($seq,$adaptLen,$scanLen);
			if ($insert =~ /N/) {
				$stats{'qc'}++;
				last;
			}
			$stats{$tag}++;
			indexSeq($tag, $insert);
			last;
		}
	}
	if ($parsed == 0) {
		$stats{'adapt'}++;
		print FAIL $seq."\n" if ($failed);
	}
	return;
}

# endParse locates the 3' adaptor sequence in the input read
# The input sequence is scanned from the min insert size to the end of the read using substr
# Adaptor matches are exact string matches to substr
# Insert sequences are stored in a hash per library to generate unique read sets with read counts
sub endParse {
  my $seq = shift;
  my $insert;
  my $parsed = 0;
  for (my $i = $scanLen; $i >= $minInsertSize; $i--) {
    foreach my $tag (@endList) {
      if (substr($seq,$i,$adaptLen) eq $tag) {
        $parsed = 1;
        $insert = substr($seq,0,$i);
        if ($insert =~ /N/) {
          $stats{'qc'}++;
          last;
        }
        $stats{$tag}++;
				indexSeq($tag, $insert);
        last;
      }
    }
    last if ($parsed == 1);
  }
  if ($parsed == 0) {
    $stats{'adapt'}++;
		print FAIL $seq."\n" if ($failed);
  }
  return;
}

# endParse2 locates the 3' index sequence and the 3' adaptor sequence in the input read
# The input sequence is scanned from the min insert size to the end of the read using substr
# Adaptor matches are exact string matches to substr
# Insert sequences are stored in a hash per library to generate unique read sets with read counts
sub endParse2 {
  my $seq = shift;
  my $insert;
  my $parsed = 0;
  for (my $i = $scanLen; $i >= $minInsertSize; $i--) {
    foreach my $tag (@indexList) {
      if (substr($seq,$i,$adaptLen) eq $tag) {
				my $amplicon = substr($seq,0,$i);
				my $ampliconLen = length($amplicon);
				my $endAdaptLen = length($endList[0]);
				my $scanLen2 = $ampliconLen - $endAdaptLen;
				for (my $j = $scanLen2; $j >= $minInsertSize; $j--) {
					if (substr($amplicon,$j,$endAdaptLen) eq $endList[0]) {
						$insert = substr($amplicon,0,$j);
						$parsed = 1;
						if ($insert =~ /N/) {
							$stats{'qc'}++;
							last;
						}
						$stats{$tag}++;
						indexSeq($tag, $insert);
						last;
					}
				}
      }
    }
    last if ($parsed == 1);
  }
  if ($parsed == 0) {
    $stats{'adapt'}++;
		print FAIL $seq."\n" if ($failed);
  }
  return;
}

# Insert sequences are stored in a hash per library to generate unique read sets with read counts
sub indexSeq {
	my $tag = shift;
  my $seq = shift;
	if (exists($seqs{$tag}->{$seq})) {
    $seqs{$tag}->{$seq}++;
  } else {
    $seqs{$tag}->{$seq} = 1;
  }
  return;
}

# var_check parses command-line options
# Activates var_error if required options missing, sets defaults
sub var_check {
  if ($opt{'h'}) {
    var_error();
  }
	# Required
  if ($opt{'f'}) {
    @files = parseFileList($opt{'f'});
  } else {
    var_error();
  }
  if ($opt{'t'}) {
    $type = $opt{'t'};
  } else {
    var_error();
  }
	if ($opt{'o'}) {
    $prefix = $opt{'o'};
  } else {
    var_error();
  }
	# Optional
	if ($opt{'s'}) {
    @startList = parseFileList($opt{'s'});
		if (!$opt{'r'}) {
			var_error();
		}
  }
  if ($opt{'e'}) {
    @endList = parseFileList($opt{'e'});
		if (!$opt{'r'}) {
			print STDERR " Read-length (-r) must be provided with the -e option\n\n";
			var_error();
		}
  }
	if ($opt{'E'}) {
		@indexList = parseFileList($opt{'E'});
		if (!$opt{'e'} || !$opt{'r'}) {
			print STDERR " Read-length (-r) and end adaptor (-e) must be provided with the -E option\n\n";
			var_error();
		}
	}
	if ($opt{'r'}) {
    $readLength = $opt{'r'};
  }
	# Optional with defaults
	if ($opt{'v'}) {
    $verbose = 1;
  } else {
    # Verbose output is off by default
    $verbose = 0;
  }
  if ($opt{'m'}) {
    $minInsertSize = $opt{'m'};
  } else {
    # Default minimum insert size
    $minInsertSize = 18;
  }
  if ($opt{'M'}) {
    $maxInsertSize = $opt{'M'};
  } else {
    $maxInsertSize = 0;
  }
	if ($opt{'a'}) {
		$failed = $opt{'a'};
	}
	if ($opt{'l'}) {
		$log = $opt{'l'};
	}
}

# var_error prints out command-line options, defaults, optional settings
sub var_error {
  print STDERR "\n\n";
  print STDERR " This script will separate reads into library sets based on adaptor sequences.\n";
  print STDERR " Usage: LibParse.pl -f <sequence file(s)> -t <type> -o <prefix>\n\n";
	print STDERR " REQUIRED:\n";
  print STDERR "   -f     Input sequence file(s). Multiple files will be aggregated into a single sample.\n\n";
  print STDERR "   -t     Type of input file. Supported: fasta, fastq, fastq_wustl, eland, seqonly, scarf\n\n";
  print STDERR "   -o     Output file prefix.\n\n";
	print STDERR " OPTIONAL:\n";
	print STDERR "   -s     Start adaptor(s). Can be a comma-separated list.\n\n";
  print STDERR "   -e     End adaptor(s). Can be comma-separated list.\n\n";
	print STDERR "   -E     End adaptor with index. Provide the end adaptor sequence to -e then list the index sequences here.\n\n";
  print STDERR "   -r     Sequencing read length. Required if -s or -e is used.\n\n";
  print STDERR "   -m     The miniumum insert length. Default = 18 nucleotides.\n\n";
  print STDERR "   -M     The maximum insert length. Default = any size > minimum.\n\n";
	print STDERR "   -a     File to save sequences that failed adaptor parsing.\n\n";
	print STDERR "   -l     File to save run statistics.\n\n";
  print STDERR "   -v     Verbose output. Prints program status to terminal. Default is quiet.\n\n";
  print STDERR "   -h     Print this menu\n\n";
  exit 1;
}
