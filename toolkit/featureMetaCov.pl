#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use POSIX qw( ceil floor);
use DBI;
use Config::Tiny;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use CommonFunctions qw(parseListToArray parseFileList);

#########################################################
# Start Variable declarations                           #
#########################################################

my (%opt, @list, $outfile, $species, $confFile, @range, $gff, $ctype, %sizes, %libs, $flankLength, %up, %in, %dw, %hits);

getopts('l:o:s:c:r:f:e:C:h', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};
my $bam = $conf->{'bam'};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

foreach my $size (@range) {
	$sizes{$size} = 1;
}
foreach my $lib (@list) {
	$libs{$lib} = 1;
}

for (my $i = 1; $i <= 100; $i++) {
	foreach my $size (@range) {
		$up{$i}->{$size} = 0;
		$in{$i}->{$size} = 0;
		$dw{$i}->{$size} = 0;
	}
}


#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

open (OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
print OUT "Position\t".join("\t", @range)."\n";

print STDERR " Calculating per small RNA hits counts... ";
open SAM, "$sam view $bam |";
while (my $line = <SAM>) {
	my @tmp = split /\t/, $line;
	my $sid = $tmp[0];
	if (exists($hits{$sid})) {
		$hits{$sid}++;
	} else {
		$hits{$sid} = 1;
	}
}
close SAM;
print STDERR "done\n";

my $count = 0;
my $sth = $dbh->prepare('SELECT * FROM `reads` WHERE `sid` = ?');
print STDERR " Processing features... $count\r";
open (GFF, $gff) or die " Cannot open $gff: $!\n\n";
while (my $line = <GFF>) {
	next if (substr($line,0,1) eq '#');
	chomp $line;
	my ($ref, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $line;
	my $fstart = $start - 500;
	my $fend = $end + 500;
	if ($fstart < 1) {
		$fstart = 1;
	}
	
	open SAM, "$sam view $bam '$ref:$fstart-$fend' |";
	while (my $line = <SAM>) {
		my @tmp = split /\t/, $line;
		my $sid = $tmp[0];
		my $length = length($tmp[9]);
		next if (!exists($sizes{$length}));
		my $reads = 0;
		$sth->execute($sid);
		while (my $row = $sth->fetchrow_hashref) {
			if (exists($libs{$row->{'library_id'}})) {
				$reads += $row->{'reads'};
			}
		}
		my $p = $tmp[3];
		if ($tmp[1] == 16) {
			$p += $length - 1;
		}
		next if ($p < $fstart || $p > $fend);
		if ($strand eq '+') {
			if ($p < $start) {
				my $rlen = $start - $fstart;
				if ($ctype eq 'reads') {
					$up{int(($p / $rlen) * 100) + 1}->{$length} += ($reads / $hits{$sid});	
				} elsif ($ctype eq 'seqs') {
					$up{int(($p / $rlen) * 100) + 1}->{$length} += (1 / $hits{$sid});	
				}
			} elsif ($p >= $start && $p <= $end) {
				my $rlen = $end - $start + 1;
				if ($ctype eq 'reads') {
					$in{int(($p / $rlen) * 100) + 1}->{$length} += ($reads / $hits{$sid});	
				} elsif ($ctype eq 'seqs') {
					$in{int(($p / $rlen) * 100) + 1}->{$length} += (1 / $hits{$sid});
				}
			} elsif ($p > $end) {
				my $rlen = $fend - $end;
				if ($ctype eq 'reads') {
					$dw{int(($p / $rlen) * 100) + 1}->{$length} += ($reads / $hits{$sid});
				} elsif ($ctype eq 'seqs') {
					$dw{int(($p / $rlen) * 100) + 1}->{$length} += (1 / $hits{$sid});
				}
			}
		} elsif ($strand eq '-') {
			if ($p > $end) {
				my $rlen = $fend - $end;
				if ($ctype eq 'reads') {
					$up{int((($fend - $p + 1) / $rlen) * 100) + 1}->{$length} += ($reads / $hits{$sid});
				} elsif ($ctype eq 'seqs') {
					$up{int((($fend - $p + 1) / $rlen) * 100) + 1}->{$length} += (1 / $hits{$sid});
				}
			} elsif ($p >= $start && $p <= $end) {
				my $rlen = $end - $start + 1;
				if ($ctype eq 'reads') {
					$in{int((($end - $p + 1) / $rlen) * 100) + 1}->{$length} += ($reads / $hits{$sid});
				} elsif ($ctype eq 'seqs') {
					$in{int((($end - $p + 1) / $rlen) * 100) + 1}->{$length} += (1 / $hits{$sid});
				}
			} elsif ($p < $start) {
				my $rlen = $start - $fstart;
				if ($ctype eq 'reads') {
					$dw{int((($start - $p + 1) / $rlen) * 100) + 1}->{$length} += ($reads / $hits{$sid});
				} elsif ($ctype eq 'seqs') {
					$dw{int((($start - $p + 1) / $rlen) * 100) + 1}->{$length} += (1 / $hits{$sid});
				}
			}
		}
	}
	close SAM;
	$count++;
	print STDERR " Processing features... $count\r";
}
close GFF;
print STDERR " Processing features... $count\n";

# upstream
for (my $i = 1; $i <= 100; $i++) {
	print OUT $i;
	foreach my $size (@range) {
		print OUT "\t".$up{$i}->{$size} / $count;
	}
	print OUT "\n";
}

# intragenic
for (my $i = 1; $i <= 100; $i++) {
	print OUT $i;
	foreach my $size (@range) {
		print OUT "\t".$in{$i}->{$size} / $count;
	}
	print OUT "\n";
}

# downstream
for (my $i = 1; $i <= 100; $i++) {
	print OUT $i;
	foreach my $size (@range) {
		print OUT "\t".$dw{$i}->{$size} / $count;
	}
	print OUT "\n";
}

close OUT;
exit;


#########################################################
# Start Subroutines                                     #
#########################################################

#########################################################
# Start of Varriable Check Subroutine "var_check"       #
#########################################################

sub var_check {
	if ($opt{'h'}) {
		var_error();
	}
	if ($opt{'l'}) {
		@list = parseListToArray($opt{'l'});
	} else {
		var_error();
	}
	if ($opt{'o'}) {
		$outfile = $opt{'o'};
	} else {
		var_error();
	}
	if ($opt{'c'}) {
		$confFile = $opt{'c'};
	} else {
		var_error();
	}
	if ($opt{'s'}) {
		$species = $opt{'s'};
	} else {
		var_error();
	}
	if ($opt{'r'}) {
		@range = parseListToArray($opt{'r'});
	} else {
		@range = parseListToArray('18-30');
	}
	if ($opt{'f'}) {
		$gff = $opt{'f'};
	} else {
		var_error();
	}
	if ($opt{'e'}) {
		$flankLength = $opt{'e'};
	} else {
		$flankLength = 500;
	}
	if ($opt{'C'}) {
		$ctype = $opt{'C'};
		unless ($ctype eq 'reads' || $ctype eq 'seqs') {
			var_error();
		}
	} else {
		$ctype = 'reads';
	}
}

#########################################################
# End of Varriable Check Subroutine "var_check"         #
#########################################################

#########################################################
# Start of Varriable error Subroutine "var_error"       #
#########################################################

sub var_error {
	print STDERR "\n  Description:\n";
	print STDERR "  This script will calculate the meta-coverage of sequencing reads for the given features.\n\n";
	print STDERR "  Usage:\n";
	print STDERR "  featureMetaCov.pl -l <library_ids> -o <output file> -c <conf_file> -s <species> -f <GFF file>\n";
	print STDERR "\n\n";
	print STDERR "  -l   The library ID's to use.\n";
	print STDERR "               example: -l '1-5'\n";
	print STDERR "               example: -l '1,7,11'\n";
	print STDERR "               example: -l '1-5,7,9-11'\n";
	print STDERR "\n";
	print STDERR "  -o   The output filename.\n";
	print STDERR "\n";
	print STDERR "  -c   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -s   The species.\n";
	print STDERR "\n";
	print STDERR "  -f   The GFF feature file.\n";
	print STDERR "\n";
	print STDERR "  -e   Extended region (upstream/downstream) length (in nucleotides). Default = 500\n";
	print STDERR "\n";
	print STDERR "  -r   The small RNA size range. Default = 18-30\n";
	print STDERR "\n";
	print STDERR "  -C   Count reads or sequences [reads | seqs]. Default = reads\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
