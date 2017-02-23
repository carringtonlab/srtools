#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use POSIX qw( ceil floor);
use DBI;
use Config::Tiny;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use CommonFunctions qw(parseListToArray);

#########################################################
# Start Variable Declaration                            #
#########################################################

my (%opt, $gffFile, $species, $confFile, @libs, %reads, %unapproved, $rpm, @sizes, %sizes, $sth);

my $phase = '.';
my $source = 'ASRP';
my $totalReads = 0;

getopts('o:S:C:l:c:hr', \%opt);
var_check();

foreach my $library (@libs) {
	unless ($library =~ /^\d+$/) {
		croak("Error: value $library is not a valid library ID (Integer, only one can be specified)\n");
	}
}

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

if (@sizes) {
	foreach my $size (@sizes) {
		$sizes{$size} = 1;
	}
}

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

open (GFF, ">$gffFile") or die " Cannot open $gffFile for writing: $!\n\n";

if ($opt{'r'}) {
	print STDERR " Getting total reads... ";
	$sth = $dbh->prepare('SELECT `total_reads` FROM `libraries` WHERE `library_id` = ?');
	foreach my $library (@libs) {
		$sth->execute($library);
		my $result = $sth->fetchrow_hashref;
		$totalReads += $result->{'total_reads'};
	}
	$rpm = 1000000 / $totalReads;
	print STDERR "done\n";
}

print STDERR " Getting reads... ";
$sth = $dbh->prepare('SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `library_id` = ?');
foreach my $library (@libs) {
	$sth->execute($library);
	while(my $row = $sth->fetchrow_hashref) {
		my $length = length($row->{'seq'});
		next if (@sizes && !exists($sizes{$length}));
		if (exists($reads{$row->{'sid'}})) {
			if ($opt{'r'}) {
				$reads{$row->{'sid'}} += sprintf("%.2f", $row->{'reads'} * $rpm);	
			} else {
				$reads{$row->{'sid'}} += $row->{'reads'};	
			}
		} else {
			if ($opt{'r'}) {
				$reads{$row->{'sid'}} = sprintf("%.2f", $row->{'reads'} * $rpm);	
			} else {
				$reads{$row->{'sid'}} = $row->{'reads'};	
			}
		}
	}
}
print STDERR "done\n";

print STDERR " Creating GFF file... ";

open SAM, "$sam view $conf->{'bam'} |";
while (my $hit = <SAM>) {
	chomp $hit;
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
	my $score = $reads{$sid};
	my $attributes = 'Name='.$sid.';';
	
	print GFF $chrom."\t";
	print GFF $source."\t";
	print GFF 'smallRNA'.join("_", @libs)."\t";
	print GFF $start."\t";
	print GFF $end."\t";
	print GFF $score."\t";
	print GFF $strand."\t";
	print GFF $phase."\t";
	print GFF $attributes."\n";
}
close SAM;

print STDERR "done\n";

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
	if ($opt{'S'}) {
		$species = $opt{'S'};
	} else {
		var_error();
	}
	if ($opt{'C'}) {
		$confFile = $opt{'C'};
	} else {
		var_error();
	}
	if ($opt{'o'}) {
		$gffFile = $opt{'o'};
	} else {
		var_error();
	}
	if ($opt{'l'}) {
		@libs = parseListToArray($opt{'l'});
	} else {
		var_error();
	}
	if ($opt{'c'}) {
		@sizes = parseListToArray($opt{'c'});
	}
}

#########################################################
# End of Varriable Check Subroutine "var_check"         #
#########################################################

#########################################################
# Start of Varriable error Subroutine "var_error"       #
#########################################################

sub var_error() {
	print STDERR "\n  Description:\n";
	print STDERR "  This script will read from the sRNAmp database and create a GFF3-formated\n";
	print STDERR "  file which lists a record for each sequence in a specified sample library. Reads will be normalized to total.\n\n";
	print STDERR "  Usage:\n";
	print STDERR "  dumpLib2gff3.pl -s <species> -c <conf file> -l <library ID(s)> -o <output file>\n";
	print STDERR "\n\n";
	print STDERR "  -l   The sample library ID(s)\n";
	print STDERR "			      example: -l 1\n";
	print STDERR "			      example: -l '1-5'\n";
	print STDERR "			      example: -l '1-3,5-10'\n";
	print STDERR "\n";
	print STDERR "  -S   The species.\n";
	print STDERR "\n";
	print STDERR "  -o   The output GFF3 file.\n";
	print STDERR "\n";
	print STDERR "  -C   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -c   The sequence sizes to use (DEFAULT = all).\n";
	print STDERR "			      example: -c '20-25'\n";
	print STDERR "			      example: -c '20,21,24'\n";
	print STDERR "			      example: -c '20-22,24-30'\n";
	print STDERR "  -r   RPM normalize reads?\n\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################

