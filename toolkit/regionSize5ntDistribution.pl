#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use POSIX qw( ceil floor);
use DBI;
use Config::Tiny;
use Env qw(HOME);
use lib "$HOME/lib/perl";
use CommonFunctions qw(parseListToArray);

#########################################################
# Start Variable Declaration                            #
#########################################################

my (%opt, $outfile, $species, $confFile, @libs, %rpm, @sizes, $sth, $region, %table, %libs);

getopts('o:S:C:l:c:R:hr', \%opt);
var_check();

my @nts = ('A', 'T', 'G', 'C');

foreach my $library (@libs) {
	unless ($library =~ /^\d+$/) {
		croak("Error: value $library is not a valid library ID (Integer, only one can be specified)\n");
	}
	$libs{$library} = 1;
}

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

if (!@sizes) {
	$sth = $dbh->prepare('SELECT DISTINCT(LENGTH(`seq`)) AS size FROM `sequences`');
	$sth->execute();
	while (my $row = $sth->fetchrow_hashref) {
		push @sizes, $row->{'size'};
	}
}

foreach my $size (@sizes) {
	foreach my $library (@libs) {
		foreach my $nt (@nts) {
			$table{$size}->{$library}->{$nt} = 0;
		}
	}
}

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

open (OUT, ">$outfile") or die " Cannot open $outfile for writing: $!\n\n";

if ($opt{'r'}) {
	print STDERR " Getting total reads... ";
	$sth = $dbh->prepare('SELECT `total_reads` FROM `libraries` WHERE `library_id` = ?');
	foreach my $library (@libs) {
		$sth->execute($library);
		my $result = $sth->fetchrow_hashref;
		$rpm{$library} = 1000000 / $result->{'total_reads'};
	}
	print STDERR "done\n";
} else {
	foreach my $library (@libs) {
		$rpm{$library} = 1;
	}
}

print STDERR " Extracting region $region... ";
$sth = $dbh->prepare('SELECT * FROM `reads` WHERE `sid` = ?');
open SAM, "$sam view $conf->{'bam'} '$region' |";
while (my $hit = <SAM>) {
	chomp $hit;
	my @fields = split /\t/, $hit;
	my $sid = $fields[0];
	my $strand = $fields[1];
	my $seq = $fields[9];
	my $size = length($seq);
	if ($strand == 16) {
		$seq = reverse $seq;
		$seq =~ tr/ATGC/TACG/;
	}
	my $first_nt = substr($seq,0,1);
	$sth->execute($sid);
	while (my $row = $sth->fetchrow_hashref) {
		if (exists($libs{$row->{'library_id'}})) {
			$table{$size}->{$row->{'library_id'}}->{$first_nt} += ($row->{'reads'} * $rpm{$row->{'library_id'}});
		}
	}
	
}
close SAM;

print STDERR "done\n";

print STDERR " Print out results... ";
foreach my $size (@sizes) {
	print OUT "Table for $size length small RNA\n";
	print OUT "Library\t".join("\t", @nts)."\n";
	foreach my $library (@libs) {
		print OUT $library;
		foreach my $nt (@nts) {
			print OUT "\t".$table{$size}->{$library}->{$nt};
		}
		print OUT "\n";
	}
	print OUT "\n";
}
print STDERR "done\n";

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
		$outfile = $opt{'o'};
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
	if ($opt{'R'}) {
		$region = $opt{'R'};
	} else {
		var_error();
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
	print STDERR "  This script will output read size and 5' nucleotide distributions for the specified region\n";
	print STDERR "  Usage:\n";
	print STDERR "  regionSize5ntDistribution.pl -s <species> -c <conf file> -l <library ID(s)> -o <output file> -R <region>\n";
	print STDERR "\n\n";
	print STDERR "  -l   The sample library ID(s)\n";
	print STDERR "			      example: -l 1\n";
	print STDERR "			      example: -l '1-5'\n";
	print STDERR "			      example: -l '1-3,5-10'\n";
	print STDERR "\n";
	print STDERR "  -S   The species.\n";
	print STDERR "\n";
	print STDERR "  -o   The output file.\n";
	print STDERR "\n";
	print STDERR "  -C   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -R   The region to extract. Format is accession:start-end\n";
	print STDERR "\n";
	print STDERR "  -c   The sequence sizes to use (DEFAULT = all).\n";
	print STDERR "			      example: -c '20-25'\n";
	print STDERR "			      example: -c '20,21,24'\n";
	print STDERR "			      example: -c '20-22,24-30'\n";
	print STDERR "\n";
	print STDERR "  -r   RPM normalize reads?\n";
	print STDERR "\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################

