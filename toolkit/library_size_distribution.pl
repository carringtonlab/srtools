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
# Start Variable declarations                           #
#########################################################

my (%opt, $libID, $outfile, $species, $confFile, @range, %read_table, %results, %nts, $sth, $rpm);

getopts('l:o:s:c:r:Rh', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");
if ($opt{'R'}) {
	$sth = $dbh->prepare("SELECT * FROM `libraries` WHERE `library_id` = $libID");
	$sth->execute();
	my $row = $sth->fetchrow_hashref;
	if ($row->{'total_reads'} > 0) {
		$rpm = 1000000 / $row->{'total_reads'};
	} else {
		print STDERR " Total reads for library $libID not set in $conf->{'db'}.libraries\n\n";
		exit 1;
	}	
} else {
	$rpm = 1;
}

foreach my $size (@range) {
	$results{$size} = 0;
	my %nt;
	$nt{'A'} = 0;
	$nt{'U'} = 0;
	$nt{'G'} = 0;
	$nt{'C'} = 0;
	$nts{$size} = \%nt;
}

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

print STDERR " Getting reads... ";

$sth = $dbh->prepare('SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `library_id` = ?');
$sth->execute($libID);
while(my $row = $sth->fetchrow_hashref){
	my $length = length($row->{'seq'});
	my $first = substr($row->{'seq'},0,1);
	$first =~ s/T/U/g;
	$results{$length} += ($row->{'reads'} * $rpm);
	$nts{$length}->{$first} += ($row->{'reads'} * $rpm);
}
print STDERR "done\n";

print STDERR " Printing out results... ";
open (OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
print OUT "Length\tTotal reads\t5'A reads\t5'U reads\t5'G reads\t5'C reads\n";
foreach my $size (@range) {
	print OUT "$size\t$results{$size}\t$nts{$size}->{'A'}\t$nts{$size}->{'U'}\t$nts{$size}->{'G'}\t$nts{$size}->{'C'}\n";
}

print STDERR "done\n";
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
		$libID = $opt{'l'};
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
		var_error();
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
	print STDERR "  This script will generate a distribution of sequencing reads by sequence length.\n\n";
	print STDERR "  Usage:\n";
	print STDERR "  library_size_distribution.pl -l <library_id> -o <output file> -c <conf_file> -s <species>\n";
	print STDERR "\n\n";
	print STDERR "  -l   The library ID.\n";
	print STDERR "\n";
	print STDERR "  -o   The output filename.\n";
	print STDERR "\n";
	print STDERR "  -c   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -s   The species.\n";
	print STDERR "\n";
	print STDERR "  -r   The small RNA size range.\n";
	print STDERR "\n";
	print STDERR "  -R   Normalize reads (Reads/Million). Default, do not normalize.\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
