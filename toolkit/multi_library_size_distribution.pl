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
# Start Variable declarations                           #
#########################################################

my (%opt, @list, $outfile, $species, $confFile, @range, %read_table, %results, %nts);

getopts('l:o:s:c:r:h', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

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

print STDERR " Getting reads... \n";

my $sth = $dbh->prepare('SELECT * FROM `reads` WHERE `library_id` = ?');
foreach my $library_id (@list) {
	print STDERR "      getting reads for library $library_id... ";
	$sth->execute($library_id);
	while(my $row = $sth->fetchrow_hashref){
		if (exists($read_table{$row->{'sid'}})) {
			$read_table{$row->{'sid'}} += $row->{'reads'};
		} else {
			$read_table{$row->{'sid'}} = $row->{'reads'};
		}
	}
	print STDERR "done\n";
}

print STDERR " Summarizing data... ";
$sth = $dbh->prepare('SELECT * FROM `sequences`');
$sth->execute();
while(my $row = $sth->fetchrow_hashref){
	next if (!exists($read_table{$row->{'sid'}}));
	my $length = length($row->{'seq'});
	my $first = substr($row->{'seq'},0,1);
	$first =~ s/T/U/g;
	$results{$length} += $read_table{$row->{'sid'}};
	$nts{$length}->{$first} += $read_table{$row->{'sid'}};
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
	print STDERR "  library_size_distribution.pl -l <library_ids> -o <output file> -c <conf_file> -s <species>\n";
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
	print STDERR "  -r   The small RNA size range.\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
