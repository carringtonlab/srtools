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

my (%opt, $libID, $species, $confFile, @range, %reads, %sizes);

getopts('l:s:c:r:h', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

foreach my $size (@range) {
	$sizes{$size} = 1;
}

my $hits = 0;

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

print STDERR " Getting reads... ";

my $sth = $dbh->prepare('SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `library_id` = ?');
$sth->execute($libID);
while(my $row = $sth->fetchrow_hashref){
	my $length = length($row->{'seq'});
	if (exists($sizes{$length})) {
		$reads{$row->{'sid'}} = $row->{'reads'};
	}
}
print STDERR "done\n";

print STDERR " Counting hits... ";
open SAM, "$sam view $conf->{'bam'} |";
while (my $line = <SAM>) {
	my @tmp = split /\t/, $line;
	if (exists($reads{$tmp[0]})) {
		$hits++;
	}
}
close SAM;
print STDERR "done\n\n";

print $hits."\n";

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
		for (18..30) {
			push @range, $_;
		}
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
	print STDERR "  This script will count the number of genome alignments for a given library.\n\n";
	print STDERR "  Usage:\n";
	print STDERR "  library_size_distribution.pl -l <library_id> -o <output file> -c <conf_file> -s <species>\n";
	print STDERR "\n\n";
	print STDERR "  -l   The library ID.\n";
	print STDERR "\n";
	print STDERR "  -c   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -s   The species. DEFAULT = all\n";
	print STDERR "\n";
	print STDERR "  -r   The small RNA size range.\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
