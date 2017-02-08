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

my (%opt, @libs, $species, $confFile, @range, %hits, %sizes, %hist, $outfile);

getopts('l:s:c:r:o:h', \%opt);
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

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

print STDERR " Getting reads... ";

my $sth = $dbh->prepare('SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `library_id` = ?');
foreach my $lib (@libs) {
	$sth->execute($lib);
	while(my $row = $sth->fetchrow_hashref){
		my $length = length($row->{'seq'});
		if (exists($sizes{$length})) {
			$hits{$row->{'sid'}} = 0;
		}
	}
}
print STDERR "done\n";

print STDERR " Counting hits... ";
open SAM, "$sam view $conf->{'bam'} |";
while (my $line = <SAM>) {
	my @tmp = split /\t/, $line;
	if (exists($hits{$tmp[0]})) {
		$hits{$tmp[0]}++;
	}
}
close SAM;
print STDERR "done\n\n";

print STDERR " Couting hit frequencies... ";
while (my ($sid, $hits) = each(%hits)) {
	if (exists($hist{$hits})) {
		$hist{$hits}++;
	} else {
		$hist{$hits} = 1;
	}
}

open (OUT, ">$outfile") or die " Cannot open $outfile: $!\n\n";
print OUT "hits\tcount\n";
foreach my $hit (sort {$a <=> $b} keys(%hist)) {
	print OUT $hit."\t".$hist{$hit}."\n";
}
close OUT;

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
		@libs = parseListToArray($opt{'l'});
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
	if ($opt{'o'}) {
		$outfile = $opt{'o'};
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
	print STDERR "  This script will count the number of genome alignments per sequence for the given libraries.\n\n";
	print STDERR "  Usage:\n";
	print STDERR "  library_size_distribution.pl -l <library_id(s)> -o <output file> -c <conf_file> -s <species>\n";
	print STDERR "\n\n";
	print STDERR "  -l   The library ID list.\n";
	print STDERR "\n";
	print STDERR "  -c   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -s   The species. DEFAULT = all\n";
	print STDERR "\n";
	print STDERR "  -o   The output file.\n";
	print STDERR "\n";
	print STDERR "  -r   The small RNA size range.\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
