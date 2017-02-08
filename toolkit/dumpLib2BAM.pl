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

my (%opt, $prefix, $species, $confFile, @libs, %reads, %unapproved, $rpm, @sizes, %sizes, $sth);

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

open (OUT, ">$prefix.sam") or die " Cannot open $prefix.sam for writing: $!\n\n";

print STDERR " Getting reads... ";
$sth = $dbh->prepare('SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `library_id` = ?');
foreach my $library (@libs) {
	$sth->execute($library);
	while(my $row = $sth->fetchrow_hashref) {
		my $length = length($row->{'seq'});
		next if (@sizes && !exists($sizes{$length}));
		if (exists($reads{$row->{'sid'}})) {
			$reads{$row->{'sid'}} += $row->{'reads'};	
		} else {
			$reads{$row->{'sid'}} = $row->{'reads'};	
		}
	}
}
print STDERR "done\n";

print STDERR " Creating SAM file... ";

open SAM, "$sam view -h $conf->{'bam'} |";
while (my $hit = <SAM>) {
	if (substr($hit,0,1) eq '@') {
		print OUT $hit;
	} else {
		my @fields = split /\t/, $hit;
		my $sid = $fields[0];
		if (exists($reads{$sid})) {
			for (my $i = 1; $i <= $reads{$sid}; $i++) {
				print OUT $hit;
			}
		}
	}
}
close SAM;
print STDERR "done\n";

print STDERR " Creating BAM file... ";
`samtools view -bS $prefix.sam -o $prefix.bam`;
`samtools sort -@ 2 -m 50G $prefix.bam $prefix.sorted`;
`samtools index $prefix.sorted.bam`;

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
		$prefix = $opt{'o'};
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
	print STDERR "  This script will read from the sRNAmp database and create a BAM-formated\n";
	print STDERR "  file which lists a record for each sequence in a specified sample library.\n\n";
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
	print STDERR "  -o   The output BAM file prefix.\n";
	print STDERR "\n";
	print STDERR "  -C   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -c   The sequence sizes to use (DEFAULT = all).\n";
	print STDERR "			      example: -c '20-25'\n";
	print STDERR "			      example: -c '20,21,24'\n";
	print STDERR "			      example: -c '20-22,24-30'\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################

