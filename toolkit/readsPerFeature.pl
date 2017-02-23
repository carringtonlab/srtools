#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use POSIX qw( ceil floor);
use DBI;
use Config::Tiny;
use HTML::Entities qw(decode_entities);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use CommonFunctions qw(parseListToArray parseFileList);

#########################################################
# Start Variable declarations                           #
#########################################################

my (%opt, @list, $outfile, $species, $confFile, @range, $gff, %sizes, %libs);

getopts('l:o:s:c:r:f:h', \%opt);
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

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

open (OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
print OUT "Feature\tSenseReads\tAntiReads\tLength\tAnnotation\n";

my $count = 0;
my $sth = $dbh->prepare('SELECT * FROM `reads` WHERE `sid` = ?');
print STDERR " Processing features... $count\r";
open (GFF, $gff) or die " Cannot open $gff: $!\n\n";
while (my $line = <GFF>) {
	next if (substr($line,0,1) eq '#');
	chomp $line;
	my ($ref, $source, $type, $start, $end, $score, $fstrand, $phase, $attributes) = split /\t/, $line;
	my $flength = $end - $start + 1;
	my @attributes = split /;/, $attributes;
	my %tags;
	foreach my $item (@attributes) {
		my ($tag, $value) = split /=/, $item;
		$tags{$tag} = $value;
	}
	my %reads;
	$reads{'sense'} = 0;
	$reads{'anti'} = 0;
	
	open SAM, "$sam view $bam '$ref:$start-$end' |";
	while (my $line = <SAM>) {
		my @tmp = split /\t/, $line;
		my $sid = $tmp[0];
		my $strand;
		if ($fstrand eq '+') {
			$strand = ($tmp[1] == 0) ? 'sense' : 'anti';	
		} elsif ($fstrand eq '-') {
			$strand = ($tmp[1] == 16) ? 'sense' : 'anti';
		}
		
		my $length = length($tmp[9]);
		next if (!exists($sizes{$length}));
		$sth->execute($sid);
		while (my $row = $sth->fetchrow_hashref) {
			if (exists($libs{$row->{'library_id'}})) {
				$reads{$strand} += $row->{'reads'};
			}
		}
	}
	close SAM;
	
	my ($feature, $annotation);
	if (exists($tags{'ID'})) {
		$feature = $tags{'ID'};
	} elsif (exists($tags{'Name'})) {
		$feature = $tags{'Name'};
	} else {
		$feature = $attributes;
	}
	if (exists($tags{'Note'})) {
		$annotation = decode_entities($tags{'Note'});
	} else {
		$annotation = 'none';
	}
	
	print OUT $feature;
	print OUT "\t".$reads{'sense'}."\t".$reads{'anti'};
	print OUT "\t".$flength."\t".$annotation."\n";
	
	$count++;
	print STDERR " Processing features... $count\r";
}
close GFF;
print STDERR " Processing features... $count\n";

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
}

#########################################################
# End of Varriable Check Subroutine "var_check"         #
#########################################################

#########################################################
# Start of Varriable error Subroutine "var_error"       #
#########################################################

sub var_error {
	print STDERR "\n  Description:\n";
	print STDERR "  This script will generate a distribution of sequencing reads by sequence length for the given features.\n\n";
	print STDERR "  Usage:\n";
	print STDERR "  readPerFeature.pl -l <library_ids> -o <output file> -c <conf_file> -s <species> -f <GFF file>\n";
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
	print STDERR "  -r   The small RNA size range. Default = 18-30\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
