#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use POSIX qw( ceil floor);
use DBI;
use Config::Tiny;
use Env qw(HOME);
use lib "$HOME/lib/perl";
use CommonFunctions qw(parseListToArray parseFileList);

#########################################################
# Start Variable declarations                           #
#########################################################

my (%opt, @list, $outfile, $species, $confFile, @range, %read_table, $gff, %results, %counts, %ftb, %match);

getopts('l:o:s:c:r:f:h', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};
my $bam = $conf->{'bam'};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

for (18..30) {
	push @range, $_;
}

foreach my $size (@range) {
	$counts{$size} = 0;
}

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

print STDERR " Initializing feature array... ";
open SAM, "$sam view -H $bam |";
while (my $line = <SAM>) {
	if (substr($line,0,3) eq '@SQ') {
		chomp $line;
		my ($tag, $refName, $length) = split /\s+/, $line;
		$refName =~ s/SN://;
		$length =~ s/LN://;
		my @p;
		for (my $i = 0; $i <= $length; $i++) {
			$p[$i] = 0;
		}
		$ftb{$refName} = \@p;
	}
}
close SAM;

print STDERR "done\n";

print STDERR " Adding features... ";
open (GFF, $gff) or die " Cannot open $gff: $!\n\n";
while (my $line = <GFF>) {
	next if (substr($line,0,1) eq '#');
	chomp $line;
	my ($ref, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $line;
	for (my $i = $start; $i <= $end; $i++) {
		$ftb{$ref}->[$i] = $type;
	}
}
close GFF;
print STDERR "done\n";

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

print STDERR " Counting reads... ";
open SAM, "$sam view $bam |";
while (my $line = <SAM>) {
	my @tmp = split /\t/, $line;
	my $sid = $tmp[0];
	next if (!exists($read_table{$sid}));
	my $ref = $tmp[2];
	my $start = $tmp[3];
	my $length = length($tmp[9]);
	my $end = $start + $length - 1;
	for (my $i = $start; $i <= $end; $i++) {
		if ($ftb{$ref}->[$i]) {
			#print $ftb{$ref}->[$i]."\n";
			$counts{$length}++;
			last;
		}
	}
}
close SAM;
print STDERR "done\n";

print STDERR " Printing out results... ";
open (OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
print OUT "Length\tHits\n";
foreach my $size (@range) {
	print OUT $size."\t".$counts{$size}."\n";
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
	#if ($opt{'r'}) {
	#	@range = parseListToArray($opt{'r'});
	#} else {
	#	var_error();
	#}
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
	print STDERR "  multi_library_featureClass_size_distribution2.pl -l <library_ids> -o <output file> -c <conf_file> -s <species> -f <GFF file>\n";
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
#	print STDERR "  -r   The small RNA size range.\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
