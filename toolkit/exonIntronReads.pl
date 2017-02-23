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

my (%opt, @list, $outfile, $species, $confFile, @range, $gff, %sizes, %libs, %genes);

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
print OUT "Feature\tExonReads\tIntronReads\n";

my $sth = $dbh->prepare('SELECT * FROM `reads` WHERE `sid` = ?');
open (GFF, $gff) or die " Cannot open $gff: $!\n\n";
while (my $line = <GFF>) {
	next if (substr($line,0,1) eq '#');
	chomp $line;
	my ($ref, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $line;
	my @attributes = split /;/, $attributes;
	my %tags;
	foreach my $item (@attributes) {
		my ($tag, $value) = split /=/, $item;
		$tags{$tag} = $value;
	}
	if ($type eq 'CDS') {
		my $gene = $tags{'Parent'};
		$gene =~ s/PITT/PITG/;
		my %hash;
		$hash{'ref'} = $ref;
		$hash{'start'} = $start;
		$hash{'end'} = $end;
		$hash{'strand'} = $strand;
		push @{$genes{$gene}}, \%hash;
	}
}	
close GFF;
	
my $count = 0;
print STDERR " Processing features... $count\r";
foreach my $gene (keys(%genes)) {
	my %reads;
	$reads{'exon'} = 0;
	$reads{'intron'} = 0;

	@{$genes{$gene}} = sort {$a->{'start'} <=> $b->{'start'}} @{$genes{$gene}};
	
	my $step = 1;
	my $last;
	foreach my $cds (@{$genes{$gene}}) {
		if ($step == 1) {
			my $reads = get_rnas($cds->{'ref'},$cds->{'start'},$cds->{'end'},$sth);
			$reads{'exon'} += $reads;
			$last = $cds->{'end'};
			$step = 2;
		} elsif ($step == 2) {
			my $in_start = $last + 1;
			my $in_end = $cds->{'start'} - 1;
			if ($in_end - $in_start > 0) {
				my $reads = get_rnas($cds->{'ref'},$last + 1,$cds->{'start'} - 1,$sth);
				$reads{'intron'} += $reads;	
			}
			my $reads = get_rnas($cds->{'ref'},$cds->{'start'},$cds->{'end'},$sth);
			$reads{'exon'} += $reads;
			$last = $cds->{'end'};
		}
	}
	
	print OUT $gene;
	print OUT "\t".$reads{'exon'}."\t".$reads{'intron'}."\n";
	$count++;
	print STDERR " Processing features... $count\r";
}
print STDERR " Processing features... $count\n";

close OUT;
exit;


#########################################################
# Start Subroutines                                     #
#########################################################

#########################################################
# Start of Varriable Check Subroutine "var_check"       #
#########################################################

sub get_rnas {
	my $ref = shift;
	my $start = shift;
	my $end = shift;
	my $sth = shift;
	my $reads = 0;
	#print $ref.':'.$start.'-'.$end."\n";
	open SAM, "$sam view $bam '$ref:$start-$end' |";
	while (my $line = <SAM>) {
		my @tmp = split /\t/, $line;
		my $sid = $tmp[0];
		my $length = length($tmp[9]);
		next if (!exists($sizes{$length}));
		$sth->execute($sid);
		while (my $row = $sth->fetchrow_hashref) {
			if (exists($libs{$row->{'library_id'}})) {
				$reads += $row->{'reads'};
			}
		}		
	}		
	close SAM;
	return $reads / ($end - $start + 1);
}

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
	print STDERR "  exonIntronReads.pl -l <library_ids> -o <output file> -c <conf_file> -s <species> -f <GFF file>\n";
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
