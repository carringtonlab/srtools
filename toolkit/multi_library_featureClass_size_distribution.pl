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

my (%opt, @list, $outfile, $species, $confFile, @range, @gffs, @types, %types, %read_table, @ftable);

getopts('l:o:s:c:r:f:t:h', \%opt);
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

foreach my $type (@types) {
	foreach my $size (@range) {
		$types{$type}->{$size} = 0;
	}
}
foreach my $size (@range) {
	$types{'none'}->{$size} = 0;
}

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

print STDERR " Getting features... ";
foreach my $file (@gffs) {
	open (GFF, $file) or die " Cannot open $file: $!\n\n";
	while (my $line = <GFF>) {
		next if (substr($line,0,1) eq '#');
		chomp $line;
		my ($ref, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $line;
		next if (!exists($types{$type}));
		my %hash;
		$hash{'ref'} = $ref;
		$hash{'start'} = $start;
		$hash{'end'} = $end;
		$hash{'strand'} = $strand;
		$hash{'type'} = $type;
		push @ftable, \%hash;
	}
	close GFF;
}
@ftable = sort {$a->{'ref'} cmp $b->{'ref'} || $a->{'start'} <=> $b->{'start'}} @ftable;
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

print STDERR " Summarizing data... ";
open SAM, "$sam view $bam |";
while (my $line = <SAM>) {
	my @tmp = split /\t/, $line;
	my $sid = $tmp[0];
	my $ref = $tmp[2];
	my $start = $tmp[3];
	my $length = length($tmp[9]);
	my $end = $start + $length - 1;
	my $a = 0;
	foreach my $f (@ftable) {
		next if ($f->{'ref'} ne $ref);
		last if ($f->{'start'} > $end);
		if ($start >= $f->{'start'} && $start <= $f->{'end'}) {
			$types{$f->{'type'}}->{$length}++;
		} elsif ($end >= $f->{'start'} && $end <= $f->{'end'}) {
			$types{$f->{'type'}}->{$length}++;
		}
	}
	if ($a == 0) {
		$types{'none'}->{$length}++;
	}
}
close SAM;
print STDERR "done\n";

print STDERR " Printing out results... ";
open (OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
print OUT "Type\t".join("\t",@range)."\n";
while (my ($type, $r) = each(%types)) {
	print OUT $type;
	foreach my $size (@range) {
		print OUT "\t".$r->{$size};
	}
	print "\n";
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
	#if ($opt{'r'}) {
	#	@range = parseListToArray($opt{'r'});
	#} else {
	#	var_error();
	#}
	if ($opt{'f'}) {
		@gffs = parseFileList($opt{'f'});
	} else {
		var_error();
	}
	if ($opt{'t'}) {
		@types = parseFileList($opt{'t'});
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
	print STDERR "  multi_library_featureClass_size_distribution.pl -l <library_ids> -o <output file> -c <conf_file> -s <species> -f <GFF file(s)> -t <type(s)>\n";
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
	print STDERR "  -f   The GFF feature file(s).\n";
	print STDERR "\n";
	print STDERR "  -t   The feature type(s).\n";
	print STDERR "\n";
#	print STDERR "  -r   The small RNA size range.\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
