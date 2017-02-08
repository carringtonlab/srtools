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

my (%opt, $ref_id, $tasi_id, $accession, $confFile, $outfile, $species, $range, $sth, %rpm, %reads, %registers, %pos, $start, $end);

getopts('r:t:c:s:o:a:p:hR', \%opt);
var_check();

open (OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

# Initialize register table
for (1..21) {
	$registers{$_} = 0;
}

# Get range values
if ($range) {
	($start, $end) = split /-/, $range;
} else {
	$start = 1;
	open SAM, "$sam view -H $conf->{'bam'} |";
	while (my $line = <SAM>) {
		if ($line =~ /$accession/) {
			chomp $line;
			my ($sq, $name, $length) = split /\t/, $line;
			$length =~ s/LN://;
			$end = $length;
		}
	}
	close SAM;
}

if (!$start || !$end) {
	print STDERR " Either the positions provided to -p were formatted incorrectly or the accession provided to -a was not found.\n\n";
	exit 1;
}

for (my $p = $start; $p <= $end; $p++) {
	$pos{$p} = 0;
}

# Calculate normalization factors
if ($opt{'R'}) {
	print STDERR " Getting total reads... ";
	$sth = $dbh->prepare('SELECT `total_reads` FROM `libraries` WHERE `library_id` = ?');
	
	my @ids;
	if ($ref_id) {
		push @ids, $ref_id;
	}
	push @ids, $tasi_id;
	
	foreach my $id (@ids) {
		$sth->execute($id);
		my $result = $sth->fetchrow_hashref;
		$rpm{$id} = 1000000 / $result->{'total_reads'};
	}
	print STDERR "done\n";
} else {
	if ($ref_id) {
		$rpm{$ref_id} = 1;	
	}
	$rpm{$tasi_id} = 1;
}

# Get syntasiRNA reads
print STDERR " Normalizing syntasiRNA reads... ";
$sth = $dbh->prepare('SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `library_id` = ?');
$sth->execute($tasi_id);
while (my $row = $sth->fetchrow_hashref) {
	if (length($row->{'seq'}) == 21) {
		$reads{$row->{'sid'}} = $row->{'reads'} * $rpm{$tasi_id};
	}
}
print STDERR "done\n";

if ($ref_id) {
	# Subtract control reads
	print STDERR " Subtracting normalized control reads... ";
	$sth = $dbh->prepare('SELECT * FROM `reads` WHERE `library_id` = ?');
	$sth->execute($ref_id);
	while (my $row = $sth->fetchrow_hashref) {
		if (exists($reads{$row->{'sid'}})) {
			$reads{$row->{'sid'}} -= $row->{'reads'} * $rpm{$ref_id};
			if ($reads{$row->{'sid'}} < 0) {
				$reads{$row->{'sid'}} = 0;
			}
		}
	}
	print STDERR "done\n";
}

# Get mapped tasiRNA
print STDERR " Calculating reads per position... ";
open SAM, "$sam view $conf->{'bam'} '$accession:$start-$end' |";
while (my $line = <SAM>) {
	my @tmp = split /\t/, $line;
	my $sid = $tmp[0];
	my $strand = $tmp[1];
	my $p = $tmp[3];
	if (exists($reads{$tmp[0]})) {
		if ($strand == 0) {
			$pos{$p} += $reads{$sid};
		} else {
			$pos{$p + 2} += $reads{$sid};
		}
	}
}
close SAM;
print STDERR "done\n";

print STDERR " Calculating register table... ";
my $register = 1;
my @ps;
for ($start..$end) {
	push @ps, $_;
}
foreach my $p (@ps) {
	$registers{$register} += $pos{$p};
	$register++;
	if ($register == 22) {
		$register = 1;
	}
}

print OUT "Register table for $accession:$start-$end\n";
print OUT "Register\tReads\n";
for (1..21) {
	print OUT "$_\t$registers{$_}\n";
}
print OUT "\n";
print OUT "Control-subtracted reads per position\n";
print OUT "Position\tReads\n";
foreach my $p (@ps) {
	print OUT $p."\t".$pos{$p}."\n";
}
close OUT;
print STDERR "done\n";
exit;

# Subroutines
sub var_check {
	if ($opt{'h'}) {
		var_error();
	}
	if ($opt{'r'}) {
		$ref_id = $opt{'r'};
	}
	if ($opt{'t'}) {
		$tasi_id = $opt{'t'};
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
	if ($opt{'a'}) {
		$accession = $opt{'a'};
	} else {
		var_error();
	}
	if ($opt{'p'}) {
		$range = $opt{'p'};
	}
}

sub var_error {
	print STDERR "\n  Description:\n";
	print STDERR "  This script will calculate a control-normalized phasing register table for a syntasiRNA construct.\n\n";
	print STDERR "  Usage:\n";
	print STDERR "  syntasiRNA_register_table.pl -r <control library ID> -t <syntasiRNA library ID> -o <output file> -c <conf_file> -s <species> -a <accession>\n";
	print STDERR "\n\n";
	print STDERR "  -r [int]      Library ID for the control (reference) sample. OPTIONAL\n\n";
	print STDERR "  -t [int]      Library ID for the syntasiRNA-transformed sample.\n\n";
	print STDERR "  -a [text]     Accession/name of the construct/reference sequence to analyze.\n\n";
	print STDERR "  -o [file]     Output filename.\n\n";
	print STDERR "  -c [file]     Configuration file.\n\n";
	print STDERR "  -s [text]     Species code.\n\n";
  print STDERR "  -R [boolean]  Normalize reads (Reads/Million). Default, do not normalize.\n\n";
	print STDERR "  -p [text]     Coordinate range to calculate register table in the format start-end. OPTIONAL\n\n";
	print STDERR "\n\n\n";
	exit 1;
}
