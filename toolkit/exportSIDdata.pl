#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use POSIX qw( ceil floor);
use DBI;
use Config::Tiny;
use Env qw(HOME);
use lib "$HOME/lib/perl";
use CommonFunctions qw(parseListToArray txt_file_check);

#########################################################
# Start Variable declarations                           #
#########################################################

my (%opt, $file, @list, $outfile, $species, $confFile, %table, $sth, %rpm);

getopts('l:o:s:c:f:hR', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");
if ($opt{'R'}) {
	$sth = $dbh->prepare('SELECT * FROM `libraries` WHERE `library_id` = ?');
	foreach my $id (@list) {
		$sth->execute($id);
		my $row = $sth->fetchrow_hashref;
		if ($row->{'total_reads'} > 0) {
			$rpm{$id} = 1000000 / $row->{'total_reads'};
		} else {
			print STDERR " Total reads for library $id not set in $conf->{'db'}.libraries\n\n";
			exit 1;
		}
	}	
}

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

print STDERR " Reading sids... ";
open (IN, $file) or die " Cannot open $file: $!\n\n";
while (my $sid = <IN>) {
  chomp $sid;
	foreach my $library_id (@list) {
		$table{$sid}->{$library_id} = 0;
	}
}
print STDERR "done\n";


print STDERR " Getting reads... \n";
$sth = $dbh->prepare('SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `sid` = ?');
while (my ($sid, $reads) = each(%table)) {
	$sth->execute($sid);
	while (my $row = $sth->fetchrow_hashref) {
		if (exists($table{$sid}->{$row->{'library_id'}})) {
			if ($opt{'R'}) {
				$row->{'reads'} *= $rpm{$row->{'library_id'}};
			}
			$table{$sid}->{$row->{'library_id'}} = $row->{'reads'};
		}
	}
}
print STDERR "done\n";

print STDERR " Finding genomic positions... ";
open(SAM, "samtools view $conf->{'bam'} |");
while (my $hit = <SAM>) {
	my @fields = split /\t/, $hit;
	my $sid = $fields[0];
	my $strand = $fields[1];
	my $chrom = $fields[2];
	my $start = $fields[3];
	my $seq = $fields[9];
	next if (!exists($table{$sid}));
	# Convert strand bitwise operator to + or - strand
	$strand = ($strand == 0) ? '+' : '-';
	# Use length of sequence to determine end position
	my $end = $start + length($seq) - 1;
	my $string = $chrom.':'.$start.'-'.$end.' '.$strand;
	if (exists($table{$sid}->{'position'})) {
		push @{$table{$sid}->{'position'}}, $string;
	} else {
		@{$table{$sid}->{'position'}} = $string;	
	}
}

close SAM;

open (OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
print OUT "sid\t".join("\t", @list)."\tPositions\n";
while (my ($sid, $data) = each(%table)) {
	print OUT $sid."\t";
	foreach my $library_id (@list) {
		print OUT $data->{$library_id}."\t";
	}
	print OUT join(';', @{$data->{'position'}})."\n";
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
  if ($opt{'f'}) {
    $file = $opt{'f'};
		txt_file_check($file);
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
	print STDERR "  This script will generate a multi-library read count table for specific sequences.\n\n";
	print STDERR "  Usage:\n";
	print STDERR "  extractSIDdata.pl -f <sid file> -l <library_ids> -o <output file> -c <conf_file> -s <species>\n";
	print STDERR "\n\n";
  print STDERR "  -f   The file containing SIDs to query. One SID per line.\n\n";
	print STDERR "  -l   The library ID's to use.\n";
	print STDERR "               example: -l '1-5'\n";
	print STDERR "               example: -l '1,7,11'\n";
	print STDERR "               example: -l '1-5,7,9-11'\n";
	print STDERR "\n";
	print STDERR "  -o   The output filename.\n\n";
	print STDERR "  -c   The configuration file.\n\n";
	print STDERR "  -s   The species.\n\n";
  print STDERR "  -R   Normalize reads (Reads/Million). Default, do not normalize.\n\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
