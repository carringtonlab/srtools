#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use POSIX qw( ceil floor);
use DBI;
use Config::Tiny;
use Statistics::R;
use Statistics::Descriptive;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use CommonFunctions qw(parseListToArray);

#########################################################
# Start Variable declarations                           #
#########################################################

my (%opt, $species, $confFile, @inputLibs, @ipLibs, $sth, %rpm, %reads, $outfile, @sizes, %sizes, $thresh, $rthresh, @data);

getopts('s:c:i:I:o:S:t:r:hR', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
our $sam = $Conf->{'PIPELINE'}->{'sam'};
our $hitDB = $Conf->{$species}->{'bam'};
our $featureDB = $conf->{'features'};

if (!$conf || !$sam || !$hitDB || !$featureDB) {
	print STDERR " Malformed configuration settings for $species!\n\n";
	exit 1;
}

foreach my $size (@sizes) {
	$sizes{$size} = 1;
}

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");
if ($opt{'R'}) {
	$sth = $dbh->prepare('SELECT * FROM `libraries` WHERE `library_id` = ?');
	foreach my $id (@inputLibs, @ipLibs) {
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

print STDERR " Getting small RNA reads... ";
$sth = $dbh->prepare('SELECT * FROM `reads` NATURAL JOIN `sequences` WHERE `library_id` = ?');
foreach my $id (@inputLibs, @ipLibs) {
	print STDERR "library $id... ";
	$sth->execute($id);
	while (my $row = $sth->fetchrow_hashref) {
		next if (!exists($sizes{length($row->{'seq'})}));
		if ($opt{'R'}) {
			$row->{'reads'} *= $rpm{$id};
		}
		$reads{$row->{'sid'}}->{$id} = $row->{'reads'};
		$reads{$row->{'sid'}}->{'seq'} = $row->{'seq'};
	}
}
print STDERR "done\n";

print STDERR " Calculating enrichment... ";
my @sids = keys(%reads);
foreach my $sid (@sids) {
	my $input = 0;
	my $ip = 0;
	my $n = 0;
	foreach my $id (@inputLibs) {
		if (exists($reads{$sid}->{$id})) {
			$input += $reads{$sid}->{$id};
			$n++;
		}
#		} else {
#			$reads{$sid}->{$id} = 0;
#		}
	}
	foreach my $id (@ipLibs) {
		if (exists($reads{$sid}->{$id})) {
			$ip += $reads{$sid}->{$id};
			$n++;
		}
#		} else {
#			$reads{$sid}->{$id} = 0;
#		}
	}
	if ($n != scalar(@inputLibs) + scalar(@ipLibs)) {
		delete($reads{$sid});
		next;
	}
	$input /= scalar(@inputLibs);
	$ip /= scalar(@ipLibs);
	if ($ip < $rthresh) {
		delete($reads{$sid});
		next;
	}
	my $ratio = ($ip + 1) / ($input + 1);
	if ($ratio >= $thresh) {
		$reads{$sid}->{'ratio'} = $ratio;	
	} else {
		delete($reads{$sid});
	}
}
print STDERR "done\n";

print STDERR " Annotating enriched small RNA... ";
open (OUT, ">$outfile") or die " Cannot open $outfile: $!\n\n";

my @header = ('sid', 'Seq', '5nt', 'Location');
foreach my $id (@inputLibs) {
	push @header, "input$id";
}
foreach my $id (@ipLibs) {
	push @header, "ip$id";
}
push @header, 'Ratio';
push @header, 'Annotation';
print OUT join("\t", @header)."\n";

open SAM, "$sam view $hitDB |";
while (my $hit = <SAM>) {
	chomp $hit;
	my @tmp = split /\t/, $hit;
	my $sid = $tmp[0];
	next if (!exists($reads{$sid}));
	my $chr = $tmp[2];
	my $start = $tmp[3];
	my $end = $start + length($tmp[9]) - 1;
	my $annotation = annotate($chr, $start, $end);
	if (!$annotation) {
		$annotation = 'intergenic';
	}
	my %hash;
	$hash{'sid'} = $sid;
	$hash{'seq'} = $reads{$sid}->{'seq'};
	$hash{'location'} = $chr.':'.$start.'-'.$end;
	foreach my $id (@inputLibs) {
		$hash{$id} = $reads{$sid}->{$id};
	}
	foreach my $id (@ipLibs) {
		$hash{$id} = $reads{$sid}->{$id};
	}
	$hash{'ratio'} = $reads{$sid}->{'ratio'};
	$hash{'annotation'} = $annotation;
	push @data, \%hash;
}
close SAM;

@data = sort {$b->{'ratio'} <=> $a->{'ratio'}} @data;

foreach my $rna (@data) {
	my @values = (
								$rna->{'sid'},
								$rna->{'seq'},
								substr($rna->{'seq'},0,1),
								$rna->{'location'}
								);
	foreach my $id (@inputLibs) {
		push @values, $rna->{$id};
	}
	foreach my $id (@ipLibs) {
		push @values, $rna->{$id};
	}
	push @values, ($rna->{'ratio'}, $rna->{'annotation'});
	print OUT join("\t", @values)."\n";
}
close OUT;
print STDERR "done\n";

exit;

#########################################################
# End Main body of Program                              #
#########################################################

#########################################################
# Start Subroutines                                     #
#########################################################

sub annotate {
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	my @features;
	open FEAT, "$sam view $featureDB '$chr:$start-$end' |";
	while (my $line = <FEAT>) {
		my @tmp = split /\t/, $line;
		push @features, $tmp[0];
	}
	close FEAT;
	return join(";", @features);
}

sub var_check {
	if ($opt{'h'}) {
		var_error();
	}
	if ($opt{'s'}) {
		$species = $opt{'s'};
	} else {
		var_error();
	}
	if ($opt{'c'}) {
		$confFile = $opt{'c'};
	} else {
		var_error();
	}
	if ($opt{'S'}) {
		@sizes = parseListToArray($opt{'S'});
	} else {
		@sizes = parseListToArray('18-30');
	}
	if ($opt{'o'}) {
		$outfile = $opt{'o'};
	} else {
		var_error();
	}
	if ($opt{'i'}) {
		@inputLibs = parseListToArray($opt{'i'});
	} else {
		var_error();
	}
	if ($opt{'I'}) {
		@ipLibs = parseListToArray($opt{'I'});
	} else {
		var_error();
	}
	if ($opt{'t'}) {
		$thresh = $opt{'t'};
	} else {
		var_error();
	}
	if ($opt{'r'}) {
		$rthresh = $opt{'r'};
	} else {
		$rthresh = 0;
	}
}

sub var_error {
	print STDERR "\n\n";
	print STDERR " This script will generate a table of IP-enriched small RNA.\n";
	print STDERR " Only small RNA that have reads in all four samples are output\n";
	print STDERR " Usage: IPenrichmentTable.pl -s <species> -c <conf file> -i <library ID(s)> -I <library ID(s)> -o <outfile>\n";
	print STDERR " REQUIRED:\n";
	print STDERR " -s     Species\n";
	print STDERR " -c     Configuration file\n";
	print STDERR " -i     Input Library ID list\n";
	print STDERR "			      example: -l 1\n";
	print STDERR "			      example: -l '1-5'\n";
	print STDERR "			      example: -l '1-3,5-10'\n";
	print STDERR " -I     IP Library ID list\n";
	print STDERR "			      example: -l 1\n";
	print STDERR "			      example: -l '1-5'\n";
	print STDERR "			      example: -l '1-3,5-10'\n";
	print STDERR " -o     Output file\n";
	print STDERR " -t     Enrichment threshold.\n";
	print STDERR " OPTIONAL:\n";
	print STDERR " -S     Small RNA size(s). Default = 18-30.\n";
	print STDERR " -R     Normalize reads (Reads/Million). Default, do not normalize.\n\n";
	print STDERR " -r     Read threshold. The minimum average reads in the IP required to output an RNA. If used with -R the threshold will be on RPM.\n\n";
	print STDERR " -h     Print this menu\n";
	print STDERR "\n\n";
	exit 1;
}

#########################################################
# End Subroutines                                       #
#########################################################
