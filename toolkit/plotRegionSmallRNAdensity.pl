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

my (%opt, $species, $confFile, @libs, $sth, %reads, $prefix, $accession, $range, $additional_features, @features, @x, @y1, @y2, $cmd);

getopts('s:c:l:p:a:r:f:hi', \%opt);
var_check();

my $cwd = getcwd;
my @colors = ('red','green','blue','yellow','orange','purple');
my $pdf = "$prefix.pdf";
#my $data_file = "$prefix.txt";

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};
my $hitDB = $Conf->{$species}->{'bam'};

my ($start, $end) = split /-/, $range;
if (!$start || !$end) {
	var_error();
}

if ($additional_features ne 'none') {
	@features = split /,/, $additional_features;
}

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

print STDERR " Getting small RNA reads... ";
my $total_reads = 0;
$sth = $dbh->prepare('SELECT * FROM `reads` WHERE `library_id` = ?');
foreach my $id (@libs) {
	print STDERR "library $id... ";
	$sth->execute($id);
	while (my $row = $sth->fetchrow_hashref) {
		if (exists($reads{$row->{'sid'}})) {
			$reads{$row->{'sid'}} += $row->{'reads'};
		} else {
			$reads{$row->{'sid'}} = $row->{'reads'};
		}
		$total_reads += $row->{'reads'};
	}
}
print STDERR "$total_reads total reads\n";

print STDERR " Getting small RNA genome hits... ";
for (my $x = $start; $x <= $end; $x++) {
	push @x, $x;
	push @y1, 0;
	push @y2, 0;
}
open SAM, "$sam view $hitDB '$accession:$range' |";
while (my $hit = <SAM>) {
	chomp $hit;
	my @tmp = split /\t/, $hit;
	next if (!exists($reads{$tmp[0]}));
	for (my $p = $tmp[3]; $p < $tmp[3] + length($tmp[9]); $p++) {
		my $i = $p - $start;
		next if ($i >= scalar(@x));
		if ($opt{'i'}) {
			$y1[$i] += $reads{$tmp[0]};
		} else {
			if ($tmp[1] == 0) {
				$y1[$i] += $reads{$tmp[0]};
			} elsif ($tmp[1] == 16) {
				$y2[$i] -= $reads{$tmp[0]};
			}
		}
	}
	#if ($tmp[1] == 0) {
	#	push @x, $tmp[3];
	#	push @y, $reads{$tmp[0]};
	#} elsif ($tmp[1] == 16) {
	#	push @x, $tmp[3] + length($tmp[9]) - 1;
	#	push @y, -$reads{$tmp[0]};
	#}
}
close SAM;
print STDERR "done\n";

print STDERR " Generating histogram... ";
my $R = Statistics::R->new();
$R->start();

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data((@y1,@y2));
my $y_max = $stat->max();
my $y_min = $stat->min();

if (!$opt{'i'}) {
	if ((abs($y_max) > abs($y_min)) || (abs($y_max) == abs($y_min) && $y_max > 0)) {
		$y_min = -$y_max;
	} elsif ((abs($y_min) > abs($y_max)) || (abs($y_max) == abs($y_min) && $y_max < 0)) {
		$y_max = abs($y_min);
	} else {
		print 'Something is wrong here:\n';
		print 'Y-max = '.$y_max."\n";
		print 'Y-max = '.$y_min."\n";
		exit 1;
	}
}

$cmd  = "x<-cbind(".join(",", @x).")\n";
$cmd .= "y1<-cbind(".join(",", @y1).")\n";
if (!$opt{'i'}) {
	$cmd .= "y2<-cbind(".join(",", @y2).")\n";
}
$cmd .= "pdf(file=\"$cwd/$pdf\",width=10,height=10,family=\"Helvetica\",paper=\"special\")\n";
$cmd .= "plot(x,y1,type='n',xaxt='n',yaxt='n',xlim=c($start,$end),ylim=c($y_min,$y_max),bty='n',xlab=NA,ylab=NA)\n";
$cmd .= "axis(side=1,pos=c($y_min,$start),lwd=1)\n";
$cmd .= "axis(side=2,pos=c($start,$y_min),lwd=1)\n";
$cmd .= "points(x=x,y=y1,type='h')\n";
if (!$opt{'i'}) {
	$cmd .= "points(x=x,y=y2,type='h')\n";
}
$cmd .= "lines(x=c($start,$end),y=c(0,0),lwd=0.5)\n";
#$cmd .= "lines(x=c($start,$end),y=c($y_min,$y_min),lwd=1)\n";


if (@features) {
	foreach my $feature (@features) {
		my ($x1, $x2) = split /-/, $feature;
		my $y1 = $y_min;
		my $y2 = $y_min + int($y_max * 0.1);
		my $color = shift(@colors);
		$cmd .= "rect(xleft=c($x1),xright=c($x2),ybottom=c($y1),ytop=c($y2),col=\'$color\')\n";
		push @colors, $color;
	}
}
$cmd .= "dev.off()\n";
$R->run($cmd);
print STDERR "done\n";
my $errors = $R->error();
if ($errors =~ /^\s*$/) {
	print STDERR " Plotting completed without any errors\n";
} else {
	print STDERR "$errors\n";
}
$R->stop();
print STDERR "done\n";

exit;

#########################################################
# End Main body of Program                              #
#########################################################

#########################################################
# Start Subroutines                                     #
#########################################################

sub var_check {
	if ($opt{'h'}) {
		var_error();
	}
	if ($opt{'s'}) {
		$species = $opt{'s'};
	} else {
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
	if ($opt{'p'}) {
		$prefix = $opt{'p'};
	} else {
		var_error();
	}
	if ($opt{'a'}) {
		$accession = $opt{'a'};
	} else {
		var_error();
	}
	if ($opt{'r'}) {
		$range = $opt{'r'};
	} else {
		var_error();
	}
	if ($opt{'f'}) {
		$additional_features = $opt{'f'};
	} else {
		$additional_features = 'none';
	}
}

sub var_error {
	print STDERR "\n\n";
	print STDERR " This script will generate a per nucleotide read density plot for a given region in any database\n";
	print STDERR " Both a data file and pdf image will be produced\n";
	print STDERR " Usage: plotRegionSmallRNAdensity.pl -s <species> -c <conf file> -l <library ID(s)> -p <outfile prefix> -a <region accession> -r <coordinate range>\n";
	print STDERR " REQUIRED:\n";
	print STDERR " -s     Species\n";
	print STDERR " -c     Configuration file\n";
	print STDERR " -l     Library ID list\n";
	print STDERR "			      example: -l 1\n";
	print STDERR "			      example: -l '1-5'\n";
	print STDERR "			      example: -l '1-3,5-10'\n";
	print STDERR " -p     Prefix for output files\n";
	print STDERR " -a     Chromosome or scaffold accession number (e.g. Scaffold_15 or NC_003070)\n";
	print STDERR " -r     Coordinate range (e.g. 1000-2000)\n";
	print STDERR " OPTIONAL:\n";
	print STDERR " -f     Additional features to be plotted. Supply as a comma separated list (e.g. '1000-2000,3000-4000')\n";
	print STDERR " -i     Ignore strand information.\n\n";
	print STDERR " -h     Print this menu\n";
	print STDERR "\n\n";
	exit 1;
}

#########################################################
# End Subroutines                                       #
#########################################################
