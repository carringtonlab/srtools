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
use GD::Simple;

#########################################################
# Start Variable Declaration                            #
#########################################################

my (%opt, $outfile, $species, $confFile, @libs, %reads, %unapproved, @sizes, %sizes, $sth, $region, $plotStrand, @hits, $color_mode);

my $totalReads = 0;
my $ntwidth = 10;

getopts('o:S:s:C:l:c:R:m:h', \%opt);
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

my ($accession, $range) = split /:/, $region;
my ($rstart, $rend) = split /-/, $range;

my $iwidth = ($rend - $rstart + 1) * $ntwidth;
my $iheight = 100;

GD::Simple->class('GD::SVG');
my $img = GD::Simple->new($iwidth,$iheight);

# Colors
my %color;
if ($color_mode eq 'Phytophthora') {
	$color{18} = $img->colorAllocate(178,224,224);
	$color{19} = $img->colorAllocate(112,205,221);
	$color{20} = $img->colorAllocate(74,151,210);
	$color{21} = $img->colorAllocate(62,91,169);
	$color{22} = $img->colorAllocate(63,83,164);
	$color{23} = $img->colorAllocate(107,82,162);
	$color{24} = $img->colorAllocate(164,81,159);
	$color{25} = $img->colorAllocate(219,68,152);
	$color{26} = $img->colorAllocate(237,28,99);
	$color{27} = $img->colorAllocate(237,34,36);
} elsif ($color_mode eq 'ASRP') {
	$color{18} = $img->colorAllocate(230,231,232);
	$color{19} = $img->colorAllocate(189,215,230);
	$color{20} = $img->colorAllocate(100,198,194);
	$color{21} = $img->colorAllocate(57,83,164);
	$color{22} = $img->colorAllocate(13,129,64);
	$color{23} = $img->colorAllocate(185,82,159);
	$color{24} = $img->colorAllocate(237,31,36);
	$color{25} = $img->colorAllocate(178,36,37);
	$color{26} = $img->colorAllocate(240,73,35);
	$color{27} = $img->colorAllocate(0,0,0);
}

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

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

print STDERR " Getting alignments... ";

open SAM, "$sam view $conf->{'bam'} '$region' |";
while (my $hit = <SAM>) {
	chomp $hit;
	my @fields = split /\t/, $hit;
	my $sid = $fields[0];
	my $strand = $fields[1];
	my $chrom = $fields[2];
	my $start = $fields[3];
	my $seq = $fields[9];
	next if (!exists($reads{$sid}));
	# Convert strand bitwise operator to + or - strand
	$strand = ($strand == 0) ? 'w' : 'c';
	next if ($strand ne $plotStrand);
	# Use length of sequence to determine end position
	my $end = $start + length($seq) - 1;
	next if ($start < $rstart || $end > $rend);
	my $score = $reads{$sid};
	$totalReads += $score;

	my %hash;
	$hash{'reads'} = $score;
	$hash{'length'} = length($seq);
	if ($strand eq 'w') {
		$hash{'x'} = $start - $rstart;
	} elsif ($strand eq 'c') {
		$hash{'x'} = $rend - $end;
	}
	push @hits, \%hash;
}
close SAM;

print STDERR "done\n";

@hits = sort {$a->{'x'} <=> $b->{'x'} || $b->{'length'} <=> $a->{'length'}} @hits;

print STDERR " Plotting... ";

my $current_y = $iheight;

$img->fgcolor(undef);
my $last_x = 0;
foreach my $hit (@hits) {
	if ($hit->{'x'} > $last_x) {
		$current_y = $iheight;
	}
	if ($hit->{'length'} >= 27) {
		$img->bgcolor($color{27});
	} else {
		$img->bgcolor($color{$hit->{'length'}});
	}
	my $x1 = $hit->{'x'} * $ntwidth;
	my $x2 = ($hit->{'x'} + $hit->{'length'}) * $ntwidth;
	my $y1 = $current_y - sprintf('%.2f', (($hit->{'reads'} / $totalReads) * 100));
	my $y2 = $current_y;
	$img->rectangle($x1,$y1,$x2,$y2);
	$last_x = $hit->{'x'} + $hit->{'length'} - 1;
	$current_y = $y1;
}

$img->fgcolor('black');
for (my $x = 0; $x <= $iwidth; $x += 10) {
	$img->moveTo($x,0);
	$img->lineTo($x,$iheight);
}

open (SVG, ">$outfile.svg") or die " Cannot open $outfile.svg for writing: $!\n\n";
print SVG $img->svg();
close SVG;
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
		$outfile = $opt{'o'};
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
	if ($opt{'R'}) {
		$region = $opt{'R'};
	} else {
		var_error();
	}
	if ($opt{'s'}) {
		$plotStrand = $opt{'s'};
	} else {
		var_error();
	}
	if ($opt{'m'}) {
		$color_mode = $opt{'m'};
		unless ($color_mode eq 'ASRP' || $color_mode eq 'Phytophthora') {
			print STDERR " Color mode $color_mode is not valid.\n\n";
			var_error();
		}
	} else {
		$color_mode = 'ASRP';
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
	print STDERR "  This script will read from the sRNAmp database and create a SVG-formatted processing plot for the specified region\n";
	print STDERR "  Usage:\n";
	print STDERR "  plotRegionProcessing.pl -S <species> -c <conf file> -l <library ID(s)> -o <output file prefix> -R <region> -s <strand>\n";
	print STDERR "\n\n";
	print STDERR "  -l   The sample library ID(s)\n";
	print STDERR "			      example: -l 1\n";
	print STDERR "			      example: -l '1-5'\n";
	print STDERR "			      example: -l '1-3,5-10'\n";
	print STDERR "\n";
	print STDERR "  -S   The species.\n";
	print STDERR "\n";
	print STDERR "  -o   The output SVG file prefix.\n";
	print STDERR "\n";
	print STDERR "  -C   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -R   The region to extract. Format is accession:start-end\n";
	print STDERR "\n";
	print STDERR "  -s   The strand to plot. (w or c)\n";
	print STDERR "\n";
	print STDERR "  -c   The sequence sizes to use (DEFAULT = all).\n";
	print STDERR "			      example: -c '20-25'\n";
	print STDERR "			      example: -c '20,21,24'\n";
	print STDERR "			      example: -c '20-22,24-30'\n";
	print STDERR "\n";
	print STDERR "  -m   The color mode. Options are [ASRP, Phytophthora] (DEFAULT = ASRP).\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################

