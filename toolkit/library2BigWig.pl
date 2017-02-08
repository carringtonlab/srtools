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

my (%opt, @list, $outfile, $species, $confFile, %read_table);

getopts('l:o:s:c:r:f:h', \%opt);
var_check();

# Get configuration settings
my $Conf = Config::Tiny->read($confFile);
my $conf = $Conf->{$species};
my $sam = $Conf->{'PIPELINE'}->{'sam'};
my $bam = $conf->{'bam'};

# Connect to the SQLite database
my $dbh = DBI->connect("dbi:SQLite:dbname=$conf->{'db'}","","");

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

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

print STDERR " Getting chromosome sizes... ";
open(OUT, ">$species.chrom.sizes") or die " Cannot open $species.chrom.sizes: $!\n\n";
open SAM, "$sam view -H $bam |";
while (my $line = <SAM>) {
	if (substr($line,0,3) eq '@SQ') {
		chomp $line;
		my ($sq, $chrom, $length) = split /\t/, $line;
		$chrom =~ s/SN://;
		$length =~ s/LN://;
		print OUT $chrom."\t".$length."\n";
	}
	
}
close SAM;
close OUT;
print STDERR "done\n";

print STDERR " Exporting SAM... ";
open(SAMOUT, ">$outfile.sam") or die "Cannot open $outfile.sam: $!\n\n";
open SAM, "$sam view -h $bam |";
while (my $line = <SAM>) {
	if (substr($line,0,1) eq '@') {
		print SAMOUT $line;
	} else {
		my @tmp = split /\t/, $line;
		my $sid = $tmp[0];
		next if (!exists($read_table{$sid}));
		for (my $r = 1; $r <= $read_table{$sid}; $r++) {
			print SAMOUT $line;
		}
	}
}
close SAM;
close SAMOUT;
print STDERR "done\n";

print STDERR " Converting SAM to sorted BAM... ";
`$sam view -@ 5 -bS $outfile.sam -o $outfile.bam`;
`$sam sort -@ 5 -m 10G $outfile.bam $outfile.sorted`;
`$sam index $outfile.sorted.bam`;
print STDERR "done\n";

print STDERR " Converting BAM to BEDGraph... ";
`bedtools genomecov -bga -ibam $outfile.sorted.bam > $outfile.sorted.bedGraph`;
print STDERR "done\n";

print STDERR " Log-transforming BEDGraph... ";
`awk '{print \$1"\t"\$2"\t"\$3"\t"log(\$4 + 1)/log(10)}' $outfile.sorted.bedGraph > $outfile.sorted.log.bedGraph`;
print STDERR "done\n";

#print STDERR " Converting Log-BEDGraph to Log-BigWig... ";
#`bedGraphToBigWig $outfile.sorted.log.bedGraph $species.chrom.sizes $outfile.log.bw`;
#print STDERR "done\n";

print STDERR " Cleaning up... ";
unlink("$species.chrom.sizes");
unlink("$outfile.sam");
unlink("$outfile.bam");
#unlink("$outfile.sorted.bam");
#unlink("$outfile.sorted.bam.bai");
#unlink("$outfile.sorted.bedGraph");
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
	print STDERR "  library2BigWig.pl -l <library_ids> -o <output file prefix> -c <conf_file> -s <species>\n";
	print STDERR "\n\n";
	print STDERR "  -l   The library ID's to use.\n";
	print STDERR "               example: -l '1-5'\n";
	print STDERR "               example: -l '1,7,11'\n";
	print STDERR "               example: -l '1-5,7,9-11'\n";
	print STDERR "\n";
	print STDERR "  -o   The output file prefix.\n";
	print STDERR "\n";
	print STDERR "  -c   The configuration file.\n";
	print STDERR "\n";
	print STDERR "  -s   The species.\n";
	print STDERR "\n";
#	print STDERR "  -r   The small RNA size range.\n";
	print STDERR "\n\n\n";
	exit 0;
}

#########################################################
# End of Varriable error Subroutine "var_error"         #
#########################################################
