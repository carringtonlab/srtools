package CommonFunctions;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(parseListToArray parseFileList formatFasta txt_file_check);


sub parseListToArray {
	my ($flatList) = @_;
	my @temp_list;
	my @final_list;
  
	# put each semicolon separated entry in an array
	if ($flatList =~ /\,/) {
		@temp_list = split (/\,/,$flatList);
	} else {
		push(@temp_list,$flatList);
	}

	# expand range entries
	foreach my $temp_element (@temp_list) {
		if ($temp_element =~ /-/) {
			my ($low_num, $high_num) = split (/\-/, $temp_element);
			if(! ($low_num < $high_num)){
				die("Error: invalid input range\n");
			}
			for my $this_num ($low_num..$high_num){
				push(@final_list,$this_num);
			}
		} else {
			push(@final_list,$temp_element);
    }
 }
 return @final_list;
}

sub parseFileList {
  my ($flatList) = @_;
  my @final_list;

  # put each semicolon separated entry in an array
  if ($flatList =~ /\,/) {
    @final_list = split (/\,/,$flatList);
  } else {
    push(@final_list,$flatList);
  }
 return @final_list;
}

sub formatFasta {
	my $seq = shift;
	my $chars = shift;
	
	my @lines;
	my $iters = int(length($seq) / $chars);
	if (length($seq) % $chars > 0) {
		$iters++;
	}
	my $p = 0;
	for (my $i = 1; $i <= $iters; $i++) {
		push @lines, substr($seq,$p,$chars);
		$p += $chars;
	}
	return join("\n", @lines);
}

sub txt_file_check {
	my $file = shift;
	
	# Does the file exist?
	unless (-e $file) {
		print STDERR "The file $file does not exists.\n\n";
		exit 1;
	}
	
	# What is the file format
	my $format = `file -b $file`;
	
	# Is the file ASCII text?
	unless ($format =~ /^ASCII text/) {
		print STDERR "The file $file format is not ASCII text:\n     $format\n";
		exit 1;
	}
	
	my ($type, $htype);
	if ($format =~ /CR line terminators/) {
		$type = 'mac';
		$htype = 'Classic Mac';
	} elsif ($format =~ /CRLF line terminators/) {
		$type = 'dos';
		$htype = 'Windows';
	} else {
		return 1;
	}
	
	my $cmd = $type.'2unix';
	
	print STDERR "The file $file is formated in $htype format, which is not readable by Perl.\n";
	print STDERR "Do you want me to run $cmd to fix the file? (y/n) ";
	while (my $response = <STDIN>) {
		chomp $response;
		if ($response =~ /^n$/i || $response =~ /^no$/i) {
			exit 1;
		} elsif ($response =~ /^y$/i || $response =~ /^yes$/i) {
			last;
		} else {
			print STDERR "$response is not a valid response. Exiting.\n";
			exit 1;
		}
	}
	
	print STDERR "OK. Do you want me to back up the file $file first? (y/n) ";
	while (my $response = <STDIN>) {
		chomp $response;
		if ($response =~ /^y$/i || $response =~ /^yes$/i) {
			print STDERR "OK. Backing up $file to $file.backup.\n";
			`cp $file $file.backup`;
			last;
		} elsif ($response =~ /^n$/i || $response =~ /^no$/i) {
			last;
		} else {
			print STDERR "$response is not a valid response. Exiting.\n";
			exit 1;
		}
	}
	
	print STDERR "Running $cmd on $file.\n";
	`$cmd $file`;
	
	return 1;
}

1;
