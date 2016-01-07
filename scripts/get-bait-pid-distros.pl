#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Data::Dumper;
use Getopt::Long;

# Lenght of baits
my $bait_length = 120;

# Minimum number of nucleotides contig must share with bait
my $min_bait_overlap = 60;

# Check PATH for required executables
my $R = check_path_for_exec("R");
my $blat = check_path_for_exec("blat");

# Parse command line options
my $pid_cutoff;
my $bait_filename;
GetOptions(
	'b|baits=s' => \$bait_filename,
	'p|pid-cutoff=f' => \$pid_cutoff,
);

# Error check input
die "You must specify the filename of the file containing utilized bait sequences.\n" if (!defined($bait_filename));
die "Could not locate file '$bait_filename', perhaps you made a typo?\n" if (!defined($bait_filename));

# Fetch filenames of input alignments
my @input = @ARGV;
die "You must specify the filename(s) of input alignment(s).\n" if (!@input);

# Parse bait sequences into memory
my %baits = parse_fasta($bait_filename);

# Calculate bait pairwise identity distribution for each alignment
my %pids;
my %final_aligns;
foreach my $file (@input) {
	print "Working on $file...\n";

	# Get name of target sequence
	(my $target = $file) =~ s/\.fasta$//;

	# Parse alignment into memory
	#my %align = parse_fasta($file);
	my %full_align = parse_fasta($file);
	$final_aligns{$target} = \%full_align;

	# Remove introns inserted while generating consensus
	#remove_introns($target, \%align);
	my %align = remove_introns($target, %full_align);

	# Extract target sequence
	my $target_seq = $align{$target};

	# Temporary check
	die "Error removing introns from '$target'.\n" if ($target_seq =~ /-+/);

	# Write target sequence to a separate, temporary file
	open(my $tmp, ">", "$target.tmp.fasta");
	print {$tmp} ">$target\n";
	print {$tmp} "$target_seq\n";
	close($tmp);

	# Run blat against alignment
	my $return = system("$blat '$bait_filename' '$target.tmp.fasta' -t=dna -q=dna -noHead '$file.psl' >/dev/null");
	die "An error occurred while running blat: '$return'." if ($return);

	# Parse blat output
	my %hits = parse_blat_output("$file.psl");

	# Calculate percent identities for each bait
	foreach my $hit (keys %hits) {

		# Extract corresponding sequence for bait
		my $bait_seq = $baits{$hit};

		# Extract corresponding contig sequences
		my %contig_seqs;
		foreach my $contig (keys %align) {
			next if ($contig eq $target);
			$contig_seqs{$contig} = substr($align{$contig}, $hits{$hit}->{START}, $hits{$hit}->{END} - $hits{$hit}->{START});

			# Temporary check
			die "Unexpected bait length: $hit (",length($contig_seqs{$contig}),").\n" if (length($contig_seqs{$contig}) != $bait_length);
		}

		# Calculate the percent identity if coverage is sufficient
		foreach my $contig (keys %contig_seqs) {
			my $seq = $contig_seqs{$contig};

			# Check that we have sufficient coverage
			(my $ungapped_seq = $seq) =~ s/-//g;
			if (length($ungapped_seq) >= $min_bait_overlap) {
				my $pid = get_pid($seq, $bait_seq);

				# Remove contig from final alignment if it doesn't meet similarity threshold
				if (defined($pid_cutoff) && $pid < $pid_cutoff) {
					delete($final_aligns{$target}->{$contig});
				}

				# Add to final array
				push(@{$pids{$target}}, $pid);
			}
		}
	}

	# Clean up temporary files 
	unlink("$file.psl");
	unlink("$target.tmp.fasta");
}

# Plot the final distributions
plot_distros(\%pids);

# Generate a file containing the new consensus sequences if a PID cutoff was specified
if (defined($pid_cutoff)) {

	# Generate consensus for each trget
	my %consensuses = consense_alignments(\%final_aligns);

	# Output consensus sequence to fasta
	open(my $out, ">", "filtered-$pid_cutoff.con.fasta");
	foreach my $consensus (keys %consensuses) {
		print {$out} ">$consensus\n";
		print {$out} "$consensuses{$consensus}\n";
	}
	close($out);
}

sub plot_distros {
	my $pids = shift;

	# Create R script which will create the plots
	my $script_name = "bait-pid-distro-plot.r";
	open(my $R_script, ">", $script_name) or die "Could not open 'bait-pid-distro-plot.r': $!.\n";
	#print {$R_script} "library(diptest);\n";
	print {$R_script} "pdf(file='bait-pid-distros.pdf');\n";

	# Extract data for each target
	my $count = 0;
	foreach my $target (sort {$a cmp $b } keys %{$pids}) {

		# Format pids so they can be written into the R script
		my $data = join(', ', @{$pids->{$target}});

		# Add data to R script
		print {$R_script} "pids$count = c($data);\n";
		print {$R_script} "hist(pids$count, freq=T, main='$target', xlab='Percent Identity');\n\n";
		#print {$R_script} "plot(density$count, main='$target', xlab='Percent Identity');\n\n";
		#print {$R_script} "density$count = density(pids$count);\n";
		#print {$R_script} "plot(density$count, main='$target');\n\n";
		#print {$R_script} "dip.test(pids$count, B=2000);\n";

		$count++;
	}
	print {$R_script} "dev.off();\n";

	# Run the R script we generated
	my $return = system("cat $script_name | $R --no-save >/dev/null");
	#my $return = system("cat $script_name | $R --no-save");
	die "An error occurred while running R: '$return'." if ($return);

	# Cleanup
	unlink($script_name);

	return;
}

sub remove_introns {
	#my ($target, $align) = @_;
	my ($target, %align) = @_;

	# Retrieve target sequence from alignment
	#my $target_seq = $align->{$target};
	my $target_seq = $align{$target};

	# Get string indices of intronic sites
	my @ends;
	my @starts;
	while ($target_seq =~ /-+/g) {
		push(@ends, $+[0]);
		push(@starts, $-[0]);
	}

	# Remove intronic sites
	foreach my $index (reverse(0 .. scalar(@starts) - 1)) {
		my $end = $ends[$index];
		my $start = $starts[$index];

		# Delete from each sequence in alignment
		#foreach my $contig (keys %{$align}) {
		foreach my $contig (keys %align) {
			#substr($align->{$contig}, $start, $end - $start) = '';
			substr($align{$contig}, $start, $end - $start) = '';
		}
	}

	return %align;
}

sub get_pid {
	my ($seq1, $seq2) = @_;

	print $seq1,"\n";
	print $seq2,"\n";

	# Final pid
	my $pid;

	# Iterate through each character and check for equality
	my $total = 0;
	my $match = 0;
	foreach my $index (0 .. length($seq1) - 1) {
		my $char1 = chop($seq1);
		my $char2 = chop($seq2);

		# Skip ambiguous sites and gaps
		next if ($char1 !~ /A|T|G|C/ || $char2 !~ /A|T|G|C/);

		# Increment counter if sites match
		if ($char1 eq $char2) {
			$match++;
		}
		$total++;
	}
	$pid = $match / $total * 100;
	print $pid,"\n\n";

	return $pid;
}

sub parse_blat_output {
	my $file = shift;

	# Stores final hits
	my %hits;

	# Open blat output and parse hits
	open(my $blat_out, "<", $file);
	while (my $line = <$blat_out>) {
		chomp($line);

		my @line = split(/\s+/, $line);

		# Associate blat headers with proper value
		my ($match, $mismatch, $rep_matches, $n_count, $q_num_inserts, $q_base_inserts, 
			$t_num_inserts, $t_base_inserts, $strand, $q_name, $q_size, $q_start, $q_end,
			$t_name, $t_size, $t_start, $t_end, $block_count, $block_sizes, $q_starts, $t_starts) = @line; 

		# Skip bait hits which aren't for the correct target
		next if ($t_name !~ /\Q$q_name\E/);

		# Baits should be an exact match
		next if ($match != $bait_length);

		# Extract start and end indices
		$hits{$t_name}->{START} = $q_start;
		$hits{$t_name}->{END} = $q_end;
	}
	close($blat_out);

	return %hits;
}

sub consense_alignments {
	my $alignment = shift;

	my %consensuses;
	foreach my $target (keys %{$alignment}) {
		my %target_alignment = %{$alignment->{$target}};

		# Determine length of alignment
		my $align_length = length((values %target_alignment)[0]);
		die "$target\n" if (!defined($align_length));

		# Loop through each site in alignment
		foreach my $index (0 .. $align_length - 1) {

			# Loop through each sequence in alignment
			my %site = ('A' => 0, 'T', => 0, 'C' => 0, 'G' => 0);
			foreach my $sequence (values %target_alignment) {
				my $char = substr($sequence, $index, 1);
				die "$target\n" if (!defined($char));
				$site{$char}++;
			}
			
			# Calculate and append consensus base
			my $consensus = consense_site(\%site);
			$consensuses{$target} .= $consensus;
		}
	}

	return %consensuses;
}

sub consense_site {
	my $site = shift;

	# Sort bases by descending order of frequency
	my @sorted_bases = sort { $site->{$b} <=> $site->{$a} } keys %{$site};

	# Return most frequent unambiguous base if possible
	foreach my $base (@sorted_bases) {
		if ($site->{$base} && $base =~ /[ACTG]/) {
			return $base;
		}
	}

	# Return most frequent ambiguity code otherwise
	foreach my $base (@sorted_bases) {
		if ($site->{$base} && $base =~ /[A-Z]/) {
			return $base;
		}
	}

	# Site was deleted by pid cutoff
	if ($site->{"-"}) {
		return '';
	}

	# Should never happen
	die "Error consensing site:\n",Dumper($site),"\n";
}

sub check_path_for_exec {
	my ($exec, $continue) = @_;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec && -x $dir.$exec && !-d $dir.$exec);
	}

	die "Could not locate: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path) && !defined($continue));
	return $exec_path;
}

sub parse_fasta {
	my $filename = shift;

	my $taxon;
	my %align;
	open(my $alignment_file, '<', $filename) 
		or die "Could not open '$filename': $!\n";

	while (my $line = <$alignment_file>) {
		$line =~ s/^\s+|\s+$//g;

		# Taxon name
		if ($line =~ /^>(\S+)/) {
			$taxon = $1;
		}
		else {
			# Taxon sequence
			#$taxon =~ s/-/_/g;
			$align{$taxon} .= uc($line);
		}
	}
	close($alignment_file);
	
	return %align;
}
