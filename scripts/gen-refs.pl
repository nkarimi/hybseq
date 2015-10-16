#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Data::Dumper;
use Getopt::Long;
use Time::HiRes qw(usleep);
use Fcntl qw(:flock SEEK_END);

# Not yet implemented
my $max_overhang_length = 300;
#my $overhang_length = 0;
my $target_perc_cov_threshold = 0;
my $contig_perc_cov_threshold = 0;

my $free_cpus = get_free_cpus();

my $cap3 = check_path_for_exec("cap3");
my $blat = check_path_for_exec("blat");
my $mafft = check_path_for_exec("mafft");

my $targets;
GetOptions(
	't|targets=s' => \$targets,
);

my $input = shift;

die "You must specify an assembly as input.\n" if (!defined($input));
die "You must specify a fasta file containing target sequences (-t).\n" if (!defined($targets));
die "You must you have specified an incorrect argument.\n" if (@ARGV);

# Output file names
#my $output_fasta = $input.".hits";
my $blat_out_name = $input.".psl";
my $blast_out_name = $input.".blast";
#(my $coverage_output = $input) =~ s/\.fa(sta)?/.cov/i;
(my $consensus_output = $input) =~ s/\.fa(sta)?/.con.fasta/i;

# Run rough preliminary blat to determine which targets and contigs align
print "Running blat against targets and contigs...\n";
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead $blat_out_name");
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead -fine $blat_out_name") if (!-e $blat_out_name);
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead -fine $blat_out_name");
my $return = system("$blat $input $targets -t=dna -q=dna -noHead -fine -repMatch=1000000 $blat_out_name");
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead -fine -repMatch=1000000 $blat_out_name") if (!-e $blat_out_name);
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead -fine -repMatch=1000000 -stepSize=2 -minMatch=1 $blat_out_name");
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead -fine -repMatch=1000000 -tileSize=6 -stepSize=1 -minMatch=1 $blat_out_name");
#my $return = system("$blat $input $targets -t=dna -q=dna -noHead -tileSize=6 -minMatch=1 -fine $blat_out_name");
die "Error running blat: '$return'.\n" if ($return);

# Parse blat output
my %targets;
open(my $blat_out, "<", $blat_out_name);
while (my $line = <$blat_out>) {
	chomp($line);

	my @line = split(/\s+/, $line);

	# Blat headers
	my ($match, $mismatch, $rep_matches, $n_count, $q_num_inserts, $q_base_inserts, 
	    $t_num_inserts, $t_base_inserts, $strand, $q_name, $q_size, $q_start, $q_end,
		$t_name, $t_size, $t_start, $t_end, $block_count, $block_sizes, $q_starts, $t_starts) = @line; 

	# Check that contig meets target coverage threshold
	if (($match + $mismatch) / $q_size >= $target_perc_cov_threshold && 
		($match + $mismatch) / $q_size >= $contig_perc_cov_threshold) {
		#push(@{$targets{$q_name}}, {'NAME' => $t_name, 'STRAND' => $strand, 'T_STARTS' => $t_starts, 'Q_STARTS' => $q_starts, 'SIZES' => $block_sizes});
		if (!exists($targets{$q_name}->{$t_name})) {
			$targets{$q_name}->{$t_name} = {'NAME' => $t_name, 'STRAND' => $strand, 'T_STARTS' => $t_starts, 'Q_STARTS' => $q_starts, 'SIZES' => $block_sizes};
		}
		else {
			my $hit = $targets{$q_name}->{$t_name};

			my $current_size = $hit->{SIZES};
			my @current_size = split(",", $current_size);

			$current_size = 0;
			foreach my $size (@current_size) {
				$current_size += $size;
			}

			my $new_size = $block_sizes;
			my @new_size = split(",", $new_size);

			$new_size = 0;
			foreach my $size (@new_size) {
				$new_size += $size;
			}

			if ($current_size > $new_size) {
				$targets{$q_name}->{$t_name} = {'NAME' => $t_name, 'STRAND' => $strand, 'T_STARTS' => $t_starts, 'Q_STARTS' => $q_starts, 'SIZES' => $block_sizes};
			}
		}
	}
}
close($blat_out);

# Load contig and target sequences
my %contig_seqs = parse_fasta($input);
my %target_seqs = parse_fasta($targets);

# Create consensus sequence for each target
my %consensuses;
my %coverage_maps;
foreach my $target (sort { $a cmp $b } keys %targets) {
	#$target = "10100|Gorai.010G161400.1|";
	#$target = "471|Gorai.006G256000.2|";
	#my @hits = @{$targets{$target}};
	my @hits = values %{$targets{$target}};

#	foreach my $hit (@hits) {
#		print "$hit->{'NAME'}\n";
#	}

	# Holds bases/introns in contigs which map to specific sites in target
	my %introns;
	#my %overhangs;
	my %alignment;

	my %prematch;
	my %postmatch;

	# Add target bases to alignment
	my $target_seq = $target_seqs{$target};
	my $target_length = length($target_seq);
	foreach my $index (0 .. length($target_seq) - 1) {
		my $base = substr($target_seq, $index, 1);
		$alignment{$index}->{"BASES"} = {};
		$alignment{$index}->{"ID"} = $index;
		$alignment{$index}->{"TARGET"} = $base;
		$alignment{$index}->{"INTRONS"} = {};
		$alignment{$index}->{"PART_INT_END"} = {};
		$alignment{$index}->{"PART_INT_START"} = {};
	}

	print "Building consensus sequence for '$target'...\n";

	my %starts;

#	open(my $hit_fasta, ">", "hits.fasta");
#	print {$hit_fasta} ">target\n";
#	print {$hit_fasta} "$target_seq\n";

	# Extract homologous sequence from each contig
	foreach my $hit (@hits) {
		print "  ",$hit->{'NAME'},"\n";

		my $seq = $contig_seqs{$hit->{'NAME'}};
		my $contig_length = length($seq);

		print "    Contig length: $contig_length\n";
		print "    Target length: $target_length\n";

		my @sizes = split(",", $hit->{'SIZES'});
		my @t_starts = split(",", $hit->{'T_STARTS'});
		my @q_starts = split(",", $hit->{'Q_STARTS'});

		$starts{$hit->{'NAME'}} = $q_starts[0];

		# Output strandedness
		if ($hit->{'STRAND'} eq "-") {
			print "    Contig will be reverse complemented (- hit).\n\n";

			$seq = rev_comp($seq);
			foreach my $index (0 .. $#t_starts) {
				$t_starts[$index] = $contig_length - $sizes[$index] - $t_starts[$index];
				$q_starts[$index] = $target_length - $sizes[$index] - $q_starts[$index];
			}
			@sizes = reverse(@sizes);
			@t_starts = reverse(@t_starts);
			@q_starts = reverse(@q_starts);
		}
		else {
			print "    Contig will not be reverse complemented (+ hit).\n\n";
		}

#		print {$hit_fasta} ">$hit->{'NAME'}\n";
#		print {$hit_fasta} $seq,"\n";
#		print {$hit_fasta} replace_ambiguities($seq),"\n";

		my $prematch = substr($seq, 0, $t_starts[0] + 1);
		my $postmatch = substr($seq, $t_starts[$#t_starts], length($seq) - $t_starts[$#t_starts]);

#		# Trim overhangs # end
#		if (length($prematch) > $max_overhang_length) {
#			#print $prematch," (",length($prematch),")\n";
#			$prematch = substr($prematch, -1 * $max_overhang_length);
#			#print $prematch," (",length($prematch),")\n";
#		}
#
#		# Trim overhangs # start
#		if (length($postmatch) > $max_overhang_length) {
#			print $postmatch," (",length($postmatch),")\n";
#			$postmatch = substr($postmatch, 0, $max_overhang_length);
#			print $postmatch," (",length($postmatch),")\n";
#		}

	#	$overhangs{$q_starts[0]}->{$prematch}++;
	#	$overhangs{$q_starts[$#q_starts] + $sizes[$#q_starts] - 1}->{$postmatch}++;
#		$prematch{$q_starts[0]}->{$prematch}++;
#		$postmatch{$q_starts[$#q_starts] + $sizes[$#q_starts] - 1}->{$postmatch}++;

#		$alignment{$q_starts[0]}->{"PART_INT_END"}->{$prematch}++;
#		$alignment{$q_starts[$#q_starts] + $sizes[$#q_starts] - 1}->{"PART_INT_START"}->{$postmatch}++;

		if ($q_starts[0] == 0) {
			$alignment{$q_starts[0]}->{"PART_INT_END"}->{$prematch}++;
		}

		if ($q_starts[$#q_starts] + $sizes[$#q_starts] - 1 == $target_length - 1) {
			$alignment{$q_starts[$#q_starts] + $sizes[$#q_starts] - 1}->{"PART_INT_START"}->{$postmatch}++;
		}

#		if (length($prematch) > 50) {
#			$alignment{$q_starts[0]}->{"PART_INT_END"}->{$prematch}++;
#		}
#		if (length($postmatch) > 50) {
#			$alignment{$q_starts[$#q_starts] + $sizes[$#q_starts] - 1}->{"PART_INT_START"}->{$postmatch}++;
#			#$alignment{$q_starts[$#q_starts] + $sizes[$#q_starts] - 1}->{"PART_INT_START"}->{$postmatch} = $hit->{'NAME'};
#		}

		# Iterate through aligned blocks of target and contig
		foreach my $index (0 .. $#t_starts) {
			my $size = $sizes[$index];
			my $q_start = $q_starts[$index];
			my $t_start = $t_starts[$index];

			# Where the intron starts 
			my $intron_start = $t_start + $size - 1;

			# If this is the last (or only) block, we set the next block 
			# start to a value which makes the intron length 0
			my $intron_length = 0;
			if ($index + 1 < scalar(@t_starts)) {
				$intron_length = $t_starts[$index + 1] - $intron_start - 1;
			}

			print "    match_start : $t_start ($q_start)\n";
			print "    match_length : $size\n";
			print "    match_end : ",($t_start + $size - 1),"\n";

			if ($intron_length > 0) {
				print "    intron_start : $intron_start (",($q_start + $intron_start - $t_start),"), intron_length : $intron_length\n\n";
			}
			else {
				print "\n\n";
			}

			# Place aligned characters into %alignment
			foreach my $offset (0 .. $size - 1) {
				my $t_index = $t_start + $offset;
				my $q_index = $q_start + $offset;

				my $char = substr($seq, $t_index, 1);
				die "Could not get character at $t_index\n" if (!$char);
				$alignment{$q_index}->{"BASES"}{$char}++;
			}

			# Place introns into %alignment
			if ($intron_length > 0) {
				my $intron = substr($seq, $intron_start + 1, $intron_length);
				die "Could not get intron at ",($intron_start + 1),"\n" if (!$intron);
				$alignment{$q_start + $intron_start - $t_start}->{"INTRONS"}{$intron}++;
				$introns{$q_start + $intron_start - $t_start}->{$intron}++;
			}
		}
	}
	#die;

#	use Data::Dumper;
	print $target,"\n";

	#my ($consensus, $coverage) = consense_alignment(\%alignment);
	#my $consensus = consense_alignment(\%alignment);
	my ($consensus, $distinct_seqs_per_intron) = consense_alignment(\%alignment);
	$consensus =~ s/^n+|n+$//g; # remove extra N's at start and end of alignment

	$consensuses{$target} = $consensus;	
	#$coverage_maps{$target} = $coverage;	
	
#	if ($distinct_seqs_per_intron > 1) {
#		$consensuses{'multi'}->{$target} = $consensus;	
#	}
#	else {
#		$consensuses{'single'}->{$target} = $consensus;	
#	}

#	die;
}

# Output consensus sequences
open(my $output, ">", $consensus_output);
foreach my $target (sort { $a cmp $b } keys %consensuses) {
	print {$output} ">$target\n";
	print {$output} "$consensuses{$target}\n";
}
close($output);

#open($output, ">", "single-copy.fasta");
#foreach my $target (sort { $a cmp $b } keys %{$consensuses{'single'}}) {
#	print {$output} ">$target\n";
#	print {$output} "$consensuses{'single'}->{$target}\n";
#}
#close($output);
#
#open($output, ">", "multi-copy.fasta");
#foreach my $target (sort { $a cmp $b } keys %{$consensuses{'multi'}}) {
#	print {$output} ">$target\n";
#	print {$output} "$consensuses{'multi'}->{$target}\n";
#}
#close($output);

## Output coverage maps
#open($output, ">", $coverage_output);
#foreach my $target (sort { $a cmp $b } keys %coverage_maps) {
#	print {$output} ">$target\n";
#	print {$output} "$coverage_maps{$target}\n\n";
#}
#close($output);

sub consense_alignment {
	my $align = shift;

	my $exon;
	my %introns;
	my %intron_seqs;

	# Consense exonic sequence first
	my %align = %{$align};
	foreach my $site (sort { $a <=> $b } keys %align) {

		my %site = %{$align{$site}};
		my %bases = %{$site{"BASES"}};
		if (scalar(keys %bases) == 0) {
			$exon .= "N";
		}
		else {
			consense_bases(\%site, \$exon);

			my %intron_sequences = (%{$site{INTRONS}}, %{$site{PART_INT_END}}, %{$site{PART_INT_START}});

			if (scalar(keys %intron_sequences)) {
				$introns{$site} = \%site;

#				foreach my $intron (keys %intron_sequences) {
#					$intron_seqs{$intron} = $site;
#				}
			}
		}
	}

	# Check that we don't have the same intron being counted multiple times due to indels
	group_all_introns(\%introns);
	#die;

	# Consense introns in reverse order to keep sites accurate
	my @intron_counts;
	my $original_length = length($exon);
	foreach my $site (sort { $b <=> $a } keys %introns) {
		#my $intron = consense_introns($introns{$site}, \$exon);
		my ($intron, $num_complete_introns) = consense_introns($introns{$site}, \$exon);

		if ($intron) {
			if ($site == 0) {
				substr($exon, $site, 0, $intron);
			}
			else {
				substr($exon, $site + 1, 0, $intron);
			}

			if ($site != 0 && $site != $original_length - 1) {
				push(@intron_counts, $num_complete_introns);
			}
		}

		#substr($exon, $site, 0) = $intron;
		#die if $site == 849;
	}
	my $distinct_seqs_per_intron = get_median(\@intron_counts);

	print "",scalar(@intron_counts)," intron(s) in target, averaging $distinct_seqs_per_intron distinct sequence(s) per site.\n";

#	my $coverage;
#	my $alignment;
#
#	my %align = %{$align};
#	SITE: foreach my $site (sort { $a <=> $b } keys %align) {
#
#		my %site = %{$align{$site}};
#		my %bases = %{$site{"BASES"}};
#
#		if (scalar(keys %bases) == 0) {
#			$coverage .= "0 ";
#			#j$alignment .= "?";
#			$alignment .= "N";
#
##			print Dumper(\%site),"\n";
##			print $site,"\n";
##			die;
#		}
#		else {
#			# An intron at site 0 isn't really an intron so it needs to be handled differently
#			if ($site == 0) {
#				consense_introns(\%site, \$coverage, \$alignment);
#				consense_bases(\%site, \$coverage, \$alignment);
#			}
#			else {
#
##				if (scalar(keys %{$site{INTRONS}}) || scalar(keys %{$site{PART_INT_END}}) || scalar(keys %{$site{PART_INT_START}})) {
##					print Dumper(\%site),"\n";
##					print $site,"\n\n";
##				}
#				consense_bases(\%site, \$coverage, \$alignment);
#				consense_introns(\%site, \$coverage, \$alignment);
#
##				if ($site == 435) {
##					print $alignment,"\n";
##					die;
##				}
#			}
#		}
#	}
#
#	return ($alignment, $coverage);
	#return $exon;
	return ($exon, $distinct_seqs_per_intron);
}

sub group_all_introns {
	my $introns = shift;
	my %introns = %{$introns};

#	print Dumper($introns{755}),"\n";
#	print Dumper($introns{756}),"\n";

	# Assign unique name to each sequence which includes type and site of occurrence
	my $count = 0;
	my %intron_seqs;
	foreach my $site (sort { $b <=> $a } keys %introns) {
		my %site = %{$introns{$site}};
		foreach my $intron (keys %{$site{INTRONS}}) {
			$intron_seqs{"full_$count"."_$site"} = $intron;
			$count++;
		}
		foreach my $intron (keys %{$site{PART_INT_END}}) {
			$intron_seqs{"end_$count"."_$site"} = $intron;
			$count++;
		}
		foreach my $intron (keys %{$site{PART_INT_START}}) {
			$intron_seqs{"start_$count"."_$site"} = $intron;
			$count++;
		}
	}
#	print Dumper(\%intron_seqs),"\n";

	# Write sequences to file
	open(my $check_fasta, ">", "all-intron-check.fasta");
	foreach my $intron (keys %intron_seqs) {
		print {$check_fasta} ">$intron\n";
		print {$check_fasta} replace_ambiguities($intron_seqs{$intron}),"\n";
		$count++;
	}
	close($check_fasta);

	# Cluster and assemble with cap3
	#my @cap3_output = `$cap3 all-intron-check.fasta -r 0 -k 0 -g 2`;
	my @cap3_output = `$cap3 all-intron-check.fasta -r 0 -k 0 -g 2 -h 90 -o 16`;
	unlink(glob("all-intron-check.fasta*"));
	
	# Parse cap3 output to determine which input introns belong to assembled contigs
	my $contig;
	my %contig_components;
	foreach my $line (@cap3_output) {
		chomp($line);

		# Extract component contigs
		if ($contig && $line =~ /(\S+_\d+_\d+)\+/) {
			while ($line =~ /(\S+_\d+_\d+)\+/g) {
				my $intron = $1;
				$contig_components{$contig}->{$1}++;
			}
		}

		# Reset in case we have multiple assembled contigs
		undef($contig) if ($line =~ /^\s*$/);

		# Start of Contig description
		if ($line =~ /\** Contig (\d+) \**/) {
			if (!exists($contig_components{"Contig$1"})) {
				$contig = "Contig$1";
			}
		}
	}

	#print Dumper(\%contig_components),"\n";

	my %site_counts;
	foreach my $contig (keys %contig_components) {
		my %components = %{$contig_components{$contig}};
		
		foreach my $component (keys %components) {
			(my $site = $component) =~ s/\S+_\d+_(\d+)/$1/;
			$site_counts{$site}++;
		}
	}
	#print Dumper(\%site_counts),"\n";
	#die;

	# Check whether any of the assembled contigs include introns from multiple sites
	foreach my $contig (keys %contig_components) {
		my %components = %{$contig_components{$contig}};
		
		# Check which sites are included
		my %included_sites;
		foreach my $component (keys %components) {
			(my $site = $component) =~ s/\S+_\d+_(\d+)/$1/;
			(my $type = $component) =~ s/(\S+)_\d+_\d+/$1/;
			push(@{$included_sites{$site}}, "$type"."_$intron_seqs{$component}");
		}
#		print "\n";

		# See if we have introns from more than one site
		if (scalar(keys %included_sites) > 1) {
			my @included_sites = sort { $a <=> $b } keys %included_sites;
			#my $base_intron = $included_sites[0];
			my $base_intron = $included_sites[$#included_sites];
			
#			my $base_intron;
#			my $most_hits = 0;
#			foreach my $site (@included_sites) {
#				if ($site_counts{$site} > $most_hits) {
#					$most_hits = $site_counts{$site};	
#					$base_intron = $site;
#				}
#			}

			# TODO: figure out better way to assign which site introns should be moved
			# TODO: add distance filter?

#			print Dumper(\%included_sites),"\n";

			# Move introns to the correct site
			print "moving intron(s) at the following sites to site $base_intron:\n";
			foreach my $index (0 .. $#included_sites - 1) {
			#foreach my $index (0 .. $#included_sites) {
				my $site = $included_sites[$index];
				my %site = %{$introns{$site}};

				#next if ($site == $base_intron);

				print "  $site\n";

				my @introns = @{$included_sites{$site}};

				foreach my $intron (@introns) {

					my $type;
					if ($intron =~ s/^(\S+)_//) {
						$type = $1;
					}

					# Delete the old location of the sequence and move it to its new location
					if ($type eq "full") {
						delete($site{INTRONS}->{$intron});
						$introns{$base_intron}{INTRONS}->{$intron}++;
					}
					elsif ($type eq "start") {
						delete($site{PART_INT_END}->{$intron});
						$introns{$base_intron}{PART_INT_END}->{$intron}++;
					}
					else {
						delete($site{PART_INT_START}->{$intron});
						$introns{$base_intron}{PART_INT_START}->{$intron}++;
					}

				}
			}
		}
	}
	print "\n";

	# Remove any sites with no introns
	foreach my $site (sort { $b <=> $a } keys %introns) {
		my %site = %{$introns{$site}};
		my %intron_sequences = (%{$site{INTRONS}}, %{$site{PART_INT_END}}, %{$site{PART_INT_START}});
		#print $site,"\n",Dumper(\%intron_sequences),"\n";
		#delete($introns{$site}) if (scalar(keys %intron_sequences) == 0);
		if (scalar(keys %intron_sequences) == 0) {
			#print join(" ",keys(%introns)),"\n";
			delete($introns->{$site});	
			#print join(" ",keys(%introns)),"\n";
		}
		#print Dumper($introns{755}),"\n";
	}
	#print Dumper($introns{755}),"\n";

	#die;

	return;
}

sub consense_bases {
	#my ($site, $coverage, $alignment) = @_;
	my ($site, $exon) = @_;

	my %bases = %{$site->{BASES}};
	my $target_base = $site->{TARGET};

	my $num_single_bases = 0;
	foreach my $base_present (keys %bases) {
		$num_single_bases += $bases{$base_present};
	}

#	foreach my $base_present (keys %bases) {
#		# Perhaps filter by proportion here? (remove base if it is present in < x% of contigs)
#
#		# Ambiguous, can't determine consensus
#		if ($base_present eq "N") {
#			#$$coverage .= "$num_single_bases ";
#			$$exon .= "N";
#			#next SITE;
#			return;
#		}
#	}

	# If only one nucleotide is found at the site output it
	if (scalar(keys %bases) == 1) {
		#$$coverage .= "$num_single_bases ";
		$$exon .= (keys %bases)[0];
	}
	# Handle ambiguous non-introns with IUPAC codes
	else {
		#$$coverage .= "$num_single_bases ";
		#$$exon .= get_IUPAC(\%bases);
		my @bases_by_frequency = sort { $bases{$b} <=> $bases{$a} } keys %bases;
		$$exon .= $bases_by_frequency[0];
	}

	return;
}

sub consense_introns {
	#my ($site, $coverage, $alignment) = @_;
	my ($site, $exon) = @_;

	my $num_complete_seqs = 0;

	#my %bases = %{$site->{BASES}};
	my $id = $site->{ID};
	my %introns = %{$site->{INTRONS}};
	my %part_int_end = %{$site->{PART_INT_END}};
	my %part_int_start = %{$site->{PART_INT_START}};

#	# Add check for differing lengths of sequences?
#	my $total_seqs = scalar(keys %introns) + scalar(keys %part_int_end) + scalar(keys %part_int_start);

#	if ($total_seqs) {

		#$$coverage .= "(";
#		if ($total_seqs == 1) {
#			my %introns = (%introns, %part_int_end, %part_int_start);
#
#			my $intron = lc((keys %introns)[0]);
#			$$exon .= $intron;
#
#			print $intron,"\n";
#			die;
#
#			#$$coverage .= "1 " x length($intron);
#		}
#		else {

			print "################CONSENSING NEW INTRON ($id)################\n";

			# Create a temporary fasta file to hold introns
			my $count = 0;
			my $longest_intron_length = 0;
			open(my $intron_fasta, ">", "introns.fa");
			foreach my $intron (keys %introns) {
				print {$intron_fasta} ">full_$count\n";	
				#print {$intron_fasta} "$intron\n";	
				print {$intron_fasta} "",replace_ambiguities($intron),"\n";	
				$longest_intron_length = length($intron) if (length($intron) > $longest_intron_length);
				$count++;
			}
			foreach my $partial (keys %part_int_start) {
				print {$intron_fasta} ">start_$count\n";	
				#print {$intron_fasta} "$partial\n";	
				print {$intron_fasta} "",replace_ambiguities($partial),"\n";	
				$count++;
			}
			foreach my $partial (keys %part_int_end) {
				print {$intron_fasta} ">end_$count\n";	
				#print {$intron_fasta} "$partial\n";	
				print {$intron_fasta} "",replace_ambiguities($partial),"\n";	
				$count++;
			}
			close($intron_fasta);

			# lower gap penalty may be needed for blat -fine
			#system("$cap3 introns.fa -r 0 -k 0 -g 2");
			#system("$cap3 introns.fa -r 0 -k 0 -g 2 > /dev/null");

			my @cap3_output = `$cap3 introns.fa -r 0 -k 0 -g 2`;
			
			# Parse cap3 output to determine which input introns belong to assembled contigs
			my $contig;
			my %contig_components;
			foreach my $line (@cap3_output) {
				chomp($line);

				# Extract component contigs
				if ($contig && $line =~ /(\S+_\d+)\+/) {
					while ($line =~ /(\S+_\d+)\+/g) {
						my $intron = $1;
						$contig_components{$contig}->{$1}++;
					}
				}

				# Reset in case we have multiple assembled contigs
				undef($contig) if ($line =~ /^\s*$/);

				# Start of Contig description
				if ($line =~ /\** Contig (\d+) \**/) {
					if (!exists($contig_components{"Contig$1"})) {
						$contig = "Contig$1";
					}
				}
			}

			# Contigs/singlets containing a 'full' intron or a joined 'start' and 'end'
			my %full_contigs;

			# Start or end of an intron
			my %end_contigs;
			my %start_contigs;

			# All contigs assembled by cap3
			my %contigs = parse_fasta("introns.fa.cap.contigs");

			# Parse sequences unassembled by cap3
			my %singlets;
			if (-s "introns.fa.cap.singlets" != 0) {
				%singlets = parse_fasta("introns.fa.cap.singlets");

				foreach my $singlet (keys %singlets) {
					#delete($singlets{$singlet}) if ($singlet !~ /full/);
					if ($singlet =~ /full/) {
						$full_contigs{$singlet} = $singlets{$singlet};
					}
					elsif ($singlet =~ /end/) {
						$end_contigs{$singlet} = $singlets{$singlet};
					}
					elsif ($singlet =~ /start/) {
						$start_contigs{$singlet} = $singlets{$singlet};
					}
				}
			}

			# Determine original components of cap3 assemblies contigs
			foreach my $contig (keys %contigs) {
				my %components = %{$contig_components{$contig}};

				my $contains_end;
				my $contains_full;
				my $contains_start;
				foreach my $component (keys %components) {
					if ($component =~ /end/) {
						$contains_end++;
					}
					if ($component =~ /full/) {
						$contains_full++;
					}
					if ($component =~ /start/) {
						$contains_start++;
					}
				}

				# Add full contigs to %full_contigs
				if ($contains_full || ($contains_end && $contains_start)) {
					$full_contigs{$contig} = $contigs{$contig};
				}
				elsif (!$contains_start && $contains_end) {
					$end_contigs{$contig} = $contigs{$contig};
				}
				elsif ($contains_start && !$contains_end) {
					$start_contigs{$contig} = $contigs{$contig};
				}
			}

			# Clean up
			unlink(glob("introns.fa*"));

			my $final_intron = '';
			if (scalar(keys %full_contigs)) {
				print "",scalar(keys %full_contigs), " complete sequence(s) at this intron.\n";

				$num_complete_seqs = scalar(keys %full_contigs);

				# TODO:
				# Rewrite to work for exonic sequence included at end of intron

				# Output longest
				my @contigs = sort { length($b) <=> length($a) } values %full_contigs;
				$final_intron = $contigs[0];

				#print $final_intron,"\n\n";

				$final_intron = remove_exon_from_intron_end($final_intron, $exon, $id);
				$final_intron = remove_exon_from_intron_start($final_intron, $exon, $id) if ($final_intron);

				#undef($final_intron) if ($final_intron =~ /^N+$/);
				$final_intron = '' if ($final_intron =~ /^N+$/);

	#			foreach my $end (keys %end_contigs) {
	#				$final_intron = remove_exon_from_intron_end($final_intron, $exon, $id);
	#			}

	#			foreach my $start (keys %start_contigs) {
	#				$final_intron = remove_exon_from_intron_start($final_intron, $exon, $id);
	#			}

#				if ($init_length != $final_length) {
#					print "init: $init_length, $final_length\n";
#					#die;
#				}
				#die "$id\n" if ($id > 755);
				#die if ($final_intron =~ /tataataggacag/i);

#				# Account for blat missing small alignments just before intron starts
#				if ($longest_intron_length > 0 && length($final_intron) > $longest_intron_length) {
#					print "Possible alignment error in blat detected:\n";
#					my $difference = length($final_intron) - $longest_intron_length;
#
#					my $contig_start = substr($final_intron, 0, $difference);
#					my $preceeding_exon = substr($$exon, -1 * $difference);
#
#					#print $preceeding_exon,"\n";
#
#					my $pid = get_pairwise_pid($preceeding_exon, $contig_start);
#					print "  First $difference bases of contig have $pid% identity to preceeding exon.\n";
#
#					if ($pid > 90) {
#						print "  Bases determined to be exon, removing from intron.\n";
#						$final_intron = substr($final_intron, $difference); 
#					}
#				}
				print "\n";
			}
			else {
				print "Full intron could not be determined.\n";

#				# Determine if we have partial intron starts/ends
#				my %end_contigs;
#				my %start_contigs;
#
#				# Check singlets
#				foreach my $singlet (keys %singlets) {
#					if ($singlet =~ /end/) {
#						$end_contigs{$singlet} = $singlets{$singlet};
#					}
#					elsif ($singlet =~ /start/) {
#						$start_contigs{$singlet} = $singlets{$singlet};
#					}
#				}

#				# Check assembled contigs
#				foreach my $contig (keys %contigs) {
#					my %components = %{$contig_components{$contig}};
#
#					my $contains_end;
#					my $contains_start;
#					foreach my $component (keys %components) {
#						if ($component =~ /end/) {
#							$contains_end++;
#						}
#						if ($component =~ /start/) {
#							$contains_start++;
#						}
#					}
#
#					if (!$contains_start && $contains_end) {
#						$end_contigs{$contig} = $contigs{$contig};
#					}
#					elsif ($contains_start && !$contains_end) {
#						$start_contigs{$contig} = $contigs{$contig};
#					}
#				}

				# Output the longest of the starts and ends present
				my @end_contigs = sort { length($b) <=> length($a) } values %end_contigs;
				my @start_contigs = sort { length($b) <=> length($a) } values %start_contigs;

				my $longest_end = $end_contigs[0];
				my $longest_start = $start_contigs[0];

				my $num_ns = 50;

				# Figure out how to account for blat alignment errors here
				if (defined($longest_end) && defined($longest_start)) {
					print "Outputting stitched intron including start and end.\n";
					#my $preceeding_exon = substr($$alignment, -1 * length($longest_start));
					$longest_end = remove_exon_from_intron_end($longest_end, $exon, $id);
					$longest_start = remove_exon_from_intron_start($longest_start, $exon, $id);

					$longest_end = '' if ($longest_end =~ /^N+$/);
					$longest_start = '' if ($longest_start =~ /^N+$/);

					if ($longest_end && $longest_start) {
						$final_intron = $longest_start.('N' x $num_ns).$longest_end;
					}
					elsif ($longest_end) {
						$final_intron = ('N' x $num_ns).$longest_end;
					}
					else {
						$final_intron = $longest_start.('N' x $num_ns);
					}
					print "\n";

					#$preceeding_exon =~ s/[A-Z]*[a-z]+//g;
				}
				elsif (!defined($longest_end)) {
					print "Outputting partial intron start.\n";
					#my $preceeding_exon = substr($$alignment, -1 * length($longest_start));
					$longest_start = remove_exon_from_intron_start($longest_start, $exon, $id);

					$longest_start = '' if ($longest_start =~ /^N+$/);

					if ($longest_start) {
						$final_intron = $longest_start.('N' x $num_ns);
					}
					print "\n";

					#$preceeding_exon =~ s/[A-Z]*[a-z]+//g;
				}
				else {
					print "Outputting partial intron end.\n\n";
					$longest_end = remove_exon_from_intron_end($longest_end, $exon, $id);

					undef($longest_end) if ($longest_end =~ /^N+$/);

					if ($longest_end) {
						$final_intron = ('N' x $num_ns).$longest_end;
					}
				}
			}

			#$$alignment .= lc($longest_contig);
			#$$exon .= lc($final_intron);
			#$$coverage .= $total_seqs; #placeholder

#		}
#
#		#chop($$coverage);
#		#$$coverage .= ") ";
#	}

	#return lc($final_intron);
	return lc($final_intron), $num_complete_seqs;
}

sub remove_exon_from_intron_end {
	#my ($intron, $alignment) = @_;
	my ($intron, $alignment, $site) = @_;

	#print "Checking for alignment errors in blat:\n";

	my $exon = substr($$alignment, $site, length($intron));
	#print "$exon\n\n";

	# Write the preceeding exon and intron we are checking to file
	open(my $check, ">", "intron-check.fasta");
	print {$check} ">intron\n";
	print {$check} $intron,"\n";
	close($check);

	open($check, ">", "target-check.fasta");
	print {$check} ">exon\n";
	print {$check} $exon,"\n";
	close($check);

	# Run blat with recommended settings for maximum sensitivity
	system("$blat target-check.fasta intron-check.fasta -fine -tileSize=6 -stepSize=1 -fine -noHead intron-check.psl >/dev/null");

	my $has_error;
	open(my $blat_out, "<", "intron-check.psl");
	while (my $line = <$blat_out>) {
		chomp($line);

		my @line = split(/\s+/, $line);

		# Blat headers
		my ($match, $mismatch, $rep_matches, $n_count, $q_num_inserts, $q_base_inserts, 
			$t_num_inserts, $t_base_inserts, $strand, $q_name, $q_size, $q_start, $q_end,
			$t_name, $t_size, $t_start, $t_end, $block_count, $block_sizes, $q_starts, $t_starts) = @line; 

		next if ($q_name eq $t_name);

		#if ($q_name eq "exon") {
		my @sizes = split(",", $block_sizes);
		my @t_starts = split(",", $t_starts); #target
		my @q_starts = split(",", $q_starts); #intron

		# Check if match exists between the intron and the current target
		#if ($t_size == $sizes[$#sizes] + $t_starts[$#t_starts]) {
		if ($t_starts[0] == 0) {
			print "Alignment error in blat detected, removing ",($q_size - $q_starts[0])," bases from intron end.\n";
			#$intron = substr($intron, 0, $q_size - $q_starts[0]);
			$intron = substr($intron, 0, $q_starts[0]);
			$has_error++;
			last;
		}
	
	}
	close($blat_out);

	# Clean up 
	unlink(glob("intron-check.*"));
	unlink(glob("target-check.fa*"));

#	if (!$has_error) {
#		print "No alignment errors detected.\n";
#	}

	return $intron;
}

sub remove_exon_from_intron_start {
	#my ($intron, $alignment) = @_;
	my ($intron, $alignment, $site) = @_;

	#print "Checking for alignment errors in blat:\n";

	#$preceeding_exon =~ s/[A-Z]*[a-z]+//g;
	#my $exon = substr($$alignment, -1 * length($intron));
	#my $exon = substr($$alignment, 0, $site);
	my $exon = substr($$alignment, 0, $site + 1);
	#my $exon = substr($$alignment, $site - length($intron) + 1, length($intron));
	#print "meow\n",$$alignment,"\n\n";

	#print "$exon\n\n";

	# Write the preceeding exon and intron we are checking to file
	open(my $check, ">", "intron-check.fasta");
	print {$check} ">intron\n";
	print {$check} $intron,"\n";
	close($check);

	open($check, ">", "target-check.fasta");
	print {$check} ">exon\n";
	print {$check} $exon,"\n";
	close($check);

	# Run blat with recommended settings for maximum sensitivity
	#system("$blat intron-check.fasta intron-check.fasta -fine -tileSize=6 -stepSize=1 -fine -noHead intron-check.psl");
	system("$blat target-check.fasta intron-check.fasta -fine -tileSize=6 -stepSize=1 -fine -noHead intron-check.psl >/dev/null");

	my $has_error;
	open(my $blat_out, "<", "intron-check.psl");
	while (my $line = <$blat_out>) {
		chomp($line);

		my @line = split(/\s+/, $line);

		# Blat headers
		my ($match, $mismatch, $rep_matches, $n_count, $q_num_inserts, $q_base_inserts, 
			$t_num_inserts, $t_base_inserts, $strand, $q_name, $q_size, $q_start, $q_end,
			$t_name, $t_size, $t_start, $t_end, $block_count, $block_sizes, $q_starts, $t_starts) = @line; 

		next if ($q_name eq $t_name);

		#if ($q_name eq "exon") {
		my @sizes = split(",", $block_sizes);
		my @t_starts = split(",", $t_starts); #target
		my @q_starts = split(",", $q_starts); #intron

		# Check if match exists between the intron and the current target
		if ($t_size == $sizes[$#sizes] + $t_starts[$#t_starts]) {
			print "Alignment error in blat detected, removing ",($q_starts[$#q_starts] + $sizes[$#sizes])," bases from intron start.\n";
#			print $$alignment,"\n\n";
#			print $intron,"\n\n";
			# Remove everything until the detected match
			$intron = substr($intron, $q_starts[$#q_starts] + $sizes[$#sizes]);
#			print $intron,"\n\n";

#			open(my $check, ">", "intron-check.fasta");
#			print {$check} ">intron\n";
#			print {$check} $intron,"\n\n";
#			close($check);
#
#			open($check, ">", "target-check.fasta");
#			print {$check} ">target\n";
#			print {$check} $$alignment,"\n";
#			close($check);
#
#			system("$blat target-check.fasta intron-check.fasta -fine -tileSize=6 -stepSize=1 -fine -noHead intron-check.psl");
			$has_error++;
			last;
		}
	
	}
	close($blat_out);

	# Clean up 
	unlink(glob("intron-check.*"));
	unlink(glob("target-check.fa*"));

#	if (!$has_error) {
#		print "  No alignment errors detected.\n";
#	}

	return $intron;
}

sub get_pairwise_pid {
	my ($seq1, $seq2) = @_;

	$seq1 = uc($seq1);
	$seq2 = uc($seq2);

	my $total = 0;
	my $match = 0;
	foreach my $index (0 .. length($seq1) - 1) {
		my $site1 = chop($seq1);
		my $site2 = chop($seq2);

		next if ($site1 !~ /[ATCG]/ || $site2 !~ /[ATCG]/);

		$total++;
		$match++ if ($site1 eq $site2);
	}
	return 0 if ($total == 0);

	return $match / $total * 100;
}

sub replace_ambiguities {
	my ($seq) = @_;

	if ($seq =~ /[MRWSYKVHDBN]/) {
		my $base;
		my $rand = rand(1);

		while ($seq =~ /R/g) {
			($rand < 0.5 ) ? $base = 'A' : $base = 'G';
			$seq =~ s/R/$base/;
		}

		while ($seq =~ /Y/g) {
			($rand < 0.5 ) ? $base = 'C' : $base = 'T';
			$seq =~ s/Y/$base/;
		}

		while ($seq =~ /S/g) {
			($rand < 0.5 ) ? $base = 'G' : $base = 'C';
			$seq =~ s/S/$base/;
		}

		while ($seq =~ /W/g) {
			($rand < 0.5 ) ? $base = 'A' : $base = 'T';
			$seq =~ s/W/$base/;
		}

		while ($seq =~ /K/g) {
			($rand < 0.5 ) ? $base = 'G' : $base = 'T';
			$seq =~ s/K/$base/;
		}

		while ($seq =~ /M/g) {
			($rand < 0.5 ) ? $base = 'A' : $base = 'C';
			$seq =~ s/M/$base/;
		}

		while ($seq =~ /B/g) {
			if ($rand < 0.3333) {
				$base = 'C';
			}
			elsif ($rand < 0.6666) {
				$base = 'G';
			}
			else {
				$base = 'T';
			}
			$seq =~ s/B/$base/;
		}

		while ($seq =~ /D/g) {
			if ($rand < 0.3333) {
				$base = 'A';
			}
			elsif ($rand < 0.6666) {
				$base = 'G';
			}
			else {
				$base = 'T';
			}
			$seq =~ s/D/$base/;
		}

		while ($seq =~ /H/g) {
			if ($rand < 0.3333) {
				$base = 'A';
			}
			elsif ($rand < 0.6666) {
				$base = 'C';
			}
			else {
				$base = 'T';
			}
			$seq =~ s/H/$base/;
		}

		while ($seq =~ /V/g) {
			if ($rand < 0.3333) {
				$base = 'A';
			}
			elsif ($rand < 0.6666) {
				$base = 'C';
			}
			else {
				$base = 'G';
			}
			$seq =~ s/V/$base/;
		}
	}
	return $seq;
}

sub get_IUPAC {
	my $bases = shift;

	my %bases = %{$bases};

	# Don't let N's confuse consensus
	delete($bases{N});

	my $code;
	if (scalar(keys %bases) == 1) {
		$code = "A" if (exists($bases{"A"}));
		$code = "T" if (exists($bases{"T"}));
		$code = "C" if (exists($bases{"C"}));
		$code = "G" if (exists($bases{"G"}));
	}
	elsif (scalar(keys %bases) == 2) {
		$code = "R" if (exists($bases{"A"}) && exists($bases{"G"}));
		$code = "Y" if (exists($bases{"C"}) && exists($bases{"T"}));
		$code = "S" if (exists($bases{"G"}) && exists($bases{"C"}));
		$code = "W" if (exists($bases{"A"}) && exists($bases{"T"}));
		$code = "K" if (exists($bases{"G"}) && exists($bases{"T"}));
		$code = "M" if (exists($bases{"A"}) && exists($bases{"C"}));
	}
	elsif (scalar(keys %bases) == 3) {
		$code = "B" if (exists($bases{"C"}) && exists($bases{"G"}) && exists($bases{"T"}));
		$code = "D" if (exists($bases{"A"}) && exists($bases{"G"}) && exists($bases{"T"}));
		$code = "H" if (exists($bases{"A"}) && exists($bases{"C"}) && exists($bases{"T"}));
		$code = "V" if (exists($bases{"A"}) && exists($bases{"C"}) && exists($bases{"G"}));
	}
#	else {
#		$code = "N";
#	}
#
#	if (!defined($code)) {
#		use Data::Dumper;
#		print Dumper(\%bases);
#		die "IUPAC code not defined\n";
#	}

	$code = "N" if (!defined($code));

	return $code;
}

sub get_partial_seq {
	my ($range, $seq) = (@_);

	my $partial_seq;

	my @range = split(",", $range);
	foreach my $segment (@range) {
		my ($start, $end) = split("-", $segment);
		$partial_seq .= substr($seq, $start, $end - $start + 1);	
	}
	
	return $partial_seq;
}

sub intersection {
	my ($new, $old) = (@_);

	my @new = split(",", $new);
	my @old = split(",", $old);

	# Populate a hash with all members of @new and @old

	my %old;
	foreach my $segment (@old) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}
	foreach my $segment (@new) {
		my ($start, $end) = split("-", $segment);
		foreach my $site ($start .. $end) {
			$old{$site}++;
		}
	}

	# Remove site if both arrays don't have it
	foreach my $site (keys %old) {
		if ($old{$site} != 2) {
			delete($old{$site});
		}
	}

	my @sites = sort { $a <=> $b } keys %old;

	return "0-0" if (scalar(@sites) == 0);

	# Convert sites to a string (i.e. 1,2,3,4,6,7,8 => 1-4,6-8)

	my @ends;
	my @starts = ($sites[0]);
	my $previous = $sites[0] - 1;
	foreach my $index (0 .. $#sites) {
		my $site = $sites[$index];
		if ($previous != $site - 1) {
			push(@ends, $sites[$index - 1]);
			push(@starts, $site);
		}
		$previous = $site;
	}
	push(@ends, $sites[$#sites]);

	my $intersection;
	foreach my $index (0 .. $#starts) {
		$intersection .= $starts[$index]."-".$ends[$index].",";
	}
	chop($intersection);

	return $intersection;
}

sub rev_comp {
	my $seq = shift;
	
	my %comp = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C');

	my $rev_comp;
	my $rev = reverse($seq);
	foreach my $index (0 .. length($rev) - 1) {
		my $char = substr($rev, $index, 1);
		#$rev_comp .= $comp{$char};
		if (exists($comp{$char})) {
			$rev_comp .= $comp{$char};
		}
		else {
			$rev_comp .= 'N';
		}
	}

	return $rev_comp;
}

sub get_median {
	my $array = shift;
	my @array = @{$array};

	@array = sort { $a <=> $b } @array;

	my @percentiles = qw(50);

	my $median;
	foreach my $percentile (@percentiles) {
		my $percentile_index = ($percentile / 100) * (scalar(@array) + 1);

		my $value;
		if ($percentile_index =~ /(\d+)(\.\d+)/) {
			my $integer_part = $1;
			my $decimal_part = $2;

			my $diff = ($array[$integer_part] - $array[$integer_part - 1]);
			$value = ($diff * $decimal_part) + $array[$integer_part - 1];
		}
		else {
			$value = $array[$percentile_index - 1];
		}
		$return = $value;
	}

	return $return;
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
		#if ($line =~ /^>(.*)/) {
		if ($line =~ /^>(\S+)/) {
			$taxon = $1;
		}
		else {
			# Taxon sequence
			$taxon =~ s/-/_/g;
			$align{$taxon} .= uc($line);
		}
	}
	close($alignment_file);
	
	return %align;
}

sub okay_to_run {
	my $pids = shift;

	# Free up a CPU by sleeping for 10 ms
	usleep(10000);

	my $current_forks = scalar(@{$pids});
	foreach my $index (reverse(0 .. $current_forks - 1)) {
		next if ($index < 0);

		my $pid = @{$pids}[$index];
		my $wait = waitpid($pid, WNOHANG);

		# Successfully reaped child
		if ($wait > 0) {
			$current_forks--;
			#splice(@pids, $index, 1);
			splice(@{$pids}, $index, 1);
		}
	}

	return ($current_forks < $free_cpus);
}

sub get_free_cpus {

	my $os_name = $^O;

	# Returns a two-member array containing CPU usage observed by top,
	# top is run twice as its first output is usually inaccurate
	my @percent_free_cpu;
	if ($os_name eq "darwin") {
		# Mac OS
		chomp(@percent_free_cpu = `top -i 1 -l 2 | grep "CPU usage"`);
	}
	else {
		# Linux
		chomp(@percent_free_cpu = `top -bn2d0.05 | grep "Cpu(s)"`);
	}

	my $percent_free_cpu = pop(@percent_free_cpu);

	if ($os_name eq "darwin") {
		# Mac OS
		$percent_free_cpu =~ s/.*?(\d+\.\d+)%\s+id.*/$1/;
	}
	else {
		# linux 
		$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
		$percent_free_cpu = eval($percent_free_cpu);
	}

	my $total_cpus;
	if ($os_name eq "darwin") {
		# Mac OS
		$total_cpus = `sysctl -n hw.ncpu`;
	}
	else {
		# linux
		$total_cpus = `grep --count 'cpu' /proc/stat` - 1;
	}

	my $free_cpus = ceil($total_cpus * $percent_free_cpu / 100);

	if ($free_cpus == 0 || $free_cpus !~ /^\d+$/) {
		$free_cpus = 1; # assume that at least one cpu can be used
	}
	
	return $free_cpus;
}
