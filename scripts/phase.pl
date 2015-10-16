#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

#TODO: account for indels (reference length > 1)

my $ref = 'A_digitata.con.fasta';
my $vcf = 'a_digitata_0.vars.vcf';
my $hapcompass_out = 'hc_MWER_solution.txt';
#my $hapcompass_out = 'hc_MWER_solution.txt.test';

my %snps = parse_vcf($vcf);
my %ref = parse_fasta($ref);
my %blocks = parse_hapcompass_out($hapcompass_out);
#print Dumper(\%snps),"\n";
#print Dumper(\%blocks),"\n";

foreach my $sequence (keys %ref) {
	my %snps = %{$snps{$sequence}};
	my %blocks = %{$blocks{$sequence}};

	# Check for SNPs that indicate error in reference, edit reference
	my @sorted_snps = sort { $b <=> $a } keys %snps;
	foreach my $index (0 .. $#sorted_snps) {
		my $snp = $sorted_snps[$index];
		my %snp = %{$snps{$snp}};
		if ($snp{GT} =~ /^(1\/?)+$/) {
			print "Replacing reference ($snp{REF} -> $snp{ALT}) at site $snp with genotype: $snp{GT} as it is not a true SNP.\n";
			my $case_sensitive_ref = substr($ref{$sequence}, $snp - 1, 1);

			# Reference was exon
			if ($case_sensitive_ref =~ /[A-Z]/) {
				substr($ref{$sequence}, $snp - 1, length($snp{REF})) = $snp{ALT};
			}
			# Reference was intron
			else {
				substr($ref{$sequence}, $snp - 1, length($snp{REF})) = lc($snp{ALT});
			}
			#substr($ref{$sequence}, $snp - 1, length($snp{REF})) = $snp{ALT};

			# Check if we introduced an offset in site indices
			my $offset_length = length($snp{REF}) - length($snp{ALT});

			# Correct for offsets
			if ($offset_length != 0) {
				#print " offset of length $offset_length detected.\n";
				foreach my $index (reverse(0 .. $index - 1)) {
					my $snp = $sorted_snps[$index];
					#print " moving $snp by ",(-1 * $offset_length),"\n";

					# Replace entry in %snps
					$snps{$snp - $offset_length} = delete($snps{$snp});

					# Replace entry in %blocks
					if (exists($blocks{$snp})) {
						$blocks{$snp - $offset_length} = delete($blocks{$snp});
					}

					# Update indices in @sorted_snps
					$sorted_snps[$index] = $snp - $offset_length;
				}
			}
			#delete($snps{$snp});
		}
	}

	my $block_id = 0;
	my $current_block;
	my @all_blocks; # includes unphased SNPs
	#foreach my $snp (sort { $a <=> $b } keys %snps) {
	foreach my $index (reverse(0 .. $#sorted_snps)) {
		my $snp = $sorted_snps[$index];
		#print $snp,"\n";
		my %snp = %{$snps{$snp}};
		#print Dumper(\%snp),"\n";
		
		# Reference does not occur, skip
		if ($snp{GT} =~ /^(1\/?)+$/) {
			#print "Skipping SNP at site $snp with genotype: $snp{GT} as it is not a true SNP.\n";
			next;
		}
		# Skip deletions
		if (length($snp{REF}) > 1) {
			print "Skipping indel at site $snp.\n";
			next;
		}

		# Check for occurrence of multiple alleles at site
		my @alts = split(",", $snp{ALT});
		if (scalar(@alts) > 1) {

			my $case_sensitive_ref = substr($ref{$sequence}, $snp - 1, 1);

			# Reference was exon
			if ($case_sensitive_ref =~ /[A-Z]/) {
				print "Skipping multiallelic variant at site $snp ($case_sensitive_ref -> N).\n";
				substr($ref{$sequence}, $snp - 1, length($snp{REF})) = 'N';
			}
			# Reference was intron
			else {
				print "Skipping multiallelic variant at site $snp ($case_sensitive_ref -> n).\n";
				substr($ref{$sequence}, $snp - 1, length($snp{REF})) = 'n';
			}
			next;
		}
		# Skip insertions
		if (length($snp{ALT}) > 1) {
			print "Skipping indel at site $snp.\n";
			next;
		}

		# Sanity check remove eventually
		my $case_sensitive_ref = substr($ref{$sequence}, $snp - 1, 1);
		if ($case_sensitive_ref =~ /[A-Z]/) {
			if ($case_sensitive_ref ne $snp{REF}) {
				die "Mismatch in reference ($snp: $ref != $snp{ALT}).\n";
			}
		}
		else {
			if ($case_sensitive_ref ne lc($snp{REF})) {
				die "Mismatch in reference ($snp: $ref != $snp{ALT}).\n";
			}
		}

		# SNP is phased
		if (exists($blocks{$snp})) {
			my %block_snp = %{$blocks{$snp}};
			#print Dumper(\%block_snp),"\n";

			# Check if we are in a new block
			$block_id++ if (defined($current_block) && $current_block != $block_snp{BLOCK});
			$current_block = $block_snp{BLOCK};

			$all_blocks[$block_id]{$snp}->{REF} = $snp{REF};
			$all_blocks[$block_id]{$snp}->{ALT} = $snp{ALT};
			$all_blocks[$block_id]{$snp}->{PHASING} = $block_snp{PHASING};
		}
		# SNP is unphased
		else {
			$block_id++;

			$all_blocks[$block_id]{$snp}->{REF} = $snp{REF};
			$all_blocks[$block_id]{$snp}->{ALT} = $snp{ALT};
			$all_blocks[$block_id]{$snp}->{PHASING} = $snp{GT};
		}
		#print "\n";
	}
	print "Reference sequence '$sequence' has ",scalar(@all_blocks)," separate block(s).\n";
	my @block_indices = get_block_indices($ref{$sequence}, \@all_blocks);
	my @phased_blocks = get_phased_blocks($ref{$sequence}, \@all_blocks, \@block_indices);

	# Output blocks to separate files
	foreach my $index (0 .. $#phased_blocks) {
		my @block_sequences = @{$phased_blocks[$index]};
		open(my $out, ">", "$sequence-block$index.fasta");
		foreach my $sequence_id (0 .. $#block_sequences) {
			my $sequence = $block_sequences[$sequence_id];
			print {$out} ">$sequence_id\n";
			print {$out} "$sequence\n";
		}
		close($out);
	}

	die;
	
	#my @phasings = get_all_possible_phasings($ref{$sequence}, \@all_blocks);
	#print Dumper(\@phasings),"\n";
	
	#die;
	
#	open(my $out, ">", "$sequence-phasings.fasta");
#	foreach my $index (0 .. $#phasings) {
#		my $phasing = $phasings[$index];	
#		print {$out} ">$index\n";
#		print {$out} "$phasing\n";
#	}
#	close($out);
}

sub get_block_indices {
	my ($reference, $all_blocks) = @_;	
	my @all_blocks = @{$all_blocks};

	my @ends;
	my @starts;
	foreach my $block (@all_blocks) {
		my @sorted_snps = sort { $a <=> $b } keys %{$block};

		# Get start and end SNP for the block
		# Subtract 1 to convert to 0-based indexing
		my $start = $sorted_snps[0] - 1;
		my $end = $sorted_snps[$#sorted_snps] - 1;

#		print "$start-$end\n";
#		print "@sorted_snps\n\n";

		push(@starts, $start);	
		push(@ends, $end);
	}

#	print Dumper(\@starts),"\n";
#	print Dumper(\@ends),"\n";

	# Set the start of the first block
	$starts[0] = 0;

	# Set the end of the last block
	$ends[$#ends] = length($reference) - 1;

	my @block_indices;
	foreach my $index (0 .. $#starts - 1) {

		my $end = $ends[$index];
		my $start = $starts[$index];
		my $next_start = $starts[$index + 1];

		# Split nonvariant sequence between the two blocks
		my $diff = $next_start - $end;
		my $offset = int($diff / 2);
		
		$ends[$index] = $end + $offset;
		$starts[$index + 1] = $next_start - $offset + 1;
	}
#	print Dumper(\@starts),"\n";
#	print Dumper(\@ends),"\n";

	foreach my $index (0 .. $#starts) {
		my $end = $ends[$index];
		my $start = $starts[$index];
#		$block_indices[$index] = {$start => $end};
#		$block_indices[$index] = {$start => $end};
		$block_indices[$index]->{END} = $end;
		$block_indices[$index]->{START} = $start;

#		$all_blocks[$index]->{END} = $end;
#		$all_blocks[$index]->{START} = $start;
	}
	
	return @block_indices;
}

sub get_phased_blocks {
	my ($reference, $all_blocks, $block_indices) = @_;	
	my @all_blocks = @{$all_blocks};
	my @block_indices = @{$block_indices};

	# Determine how many phasings we have
	my $block = $all_blocks[0];
	my $snp = (keys %{$block})[0];
	my $phasing = $block->{$snp}->{PHASING};
	my @states = split("/", $phasing);

	my @phased_blocks;
	foreach my $index (0 .. $#all_blocks) {
		my $block = $all_blocks[$index];
		my $end = $block_indices[$index]->{END};
		my $start = $block_indices[$index]->{START};

		# Create copies of the reference for modification
		my @block_sequences = (($reference) x scalar(@states));
		#my $block = substr($reference, $start, $end - $start + 1);

		# Look at SNPs in each block to determine what the nucleotide should be
		foreach my $snp (sort { $a <=> $b } keys %{$block}) {
			my $ref = $block->{$snp}->{REF};
			my $alt = $block->{$snp}->{ALT};
			my $phasing = $block->{$snp}->{PHASING};

			my @states = split("/", $phasing);
			#foreach my $state (@states) {
			foreach my $index (0 .. $#states) {
				my $state = $states[$index];

				# Store possible alleles with corresponding ID
				my @nucleotides = ($ref);
				push(@nucleotides, split(",", $alt));

				# Determine if reference was exon or intron
				#my $case_sensitive_ref = substr($sequence, $snp - 1, 1);
				my $case_sensitive_ref = substr($block_sequences[$index], $snp - 1, 1);

				# Reference was exon
				if ($case_sensitive_ref =~ /[A-Z]/) {
					#substr($sequence, $snp - 1, 1) = $nucleotides[$state];
					substr($block_sequences[$index], $snp - 1, 1) = $nucleotides[$state];
				}
				# Reference was intron
				else {
					#substr($sequence, $snp - 1, 1) = lc($nucleotides[$state]);
					substr($block_sequences[$index], $snp - 1, 1) = lc($nucleotides[$state]);
				}
			}
		}

		# Trim to start and end indices
		foreach my $index (0 .. $#block_sequences) {
			$block_sequences[$index] = substr($block_sequences[$index], $start, $end - $start + 1);
		}

		push(@phased_blocks, \@block_sequences);
	}
	#print Dumper(@phased_blocks),"\n";
	#die;

	return @phased_blocks;
}

sub get_all_possible_phasings {
	my ($reference, $all_blocks) = @_;	
	my @all_blocks = @{$all_blocks};

	my @phasings;

	# Determine number of haplotypes present in each block
	my $block = $all_blocks[0];
	my $snp = (keys %{$block})[0];
	my $phasing = $block->{$snp}->{PHASING};
	my @states = split("/", $phasing);

	# Add in separator
	my @digits = '0' .. $#states;
	foreach my $index (0 .. $#digits) {
		$digits[$index] .= '/';
	}

	# Generate possible combinations
	my $digits = do { local $" = ","; "{@digits}" };
	my @possible_phasings = @{[glob $digits x scalar(@all_blocks)]};
	chop(@possible_phasings); # Remove trailing '/'

	# Convert blocks and their haplotype to actual sequences
	foreach my $possible_phasing (@possible_phasings) {
		my $sequence = $reference;

		#print "phase info: $possible_phasing\n";

		# Stores which haplotype each block has
		my @hap_indices = split("/", $possible_phasing);
		foreach my $index (0 .. $#hap_indices) {
			my $block = $all_blocks[$index];
			my $hap_index = $hap_indices[$index];

			#print " block: $index\n";
			#print " haplotype: $hap_index\n";

			# Look at SNPs in each block to determine what the nucleotide should be
			foreach my $snp (sort { $a <=> $b } keys %{$block}) {
				my $ref = $block->{$snp}->{REF};
				my $alt = $block->{$snp}->{ALT};
				my $phasing = $block->{$snp}->{PHASING};
				my @states = split("/", $phasing);

				# The value of the SNP for this haplotype
				my $state = $states[$hap_index];

				# Store possible alleles with corresponding ID
				my @nucleotides = ($ref);
				push(@nucleotides, split(",", $alt));

				# Determine if reference was exon or intron
				my $case_sensitive_ref = substr($sequence, $snp - 1, 1);

				# Reference was exon
				if ($case_sensitive_ref =~ /[A-Z]/) {
					substr($sequence, $snp - 1, 1) = $nucleotides[$state];
				}
				# Reference was intron
				else {
					substr($sequence, $snp - 1, 1) = lc($nucleotides[$state]);
				}

				#print "$snp: $ref ",substr($sequence, $snp - 1, 1),"\n";
				#substr($sequence, $snp - 1, 1) = $nucleotides[$state];
			}
		}
		#die;
		push(@phasings, $sequence);
	}
	#print scalar(@phasings),"\n";

	return @phasings;
}

sub parse_vcf {
	my $vcf = shift;

	# Open vcf file
	my %snps;
	open(my $vcf_file, "<", $vcf);
	while (my $line = <$vcf_file>) {
		next if ($line =~ /^\s*#/);

		# Extract required information about SNPs
		my @line = split(/\s+/, $line);	
		my ($seq, $site, $ref, $alt, $info, $genotype) = ($line[0], $line[1], $line[3], $line[4], $line[7], $line[9]);
		$genotype =~ s/(\S+?):.*/$1/;

		$snps{$seq}{$site}->{REF} = $ref;
		$snps{$seq}{$site}->{ALT} = $alt;
		$snps{$seq}{$site}->{GT} = $genotype;
	}
	close($vcf_file);

	return %snps;
}

sub parse_hapcompass_out {
	my $hapcompass_out = shift;

	my $ref;
	my %blocks;
	my $block_id = -1;
	open(my $hapcompass, "<", $hapcompass_out);
	while(my $line = <$hapcompass>) {
		chomp($line);

		my @line = split(/\s+/, $line);

		# Start of new phased block
		if ($line =~ /^BLOCK.*\s+(\S+)$/) {
			$ref = $1;	
			$block_id++;
			next;
		}

		# Block phasing
		if ($block_id >= 0 && $line) {
			my $site = $line[1];	

			my @phasing = @line[3 .. $#line];
			#my $phasing = join("|", @phasing);
			my $phasing = join("/", @phasing);

			#print "$ref $site @phasing $phasing\n";

			$blocks{$ref}{$site}->{BLOCK} = $block_id;
			$blocks{$ref}{$site}->{PHASING} = $phasing;
		}
	}
	close($hapcompass);

	return %blocks;
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
			$align{$taxon} .= $line;
		}
	}
	close($alignment_file);
	
	return %align;
}
