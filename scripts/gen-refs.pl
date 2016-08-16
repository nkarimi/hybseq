#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Data::Dumper;
use Getopt::Long;
use Fcntl qw(:flock SEEK_END);

# Turn on autoflush
$|++;

# Lengths in base pairs, ideal values depend on sequencing strategy
my $min_exon_length = 10;
#my $min_exon_length = 15;
#my $min_exon_length = 20;
#my $min_misaligned_intron_length = 100;
my $min_misaligned_intron_length = 30;
#my $min_overhang_length = 40;
#my $min_overhang_length = 50;
my $min_overhang_length = 30;
#my $max_overhang_length = 500; # not used yet
my $max_overhang_length = 600;
#my $max_overhang_length = 200; # not used yet
my $min_overhang_length_discrepancy = 10; # not used yet

# We want to multithread when we can
my $free_cpus = &get_free_cpus;

# Check that required executables are in path
my $blat = check_path_for_exec("blat");
my $mafft = check_path_for_exec("mafft");
my $spaln = check_path_for_exec("spaln");

# Parse command line invocation
my $targets;
#my $allow_N = 0;
#my $allow_N = 1;
my $exon_only;
my $contig_coverage;
my $exclude_unmapped;
my $include_unspanned;
GetOptions(
	't|targets=s' => \$targets,
	'e|exon-only' => \$exon_only,
	'c|contig-coverage' => \$contig_coverage,
	'x|exclude-unmapped' => \$exclude_unmapped,
	'u|include-unspanned' => \$include_unspanned,
#	'N|allow-N' => \$allow_N,
);

my $input = shift;

# Input error checking
die "You must specify an assembly as input.\n" if (!defined($input));
die "You must specify a fasta file containing target sequences (-t).\n" if (!defined($targets));
die "You must you have specified an incorrect argument.\n" if (@ARGV);

# Output file names
my $blat_out_name = $input.".psl";
#(my $consensus_output = $input) =~ s/\.fa(sta)?/.con.fasta/i;
(my $consensus_output = $input) =~ s/(\.con)?\.fa(sta)?/.con.fasta/i;

print "$blat_out_name\n\n";

# Run preliminary blat to determine which targets and contigs have alignable sequence
print "Running blat against targets and contigs...\n";
my $return = system("$blat '$input' '$targets' -t=dna -q=dna -noHead -repMatch=1000000 '$blat_out_name'") if (!-e $blat_out_name);
#my $return = system("$blat '$input' '$targets' -t=dna -q=dna -noHead -repMatch=1000000 '$blat_out_name'");
die "Error running blat: '$return'.\n" if ($return);
print "Blat completed, parsing hits...\n";

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

	# Add hit to target
	$targets{$q_name}->{$t_name}++;
}
close($blat_out);
print "Hits successfully parsed.\n";

# Load contig and target sequences
my %contig_seqs = parse_fasta($input);
my %target_seqs = parse_fasta($targets);

my %ignore;
#if (!$allow_N) {
#	foreach my $contig (keys %contig_seqs) {
#		$ignore{$contig}++ if ($contig_seqs{$contig} =~ /N/);
#	}
#}

# Create consensus sequence for each target
my $count = 0;
my %troublesome;
my %final_alignments;
print("list of targets after blat:\n");
for (keys %targets){ print("$_\n"); }

foreach my $target (keys %targets) {
	# list of defeated troublesome contigs
	#$target = '161|Gorai.012G146100.2|';
	#$target = '596|Gorai.006G215300.1|';
	#$target = '959|Gorai.001G017000.1|';
	#$target = '11792|Gorai.007G202700.5|';
	#$target = '13902|Gorai.011G084700.1|';
	#$target = '1984|Gorai.001G079800.1|';
	#$target = '7115|Gorai.005G200400.1|';
	#$target = '20778|Gorai.005G183500.2|';
	#$target = '20311|Gorai.011G239400.1|';
	#$target = '2996|Gorai.004G267500.2|';
	#$target = '6227|Gorai.010G034900.1|';
	#$target = '13870|Gorai.008G045800.3|';
	#$target = '3168|Gorai.013G033700.2|';
	#$target = '9480|Gorai.002G201200.1|';	
	#$target = '1127|Gorai.011G245700.4|';	
	#$target = '471|Gorai.006G256000.2|';
	#$target = '8547|Gorai.003G086700.1|';
	#$target = '6801|Gorai.008G191000.1|';
	#$target = '4196|Gorai.005G137100.2|';

	my @hits = keys %{$targets{$target}};
	my $target_seq = $target_seqs{$target};

	print "\nAligning the following contigs to target '$target':\n";

	# Stores locations and sequences of introns
	my %introns;
	my %potential_introns; 
	$potential_introns{ENDS} = {};
	$potential_introns{STARTS} = {};
	$final_alignments{$target}->{$target} = $target_seq;

	# Use spaln to align exons of each contig with the reference and identify introns
	foreach my $seqname (@hits) {
		my $seq = $contig_seqs{$seqname};

		next if ($ignore{$seqname});

		# Perform alignment with spaln
		align_contig({CONTIG_NAME => $seqname, 
		              CONTIG_SEQ => $seq, 
					  TARGET_NAME => $target, 
					  TARGET_SEQ => $target_seq, 
					  INTRONS => \%introns,
					  POT_INTRONS => \%potential_introns});
	}

	# Think remove_intron_from_exon needs to be run first?

	# Check for and correct misalignments
	move_misaligned_introns(\%introns, \%potential_introns, $target, $target_seq);

	# Insert unspanned introns into the alignment
	if ($include_unspanned) {
		insert_potential_introns(\%introns, \%potential_introns, $target, $target_seq);
	}

	print "  Aligning and inserting introns (",(scalar(keys %introns)),")...\n";

	# Insert intronic sequence into exon
	foreach my $site (sort { $b <=> $a } keys %introns) {
		print "    Consensing introns at site: $site...\t";

		# Determine which contigs don't have an intron at this site
		my %intronless_seqs = map { ($_ ne $target) ? ($_ => 1) : () } keys %{$final_alignments{$target}};
		foreach my $intron (@{$introns{$site}}) {
			#my $seq = $intron->{SEQ};
			my $seqname = $intron->{SEQNAME};
			delete($intronless_seqs{$seqname});
		}

		# Correct for minor misalignments which can occur with spaln
		foreach my $seq (keys %intronless_seqs) {
			remove_intron_from_exon($site, $seq, $target, \%introns);
		}

		# Align all introns present at this site
		my %aligned;
		if (scalar(@{$introns{$site}}) > 1) {

			# Output introns to fasta file for alignment
			open(my $tmp, ">", "intron-$site.fasta");
			foreach my $intron (@{$introns{$site}}) {
				my $seq = $intron->{SEQ};
				my $seqname = $intron->{SEQNAME};

				# Correct gaps put into sequence by having ambiguous sites at intron end
				if ($seq =~ /(N+$)/) {
					my $num_ns = length($1);
					#my $following_exon = substr($final_alignments{$target}->{$seqname}, $site, $min_exon_length);
					my $following_exon = substr($final_alignments{$target}->{$seqname}, $site, length($final_alignments{$target}->{$seqname}) - $site);
					$following_exon =~ s/-+$//;

					print "$seqname\n";
					print "$seq\n";
					print "$following_exon\n";
					print "$num_ns\n";

					# Check if there are misplaced gaps
					if ($following_exon =~ /(^-+)/) {

						# Move N's out of intron and into exon
						my $mismatch_length = length($1);
						if ($mismatch_length && $mismatch_length <= $num_ns) {
							#print "\nadjusting $seqname end ($site)...\n";
							substr($intron->{SEQ}, length($intron->{SEQ}) - $mismatch_length, $mismatch_length) = '';
							substr($final_alignments{$target}->{$seqname}, $site, $mismatch_length) = 'N' x $mismatch_length;
							$intron->{T_END} -= $mismatch_length;
						}
					}
					print "\n";
				}
				# Correct gaps put into sequence by having ambiguous sites at intron start
				if ($seq =~ /(^N+)/) {
					my $num_ns = length($1);
					#my $preceeding_exon = substr($final_alignments{$target}->{$seqname}, $site - $min_exon_length, $min_exon_length);
					my $preceeding_exon = substr($final_alignments{$target}->{$seqname}, 0, $site);
					$preceeding_exon =~ s/^-+//;

					# Check if there are misplaced gaps
					if ($preceeding_exon =~ /(-+$)/) {

						# Move N's out of intron and into exon
						my $mismatch_length = length($1);
						if ($mismatch_length && $mismatch_length <= $num_ns) {
							#print "\nadjusting $seqname start ($site)...\n";
							substr($intron->{SEQ}, 0, $mismatch_length) = '';
							substr($final_alignments{$target}->{$seqname}, $site - $mismatch_length, $mismatch_length) = 'N' x $mismatch_length;
							$intron->{T_START} += $mismatch_length;
							$intron->{REF_START} += $mismatch_length;
						}
					}
				}
				$seq = $intron->{SEQ};

				# Pad with exon to facilitate alignment
				my $start_exon = substr($contig_seqs{$seqname}, $intron->{T_START} - $min_exon_length, $min_exon_length);
				my $end_exon = substr($contig_seqs{$seqname}, $intron->{T_END}, $min_exon_length);
				my $padded_seq = $start_exon.$seq.$end_exon;

				print {$tmp} ">$seqname\n";
				print {$tmp} "$seq\n";
				#print {$tmp} "$padded_seq\n";
			}
			close($tmp);

			# Align with mafft, parse results, and clean up
			my $return = system("$mafft --maxiterate 1000 --localpair --thread $free_cpus intron-$site.fasta > intron-$site.aln.fasta 2> /dev/null");
			die "Error occured while running mafft: $!.\n" if ($return);
			%aligned = parse_fasta("intron-$site.aln.fasta");
			unlink("intron-$site.aln.fasta");
		}
		else {

			die "MEOWMEOMWEOMWEOWMEOEMW\n" if (scalar(@{$introns{$site}} == 0));

			# Only one contig has an intron here, alignment isn't needed
			my $intron = @{$introns{$site}}[0];
			my $seq = $intron->{SEQ};
			my $seqname = $intron->{SEQNAME};

			# Pad with exon so we can handle insertion the same way as when there are multiple sequences
			my $start_exon = substr($contig_seqs{$seqname}, $intron->{T_START} - $min_exon_length, $min_exon_length);
			my $end_exon = substr($contig_seqs{$seqname}, $intron->{T_END}, $min_exon_length);
			my $padded_seq = $start_exon.$seq.$end_exon;

			%aligned = ($seqname => $seq);
			#%aligned = ($seqname => $padded_seq);
		}

		# How many nucleotides are in the final alignment
		my $align_length = length((values %aligned)[0]);
		#my $align_length = length((values %aligned)[0]) - 2 * $min_exon_length;
		print "finished intron consensus, length = $align_length.\n";

		# Insert aligned sequence or required number of gaps
		foreach my $seq (keys %{$final_alignments{$target}}) {

			# Intronic sequence identified for this contig
			if (exists($aligned{$seq})) {
				substr($final_alignments{$target}->{$seq}, $site, 0) = $aligned{$seq};	
				#substr($final_alignments{$target}->{$seq}, $site - $min_exon_length, 2 * $min_exon_length) = $aligned{$seq};	
			}
			# Intronic sequence was not identified for this contig, add gaps
			else {
				substr($final_alignments{$target}->{$seq}, $site, 0) = "-" x $align_length;	
			}
		}
		unlink("intron-$site.fasta");
	}
	print "\n";

	# Output contigs aligned to target sorted by starting index
	open(my $out, ">", "$target.fasta");
#	print {$out} ">$target\n";
#	print {$out} "$final_alignments{$target}->{$target}\n";
	foreach my $contig (sort { local $a = $a;
							   local $b = $b;
							   ($a = $final_alignments{$target}->{$a}) =~ s/^(-*).*/$1/;
							   ($b = $final_alignments{$target}->{$b}) =~ s/^(-*).*/$1/;
							   length($a) <=> length($b) } keys %{$final_alignments{$target}}) {
		next if ($contig eq $target);

		# Add coverage to contig name
		my $start = 0;
		if ($final_alignments{$target}->{$contig} =~ /(^-+)/) {
			$start = length($1);	
		}

		my $end = length($final_alignments{$target}->{$contig}) - 1;
		if ($final_alignments{$target}->{$contig} =~ /(-+$)/) {
			$end = $end - length($1);	
		}

		# Output coverage if specified by user
		if ($contig_coverage) {
			print {$out} ">$contig|$start-$end|\n";
		}
		else {
			print {$out} ">$contig\n";
		}

		print {$out} "$final_alignments{$target}->{$contig}\n";
		#print "$contig: ".length($final_alignments{$target}->{$contig})."\n";
	}
	close($out);

	$count++;
	#last if ($count >= 5);
	#last;
	#print Dumper($final_alignments{$target}),"\n";
	#die;
}

# Create a consensus for each target which can be used for SNP/haplotype calling
my %consensuses = consense_alignments(\%final_alignments);

# Output consensus sequence to fasta
open(my $out, ">", $consensus_output);
foreach my $consensus (keys %consensuses) {
	print {$out} ">$consensus\n";
	print {$out} "$consensuses{$consensus}\n";
}
close($out);

#if (!$allow_N) {
#	print "perl $0 $input -t $consensus_output -N\n";
#	system("perl $0 $input -t $consensus_output -N");
#}

## Output troublesome targets
#open($out, ">", "troublesome.txt");
#foreach my $target (keys %troublesome) {
#	print {$out} "$target\n";
#}
#close($out);

sub align_contig {
	my $opts = shift;
	my $seq = $opts->{CONTIG_SEQ};
	my $name = $opts->{CONTIG_NAME};
	my $target_seq = $opts->{TARGET_SEQ};
	my $target_name = $opts->{TARGET_NAME};
	my $introns = $opts->{INTRONS};
	my $potential_introns = $opts->{POT_INTRONS};
	my $ignore_introns = $opts->{IGNORE_INTRONS};

	$ignore_introns++ if ($exon_only);

	# Output contig to file
	open(my $tmp, ">", "$name.fasta");
	print {$tmp} ">$name\n";
	print {$tmp} "$seq\n";
	close($tmp);

	my $contig_length = length($seq);

	# Output reference to file
	open($tmp, ">", "target.fasta");
	print {$tmp} ">$target_name\n";
	print {$tmp} "$target_seq\n";
	close($tmp);

	# For more human readable output 
#	system("$spaln -Q1 -S1 -LS -O4 $name.fasta target.fasta");
#	print "\n\n";
#	system("$spaln -Q1 -S1 -LS -O5 $name.fasta target.fasta");
#	print "\n\n";

	# Align with spaln -08 output is cigar format
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -t$free_cpus -ya2 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -M1 -t$free_cpus -ya2 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -t$free_cpus -ya2 -yy30 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE10 -yi50 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE10 -yi6 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE10 -yi11,8,11 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -M1 -yE$min_exon_length -yi11,8,11 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -M1 -yE$min_exon_length -yi20,8,20 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE$min_exon_length -yi20,8,20 -vN7,5,8 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE$min_exon_length -yi20,8,20 -vN8,6,9 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE$min_exon_length -yi20,8,20 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE$min_exon_length -t$free_cpus -ya3 $name.fasta target.fasta`);
	
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE$min_exon_length -yi20,8,20 -t$free_cpus -ya3 $name.fasta target.fasta`);
	#chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE$min_exon_length -yi20,8,20 -t$free_cpus -ya2 $name.fasta target.fasta`);
	chomp(my @spaln_out = `$spaln -Q3 -S1 -LS -O8 -yE$min_exon_length -yi40,8,40 -t$free_cpus -ya2 '$name.fasta' target.fasta`);
	die "Error occured while aligning with spaln: $!.\n" if ($?);
	#exit(0) if $ignore_introns;
	unlink("$name.fasta");
	unlink("target.fasta");

#	print Dumper(\@spaln_out),"\n";

	my $line = shift(@spaln_out);
	return 0 if (!defined($line));

	# Remove cigar identifier
	$line =~ s/cigar: //;

	my @line = split(/\s+/, $line);

	# Extract general alignment info from cigar
	my ($qid, $qstart, $qend, $qstrand, $tid, $tstart, $tend, $tstrand, $score) = splice(@line, 0, 9);

	# If spaln prefers reverse complement hit (even though -S1 should make only forward hits...)
	if ($tstrand eq '-') {
		$seq = rev_comp($seq);
		$contig_seqs{$name} = $seq;

		$tstrand = "+";
		$tstart = length($seq) - $tstart - 1;
		$tend = length($seq) - $tend - 1;
	}

	# Extract alignment info
	my $cigar = join("", @line);

	# Check for and remove exons which are too small, assume they are part of the flanking intron
	while ($cigar =~ /(N(\d+)(.*?)N(\d+))/g) {
	#while ($cigar =~ /(N(\d+)(.*?)N(\d+)(\S)(\d+))/g) {
		my $full = $1;
		my $intron1_length = $2;
		my $intron2_length = $4;
		my $exon = $3;
#		my $following_type = $5;
#		my $following_length = $6;

		# Calculate length of exon between introns
		my $exon_length = 0;
		while ($exon =~ /(\d+)/g) {
			$exon_length += $1;
		}

		# Remove exon and join introns if it's too short
		if ($exon_length <= $min_exon_length) {
			print "  $cigar\n";
			print "  $1\n";
			
##			my $total_length = $intron1_length + $intron2_length + $exon_length;
##			my $sub = "N$total_length"."I$exon_length"."$following_type$following_length";
##			$cigar =~ s/\Q$full\E/$sub/;
#
#			my $total_length = $intron1_length + $intron2_length + $exon_length;
#			$cigar =~ s/\Q$full\E/N$total_length/;
			
			die "Exon with length less than min_exon_length input by spaln.\n";
		}
	}

	print "  $qid $qstart $qend $tid $tstart $tend\n";
	print "  $cigar\n\n";

	# Handle potential introns occuring at alignment overhangs
	if ($potential_introns && !$ignore_introns) {
		# Might allow overhangs, especially if $qstart == 0 or $qend == length of reference sequence?
		# Set maximum length, 300 bp? truncate or completely remove?
	
		# Extract unaligned overhangs from start and end of contigs
		my $end_overhang = substr($contig_seqs{$tid}, 0, $tstart);
		my $start_overhang = substr($contig_seqs{$tid}, $tend, length($contig_seqs{$tid}) - $tend);

		# Add end overhang to list of potential introns
		if ($end_overhang && length($end_overhang) >= $min_overhang_length && length($end_overhang) <= $max_overhang_length) {
			push(@{$potential_introns->{ENDS}->{$qstart}}, {SEQNAME => $tid,
												        	REF_START => $qstart,
												        	T_START => 0,
												        	T_END => $tstart,
												        	SEQ => $end_overhang});	
		}
#		elsif ($end_overhang && length($end_overhang) > $max_overhang_length) {
#			substr($end_overhang, length($end_overhang) - $max_overhang_length, $max_overhang_length) = '';
#			push(@{$potential_introns->{ENDS}->{$qstart}}, {SEQNAME => $tid,
#												        	REF_START => $qstart,
#												        	#T_START => 0,
#															T_START => $tstart - $max_overhang_length,
#												        	T_END => $tstart,
#												        	SEQ => $end_overhang});	
#		}
		# Add start overhang to list of potential introns
		if ($start_overhang && length($start_overhang) >= $min_overhang_length && length($start_overhang) <= $max_overhang_length) {
			push(@{$potential_introns->{STARTS}->{$qend}}, {SEQNAME => $tid,
												            REF_START => $qend,
												            T_START => $tend,
												            T_END => length($contig_seqs{$tid}),
												            SEQ => $start_overhang});	
		}
#		elsif ($start_overhang && length($start_overhang) > $max_overhang_length) {
#			substr($start_overhang, $max_overhang_length, length($start_overhang) - $max_overhang_length) = '';
#			push(@{$potential_introns->{STARTS}->{$qend}}, {SEQNAME => $tid,
#												            REF_START => $qend,
#												            T_START => $tend,
#												            T_END => $tend + $max_overhang_length,
#												            SEQ => $start_overhang});	
#		}
	}

	# Add gaps before alignment begins
	$final_alignments{$target_name}->{$tid} .= "-" x ($qstart);

	# Loop through cigar alignment
	my $query_index = 0;
	my $target_index = 0;
	while ($cigar =~ /([A-Z])(\d+)/g) {
		my $type = $1;
		my $length = $2;

		# Deletion is an intron
		if ($type eq "N") {
			if (!$ignore_introns) {
				push(@{$introns->{$qstart + $query_index}}, {SEQNAME => $tid, REF_START => $qstart + $query_index, T_START => $tstart + $target_index,
															 T_END => $tstart + $target_index + $length, SEQ => substr($seq, $tstart + $target_index, $length)});
			}
			$target_index += $length;
		}
		# Deletion, missing base in reference compared to contig
		elsif ($type eq "D") {
			$target_index += $length;
		}
		# Insertion, extra base in reference compared to contig
		elsif ($type eq "I") {
			$final_alignments{$target_name}->{$tid} .= "-" x $length;
			$query_index += $length;
		}
		# Alignable sequence
		else {
			$final_alignments{$target_name}->{$tid} .= substr($seq, $tstart + $target_index, $length);
			$target_index += $length;
			$query_index += $length;
		}
	}
#	print "\n";

	# Add gaps until alignment ends
	$final_alignments{$target_name}->{$tid} .= "-" x (length($target_seq) - $qend);

#	if ($tid eq "s_micrantha_rep_c58155") {
#		print $final_alignments{$target_name}->{$tid},"\n";
#	}

	return;
}

sub adjust_ambiguous_intron_exon_boundaries {
	my ($introns, $target) = @_;

	#foreach my $site (sort { $b <=> $a } keys %introns) {
	foreach my $site (sort { $b <=> $a } keys %{$introns}) {
		#foreach my $intron (@{$introns{$site}}) {
		foreach my $intron (@{$introns->{$site}}) {
			my $seq = $intron->{SEQ};
			my $seqname = $intron->{SEQNAME};

			# Correct gaps put into sequence by having ambiguous sites at intron end
			if ($seq =~ /(N+$)/) {
				my $num_ns = length($1);
				my $following_exon = substr($final_alignments{$target}->{$seqname}, $site, $min_exon_length);
				#my $following_exon = substr($final_alignments{$target}->{$seqname}, $site, length($final_alignments{$target}->{$seqname}) - $site);
				#$following_exon =~ s/-+$//;
				print "$seqname\n";
				print "$seq\n";
				print "$following_exon\n";
				print "$num_ns\n";

				# Check if there are misplaced gaps
				if ($following_exon =~ /(^-+)/) {

					# Move N's out of intron and into exon
					my $mismatch_length = length($1);
					print "$mismatch_length\n";

					if ($mismatch_length && $mismatch_length <= $num_ns) {
						print "\nadjusting $seqname end ($site)...\n";
						substr($intron->{SEQ}, length($intron->{SEQ}) - $mismatch_length, $mismatch_length) = '';
						substr($final_alignments{$target}->{$seqname}, $site, $mismatch_length) = 'N' x $mismatch_length;
						$intron->{T_END} -= $mismatch_length;
						#$intron->{T_START} -= $mismatch_length;
					}
				}
				print "\n";
			}
			# Correct gaps put into sequence by having ambiguous sites at intron start
			if ($seq =~ /(^N+)/) {
				my $num_ns = length($1);
				#my $preceeding_exon = substr($final_alignments{$target}->{$seqname}, $site - $min_exon_length, $min_exon_length);
				my $preceeding_exon = substr($final_alignments{$target}->{$seqname}, 0, $site);
				$preceeding_exon =~ s/^-+//;

				# Check if there are misplaced gaps
				if ($preceeding_exon =~ /(-+$)/) {

					# Move N's out of intron and into exon
					my $mismatch_length = length($1);
					if ($mismatch_length && $mismatch_length <= $num_ns) {
						print "\nadjusting $seqname start ($site)...\n";
						substr($intron->{SEQ}, 0, $mismatch_length) = '';
						substr($final_alignments{$target}->{$seqname}, $site - $mismatch_length, $mismatch_length) = 'N' x $mismatch_length;
						$intron->{T_START} += $mismatch_length;
						$intron->{REF_START} += $mismatch_length;
					}
				}
			}
		}
	}

	return;
}

sub remove_intron_from_exon {
	my ($site, $seqname, $target, $introns) = @_;
	my %introns = %{$introns};

	my $seq = $final_alignments{$target}->{$seqname};

	# Extract sequence before and after intron occurence
	my $pre_intron = substr($seq, 0, $site);
	my $post_intron = substr($seq, $site, length($seq) - $site);

	# Remove gaps from sequences
	(my $pre_intron_ungapped = $pre_intron) =~ s/-//g;
	(my $post_intron_ungapped = $post_intron) =~ s/-//g;

	# Check that the intron actually occurs in a region covered by this sequence
	if ($pre_intron_ungapped && $post_intron_ungapped) {

		# Check length of sequence is under threshold
		if (length($pre_intron_ungapped) <= $min_misaligned_intron_length || length($post_intron_ungapped) <= $min_misaligned_intron_length) {
			
			# Sequence before intron is most likely intron misaligned as exon
			if (length($pre_intron_ungapped) <= $min_misaligned_intron_length) {
				$final_alignments{$target}->{$seqname} = "-" x length($pre_intron) . $post_intron;
			}
			# Sequence after intron is most likely intron misaligned as exon
			else {
				$final_alignments{$target}->{$seqname} = $pre_intron . "-" x length($post_intron);
			}
		}
	}
	return;
}

sub insert_potential_introns {
	#my ($potential_introns, $target, $target_seq) = @_;
	my ($introns, $potential_introns, $target, $target_seq) = @_;

	# Extract potential intron ends
	my %ends = %{$potential_introns->{ENDS}};
	my @end_sites = sort { $a <=> $b } keys %ends;

	# Extract potential intron starts
	my %starts = %{$potential_introns->{STARTS}};
	my @start_sites = sort { $a <=> $b } keys %starts;

	# Join locations of potential introns
	my @all_sites = sort { $a <=> $b } (@end_sites, @start_sites);

	# Determine number of contigs which contain each intron
	my %counts;
	my %ambiguous_counts;
	foreach my $site (keys %starts) {
		foreach my $intron (@{$starts{$site}}){
			$counts{$site}++;

#			if ($intron->{SEQ} =~ /(^N+)/) {
#				$ambiguous_counts{START}->{$site} = length($1);	
#			}
#			elsif ($intron->{SEQ} =~ /(N+$)/) {
#				$ambiguous_counts{START}->{$site} = length($1);	
#			}
		}
	}
	foreach my $site (keys %ends) {
		foreach my $intron (@{$ends{$site}}){
			$counts{$site}++;

#			if ($intron->{SEQ} =~ /(^N+)/) {
#				$ambiguous_counts{START}->{$site} = length($1);	
#			}
#			elsif ($intron->{SEQ} =~ /(N+$)/) {
#				$ambiguous_counts{START}->{$site} = length($1);	
#			}
		}
	}

#	my %clusters;
#	my %in_cluster;
#	my @sites = sort { $a <=> $b } keys %{$introns};
#	#foreach my $index (0 .. $#sites - 1) {
#	foreach my $index (0 .. $#all_sites - 1) {
#		#my $site1 = $sites[$index];
#		my $site1 = $all_sites[$index];
#
#		my $site1_end_offset = 0; 
#		my $site1_start_offset = 0;
#
#		if (exists($ambiguous_counts{START}->{$site1})) {
#			$site1_start_offset = $ambiguous_counts{START}->{$site1};
#		}
#		if (exists($ambiguous_counts{END}->{$site1})) {
#			$site1_end_offset = $ambiguous_counts{END}->{$site1};
#		}
#		my $site1_offset = $site1_start_offset + $site1_end_offset;
#
#		next if ($in_cluster{$site1});
##		foreach my $index2 ($index + 1 .. $#sites) {
#		foreach my $index2 ($index + 1 .. $#all_sites) {
#			#my $site2 = $sites[$index2];
#			my $site2 = $all_sites[$index2];
#
#			my $site2_end_offset = 0; 
#			my $site2_start_offset = 0;
#
#			if (exists($ambiguous_counts{START}->{$site2})) {
#				$site2_start_offset = $ambiguous_counts{START}->{$site2};
#			}
#			if (exists($ambiguous_counts{END}->{$site2})) {
#				$site2_end_offset = $ambiguous_counts{END}->{$site2};
#			}
#			my $site2_offset = $site2_start_offset + $site2_end_offset;
#
#			if ($site2 - $site2_offset - $site1 - $site1_offset <= $min_exon_length) {
#				if (!exists($in_cluster{$site1})) {
#					push(@{$clusters{$site1}}, $site2);
#					$in_cluster{$site2} = $site1;
#				}
#				else {
#					push(@{$clusters{$in_cluster{$site1}}}, $site2);
#					$in_cluster{$site2} = $site1;
#				}
#			}
#			else {
#				last;
#			}
#		}
#	}
#
#	foreach my $site (@all_sites) {
#		if (!$in_cluster{$site}) {
#			push(@{$clusters{$site}}, $site);
#		}
#	}

	# Stores what's already been clustered
	my %clusters;
	my %in_cluster;
	foreach my $index (0 .. $#all_sites - 1) {

		my $site1 = $all_sites[$index];
		next if ($in_cluster{$site1});
		foreach my $index2 ($index + 1 .. $#all_sites) {

			my $site2 = $all_sites[$index2];
			if ($site2 - $site1 <= $min_overhang_length_discrepancy) {
				if (!exists($in_cluster{$site1})) {
					push(@{$clusters{$site1}}, $site2);
					$in_cluster{$site2} = $site1;
				}
				else {
					push(@{$clusters{$in_cluster{$site1}}}, $site2);
					$in_cluster{$site2} = $site1;
				}
			}
			else {
				last;
			}
		}
	}
##	print Dumper(\%clusters),"\n";
##	print Dumper(\%ends),"\n";
##	print Dumper(\%starts),"\n";
##	die;

	# Remove clusters which cluster to a spanned intron, add sequences to that intron
	CLUSTER: foreach my $cluster (keys %clusters) {
		unshift(@{$clusters{$cluster}}, $cluster);
			
		# Iterate through each spanned intron location
		foreach my $spanned_intron_site (keys %{$introns}) {

			# Check each site in cluster against known intron site
			foreach my $incomplete_intron_site (@{$clusters{$cluster}}) {

				# Check if cluster and intron fall within minimum distance required for joining
				if (abs($incomplete_intron_site - $spanned_intron_site) <= $min_exon_length) {

					# Adding partials to spanned introns, lot of work to get right, not much benefit
					#
					# Add each partial sequence to the complete intron
#					foreach my $site (@{$clusters{$cluster}}) {
#						#print Dumper($introns->{$spanned_intron_site}),"\n";
#						#print Dumper($ends{$site}),"\n";
#						#print Dumper($starts{$site}),"\n";
##						foreach my $end (@{$ends{$site}}) {
##							#my %meow = %{$end};
##							#push(@{$introns->{$spanned_intron_site}}, $end);
##							push(@{$introns->{$spanned_intron_site}}, $end);
##							print "adding partial to $site\n";
##						}
##						foreach my $start (@{$starts{$site}}) {
##							#my %meow = %{$start};
##							#push(@{$introns->{$spanned_intron_site}}, $start);
##							push(@{$introns->{$spanned_intron_site}}, $start);
##							print "adding partial to $site\n";
##						}
#						print "adding partial to $site\n";
#						my $starts = delete($starts{$site});	
#						my $ends = delete($ends{$site});	
#						#print Dumper($starts);
#						#print Dumper($ends);
#						#print Dumper(\@{$introns->{$spanned_intron_site}}),"\n";
#						push(@{$introns->{$spanned_intron_site}}, @{$ends}) if (defined($ends));
#						push(@{$introns->{$spanned_intron_site}}, @{$starts}) if (defined($starts));
#						#print Dumper(\@{$introns->{$spanned_intron_site}}),"\n";
#
#					#	print @{$ends{$site}},"\n" if (exists($ends{$site}));
#					#	print @{$starts{$site}},"\n" if (exists($starts{$site}));
#				#		print $ends{$site},"\n" if (exists($ends{$site}));
#				#		print $starts{$site},"\n" if (exists($starts{$site}));
#						#push(@{$introns->{$spanned_intron_site}}, @{$ends{$site}}) if (exists($ends{$site}));
#						#push(@{$introns->{$spanned_intron_site}}, @{$starts{$site}}) if (exists($starts{$site}));
#					#	push(@{$introns->{$spanned_intron_site}}, @{$ends{$site}}) if (exists($ends{$site}));
#					#	push(@{$introns->{$spanned_intron_site}}, @{$starts{$site}}) if (exists($starts{$site}));
#					#	die;
#						#print Dumper($introns->{$spanned_intron_site}),"\n";
#						#die;
#					}
					
					# Remove this cluster
					#splice(@{$clusters{$cluster}}, $index, 1);
					delete($clusters{$cluster});
					next CLUSTER;
				}
#				elsif (scalar(@{$clusters{$cluster}} == 1)) {
#					delete($clusters{$cluster});
#					next CLUSTER;
#				}
			}
		}
	}
#	print Dumper($introns->{345}),"\n";
	#print Dumper($introns),"\n";
	#die;
	print "  Aligning and inserting unspanned introns (",(scalar(keys %clusters)),")...\n";

	#my $min_offset = 0;
	
	# Iterate through each cluster
	my %cluster_lengths;
	foreach my $cluster (sort { $b <=> $a } keys %clusters) {
	#foreach my $cluster (sort { $a <=> $b } keys %clusters) {
		#my $cluster = 1186;
		#my $cluster = 470;

		my $cluster_id = join("-", @{$clusters{$cluster}});

		# Check that this cluster contains both the start and end of this intron
		my $has_end;
		my $has_start;
		foreach my $site (@{$clusters{$cluster}}) {
			$has_end++ if (exists($ends{$site}));
			$has_start++ if (exists($starts{$site}));
		}

		# Only do something if we have the start and end of this intron
		if ($has_start && $has_end) {
			
			#my $exon_length = 20;
			my $exon_length = 10;

			# Determine maximum start and minimum end
			my %sites;
			my %contigs;
			my ($max_start, $min_end);
			foreach my $site (@{$clusters{$cluster}}) {
				next if ($sites{$site});

				# Calculate minimum end
				if (exists($ends{$site})) {
					foreach my $intron (@{$ends{$site}}) {
						my $seqname = $intron->{SEQNAME};
						$contigs{$seqname} = $site;
						$min_end = $site if (!defined($min_end) || $site < $min_end);
						$sites{$site}++;
					}
				}

				# Calculate maximum start
				if (exists($starts{$site})) {
					foreach my $intron (@{$starts{$site}}) {
						my $seqname = $intron->{SEQNAME};
						$contigs{$seqname} = $site;
						$max_start = $site if (!defined($max_start) || $site > $max_start);
						$sites{$site}++;
					}
				}
			}
			undef(%sites);
			my $max_offset = $max_start - $min_end;	

			# Account for nonoverlapping intron start and ends
			if ($max_offset < 0) {
				($max_start, $min_end) = fix_nonoverlapping_unspanned_introns($target, $target_seq, \%starts, \%ends, $max_start, $min_end, \%contigs, $clusters{$cluster});
				$max_offset = $max_start - $min_end;	
			}

			# Add spacer exon to facilitate alignment
			foreach my $site (@{$clusters{$cluster}}) {
				next if ($sites{$site});

				# Add exon to alignments containing intron ends
				if (exists($ends{$site})) {
					foreach my $intron (@{$ends{$site}}) {
						my $seqname = $intron->{SEQNAME};
						my $tend = $intron->{T_END};
						my $diff = $contigs{$seqname} - $min_end;
						my $exon = substr($contig_seqs{$seqname}, $tend, $exon_length - $diff);

						# Shouldn't happen with an exon length of 10
						die "Could not obtain exon for alignment of $cluster, $seqname ($tend, ",($exon_length - $diff),")\n" if (!defined($exon));
						$intron->{SEQ} .= $exon;
						$sites{$site}++;
					}
				}

				# Add exon to alignments containing intron starts
				if (exists($starts{$site})) {
					foreach my $intron (@{$starts{$site}}) {
						my $seqname = $intron->{SEQNAME};
						my $tstart = $intron->{T_START};
						my $diff = $max_start - $contigs{$seqname};
						my $exon = substr($contig_seqs{$seqname}, $tstart - $exon_length + $diff, $exon_length - $diff);

						# Shouldn't happen with an exon length of 10
						die "Could not obtain exon for alignment of $cluster, $seqname (",($tstart - $exon_length + $diff),", ",($exon_length - $diff),")\n" if (!defined($exon));
						$intron->{SEQ} = $exon.$intron->{SEQ};
						$sites{$site}++;
					}
				}
			}
			undef(%sites);

			print "    Inserting partial intron near site: $cluster...\t";

			# Outputting ending sequences of intron
			my $end_count = 0;
			open(my $tmp, ">", "intron-$cluster_id-ends.fasta");
			foreach my $site (@{$clusters{$cluster}}) {
				next if ($sites{$site});

				# Output intron ends
				if (exists($ends{$site})) {
					foreach my $intron (@{$ends{$site}}) {
						my $seq = $intron->{SEQ};
						my $seqname = $intron->{SEQNAME};
						print {$tmp} ">$seqname\n";
						print {$tmp} "$seq\n";
						#$contigs{$seqname} = $site;
						$end_count++;
						#$min_end = $site if (!defined($min_end) || $site < $min_end);
						$sites{$site}++;
					}
				}
			}
			close($tmp);
			undef(%sites);

			# Align sequences if needed
			my %aligned_ends;
			if ($end_count > 1) {
				# Align with mafft, parse results, and clean up
				my $return = system("$mafft --maxiterate 1000 --localpair --thread $free_cpus intron-$cluster_id-ends.fasta > intron-$cluster_id-ends.aln.fasta 2> /dev/null");
				die "Error occured while running mafft: $!.\n" if ($return);
				%aligned_ends = parse_fasta("intron-$cluster_id-ends.aln.fasta");
			}
			else {
				%aligned_ends = parse_fasta("intron-$cluster_id-ends.fasta");
			}
			unlink("intron-$cluster_id-ends.fasta");
			unlink("intron-$cluster_id-ends.aln.fasta");
			
			# Outputting starting sequences of intron
			my $start_count = 0;
			open($tmp, ">", "intron-$cluster_id-starts.fasta");
			foreach my $site (@{$clusters{$cluster}}) {
				next if ($sites{$site});

				# Output intron starts
				if (exists($starts{$site})) {
					foreach my $intron (@{$starts{$site}}) {
						my $seq = $intron->{SEQ};
						my $seqname = $intron->{SEQNAME};
						print {$tmp} ">$seqname\n";
						print {$tmp} "$seq\n";
						#$contigs{$seqname} = $site;
						$start_count++;
						#$max_start = $site if ($site > $max_start);
						$sites{$site}++;
					}
				}
			}
			close($tmp);

			# Align sequences if needed
			my %aligned_starts;
			if ($start_count > 1) {
				# Align with mafft, parse results, and clean up
				$return = system("$mafft --maxiterate 1000 --localpair --thread $free_cpus intron-$cluster_id-starts.fasta > intron-$cluster_id-starts.aln.fasta 2> /dev/null");
				#$return = system("$mafft --maxiterate 1000 --genafpair --thread $free_cpus intron-$cluster_id-starts.fasta > intron-$cluster_id-starts.aln.fasta 2> /dev/null");
				die "Error occured while running mafft: $!.\n" if ($return);
				%aligned_starts = parse_fasta("intron-$cluster_id-starts.aln.fasta");
			}
			else {
				%aligned_starts = parse_fasta("intron-$cluster_id-starts.fasta");
			}
			unlink("intron-$cluster_id-starts.fasta");
			unlink("intron-$cluster_id-starts.aln.fasta");

			# Determine lengths of start and end alignments
			my $aligned_end_length = length((values %aligned_ends)[0]) - $exon_length;
			my $aligned_start_length = length((values %aligned_starts)[0]) - $exon_length;

			# For correcting spanned intron indices after we insert the partial introns
			$cluster_lengths{$cluster} = $aligned_end_length + $aligned_start_length + $max_offset;
			
			# Insert aligned intron ends into alignment
			foreach my $seq (keys %{$final_alignments{$target}}) {
				#next if (exists($aligned_starts{$seq}));

				# Intronic sequence identified for this contig
				if (exists($aligned_ends{$seq})) {
					#print "inserting end at $contigs{$seq}...\n";
					#my $diff = $contigs{$seq} - $min_end;
					
					my $offset = $max_offset;
					substr($final_alignments{$target}->{$seq}, $min_end, $exon_length) = "-" x ($max_offset).$aligned_ends{$seq};
				}
				# Intronic sequence was not identified for this contig, add gaps
				else {
					substr($final_alignments{$target}->{$seq}, $max_start, 0) = "-" x ($aligned_end_length + $max_offset);
				}
			}

			# Insert aligned intron starts into alignment
			foreach my $seq (keys %{$final_alignments{$target}}) {
				#next if (exists($aligned_ends{$seq}));

				# Intronic sequence identified for this contig
				if (exists($aligned_starts{$seq})) {
					#my $diff = $max_start - $contigs{$seq};
					#print "inserting start at $contigs{$seq}...\n";
					substr($final_alignments{$target}->{$seq}, $max_start - $exon_length, $exon_length) = $aligned_starts{$seq};
				}
				# Intronic sequence was not identified for this contig, add gaps
				else {
					substr($final_alignments{$target}->{$seq}, $max_start, 0) = "-" x $aligned_start_length;	
				}
			}
#			foreach my $seq (keys %{$final_alignments{$target}}) {
#				if (exists($aligned_starts{$seq})) {
#					print "$seq: ",length($final_alignments{$target}->{$seq})," (start)\n";
#				}
#				elsif (exists($aligned_ends{$seq})) {
#					print "$seq: ",length($final_alignments{$target}->{$seq})," (end)\n";
#				}
#				else {
#					print "$seq: ",length($final_alignments{$target}->{$seq}),"\n";
#				}
#			}
			print "finished intron insertion, length of inserted sequence = ",($aligned_start_length + $aligned_end_length + $max_offset),".\n";
		}
		else {
			print "    $cluster == cluster? more like laCkLUSTER\n";
			delete($clusters{$cluster});
		}
		#die;
	}
	#die;

	# Correct intron indices being off now
	print "\n  Updating indices...\n";
	foreach my $site (sort { $b <=> $a } keys %{$introns}) {
		my $new_index = $site;
		foreach my $cluster (sort { $a <=> $b } keys %clusters) {
			if ($site > $cluster) {
				$new_index += $cluster_lengths{$cluster};
			}
		}
		$introns->{$new_index} = delete($introns->{$site});
		print "    Intron at $site is now at $new_index.\n";
	}
	print "\n";
	
	#die;

	# Special cases for end at 0 and length of target sequence

	return;
	#return 1 if ($min_offset < 0 && $target ne "20778|Gorai.005G183500.2|");
}

sub fix_nonoverlapping_unspanned_introns {
	my ($target, $target_seq, $starts, $ends, $max_start, $min_end, $contigs, $cluster) = @_;

	my @cluster = @{$cluster};

	my $exon_length = 10;
	my $max_offset = $max_start - $min_end;

	my $unaligned_exon = substr($target_seq, $max_start, abs($max_offset));
	my $unaligned_exon_padded = substr($target_seq, $max_start - 10, abs($max_offset) + 20);

	# Create a fasta file to facilitate placement of unaligned exon
	open(my $tmp, ">", "intron-overlap.fasta");
	print {$tmp} ">unaligned_exon\n";
	print {$tmp} "$unaligned_exon_padded\n";

	# Get the intron end sequence as well its following exon
	my $end_intron = @{$ends->{$min_end}}[0];
	my $end_seqname = $end_intron->{SEQNAME};
	my $tend = $end_intron->{T_END};
	my $end_diff = $contigs->{$end_seqname} - $min_end;
	my $end_exon = substr($contig_seqs{$end_seqname}, $tend, $exon_length - $end_diff);
	$end_intron = substr($end_intron->{SEQ}, length($end_intron->{SEQ}) - $exon_length, $exon_length);

	# Get the intron start sequence as well its preceeding exon
	my $start_intron = @{$starts->{$max_start}}[0];
	my $start_seqname = $start_intron->{SEQNAME};
	my $tstart = $start_intron->{T_START};
	my $start_diff = $max_start - $contigs->{$start_seqname};
	my $start_exon = substr($contig_seqs{$start_seqname}, $tstart - $exon_length + $start_diff, $exon_length - $start_diff);
	$start_intron = substr($start_intron->{SEQ}, 0, $exon_length);

	my $junction = "$start_exon$start_intron".("n" x 20)."$end_intron$end_exon";
	print {$tmp} ">junction\n";
	print {$tmp} "$junction\n";
	close($tmp);

	#my $return = system("$mafft --maxiterate 1000 --localpair --lop -20 --lep 0 --thread $free_cpus intron-$cluster_id-overlap.fasta > intron-$cluster_id-overlap.aln.fasta 2> /dev/null");
	my $return = system("$mafft --maxiterate 1000 --genafpair --ep 0 --op 10 --thread $free_cpus intron-overlap.fasta > intron-overlap.aln.fasta 2> /dev/null");
	die "Error occured while running mafft: $!.\n" if ($return);

	my %align = parse_fasta("intron-overlap.aln.fasta");
	# Testing
	#%align = ('unaligned_exon' => 'ACAAGTACAG-------------------------------CACCTCCCCTGAAATTAAT', 'junction' => 'ACAAGTACAGGTTGGTGCTGNNNNNNNNNNNNNNNNNNNNGCTCCACCTTTGAAATTAAT');
	#%align = ('unaligned_exon' => 'ACAAGTACAGCACCTCCCC-------------------------------TGAAATTAAT', 'junction' => 'ACAAGTACAGGTTGGTGCTGNNNNNNNNNNNNNNNNNNNNGCTCCACCTTTGAAATTAAT');
	#%align = ('unaligned_exon' => 'ACAAGTACAGCACC-------------------------------TCCCCTGAAATTAAT', 'junction' => 'ACAAGTACAGGTTGGTGCTGNNNNNNNNNNNNNNNNNNNNGCTCCACCTTTGAAATTAAT');

	my $aligned_exon = $align{unaligned_exon};

	my %sites;
	my $new_min_end = $min_end;
	my $new_max_start = $max_start;

	# Check if ungapped exon is present in alignment or if it was split by gaps
	if ($aligned_exon =~ /\Q$unaligned_exon\E/) {
		my $start = $-[0];

		# Unaligned exon matches with 'intron' start
		if ($start < length($junction) / 2) {

			# Remove the sequence from intron and add to exon
			foreach my $site (@cluster) {
				next if ($sites{$site});

				# Fix intron starts
				if (exists($starts->{$site})) {
					foreach my $intron (@{$starts->{$site}}) {
						my $seqname = $intron->{SEQNAME};
						my $tstart = $intron->{T_START};

						# Extract the misaligned exon
						my $start_diff = abs($max_offset);
						$start_intron = substr($intron->{SEQ}, 0, $start_diff);

						# Remove misaligned exon and add to exon alignment
						substr($intron->{SEQ}, 0, $start_diff) = '';
						substr($final_alignments{$target}->{$seqname}, $contigs->{$seqname}, $start_diff) = $start_intron;

						# Update indices
						$intron->{T_START} += $start_diff;
						$contigs->{$seqname} += $start_diff;
						$intron->{REF_START} += $start_diff;

						$new_max_start = $contigs->{$seqname} if (!defined($new_max_start) || $contigs->{$seqname} > $new_max_start);
					}
				}
			}
			undef(%sites);
		}
		# Unaligned exon matches with 'intron' end
		else {

			# Remove the sequence from intron and add to exon
			foreach my $site (@cluster) {
				next if ($sites{$site});

				# Fix intron ends
				if (exists($ends->{$site})) {
					foreach my $intron (@{$ends->{$site}}) {
						my $seqname = $intron->{SEQNAME};
						my $tend = $intron->{T_END};

						# Extract the misaligned exon
						my $end_diff = abs($max_offset);
						$end_intron = substr($intron->{SEQ}, length($intron->{SEQ}) - $end_diff, $end_diff);

						# Remove misaligned exon and add to exon alignment
						substr($intron->{SEQ}, length($intron->{SEQ}) - $end_diff, $end_diff) = '';
						substr($final_alignments{$target}->{$seqname}, $contigs->{$seqname} - $end_diff, $end_diff) = $end_intron;

						# Update indices
						$intron->{T_END} -= $end_diff;
						$contigs->{$seqname} -= $end_diff;

						$new_min_end = $contigs->{$seqname} if (!defined($new_min_end) || $contigs->{$seqname} < $new_min_end);
					}
				}
			}
			undef(%sites);
		}
	}
	else {

		# Determine where gap was placed and exon was split
		if ($aligned_exon =~ /-+/) {

			# Where gap was placed
			my $start = $-[0];

			# Where unaligned exon was split
			my $exon_end = $start - $exon_length;

			# Fix the sequences in this cluster
			foreach my $site (@cluster) {
				next if ($sites{$site});

				# Fix intron starts
				if (exists($starts->{$site})) {
					foreach my $intron (@{$starts->{$site}}) {
						my $seqname = $intron->{SEQNAME};
						my $tstart = $intron->{T_START};

						# Extract the misaligned exon
						my $start_diff = $exon_end;
						$start_intron = substr($intron->{SEQ}, 0, $start_diff);

						# Remove misaligned exon and and to exon alignment
						substr($intron->{SEQ}, 0, $start_diff) = '';
						substr($final_alignments{$target}->{$seqname}, $contigs->{$seqname}, $start_diff) = $start_intron;

						# Update indices
						$intron->{T_START} += $start_diff;
						$contigs->{$seqname} += $start_diff;
						$intron->{REF_START} += $start_diff;

						$new_max_start = $contigs->{$seqname} if (!defined($new_max_start) || $contigs->{$seqname} > $new_max_start);
					}
				}

				# Fix intron ends
				if (exists($ends->{$site})) {
					foreach my $intron (@{$ends->{$site}}) {
						my $seqname = $intron->{SEQNAME};
						my $tend = $intron->{T_END};

						# Extract the misaligned exon
						my $end_diff = length($unaligned_exon) - $exon_end;
						$end_intron = substr($intron->{SEQ}, length($intron->{SEQ}) - $end_diff, $end_diff);

						# Remove misaligned exon and and to exon alignment
						substr($intron->{SEQ}, length($intron->{SEQ}) - $end_diff, $end_diff) = '';
						substr($final_alignments{$target}->{$seqname}, $contigs->{$seqname} - $end_diff, $end_diff) = $end_intron;

						# Update indices
						$intron->{T_END} -= $end_diff;
						$contigs->{$seqname} -= $end_diff;

						$new_min_end = $contigs->{$seqname} if (!defined($new_min_end) || $contigs->{$seqname} < $new_min_end);
					}
				}
			}
			undef(%sites);
		}
	}
	unlink("intron-overlap.fasta");
	unlink("intron-overlap.aln.fasta");

	return ($new_max_start, $new_min_end);
}

sub move_misaligned_introns {
	#my ($introns, $target, $target_seq) = @_;
	my ($introns, $potential_introns, $target, $target_seq) = @_;

#	print Dumper($introns->{3724}),"\n";
#	print Dumper($introns->{3727}),"\n";
#	print Dumper($introns->{3798}),"\n";
#	print Dumper($introns->{3801}),"\n";

	# Determine number of contigs which contain each intron
	my %counts;
	my %ambiguous_counts;
	foreach my $site (keys %{$introns}) {
		foreach my $intron (@{$introns->{$site}}){
			$counts{$site}++;

#			if ($intron->{SEQ} =~ /(^N+|N+$)/) {
#				$ambiguous_counts{$site} = length($1);	
#				print "$site: $ambiguous_counts{$site}\n";
#			}
			if ($intron->{SEQ} =~ /(^N+)/) {
				$ambiguous_counts{START}->{$site} = length($1);	
				#$ambiguous_counts{END}->{$site} = -1 * length($1);	
			}
			elsif ($intron->{SEQ} =~ /(N+$)/) {
				$ambiguous_counts{START}->{$site} = length($1);	
				#$ambiguous_counts{END}->{$site} = -1 * length($1);	
			}
		}
	}

	my %clusters;
	my %in_cluster;
	my @sites = sort { $a <=> $b } keys %{$introns};
	foreach my $index (0 .. $#sites - 1) {
		my $site1 = $sites[$index];

		my $site1_end_offset = 0; 
		my $site1_start_offset = 0;

		if (exists($ambiguous_counts{START}->{$site1})) {
			$site1_start_offset = $ambiguous_counts{START}->{$site1};
		}
		if (exists($ambiguous_counts{END}->{$site1})) {
			$site1_end_offset = $ambiguous_counts{END}->{$site1};
		}
		my $site1_offset = $site1_start_offset + $site1_end_offset;

		next if ($in_cluster{$site1});
		foreach my $index2 ($index + 1 .. $#sites) {
			my $site2 = $sites[$index2];

			my $site2_end_offset = 0; 
			my $site2_start_offset = 0;

			if (exists($ambiguous_counts{START}->{$site2})) {
				$site2_start_offset = $ambiguous_counts{START}->{$site2};
			}
			if (exists($ambiguous_counts{END}->{$site2})) {
				$site2_end_offset = $ambiguous_counts{END}->{$site2};
			}
			my $site2_offset = $site2_start_offset + $site2_end_offset;

			#if ($site2 - $site1 - $ambig_offset <= $min_exon_length) {
#			print "$site2 - $site2_offset - $site1 - $site1_offset <= $min_exon_length\n";
#			print "",($site2 - $site2_offset - $site1 - $site1_offset),"\n";
			if ($site2 - $site2_offset - $site1 - $site1_offset <= $min_exon_length) {
				if (!exists($in_cluster{$site1})) {
					push(@{$clusters{$site1}}, $site2);
					$in_cluster{$site2} = $site1;
				}
				else {
					push(@{$clusters{$in_cluster{$site1}}}, $site2);
					$in_cluster{$site2} = $site1;
				}
			}
			else {
				last;
			}
		}
	}

#	die;
#
#	# Cluster contigs based on distance
#	my %clusters;
#	my %in_cluster;
#	my @sites = sort { $a <=> $b } keys %{$introns};
#	foreach my $index (0 .. $#sites - 1) {
#
#		my $ambig_offset = 0; # TODO: make this better
#		my $site1 = $sites[$index];
#		next if ($in_cluster{$site1});
#		foreach my $index2 ($index + 1 .. $#sites) {
#
#			my $site2 = $sites[$index2];
#
#			#my $ambig_offset = 0;
#			if (exists($ambiguous_counts{$site2})) {
#				#$ambig_offset = $ambiguous_counts{$site2};
#				$ambig_offset += $ambiguous_counts{$site2};
#			}
#
#			#if ($site2 - $site1 <= $min_exon_length) {
#			if ($site2 - $site1 - $ambig_offset <= $min_exon_length) {
#				if (!exists($in_cluster{$site1})) {
#					push(@{$clusters{$site1}}, $site2);
#					$in_cluster{$site2} = $site1;
#				}
#				else {
#					push(@{$clusters{$in_cluster{$site1}}}, $site2);
#					$in_cluster{$site2} = $site1;
#				}
#			}
#			else {
#				last;
#			}
#		}
#	}
#	print Dumper(\%counts),"\n";
#	print Dumper(\%clusters),"\n";
#	die;

	# Move introns which occur too closely to the site with the highest frequency
	foreach my $cluster (keys %clusters) {
		unshift(@{$clusters{$cluster}}, $cluster);

		my $most_frequent;
		my $highest_frequency = 0;
		foreach my $site (@{$clusters{$cluster}}) {
			#print "$site: $counts{$site}\n";
			if ($counts{$site} > $highest_frequency) {
				$highest_frequency = $counts{$site};
				$most_frequent = $site;
			}
			elsif ($counts{$site} == $highest_frequency) {

				# Calculate average difference from canonical splice pattern for new site
				my $new_score_sum = 0;
				foreach my $intron (@{$introns->{$site}}) {
					my $seq = $intron->{SEQ};
					my $splice = substr($seq, 0, 2).substr($seq, length($seq) - 2, 2);
					my $score = get_canonical_splice_diff($splice);
					$new_score_sum += $score;
				}

				# Calculate average difference from canonical splice pattern for old site
				my $old_score_sum = 0;
				foreach my $intron (@{$introns->{$most_frequent}}) {
					my $seq = $intron->{SEQ};
					my $splice = substr($seq, 0, 2).substr($seq, length($seq) - 2, 2);
					my $score = get_canonical_splice_diff($splice);
					$old_score_sum += $score;
				}

				my $new_average = $new_score_sum / $highest_frequency;
				my $old_average = $old_score_sum / $highest_frequency;

				# Opt for site closest to canonical splice pattern
				if ($new_average < $old_average) {
					$most_frequent = $site;
				}
				elsif ($new_average == $old_average) {
					# TODO: figure out what to do when equal???
				}
			}
		}
		#print "most frequent for cluster $cluster: $most_frequent ($highest_frequency)\n";

		# Set cluster id to most frequently occuring site
		if ($most_frequent != $cluster) {
			$clusters{$most_frequent} = delete($clusters{$cluster});
		}
	}
#	print "clusters:\n";
#	print Dumper(\%clusters),"\n";
#	print Dumper(\%counts),"\n";
#	die;

	my %realign;
	# Move introns from cluster to correct site
	#foreach my $cluster (keys %clusters) {
	foreach my $cluster (sort { $b <=> $a } keys %clusters) {
		$troublesome{$target}++;
		print "\nIntron(s) occurred too closely to site $cluster\n";

		#foreach my $site (reverse(@{$clusters{$cluster}})) {
		foreach my $site (@{$clusters{$cluster}}) {
			#last if ($site == $cluster);
			next if ($site == $cluster);

			print "  moving introns at $site...\n";

			foreach my $intron (@{$introns->{$site}}) {
#				print Dumper(\$intron),"\n";

				my $diff = $site - $cluster;
				my $t_start = $intron->{T_START};
				my $t_end = $intron->{T_END};
				my $seqname = $intron->{SEQNAME};

				# Add sequence to realignment hash
				if (!exists($realign{$seqname})) {
					$realign{$seqname} = $contig_seqs{$seqname};
				}

				# Moving intron forwards
				if ($diff < 0) {
					my $misaligned_tail = substr($contig_seqs{$seqname}, $t_start, -1 * $diff);

					# Remove misaligned sequence from intron
					$intron->{SEQ} =~ s/^\Q$misaligned_tail\E//;
				}
				# Moving intron backwards
				else {
					my $misaligned_tail = substr($contig_seqs{$seqname}, $t_start - $diff, $diff);
					print "diff: ",$diff,"\n";
					print "tail: ",$misaligned_tail,"\n";
					print Dumper($intron),"\n";

					# Append misaligned sequence which occurred between the two sites to beginning
					$intron->{SEQ} = $misaligned_tail.$intron->{SEQ};
				}

				# Correct intron's start indices
				$intron->{T_START} -= $diff;
				$intron->{REF_START} -= $diff;

				# Move this intron to its proper location
				push(@{$introns->{$cluster}}, $intron);

#				# Remove intron from sequence in realignment hash
#				substr($realign{$seqname}, $intron->{T_START}, $intron->{T_END} - $intron->{T_START}) = '';
			}
			delete($introns->{$site});
		}
	}

	# Join introns which may have been split into multiple introns
	foreach my $site (sort { $b <=> $a } keys %{$introns}) {

		# Determine if the scenario exists for certain contigs
		my %contigs;
		foreach my $intron (sort { $b->{T_START} <=> $a->{T_START} } @{$introns->{$site}}) {
			my $seqname = $intron->{SEQNAME};
			$contigs{$seqname}++;
		}

		# Apply fix to effected contigs
		foreach my $contig (keys %contigs) {

			if ($contigs{$contig} > 1) {

				# Information we need for new intron
				my $t_end;
				my $t_start;
				my $ref_start;
				my $full_intron;

				# Determine information of true intron, delete old pieces
				foreach my $index (reverse(0 .. scalar(@{$introns->{$site}}) - 1)) {
					my $intron = @{$introns->{$site}}[$index];
					if ($intron->{SEQNAME} eq $contig) {
						$t_start = $intron->{T_START} if (!defined($t_start) || $intron->{T_START} < $t_start);	
						$t_end = $intron->{T_END} if (!defined($t_end) || $intron->{T_END} < $t_end);	
						splice(@{$introns->{$site}}, $index, 1);
					}
				}

				# Add corrected intron
				$full_intron = {SEQNAME => $contig, REF_START => $site, T_END => $t_end, T_START => $t_start, SEQ => substr($contig_seqs{$contig}, $t_start, $t_end - $t_start)};
				push(@{$introns->{$site}}, $full_intron);
			}
		}
	}

	# Delete all introns from sequences we are realigning to minimize chance of misalignment
	foreach my $site (sort { $b <=> $a } keys %{$introns}) {
		foreach my $intron (@{$introns->{$site}}) {
		#foreach my $intron (sort { $b->{T_START} <=> $a->{T_START} } @{$introns->{$site}}) {
			my $seqname = $intron->{SEQNAME};
			if (exists($realign{$seqname})) {

#				if ($seqname eq "a_digitata_rep_c8779") {
#					print Dumper($intron),"\n";
#					#print "removing ",$intron->{T_START}," to ",$intron->{T_END} - $intron->{T_START}, " length: ",length($realign{$seqname}),"\n";
#					print "removing ",$intron->{T_START}," to ",$intron->{T_END}, " length: ",length($realign{$seqname}),"\n";
#				}

				die "Error removing intron from $seqname, sequence most likely has incorrect second intron at site.\n" if (!defined(substr($realign{$seqname}, $intron->{T_START}, $intron->{T_END} - $intron->{T_START})));
				substr($realign{$seqname}, $intron->{T_START}, $intron->{T_END} - $intron->{T_START}) = '';
			}
		}
	}

	# Realign exons of misaligned sequences
	foreach my $seqname (keys %realign) {
		my $seq = $realign{$seqname};

		# Delete incorrect alignment we currently have
		delete($final_alignments{$target}->{$seqname});

		# Realign exonic sequence after alignment error
		# TODO: add large penalty to prevent additional introns?
		print "  Realigning exon...\n";
		align_contig({CONTIG_NAME => $seqname, CONTIG_SEQ => $seq, TARGET_NAME => $target, TARGET_SEQ => $target_seq, IGNORE_INTRONS => 1});

		# Delete contig's introns and potential introns if we could not successfully realign
		if (!exists($final_alignments{$target}->{$seqname})) {

			# Check all introns at all sites
			foreach my $site (keys %{$introns}) {
				foreach my $index (reverse(0 .. scalar(@{$introns->{$site}}) - 1))  {
					my $intron = (@{$introns->{$site}})[$index];
					my $intron_seqname = $intron->{SEQNAME};

					# Remove the intron associated with this sequence
					if ($intron_seqname eq $seqname) {
						splice(@{$introns->{$site}}, $index, 1);
					}
				}
				
				# Remove if we no longer have sequences
				if (scalar(@{$introns->{$site}}) == 0) {
					delete($introns->{$site});
				}
			}

			# Check all intron starts
			foreach my $site (keys %{$potential_introns->{STARTS}}) {
				foreach my $index (reverse(0 .. scalar(@{$potential_introns->{STARTS}->{$site}}) - 1))  {
					my $intron = (@{$potential_introns->{STARTS}->{$site}})[$index];
					my $intron_seqname = $intron->{SEQNAME};

					# Remove the intron associated with this sequence
					if ($intron_seqname eq $seqname) {
						splice(@{$potential_introns->{STARTS}->{$site}}, $index, 1);
					}
				}

				# Remove if we no longer have sequences
				if (scalar(@{$potential_introns->{STARTS}->{$site}}) == 0) {
					delete($potential_introns->{STARTS}->{$site});
				}
			}

			# Check all intron ends
			foreach my $site (keys %{$potential_introns->{ENDS}}) {
				foreach my $index (reverse(0 .. scalar(@{$potential_introns->{ENDS}->{$site}}) - 1))  {
					my $intron = (@{$potential_introns->{ENDS}->{$site}})[$index];
					my $intron_seqname = $intron->{SEQNAME};

					# Remove the intron associated with this sequence
					if ($intron_seqname eq $seqname) {
						splice(@{$potential_introns->{ENDS}->{$site}}, $index, 1);
					}
				}

				# Remove if we no longer have sequences
				if (scalar(@{$potential_introns->{ENDS}->{$site}}) == 0) {
					delete($potential_introns->{ENDS}->{$site});
				}
			}
		}
	}

	return;
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
			$taxon =~ s/-/_/g;
			$align{$taxon} .= uc($line);
		}
	}
	close($alignment_file);
	
	return %align;
}

sub rev_comp {
	my $seq = shift;
	
	#my %comp = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C');
	my %comp = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', '-' => '-');

	my $rev_comp;
	my $rev = reverse($seq);
	foreach my $index (0 .. length($rev) - 1) {
		my $char = substr($rev, $index, 1);

		if (exists($comp{$char})) {
			$rev_comp .= $comp{$char};
		}
		else {
			$rev_comp .= 'N';
		}
	}

	return $rev_comp;
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
			#foreach my $sequence (values %target_alignment) {
			foreach my $seqname (keys %target_alignment) {
				my $sequence = $target_alignment{$seqname};

				# Exclude target from consensus
				next if ($seqname eq $target); 

				# Get sequence's base at this site
				my $char = substr($sequence, $index, 1);
				die "$target\n" if (!defined($char));
				$site{$char}++;
			}
			
			# Calculate and append consensus base
			my $consensus = consense_site(\%site);
			die "Error consensing site $index in $target\n".Dumper(\%site)."\n" if ($consensus =~ /\d/ && $consensus == -1);
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
		if ($site->{$base} && $base =~ /[ACTG]/) { ## TODO: check affect
			return $base;
		}
	}

	# Return most frequent ambiguity code otherwise
	foreach my $base (@sorted_bases) {
		if ($site->{$base} && $base =~ /[A-Z]/) { ## TODO: check affect
			return $base;
		}
	}

#	# No site from a contig mapped to target here, use ambiguity code
#	return 'N';

	# Depending on whether or not we want to include unmapped bases in target...
	# Return nothing
	if ($exclude_unmapped) {
		return '';
	}
	# No site from a contig mapped to target here, use ambiguity code
	else {
		return 'N';
	}

	#return '';

	# Should never happen if target is included
	#die "Error consensing site:\n",Dumper($site),"\n";
	return -1;
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
		chomp(@percent_free_cpu = `top -b -n2 -d0.05 | grep "Cpu(s)"`);
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

sub get_canonical_splice_diff {
	my $query = shift;
	my $canonical = "GTAG";

	# Calculate total number of differences from canonical "GTAG"
	my $diff = 0;
	foreach my $char (0 .. length($query) - 1) {
		my $char1 = chop($query);
		my $char2 = chop($canonical);

		$diff++ if ($char1 ne $char2);
	}

	return $diff;
}
