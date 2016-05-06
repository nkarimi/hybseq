#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Cwd "abs_path";
use File::Copy;
use File::Path "remove_tree";
use Data::Dumper;

# TODO:
# add SIGINT handlers to delete in progress files
# allow additional bwa options?
# add time stamps to output
# figure out terminal output
# remove loops around most method calls, make %config only argument
# use nct to multithread haplotype caller

my $targets;
my %targets;
my $project_name;
my %superassembly;
my $superassembly;
my $membership_groups;
my $superassembly_dir;
my $superassembly_consensus;
my %superassembly_consensus;
my $free_cpus = get_free_cpus();
(my $script_dir = abs_path($0)) =~ s/(.*\/).*/$1/;

# Required executables
my $bwa = check_path_for_exec("bwa");
my $java = check_path_for_exec("java");
my $mira = check_path_for_exec("mira");
my $blat = check_path_for_exec("blat"); # for consensus sequences
my $raxml = check_path_for_exec("raxmlHPC");
my $samtools = check_path_for_exec("samtools");
my $bam2consensus = check_path_for_exec("bam2consensus");

# Required jar files
my $gatk = check_path_for_jar("GenomeAnalysisTK.jar");
my $hapcompass = check_path_for_jar("hapcompass.jar");
my $mark_dups = check_path_for_jar("MarkDuplicates.jar");
my $create_seq_dict = check_path_for_jar("CreateSequenceDictionary.jar");

# Check for required input
my $config = shift(@ARGV);
die "You must specify a configuration file.\n" if (!defined($config));
die "Could not locate '$config', perhaps you made a typo?\n" if (!-e $config);

# Parse config file
my %config = parse_config($config);

initialize_directory_structure(\%config);

# Run mira once per species
foreach my $species_name (keys %config) {
	my $species = $config{$species_name};
	my $assembly = $species->{ASSEMBLY};
	my $reference = $species->{REFERENCE};

	chdir("$project_name/$species_name");

	# Only run if user didn't specify an assembly or a reference
	if (!-e "$species_name.fasta" && !-e "$species_name.con.fasta") {
		$species->{ASSEMBLY} = run_mira($species_name, $species);	
	}
	else {
		$species->{ASSEMBLY} = "$project_name/$species_name/$species_name.fasta";	
	}
}

# Old script run behavior
## Create consensus sequences from each mira assembly
#foreach my $species_name (keys %config) {
#	my $species = $config{$species_name};
#	my $assembly = $species->{ASSEMBLY};
#	my $reference = $species->{REFERENCE};
#
#	chdir("$project_name/$species_name");
#
#	# Only run if user didn't specify a reference
#	#if ($reference eq '') {
#	if (!-e "$species_name.con.fasta") {
#		#print "$species_name\n";
#		print "Generating reference sequences for targets in '$targets' using contigs in '$assembly'...\n";
#		#my $return = system("$script_dir/gen-refs.pl $assembly -t $targets");
#		#my $return = system("$script_dir/gen-refs.pl.2 $assembly -t $targets");
#		my $return = system("$script_dir/gen-refs.pl $assembly -t $targets");
#		die "Error generating references '$return' for $assembly.\n" if ($return); 
#
#		print "Reference sequence generation complete.\n";
#
#		# return output filename to avoid this?
#		($species->{REFERENCE} = $assembly) =~ s/\.fa(sta)?$/.con.fasta/;
#		#($species->{REFERENCE} = "$project_name/$species_name/$assembly") =~ s/\.fa(sta)?$/.con.fasta/;
#	}
#	else {
#		($species->{REFERENCE} = $assembly) =~ s/\.fa(sta)?$/.con.fasta/;
#		#($species->{REFERENCE} = "$project_name/$species_name/$assembly") =~ s/\.fa(sta)?$/.con.fasta/;
#	}
#}

# Create mira superassembly
if (!-e $superassembly) {

	# Concat the mira assembly from each species
	foreach my $species_name (keys %config) {
		my $species = $config{$species_name};
		my $assembly = $species->{ASSEMBLY};
		my $reference = $species->{REFERENCE};

		# Add species' assembly to superassembly
		system("cat '$assembly' >> '$superassembly'");
		#($species->{REFERENCE} = $assembly) =~ s/\.fa(sta)?$/.con.fasta/;
	}
}
else {
	foreach my $species_name (keys %config) {
		my $species = $config{$species_name};
		my $assembly = $species->{ASSEMBLY};
		my $reference = $species->{REFERENCE};
		#($species->{REFERENCE} = $assembly) =~ s/\.fa(sta)?$/.con.fasta/;
	}
}

# ID target paralogs and their consensus
identify_paralogs_and_create_consensuses(2);

sub identify_paralogs_and_create_consensuses {
	
	# Number of times we will generate consensus sequences
	# TODO: Setup way to run this until we no longer have nonmonophyletic sequences
	my $max_iterations = shift;

	# Move into superassembly directory
	chdir("$superassembly_dir");

	# Determine which iteration we are on
	my $current_iteration = 1;
	while (-e "run$current_iteration") {
		$current_iteration++;
	}

	# This will need to be reworked
	return if ($current_iteration > $max_iterations);

	# Clean up any previous output so we can run this iteration
	foreach my $species_name (keys %config) {
		my $species = $config{$species_name};

		# Remove old consensus files
		print glob("$project_name/$species_name/$species_name.con.*"),"\n";
		unlink(glob("$project_name/$species_name/$species_name.con.*"));

		# Delete gmap/GATK output for rerun
		my %accessions = %{$species->{ACCESSIONS}};
		foreach my $accession_name (keys %accessions) {
			unlink(glob("$project_name/$species_name/$accession_name/*"));
			remove_tree("$project_name/$species_name/$accession_name/gmap");
		}
	}
	#die;

	# Check if an alignment from the superassembly exists for all targets
	my $generate_consensus;
	foreach my $target (keys %targets) {
		$generate_consensus++ if (!-e "$target.fasta");	
	}

	# Generate references for targets using the superassembly
	if ($generate_consensus) {

		# Run once to determine intron locations
		my $return = system("$script_dir/gen-refs.pl '$superassembly' -t '$targets' -c");
		die "Error generating references '$return' for '$superassembly'.\n" if ($return); 

		# Run again using previous consensus as targets improve contig coverage, don't allow introns this time
		$return = system("$script_dir/gen-refs.pl '$superassembly' -t '$superassembly_consensus' -c -e");
		die "Error generating references '$return' for '$superassembly'.\n" if ($return); 
	}

	# Parse superassembly into memory
	%superassembly = parse_fasta($superassembly);
	
	# Parse superassembly consensuses into memory
	%superassembly_consensus = parse_fasta($superassembly_consensus);
	
	# TODO: Create new %targets with these consensuses?

	# Run RAxML on each target alignment
	my %consensuses;
	foreach my $target (keys %targets) {
		my $target_alignment = $target.".fasta";

		#next if (!-e "$target_alignment"); # For debugging
		
		# Write target consensus to separate file
		open(my $out_fasta, ">", "$target.con.fasta");
		print {$out_fasta} ">$target\n";	
		print {$out_fasta} "$superassembly_consensus{$target}\n";	
		close($out_fasta);

		# Run RAxML rapid bootstrap
		#if (!-e "RAxML_bipartitions.$target") {
		if (!-e "RAxML_bestTree.$target") {

			# RAxML doesn't use hyperthreaded cores well
			my $raxml_cpus = int($free_cpus / 2);
			$raxml_cpus = 1 if ($raxml_cpus < 1);

			# TODO: add conversion to phylip to avoid issues with RAxML version < 8?
			# TODO: add raxml settings to global options, add to contig file
			# Bootstrapped RAxML
			#my $return = system("$raxml -f a -m GTRGAMMA -s '$target_alignment' -n '$target' -p 123 -x 321 -# 100 -T $raxml_cpus");
			my $return = system("$raxml -m GTRGAMMA -s '$target_alignment' -n '$target' -p 123 -T $raxml_cpus");
			die "Error running RAxML on '$target_alignment'.\n" if ($return); 

			# Clean up unneeded RAxML output files
			# Files to remove if a bootstrap is run
	#		unlink("RAxML_bestTree.$target", "RAxML_bipartitionsBranchLabels.$target", 
	#			   "RAxML_bootstrap.$target", "RAxML_info.$target");
			unlink("RAxML_info.$target", "RAxML_log.$target", 
				   "RAxML_parsimonyTree.$target", "RAxML_result.$target", "$target.reduced");
		}

		# Check if paralogs already exist for this target, otherwise use separate-contigs.r to split
		my @paralogs = glob($target."_paralog*.tre");
		if (!@paralogs) {

			# Cluster contigs into subtrees based on tree topology and membership groups
			my $return = system("$script_dir/separate-contigs.r 'RAxML_bestTree.$target' '$membership_groups' >/dev/null");
			die "Error separating contigs in 'RAxML_bestTree.$target'.\n" if ($return); 

			# Reassign values
			@paralogs = glob($target."_paralog*.tre");
		}

		# Determine contigs present in each paralog, and write their sequences to two separate files
		my %used_contigs;
		foreach my $paralog (@paralogs) {

			# Extract paralog name
			(my $paralog_name = $paralog) =~ s/\.tre$//;
			my $consensus = $paralog_name.".con.fasta";

			print $paralog_name," = paralog name\n";
			print $consensus," = consensus\n";
			print "$target.con.fasta = target\n";

			# Read tree into memory
			open(my $tree, "<", $paralog) || die "Could not open '$paralog': $!.\n";
			chomp(my $lines = <$tree>);
			close($tree);
	
			# Extract contigs present in tree
			my @contigs;
			while ($lines =~ /[\(,]([^\(\),]+?)(:[^\(\),]+)?(?=[\),])/g) {
				push(@contigs, $1) if ($1 !~ /^\d+$/);
			}
	
			# Write sequences of contigs present in paralog to separate file
			open(my $out_fasta, ">", "$paralog_name.fasta");
			foreach my $contig (@contigs) {
				print {$out_fasta} ">$contig\n";	
				print {$out_fasta} "$superassembly{$contig}\n";	
				$used_contigs{$contig}++;
			}
			close($out_fasta);
	
			# Run reference generation script
			# TODO: run using previous consensus?
			#$return = system("$script_dir/gen-refs.pl.2 '$paralog_name.fasta' -t '$targets'");
	#		$return = system("$script_dir/gen-refs.pl.2 '$paralog_name.fasta' -t '$target.con.fasta' -e >/dev/null");
			my $return = system("$script_dir/gen-refs.pl '$paralog_name.fasta' -t '$target.con.fasta' -e >/dev/null");
			die "Error generating references '$return' for $superassembly.\n" if ($return); 
			
			# Parse paralog consensus and add to hash of final consensuses
			# TODO: delete invalid consensuses from hash
			my %paralog_consensus = parse_fasta($consensus);
			foreach my $consensus (keys %paralog_consensus) {
				$consensuses{$paralog_name} = $paralog_consensus{$consensus};
			}

			# TODO: further clean up?
			#unlink($consensus);
			unlink("$paralog_name.fasta.psl");
		}
		#unlink("$target.con.fasta");

		# Output contigs which weren't grouped with a paralog
		# Not currently used so potentially can be removed along with %used_contigs
		open(my $out, ">", $target."_unused.fasta");
		my %contigs = parse_fasta($target_alignment);
		foreach my $contig (keys %contigs) {
			if (!exists($used_contigs{$contig})) {
				print {$out} ">$contig\n";
				print {$out} "$contigs{$contig}\n";
			}
		}
		close($out);
	}

	# Write consensus for each paralog to file
	foreach my $species_name (keys %config) {
		my $species = $config{$species_name};
		my $assembly = $species->{ASSEMBLY};
		(my $reference = $assembly) =~ s/\.fa(sta)?$/.con.fasta/;
		$species->{REFERENCE} = $reference;
	
		# Write consensuses to file
		if (!-e $reference) {
	
			# Complicated sorting of target names is required for gsnap to work
			open(my $reference_file, ">", $reference) || die "Could not open '$reference': $!.\n";
			foreach my $consensus (sort { my $new_a = $a;
			                              my $new_b = $b;
										  $new_a =~ s/^(\d+).*/$1/;
	                                      $new_b =~ s/^(\d+).*/$1/;
	
										  if ($new_a == $new_b) {
	                                      	$a cmp $b;
										  }
										  else {
	                                      	$new_a <=> $new_b; 
										  } 
	
										  } keys %consensuses) {
				print {$reference_file} ">$consensus\n";
				print {$reference_file} "$consensuses{$consensus}\n";
			}
			close($reference_file);
		}
	}

	# Map raw reads back to consensus sequences with (bwa, gsnap?)
	foreach my $species_name (keys %config) {
		my $species = $config{$species_name};
		my $reference = $species->{REFERENCE};

		# Map each accessions reads back to the references
		my %accessions = %{$species->{ACCESSIONS}};
		map_accessions_to_references(\%accessions, $species_name, $reference);
	}

	# For faster debugging, uncomment for actual runs
	# Run GATK best practices on each mapping
	foreach my $species_name (keys %config) {
		my $species = $config{$species_name};
		preprocess_species_bam_files($species, $species_name);
	}
	# For faster debugging, remove for actual runs
#	foreach my $species_name (keys %config) {
#		my $species = $config{$species_name};
#
#		my %accessions = %{$species->{ACCESSIONS}};
#		foreach my $accession (keys %accessions) {
#			my $bam = $accessions{$accession}->{BAM_ALIGN};
#
#			$accessions{$accession}->{BAM_ALIGN_RECAL} = abs_path($bam);
#			next;
#		}
#	}

	# Determine our current directory
	my $init_dir = abs_path(getcwd());

	# Create consensus for each paralog using bam2consensus
	my %final_align;
	my $total_accessions = 0;
	foreach my $species_name (keys %config) {
		my $species = $config{$species_name};

		# Get accession names
		my %accessions = %{$species->{ACCESSIONS}};

		# Create consensus for each accession
		foreach my $accession (keys %accessions) {
			$total_accessions++;

			# Go into accession's output directory	
			chdir("$project_name/$species_name/$accession");

			# For faster debugging, uncomment for actual runs
			# Get name of bam file
			my $recal_bam = $accession.".recal.bam";
			# For faster debugging, remove for actual runs
			#my $recal_bam = $accession.".bam";
			die "Recalibrated bam file does not exist for '$accession'.\n" if (!-e $recal_bam);
				
			# Create consensus
			print "Generating consensus with bam2consensus...\n";
			my @consensus = `$bam2consensus '$recal_bam'`;

			# Parse fasta formatted output
			my $taxon;
			my %align;
			foreach my $line (@consensus) {
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

			# Add sequences to %final_align, with accession id
			foreach my $seq (keys %align) {
				(my $target = $seq) =~ s/_paralog\d+(_cov.*)?$//;

				# Remove -'s from name for file compatibility issues
				(my $formatted_name = $accession) =~ s/-/_/g;
				#$final_align{$target}->{"$seq"."_$accession"} = $align{$seq};
				
				# Only include the sequence if it isn't entirely N's
				# TODO: how this is handled really needs to be looked since it removes tips from the tree
				if ($align{$seq} !~ /^N+$/) {
					$final_align{$target}->{"$seq"."_$formatted_name"} = $align{$seq};
				}
			}

			# Go back to initial working directory
			chdir($init_dir);
		}
	}

	# Create and move into directory which will store output for this iteration
	mkdir("run$current_iteration");
	chdir("run$current_iteration");

	# Output final alignments
	foreach my $target (keys %final_align) {
		open(my $out, ">", "$target.fasta");
		foreach my $seq (keys %{$final_align{$target}}) {
			print {$out} ">$seq\n";	
			print {$out} "$final_align{$target}->{$seq}\n";	
		}
		close($out);
	}

	# Potential binning directories
	mkdir("single");
	mkdir("multi-mono");
	mkdir("multi-non-mono"); # I want to make this "multi-moNO" so badly...

	# Check for monophyly of paralog tips (bin.pl)
	my @good_targets;
	my @aligns = glob("*.fasta");
	foreach my $align (@aligns) {
		my %align = parse_fasta($align);

		# Determine how many paralogs we have
		my %paralogs = map { local $_ = $_; $_ =~ s/.*(_paralog\d+_).*/$1/; $_ => 1 } keys %align;

		# Single-copy
		if (scalar(keys %paralogs) == 1) {
			system("mv '$align' single/");
		}
		# Multi-copy
		else {
			(my $target = $align) =~ s/\.fasta$//;

			if (!-e "RAxML_bestTree.$target") { # Probably not needed

				# Input file names
				my $align_file = "$target"."_paralog_check.fasta";
				my $align_file_reduced = "$target"."_paralog_check.fasta.reduced"; # RAxML may create this

				# Write the alignment to file
				open(my $fasta, ">", $align_file);
				foreach my $paralog (keys %align) {
					print {$fasta} ">$paralog\n";
					print {$fasta} "$align{$paralog}\n";
				}
				close($fasta);

				# RAxML doesn't use hyperthreaded cores well
				my $raxml_cpus = int($free_cpus / 2);
				$raxml_cpus = 1 if ($raxml_cpus < 1);

				# Run RAxML
				my $return = system("$raxml -m GTRGAMMA -s '$align_file' -n '$target' -p 123 -T $raxml_cpus >/dev/null");
				die "Error running RAxML on '$align'.\n" if ($return); 

				# Clean up
				unlink("RAxML_info.$target", "RAxML_log.$target", 
					   "RAxML_parsimonyTree.$target", "RAxML_result.$target", "$target.reduced");

				unlink("$align_file");
				unlink("$align_file_reduced");
			}

			# Check monophyly with R
			my $is_mono = system("$script_dir/check-paralog-monophyly.r 'RAxML_bestTree.$target' >/dev/null");

			# Move alignment into its corresponding directory
			if ($is_mono) {
				system("mv '$align' multi-mono/");	
			}
			else {
				system("mv '$align' multi-non-mono/");	
			}
		}
	}
	
	# Join unsupported paralogs into single paralog
	
	# Parse CSV containing which paralogs are monophyletic
	open(my $csv, "<", "counts.csv");
	chomp(my @csv = <$csv>);
	close($csv);

	# Move back into superassembly's directory
	chdir($init_dir);

	# Iterate through each line of the CSV
	my $bad_target_count = 0;
	foreach my $line (@csv) {
		my @line = split(/,/, $line);

		# Extract CSV values from line
		my ($target, $total, $monophyletic, $ids) = @line;
		$ids = '' if (!defined($ids));
		my @ids = split(" ", $ids);

		# Get the filenames of all paralogs
		my @all_paralogs = glob("$target"."_paralog*.fasta");
		@all_paralogs = grep { !/\.con\.fasta/ } @all_paralogs;

		# Remove supported paralogs from @bad_paralogs
		my @bad_paralogs = @all_paralogs;
		foreach my $index (reverse(0 .. $#bad_paralogs)) {
			my $paralog = $bad_paralogs[$index];

			# Remove from @bad_paralogs if $id matches a monophyletic id
			foreach my $id (@ids) {
				if ($paralog =~ /\Q$id\E/) {
					splice(@bad_paralogs, $index, 1);
				}
			}
		}

		print "@bad_paralogs\n";
		
		# Join contigs used to form bad paralogs into a single file
		my %bad_contigs;
		foreach my $paralog (@bad_paralogs) {
			#(my $paralog_id = $paralog) =~ s/.*(_paralog\d+)_?.*/$1/;
			(my $paralog_id = $paralog) =~ s/.*(_paralog\d+).*/$1/;

			my %align = parse_fasta($paralog);
			%bad_contigs = (%align, %bad_contigs);

			# Remove the old output files for this paralog
			unlink(glob("$target$paralog_id".".*"));
		}

		# Join the bad paralogs if we have any
		if (@bad_paralogs) {
			$bad_target_count++;

			# Output the joined contigs to a new file
			my $new_id;
			foreach my $paralog (sort { $a cmp $b } @all_paralogs) {
				if (!-e $paralog) {
					#($new_id = $paralog) =~ s/.*(_paralog\d+)_?.*/$1/;
					($new_id = $paralog) =~ s/.*(_paralog\d+).*/$1/;

					# Output sequences as FASTA format
					open(my $out, ">", "$target$new_id.fasta");
					foreach my $contig (keys %bad_contigs) {
						print {$out} ">$contig\n";
						print {$out} "$bad_contigs{$contig}\n";
					}
					close($out);

					# For integration into current pipeline, create dummy tree
					open($out, ">", "$target$new_id.tre");
					print {$out} "(".join(",", keys %bad_contigs).");\n";
					close($out);
					last;
				}
			}
			###die if (!$new_id);

			# Create a new consensus for these sequences
			my $return = system("$script_dir/gen-refs.pl '$target$new_id.fasta' -t '$target.con.fasta' -e >/dev/null");
		}
	}
	
	# Rerun if needed
	if ($bad_target_count == 0) {
		return;	
	}
	else {
		identify_paralogs_and_create_consensuses($max_iterations);
	}
}

## Map raw reads back to consensus sequences with (bwa, gsnap?)
#foreach my $species_name (keys %config) {
#	my $species = $config{$species_name};
#	my $reference = $species->{REFERENCE};
#
#	# Map each accessions reads back to the references
#	my %accessions = %{$species->{ACCESSIONS}};
#	map_accessions_to_references(\%accessions, $species_name, $reference);
#}
#
## Run GATK best practices on each mapping
#foreach my $species_name (keys %config) {
#	my $species = $config{$species_name};
#	preprocess_species_bam_files($species, $species_name);
#}

# Run GATK's haplotype caller
foreach my $species_name (keys %config) {
	my $species = $config{$species_name};
	call_species_haplotypes($species, $species_name);
}

# Haven't figured out stuff yet after here
die;

# Phase haplotypes with HapCompass
foreach my $species_name (keys %config) {
	my $species = $config{$species_name};
	phase_species_haplotypes($species, $species_name);
}

sub initialize_directory_structure {
	my $config = shift;

	my $init_dir = abs_path(getcwd());
	chdir($project_name);

	my %config = %{$config};

	# Iterate through each species
	foreach my $species_name (keys %config) {
		my $species = $config{$species_name};

		# Create an output directory for each species
		if (!-e $species_name) {
			mkdir($species_name) || die "Could not create directory for species '$species_name'.\n";
		}
		elsif (!-d $species_name) {
			die "File named '$species_name' already exists and is not a directory.\n";
		}

		# Create symlink for user specified assemblies
		if ($species->{ASSEMBLY} && !-e "$species_name/$species_name.fasta") {
			symlink($species->{ASSEMBLY}, "$species_name/$species_name.fasta") || die "Failed to create symlink: $!.\n";
			$species->{ASSEMBLY} = "$project_name/$species_name/$species_name.fasta";
		}
		elsif ($species->{ASSEMBLY} && -l "$species_name/$species_name.fasta") {
			$species->{ASSEMBLY} = "$project_name/$species_name/$species_name.fasta";
		}

		# Create symlink for user specified references
		if ($species->{REFERENCE} && !-e "$species_name/$species_name.con.fasta") {
			symlink($species->{REFERENCE}, "$species_name/$species_name.con.fasta") || die "Failed to create symlink: $!.\n";
			$species->{REFERENCE} = "$project_name/$species_name/$species_name.con.fasta";
		}
		elsif ($species->{REFERENCE} && -l "$species_name/$species_name.con.fasta") {
			$species->{REFERENCE} = "$project_name/$species_name/$species_name.con.fasta";
		}

		# Iterate through each accession of the species
		my %accessions = %{$species->{ACCESSIONS}};
		foreach my $accession_name (keys %accessions) {

			# Create an output directory for each accession
			if (!-e "$species_name/$accession_name") {
				mkdir("$species_name/$accession_name") || die "Could not create directory for accession '$accession_name'.\n";
			}
			elsif (!-d "$species_name/$accession_name") {
				die "File named '$species_name/$accession_name' already exists and is not a directory.\n";
			}
		}
	}

	# Return to starting directory
	chdir($init_dir);

	print "Output directory structure successfully initiated.\n\n";

	return;
}

sub run_mira {
	my ($species_name, $species) = @_;

	my %accessions = %{$species->{ACCESSIONS}};
	my $base_accession_name = (sort { $a cmp $b } keys %accessions)[0];
	my $base_accession = $accessions{$base_accession_name};

	# Get associated reads for base accession
	my %pe_reads = %{$base_accession->{PE_READS}};
	#my @se_reads = @{$base_accession->{SE_READS}};

	# Write configuration file for mira
	my $conf_file = "$base_accession_name.mira.conf";
	my $mira_output = "$base_accession_name.mira.out";
	open(my $mira_conf, ">", $conf_file);
	#print {$mira_conf} "project = $base_accession_name\n";
	print {$mira_conf} "project = $species_name\n";
	print {$mira_conf} "job = est,denovo,accurate\n";
	print {$mira_conf} "parameters = -GENERAL:number_of_threads=$free_cpus \\\n";
	print {$mira_conf} "             COMMON_SETTINGS -NW:cmrnl=warn -SK:mmhr=1\\\n";
	#print {$mira_conf} "             COMMON_SETTINGS -NW:cmrnl=no -SK:mmhr=1\\\n";
	#print {$mira_conf} "             COMMON_SETTINGS -NW:cmrnl=warn \\\n";
	print {$mira_conf} "             SOLEXA_SETTINGS -CL:pec \\\n";
	print {$mira_conf} "             -CO:fnicpst=yes \\\n\n";

	# PE reads
	my $read_group_count = 0;
	foreach my $read_pair (keys %pe_reads) {
		my $insert_size = $pe_reads{$read_pair};

		print {$mira_conf} "readgroup = pe_reads_$read_group_count\n";
		print {$mira_conf} "data = $read_pair\n";
		#print {$mira_conf} "template = 100 1000 autorefine\n";
		#print {$mira_conf} "template = $insert_size $insert_size autorefine\n";
		#print {$mira_conf} "segment_placement = FR\n";
		print {$mira_conf} "autopairing\n";
		print {$mira_conf} "technology = solexa\n\n";

		$read_group_count++;
	}

#	# SE reads
#	$read_group_count = 0;
#	foreach my $single_end (@se_reads) {
#		print {$mira_conf} "readgroup = se_reads_$read_group_count\n";
#		print {$mira_conf} "data = $single_end\n";
#		#print {$mira_conf} "template = 100 1000 autorefine\n";
#		print {$mira_conf} "technology = solexa\n\n";
#
#		$read_group_count++;
#	}
	close($mira_conf);

	# Run mira
	print "Running mira for '$species_name' using reads from accession '$base_accession_name'...\n";
	my $return = system("$mira $conf_file > $mira_output");
	die "ERROR: Mira returned '$return', check output in '$mira_output' for further details.\n" if ($return);
	print "Finished mira.\n";
	
	my $assembly = $species_name."_assembly/$species_name"."_d_results/$species_name"."_out.unpadded.fasta";
	move($assembly, getcwd()."/$species_name.fasta") or die "Moving of '$assembly' failed: $!.\n";
	$assembly = getcwd()."/$species_name.fasta";

	#my $assembly = getcwd()."/$species_name".".fasta";

	# Remove unneeded mira output
	remove_tree($species_name."_assembly");

	return $assembly;
}

sub map_accessions_to_references {
	my ($accessions, $species_name, $reference) = @_;
	#$reference = "../$reference";

	my $init_dir = abs_path(getcwd());

	my %accessions = %{$accessions};
	foreach my $accession_name (keys %accessions) {
		my $output = $accession_name;

		chdir("$project_name/$species_name/$accession_name");

		# Skip if this has been completed in a previous run
		if (-e "$accession_name.bam") {
			$accessions->{$accession_name}{BAM_ALIGN} = abs_path("$output.bam");
			next;
		}

		my $accession = $accessions{$accession_name};
		my %pe_reads = %{$accession->{PE_READS}};

		# If we have multiple paired end libraries, merge them into 2 files
		my $return;
		my @f_reads;
		my @r_reads;
		foreach my $read_pair (keys %pe_reads) {
			my @reads = split(/\s+/, $read_pair);
			push(@f_reads, $reads[0]);
			push(@r_reads, $reads[1]);

			# Join reads, unzip if needed
			if ($reads[0] =~ /\.gz$/) {
				my $return = system("zcat $reads[0] > bwa.R1.fastq");
				die "Error concatenating reads '$return'" if ($return);
			}
			else {
				$return = system("cat $reads[0] > bwa.R1.fastq");
				die "Error concatenating reads '$return'" if ($return);
			}

			# Join reads, unzip if needed
			if ($reads[1] =~ /\.gz$/) {
				$return = system("zcat $reads[1] > bwa.R2.fastq");
				die "Error concatenating reads '$return'" if ($return);
			}
			else {
				$return = system("cat $reads[1] > bwa.R2.fastq");
				die "Error concatenating reads '$return'" if ($return);
			}
		}

#		my $return = system("cat @f_reads > bwa.R1.fastq");
#		die "Error concatenating reads '$return'" if ($return);
#		$return = system("cat @r_reads > bwa.R2.fastq");
#		die "Error concatenating reads '$return'" if ($return);

		# gmap_build -d meow -D test GroverBaumWendelbaits.fa
		mkdir("gmap") if (!-e "gmap");
		$return = system("gmap_build -d '$accession_name' -D gmap '$reference'");
		die "Error running gmap_build for '$accession_name': '$return'" if ($return);
		
		# gsnap -A sam -N 1 --gunzip -d meow -D test -B 5 -t 40 /scratch/nstenz/hyb-seq-analysis/reads/assembly_with_transcript/MiSeq_AD1R4/E_S6_L001_R*.gz > gsnap.sam
		#   --read-group-id=STRING         Value to put into read-group id (RG-ID) field
		#     --read-group-name=STRING       Value to put into read-group name (RG-SM) field
		#       --read-group-library=STRING    Value to put into read-group library (RG-LB) field
		#         --read-group-platform=STRING   Value to put into read-group library (RG-PL) field
		$return = system("gsnap -A sam -d '$accession_name' -D gmap -B 5 -t $free_cpus --read-group-id=test --read-group-name=test --read-group-library=lib1 --read-group-platform=illumina bwa.R1.fastq bwa.R2.fastq > '$output.sam'");
		die "Error running gsnap for '$accession_name': '$return'" if ($return);


#		# Index reference file
#		$return = system("$bwa index '$reference'");
#		print "$bwa index $reference\n";
#		die "Error indexing reference fasta: '$return'.\n" if ($return);
#
#		#my $free_cpus = get_free_cpus();
#
#		# Run bwa mem
#		#$return = system("$bwa mem $reference bwa.R1.fastq bwa.R2.fastq $bwa_opts -M -t $threads -R '\@RG\\tID:test\\tSM:test\\tPL:illumina\\tLIB:lib1\\tPU:unit1' > $output.sam");
#		$return = system("$bwa mem '$reference' bwa.R1.fastq bwa.R2.fastq -M -t $free_cpus -R '\@RG\\tID:test\\tSM:test\\tPL:illumina\\tLIB:lib1\\tPU:unit1' > '$output.sam'");
#		die "Error running bwa mem: '$return'.\n" if ($return);


		# Convert bwa output to bam
		# Remove unmapped reads
		# TODO: include -@ number of threads
		#$return = system("$samtools view $output.sam -b -S > $output.bam");
		$return = system("$samtools view -@ $free_cpus -F4 '$output.sam' -b -S > '$output.bam'");
		die "Error converting bwa's sam output to bam: '$return'.\n" if ($return);

		# Sort bam output
		# TODO: include -@ number of threads, -m memory per thread
		$return = system("$samtools sort -@ $free_cpus '$output.bam' '$output'");
		die "Error sorting bam output: '$return'.\n" if ($return);

		# Index bam output
		$return = system("$samtools index '$output.bam'");
		die "Error error indexing bam output: '$return'.\n" if ($return);

		# Clean up unneeded intermediate files
		unlink("$output.sam", "$reference.amb", "$reference.ann", "$reference.pac", "$reference.bwt",
			   "$reference.sa", "bwa.R1.fastq", "bwa.R2.fastq");

		#$accessions->{$accession_name}{BAM_ALIGN} = "$output.bam";
		$accessions->{$accession_name}{BAM_ALIGN} = abs_path("$output.bam");
	}
	chdir($init_dir);

	return;
}

# bwa version
#sub map_accessions_to_references {
#	my ($accessions, $species_name, $reference) = @_;
#	#$reference = "../$reference";
#
#	my $init_dir = abs_path(getcwd());
#
#	my %accessions = %{$accessions};
#	foreach my $accession_name (keys %accessions) {
#		my $output = $accession_name;
#
#		chdir("$project_name/$species_name/$accession_name");
#
#		# Skip if this has been completed in a previous run
#		if (-e "$accession_name.bam") {
#			$accessions->{$accession_name}{BAM_ALIGN} = abs_path("$output.bam");
#			next;
#		}
#
#		my $accession = $accessions{$accession_name};
#		my %pe_reads = %{$accession->{PE_READS}};
#
#		# If we have multiple paired end libraries, merge them into 2 files
#		my @f_reads;
#		my @r_reads;
#		foreach my $read_pair (keys %pe_reads) {
#			my @reads = split(/\s+/, $read_pair);
#			push(@f_reads, $reads[0]);
#			push(@r_reads, $reads[1]);
#		}
#
#		my $return = system("cat @f_reads > bwa.R1.fastq");
#		die "Error concatenating reads '$return'" if ($return);
#		$return = system("cat @r_reads > bwa.R2.fastq");
#		die "Error concatenating reads '$return'" if ($return);
#
#		# Index reference file
#		$return = system("$bwa index $reference");
#		print "$bwa index $reference\n";
#		die "Error indexing reference fasta: '$return'.\n" if ($return);
#
#		#my $free_cpus = get_free_cpus();
#
#		# Run bwa mem
#		#$return = system("$bwa mem $reference bwa.R1.fastq bwa.R2.fastq $bwa_opts -M -t $threads -R '\@RG\\tID:test\\tSM:test\\tPL:illumina\\tLIB:lib1\\tPU:unit1' > $output.sam");
#		$return = system("$bwa mem '$reference' bwa.R1.fastq bwa.R2.fastq -M -t $free_cpus -R '\@RG\\tID:test\\tSM:test\\tPL:illumina\\tLIB:lib1\\tPU:unit1' > '$output.sam'");
#		die "Error running bwa mem: '$return'.\n" if ($return);
#
#
#		# Convert bwa output to bam
#		# Remove unmapped reads
#		# TODO: include -@ number of threads
#		#$return = system("$samtools view $output.sam -b -S > $output.bam");
#		$return = system("$samtools view -@ $free_cpus -F4 '$output.sam' -b -S > '$output.bam'");
#		die "Error converting bwa's sam output to bam: '$return'.\n" if ($return);
#
#		# Sort bam output
#		# TODO: include -@ number of threads, -m memory per thread
#		$return = system("$samtools sort -@ $free_cpus '$output.bam' '$output'");
#		die "Error sorting bam output: '$return'.\n" if ($return);
#
#		# Index bam output
#		$return = system("$samtools index '$output.bam'");
#		die "Error error indexing bam output: '$return'.\n" if ($return);
#
#		# Clean up unneeded intermediate files
#		unlink("$output.sam", "$reference.amb", "$reference.ann", "$reference.pac", "$reference.bwt",
#			   "$reference.sa", "bwa.R1.fastq", "bwa.R2.fastq");
#
#		#$accessions->{$accession_name}{BAM_ALIGN} = "$output.bam";
#		$accessions->{$accession_name}{BAM_ALIGN} = abs_path("$output.bam");
#	}
#	chdir($init_dir);
#
#	return;
#}

sub preprocess_species_bam_files {
	#my $species = shift;	
	my ($species, $species_name) = @_;

	my $ploidy = $species->{PLOIDY};
	my $reference = $species->{REFERENCE};
	my @known_variants = @{$species->{KNOWN_VARIANTS}};
	(my $dict = $reference) =~ s/\.fa(sta)?$/.dict/i;
	#(my $dict = $reference) =~ s/\.fa(sta)?$/.dict.fasta/i;

	my $init_dir = abs_path(getcwd());

	# Move into species' directory
	chdir("$project_name/$species_name/");

	# Create dictionary with Picard tools for this species
	my $return = system("$java -jar $create_seq_dict R=$reference O=$dict");
	die "Error creating dictionary for '$reference'.\n" if ($return);

	# Index the reference sequences for this species
	$return = system("$samtools faidx $reference");
	die "Error indexing reference '$reference'.\n" if ($return);

	my %accessions = %{$species->{ACCESSIONS}};
	foreach my $accession (keys %accessions) {
		my $bam = $accessions{$accession}->{BAM_ALIGN};

		chdir("$project_name/$species_name/$accession/");

		# Skip if this has been completed in a previous run
		if (-e "$accession.recal.bam") {
			$accessions{$accession}->{BAM_ALIGN_RECAL} = abs_path("$accession.recal.bam");
			next;
		}

		#TODO: clean up intermediate output

#		my $return;

		# Create dictionary for references sequences
##		if (!-e $dict) {

#			# Create dictionary with Picard tools
#			$return = system("$java -jar $create_seq_dict R=$reference O=$dict");
#			die "Error creating dictionary for '$reference'.\n" if ($return);
#
#			# Index the reference
#			$return = system("$samtools faidx $reference");
#			die "Error indexing reference '$reference'.\n" if ($return);
##		}

		# MarkDuplicates
		#system("java -jar ~/private/phyloPrograms/picard-tools-1.119/MarkDuplicates.jar REMOVE_DUPLICATES=true METRICS_FILE=$accession.dedup.txt INPUT=$bam OUTPUT=$accession.dedup.bam");
		system("$java -jar $mark_dups REMOVE_DUPLICATES=true METRICS_FILE=$accession.dedup.txt INPUT=$bam OUTPUT=$accession.dedup.bam");
		$return = system("$samtools index $accession.dedup.bam");
		die "Error indexing deduplicated bam file '$accession.dedup.bam'.\n" if ($return);

		# Realign introns first step (identify intervals)
		#$return = system("java -jar ~/private/phyloPrograms/GenomeAnalysisTK.jar \\
		$return = system("$java -jar $gatk \\
				-T RealignerTargetCreator \\
				-R $reference \\
				-I $accession.dedup.bam \\
				-o forIndelRealigner.intervals");
		die "Error identifying intervals for intron realignment for '$accession.dedup.bam'.\n" if ($return);

		# Realign introns second step, application
		#$return = system("java -jar ~/private/phyloPrograms/GenomeAnalysisTK.jar \\
		$return = system("$java -jar $gatk \\
				-T IndelRealigner \\
				-R $reference \\
				-I $accession.dedup.bam \\
				-targetIntervals forIndelRealigner.intervals \\
				-o $accession.realigned.bam");
		die "Error realigning introns for '$accession.dedup.bam'.\n" if ($return);

		$return = system("$samtools index $accession.realigned.bam");
		die "Error indexing intron-realigned bam file '$accession.dedup.bam'.\n" if ($return);

		# Recalibrate base qualities
		if (!@known_variants) {

			# Since we don't have known variants we will estimate by running HaplotypeCaller
			#$return = system("java -jar ~/private/phyloPrograms/GenomeAnalysisTK.jar \\
			$return = system("$java -jar $gatk \\
					-T HaplotypeCaller \\
					-R $reference \\
					-I $accession.realigned.bam \\
					-ploidy $ploidy \\
					--genotyping_mode DISCOVERY \\
					-stand_emit_conf 10 \\
					-stand_call_conf 30 \\
					-o $accession.uncal.vars.vcf");
			die "Error calling haplotypes in bam file '$accession.realigned.bam'.\n" if ($return);

			# Recalibrate base scores based on HaplotypeCaller output
			$return = system("$java -jar $gatk \\
					-T BaseRecalibrator \\
					-R $reference \\
					-I $accession.realigned.bam \\
					-knownSites $accession.uncal.vars.vcf \\
					-o recal_data.table");
			die "Error calculating base quality recalibrations for bam file '$accession.realigned.bam'.\n" if ($return);

			# Rewrite bam file to include recalibrated scores
			$return = system("$java -jar $gatk \\
					-T PrintReads \\
					-R $reference \\
					-I $accession.realigned.bam \\
					-BQSR recal_data.table \\
					-o $accession.recal.bam");
			die "Error recalibrating base qualities in bam file '$accession.realigned.bam'.\n" if ($return);

			# Index output
			$return = system("$samtools index $accession.recal.bam");
			die "Error indexing base quality recalibrated bam file '$accession.recal.bam'.\n" if ($return);
		}
		else {

			# Join file names into a single command for the program
			my $known_sites_opt;
			foreach my $known_sites (@known_variants) {
				$known_sites_opt .= "-knownSites $known_sites \\\n";
			}

			# Recalibrate base scores
			$return = system("$java -jar $gatk \\
					-T BaseRecalibrator \\
					-R $reference \\
					-I $accession.realigned.bam \\
					$known_sites_opt
					-o recal_data.table");
			die "Error calculating base quality recalibrations for bam file '$accession.realigned.bam'.\n" if ($return);

			# Rewrite bam file to include recalibrated scores
			$return = system("$java -jar $gatk \\
					-T PrintReads \\
					-R $reference \\
					-I $accession.realigned.bam \\
					-BQSR recal_data.table \\
					-o $accession.recal.bam");
			die "Error recalibrating base qualities in bam file '$accession.realigned.bam'.\n" if ($return);

			$return = system("$samtools index $accession.recal.bam");
			die "Error indexing base quality recalibrated bam file '$accession.recal.bam'.\n" if ($return);
		}
		$accessions{$accession}->{BAM_ALIGN_RECAL} = abs_path("$accession.recal.bam");
	}
	chdir($init_dir);

	return;
}

sub call_species_haplotypes {
	my ($species, $species_name) = @_;

	my $ploidy = $species->{PLOIDY};
	my $reference = $species->{REFERENCE};

	my $init_dir = abs_path(getcwd());

	my %accessions = %{$species->{ACCESSIONS}};
	foreach my $accession (keys %accessions) {
		my $bam = $accessions{$accession}->{BAM_ALIGN_RECAL};

		chdir("$project_name/$species_name/$accession/");

		# Skip if this has been completed in a previous run
		#next if (-e "$accession.vars.vcf");
		if (-e "$accession.vars.vcf") {
			$accessions{$accession}->{VARIANTS} = "$accession.vars.vcf";
			next;
		}

		# Run HaplotypeCaller on recalibrated bam file
		my $return = system("$java -jar $gatk \\
				-T HaplotypeCaller \\
				-R $reference \\
				-I $bam \\
				-ploidy $ploidy \\
				--genotyping_mode DISCOVERY \\
				-stand_emit_conf 10 \\
				-stand_call_conf 30 \\
				-o $accession.vars.vcf");
		$accessions{$accession}->{VARIANTS} = "$accession.vars.vcf";
		die "Error calling haplotypes in bam file '$bam'.\n" if ($return);
	}
	chdir($init_dir);

	return;
}

sub phase_species_haplotypes {
	my ($species, $species_name) = @_;

	my $ploidy = $species->{PLOIDY};
	my $reference = $species->{REFERENCE};

	my $init_dir = abs_path(getcwd());

	my %accessions = %{$species->{ACCESSIONS}};
	foreach my $accession (keys %accessions) {
		my $bam = $accessions{$accession}->{BAM_ALIGN_RECAL};
		my $vcf = $accessions{$accession}->{VARIANTS};

		chdir("$project_name/$species_name/$accession/");

		# Skip if this has been completed in a previous run
		#next if (-e "$accession.vars.vcf");
		#next if (-e "hc_MWER_solution.txt");
		if (-e "hc_MWER_solution.txt") {
			$accessions{$accession}->{PHASED_HAPS} = "hc_MWER_solution";
			next;
		}

		# Run HaplotypeCaller on recalibrated bam file
		# java -Xmx200g -jar ~/private/phyloPrograms/hapcompass_v0.7.7/hapcompass.jar  --bam a_digitata_0.recal.bam --vcf a_digitata_0.vars.vcf -o ./ -p 4
		my $return = system("$java -Xmx200g -jar $hapcompass --bam $bam --vcf $vcf -o ./ -p $ploidy");
		$accessions{$accession}->{PHASED_HAPS} = "hc_MWER_solution";
		die "Error phasing haplotypes in vcf file '$vcf'.\n" if ($return);
	}
	chdir($init_dir);

	return;
}

sub parse_config {
	my $file = shift;

	print "Reading config file '$file'...\n";

	# Slurp config file
	open(my $config, "<", $file) || die "Could not open '$file': $!.\n";
	chomp(my @lines = <$config>);
	close($config);

	# Preliminary parse to set variables and remove comments/whitespace
	my %vars;
	foreach my $line (@lines) {
		chomp($line);

		$line =~ s/#.*//; # remove comments
		$line =~ s/^\s+|\s+$//g; # remove trailing/leading whitespace

		# Define user variables
		if ($line =~ /^set\s+(\$\S+)\s*=\s*(.*)/i) {
			my ($name, $value) = ($1, $2);
			$vars{$name} = $value;
		}
	}

	# Replace variable names with their values
	foreach my $var (keys %vars) {
		@lines = map { s/\Q$var\E/$vars{$var}/gi; $_ } @lines;
	}

	# Fully parse input
	my %config;

	# TODO
	# Stores which index in the ACCESSIONS array a particular name is
	my %accession_index_map;

	my $species;
	my $accession;
	foreach my $index (0 .. $#lines) {
		my $line = $lines[$index];
		$. = $index + 1; # Manually set line number for helping with error output

		if ($line =~ /^targets\s*=\s*(.*)/i) {
			$targets = $1;

			if (!-e $targets) {
				die "ERROR: could not locate target file '$targets' specified at line $., perhaps you made a typo?\n";
				#print "ERROR: could not locate target file '$targets' specified at line $., perhaps you made a typo?\n"; # warning for debug
			}
			$targets = abs_path($targets);
			%targets = parse_fasta($targets);
			next;
		}

		if ($line =~ /^membership\s*=\s*(.*)/i) {
			$membership_groups = $1;

			# Quick check for proper format
			my $left_count = 0;
			while ($membership_groups =~ /\(/g) {
				$left_count++;
			}

			my $right_count = 0;
			while ($membership_groups =~ /\)/g) {
				$right_count++;
			}

			die "ERROR: incorrect membership group formatting, mismatching parentheses at line $.?\n" if ($left_count != $right_count);

			next;
		}

		if ($line =~ /^project\s*=\s*(.*)/i) {
			$project_name = $1;

			if (!-e $project_name) {
				mkdir($project_name) || die "ERROR: could not create directory for project name specified at line $.?\n";
			}
			elsif (!-d $project_name) {
				die "ERROR: project name specified at line $. already exists and is not a directory.\n";
			}

			$project_name = abs_path($project_name);

			$superassembly_dir = "$project_name/_superassembly";
			$superassembly = "$project_name/_superassembly/superassembly.fasta";
			$superassembly_consensus = "$project_name/_superassembly/superassembly.con.fasta";

			if (!-e $superassembly_dir) {
				mkdir($superassembly_dir) || die "ERROR: could not create superassembly directory for project name specified at line $..\n";
			}

			next;
		}

		# Skip blank and user variable defining lines
		next if ($line eq '' || $line =~ /^set\s+/i);

		# Characterize the accession
		if ($accession) {
			if ($line =~ /^pe_reads\s*=\s*(.*)/i) {
				my $reads = $1;

				if (!defined($accession)) {
					die "ERROR: no accession defined for pe_reads '$reads' at line $..\n";
				}

				my $insert_size;
				if ($reads =~ s/\[(\d+)\]$//) {
					$insert_size = $1;
				}
				if (!defined($insert_size)) {
					die "ERROR: no insert size or non-integer insert size specified at line $..\n";
				}

				my @reads = split(/\s+/, $reads);

				# Check that the proper number of files was specified, don't allow interleaved paired end reads
				if (scalar(@reads) > 2) {
					die "ERROR: more than two sets of paired end reads specified at line $..\n";
				}
				elsif (scalar(@reads) == 1) {
					die "ERROR: only one file specified for the reads of '$accession', if this file contains interleaved reads separate them.\n";
				}

				undef($reads);
				foreach my $read (@reads) {
					if (!-e $read) {
						die "Could not locate read file '$read' specified at line $., perhaps you made a typo?\n";
						#print "ERROR: could not locate read file '$read' specified at line $., perhaps you made a typo?\n"; # warning for debug
					}
#					elsif ($read =~ /(.*)\.gz$/) {
#						my $unzipped_read = $1;
#
#						print "Unzipping read '$read'...\n";
#						my $return = system("gunzip -c $read > $unzipped_read");
#						die "Error unzipping '$read': gunzip returned '$return'.\n" if ($return);
#
#						$read = $unzipped_read;
#					}

					$reads .= abs_path($read)." ";
				}
				chop($reads);

				#push(@{$config{$species}{ACCESSIONS}{$accession}->{PE_READS}}, $reads);
				$config{$species}{ACCESSIONS}{$accession}{PE_READS}{$reads} = $insert_size;
				next;
			}
#			elsif ($line =~ /^se_reads\s*=\s*(.*)/i) {
#				my $reads = $1;
#				my @reads = split(/\s+/, $reads);
#
#				if (!defined($accession)) {
#					die "ERROR: no accession defined for se_reads '$reads' at line $..\n";
#				}
#				if ($reads =~ /\[(\d+)\]$/) {
#					die "ERROR: insert size specified for se_reads '$reads' at line $..\n";
#				}
#
#				undef($reads);
#				foreach my $read (@reads) {
#					if (!-e $read) {
#						#die "Could not locate read file '$read' specified at line $., perhaps you made a typo?\n";
#						print "ERROR: could not locate read file '$read' specified at line $., perhaps you made a typo?\n"; # warning for debug
#					}
#					$reads .= abs_path($read)." ";
#				}
#				chop($reads);
#
#				$config{$species}{ACCESSIONS}{$accession}->{SE_READS} = \@reads;
#				next;
#			}
#			elsif ($line =~ /^insert_size\s*=\s*(.*)/i) {
#				my $insert_size = $1;
#				if ($insert_size !~ /^\d+$/) {
#					print $line,"\n";
#					print "'$insert_size'\n";
#					die "Error in config file line $.. Specified insert size is not an integer.\n";
#				}
#
#				$config{$species}{ACCESSIONS}{$accession}->{INSERT_SIZE} = $insert_size;
#			}
		}

		# Characterize the species
		if ($species) {
			if ($line =~ /^ploidy\s*=\s*(.*)/i) {
				my $ploidy = $1;
				if ($ploidy !~ /^\d+$/) {
					#die "Error in config file line $.. Specified ploidy '$ploidy' is not an integer.\n";
					die "ERROR: specified ploidy '$ploidy' is not an integer at line $..\n";
				}
				$config{$species}{PLOIDY} = $ploidy;
				next;
			}
			elsif ($line =~ /^assembly\s*=\s*(.*)/i) {
				my $assembly = $1;

				if (!-e $assembly) {
					die "ERROR: could not locate species assembly '$assembly' specified at line $., perhaps you made a typo?\n";
					#print "ERROR: could not locate species assembly '$assembly' specified at line $., perhaps you made a typo?\n"; # warning for debugging
				}
				$assembly = abs_path($assembly);

				$config{$species}{ASSEMBLY} = $assembly;
				next;
			}
			elsif ($line =~ /^reference\s*=\s*(.*)/i) {
				my $reference = $1;

				if (!-e $reference) {
					die "ERROR: could not locate reference sequences '$reference' specified at line $., perhaps you made a typo?\n";
					#print "ERROR: could not locate reference sequences '$reference' specified at line $., perhaps you made a typo?\n"; # warning for debugging
				}
				$reference = abs_path($reference);

				$config{$species}{REFERENCE} = $reference;
				next;
			}
			elsif ($line =~ /^known_variants\s*=\s*(.*)/i) {
				my $variants = $1;
				my @variants = split(/\s+/, $variants);

				if (!defined($species)) {
					die "ERROR: no species defined for known_variants '$variants' at line $..\n";
				}
				if (scalar(@variants) > 2) { # Not sure if this is actually a requirement
					die "ERROR: no more than two known variant files can be specified at line $..\n";
				}

				my @abs_path_variants;
				foreach my $variant (@variants) {
					if (!-e $variant) {
						die "Could not locate variant file '$variant' specified at line $., perhaps you made a typo?\n";
						#print "ERROR: could not locate variant file '$variant' specified at line $., perhaps you made a typo?\n"; # warning for debugging
					}

					push(@abs_path_variants, abs_path($variant));
				}
				#$config{$species}{KNOWN_VARIANTS} = \@variants;
				$config{$species}{KNOWN_VARIANTS} = \@abs_path_variants;
				next;
			}
		}

		# Accession definition
		if ($line =~ /^accession\s*=\s*(.*)/i) {
			$accession = $1;

			if (!defined($species)) {
				die "ERROR: no species defined for accession '$accession' at line $..\n";
			}
			if (exists($config{$species}{ACCESSIONS}->{$accession})) {
				die "ERROR: duplicate accession name '$accession' for species '$species' declared at line $..\n";
			}

			$config{$species}{ACCESSIONS}->{$accession} = {};
			#$config{$species}{ACCESSIONS}{$accession}->{SE_READS} = [];
			$config{$species}{ACCESSIONS}{$accession}->{PE_READS} = {};
			next;
		}

		# Species definition
		if ($line =~ /^species\s*=\s*(.*)/i) {
			$species = $1;

			if ($config{$species}) {
				die "ERROR: duplicate species name '$species' declared at line $..\n";
			}

			$config{$species}{PLOIDY} = 2;
			$config{$species}{ASSEMBLY} = '';
			$config{$species}{REFERENCE} = '';
			$config{$species}{ACCESSIONS} = {};
			$config{$species}{KNOWN_VARIANTS} = [];
			next;
		}

		# If we end up here user messed something up
		die "ERROR: unknown command setting '$line' at line $..\n";
	}

	# Check that our target file was defined
	if (!defined($targets)) {
		die "ERROR: target file containing hyb-seq baits was not defined.\n";
	}

	# Check that our project name was defined
	if (!defined($project_name)) {
		die "ERROR: project name was not defined.\n";
	}

	# Check for additional input errors
	foreach my $species (keys %config) {
		my %species = %{$config{$species}};

		my %accessions = %{$species{ACCESSIONS}};
		if (scalar(keys %accessions) == 0) {
			die "ERROR: species '$species' has no members.\n";
		}

		foreach my $accession (keys %accessions) {
			my %pe_reads = %{$accessions{$accession}->{PE_READS}};
			#my @se_reads = @{$accessions{$accession}->{SE_READS}};

			#my @all_reads = push(@se_reads, keys %pe_reads);

			#if (scalar(@se_reads) == 0 && scalar(keys %pe_reads) == 0) {
			if (scalar(keys %pe_reads) == 0) {
				die "ERROR: no reads were specified for '$accession'.\n";
			}
		}
	}

	print "Config file successfully parsed with no errors.\n\n";

	return %config;
}

sub check_path_for_exec {
	my $exec = shift;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec && -x $dir.$exec && !-d $dir.$exec);
	}

	die "Could not find the following executable: '$exec'. This script requires this program in your path.\n" if (!defined($exec_path));
	return $exec_path;
}

sub check_path_for_jar {
	my $jar = shift;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $jar_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$jar_path = $dir.$jar if (-e $dir.$jar);
	}

	die "Could not find the following jar file: '$jar'. This script requires this program in your path.\n" if (!defined($jar_path));
	return $jar_path;
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

	my $free_cpus = int($total_cpus * $percent_free_cpu / 100) + 1;

	if ($free_cpus == 0 || $free_cpus !~ /^\d+$/) {
		$free_cpus = 1; # assume that at least one cpu can be used
	}
	
	return $free_cpus;
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
			die "ERROR: Invalid symbol '-', present in taxon name '$taxon', file $filename, line $..\n" if ($taxon =~ /-/);
		}
		else {
			# Taxon sequence
			$align{$taxon} .= uc($line);
		}
	}
	close($alignment_file);
	
	return %align;
}
