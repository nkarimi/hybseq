#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Cwd "abs_path";
use File::Copy;
use File::Path "remove_tree";
use Data::Dumper;

# TODO:
# allow additional bwa options?
# add time stamps to output
# figure out terminal output
# remove loops around most method calls, make %config only argument
# use nct to multithread haplotype caller

my $targets;
my $project_name;
my $script_dir = abs_path($0);
my $free_cpus = get_free_cpus();

# Required executables
my $bwa = check_path_for_exec("bwa");
my $java = check_path_for_exec("java");
my $mira = check_path_for_exec("mira");
my $cap3 = check_path_for_exec("cap3"); # for consensus sequences
my $blat = check_path_for_exec("blat"); # for consensus sequences
my $samtools = check_path_for_exec("samtools");

# Required jar files
my $gatk = check_path_for_jar("GenomeAnalysisTK.jar");
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
	#if ($assembly eq '' && $reference eq '') {
	if (!-e "$species_name.fasta" && !-e "$species_name.con.fasta") {
		$species->{ASSEMBLY} = run_mira($species_name, $species);	
	}
	else {
		$species->{ASSEMBLY} = "$project_name/$species_name/$species_name.fasta";	
	}
}

# Create consensus sequences from each mira assembly
foreach my $species_name (keys %config) {
	my $species = $config{$species_name};
	my $assembly = $species->{ASSEMBLY};
	my $reference = $species->{REFERENCE};

	chdir("$project_name/$species_name");

	# Only run if user didn't specify a reference
	#if ($reference eq '') {
	if (!-e "$species_name.con.fasta") {
		#print "$species_name\n";
		print "Generating reference sequences for targets in '$targets' using contigs in '$assembly'...\n";
		my $return = system("$script_dir/gen-refs.pl $assembly -t $targets");
		die "Error generating references '$return' for $assembly.\n" if ($return); 

		print "Reference sequence generation complete.\n";

		# return output filename to avoid this?
		($species->{REFERENCE} = $assembly) =~ s/\.fa(sta)?$/.con.fasta/;
		#($species->{REFERENCE} = "$project_name/$species_name/$assembly") =~ s/\.fa(sta)?$/.con.fasta/;
	}
	else {
		($species->{REFERENCE} = $assembly) =~ s/\.fa(sta)?$/.con.fasta/;
		#($species->{REFERENCE} = "$project_name/$species_name/$assembly") =~ s/\.fa(sta)?$/.con.fasta/;
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

# Run GATK best practices on each mapping
foreach my $species_name (keys %config) {
	my $species = $config{$species_name};
	preprocess_species_bam_files($species, $species_name);
}

# Run GATK's haplotype caller
foreach my $species_name (keys %config) {
	my $species = $config{$species_name};
	call_species_haplotypes($species, $species_name);
}

# Phase haplotypes with (HapHunt?)


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
	#print {$mira_conf} "             COMMON_SETTINGS -NW:cmrnl=warn -SK:mmhr=1\\\n";
	#print {$mira_conf} "             COMMON_SETTINGS -NW:cmrnl=no -SK:mmhr=1\\\n";
	print {$mira_conf} "             COMMON_SETTINGS -NW:cmrnl=warn \\\n";
	print {$mira_conf} "             SOLEXA_SETTINGS -CL:pec \\\n\n";

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

#	# Run mira
#	print "Running mira for '$species_name' using reads from accession '$base_accession_name'...\n";
#	my $return = system("$mira $conf_file > $mira_output");
#	die "ERROR: Mira returned '$return', check output in '$mira_output' for further details.\n" if ($return);
#	print "Finished mira.\n";
#	
#	my $assembly = $species_name."_assembly/$species_name"."_d_results/$species_name"."_out.unpadded.fasta";
#	move($assembly, getcwd()."/$species_name.fasta") or die "Moving of '$assembly' failed: $!.\n";
#	$assembly = getcwd()."/$species_name.fasta";

	my $assembly = getcwd()."/$species_name".".fasta";

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
		my @f_reads;
		my @r_reads;
		foreach my $read_pair (keys %pe_reads) {
			my @reads = split(/\s+/, $read_pair);
			push(@f_reads, $reads[0]);
			push(@r_reads, $reads[1]);
		}

		my $return = system("cat @f_reads > bwa.R1.fastq");
		die "Error concatenating reads '$return'" if ($return);
		$return = system("cat @r_reads > bwa.R2.fastq");
		die "Error concatenating reads '$return'" if ($return);

		# Index reference file
		$return = system("$bwa index $reference");
		print "$bwa index $reference\n";
		die "Error indexing reference fasta: '$return'.\n" if ($return);

		#my $free_cpus = get_free_cpus();

		# Run bwa mem
		#$return = system("$bwa mem $reference bwa.R1.fastq bwa.R2.fastq $bwa_opts -M -t $threads -R '\@RG\\tID:test\\tSM:test\\tPL:illumina\\tLIB:lib1\\tPU:unit1' > $output.sam");
		$return = system("$bwa mem $reference bwa.R1.fastq bwa.R2.fastq -M -t $free_cpus -R '\@RG\\tID:test\\tSM:test\\tPL:illumina\\tLIB:lib1\\tPU:unit1' > $output.sam");
		die "Error running bwa mem: '$return'.\n" if ($return);

		# Convert bwa output to bam
		# TODO: include -@ number of threads
		$return = system("$samtools view $output.sam -b -S > $output.bam");
		die "Error converting bwa's sam output to bam: '$return'.\n" if ($return);

		# Sort bam output
		# TODO: include -@ number of threads, -m memory per thread
		$return = system("$samtools sort $output.bam $output");
		die "Error sorting bam output: '$return'.\n" if ($return);

		# Index bam output
		$return = system("$samtools index $output.bam");
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

sub preprocess_species_bam_files {
	#my $species = shift;	
	my ($species, $species_name) = @_;

	my $ploidy = $species->{PLOIDY};
	my $reference = $species->{REFERENCE};
	my @known_variants = @{$species->{KNOWN_VARIANTS}};
	(my $dict = $reference) =~ s/\.fa(sta)?$/.dict/i;

	my $init_dir = abs_path(getcwd());

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

		my $return;

		# Create dictionary for references sequences
		if (!-e $dict) {

			# Create dictionary with Picard tools
			#$return = system("java -jar ~/private/phyloPrograms/picard-tools-1.119/CreateSequenceDictionary.jar R=$reference O=$dict");
			$return = system("$java -jar $create_seq_dict R=$reference O=$dict");
			die "Error creating dictionary for '$reference'.\n" if ($return);

			# Index the reference
			$return = system("$samtools faidx $reference");
			die "Error indexing reference '$reference'.\n" if ($return);
		}

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
		next if (-e "$accession.vars.vcf");

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

	my $free_cpus = int($total_cpus * $percent_free_cpu / 100) + 1;

	if ($free_cpus == 0 || $free_cpus !~ /^\d+$/) {
		$free_cpus = 1; # assume that at least one cpu can be used
	}
	
	return $free_cpus;
}
