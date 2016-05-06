#!/usr/bin/perl
use strict;
use warnings;

# Target sequence file
my $targets = '/scratch/nstenz/hyb-seq-analysis/GroverBaumWendelbaits.fa';
my %targets = parse_fasta($targets);

# Where reads are stored
my $read_dir = '/scratch/nstenz/hyb-seq-analysis/latest/subset/MiSeq_AD1R4';

# Read filenames
my @f_reads = glob("$read_dir/*_R1_*.fastq");
my @r_reads = glob("$read_dir/*_R2_*.fastq");

# Selected targets for testing
my @test_targets = (
	"2708|Gorai.011G225900.1|", #single copy after 2nd round
	"16415|Gorai.010G142300.1|", #failed 2nd round
	"15617|Gorai.012G082900.5|", #passed 2nd round
); 

# Create target subset sequence file
open(my $out, ">", "subset-targets.fasta");
foreach my $test_target (@test_targets) {
	print {$out} ">$test_target\n";
	print {$out} "$targets{$test_target}\n";

}
close($out);

# Map reads to targets
my @bam_files;
foreach my $f_read (@f_reads) {
	(my $r_read = $f_read) =~ s/_R1_/_R2_/;

	# Map reads to our chosen targets
	if (!-e "$f_read.bam") {
		system("bwa index subset-targets.fasta");
		system("bwa mem subset-targets.fasta '$f_read' '$r_read' -M -t 47 -R '\@RG\tID:test\tSM:test\tPL:illumina\tLIB:lib1\tPU:unit1' > '$f_read.sam'");
		system("samtools view -@ 47 '$f_read.sam' -b -S > '$f_read.bam'");
		system("samtools sort -@ 47 '$f_read.bam' '$f_read'");
		system("samtools index '$f_read.bam'");
	}
	push(@bam_files, "$f_read.bam");
}

# Get the reads we want
my %reduced_fastqs;
foreach my $test_target (@test_targets) {
	(my $id = $test_target) =~ s/(^\d+).*/$1/;
	print $test_target,"\n";

	# Extract reads which were mapped to each target
	foreach my $bam_file (@bam_files) {

		# Get the names of the original reads used
		(my $original_read1 = $bam_file) =~ s/\.bam$//;
		(my $original_read2 = $original_read1) =~ s/_R1_/_R2_/g;

		# Create names of the new output reads
		(my $filtered_read1 = $bam_file) =~ s/.*\/(.*)\.fastq\.bam$/$1.$id.fastq/;
		(my $filtered_read2 = $filtered_read1) =~ s/_R1_/_R2_/g;

		# Use samtools and filterbyname.sh to pick only the reads which mapped to our chosen targets
		system("samtools view '$bam_file' '$test_target' | cut -f1 > $id.reads.txt");
		system("filterbyname.sh in='$original_read1' out='$filtered_read1' names=$id.reads.txt include=t overwrite=t");
		system("filterbyname.sh in='$original_read2' out='$filtered_read2' names=$id.reads.txt include=t overwrite=t");

		# Store output filename
		push(@{$reduced_fastqs{$bam_file}->{R1}}, $filtered_read1);
		push(@{$reduced_fastqs{$bam_file}->{R2}}, $filtered_read2);

		# Clean up
		unlink("$id.reads.txt");
	}
}

# Join into a single fastq
foreach my $bam_file (keys %reduced_fastqs) {

	# Join each read direction (R1 and R2) separately
	foreach my $read_direction (keys %{$reduced_fastqs{$bam_file}}) {

		# Create name of the read file which contains reads for our chosen targets
		(my $new_read = $bam_file) =~ s/.*\/(.*)\.fastq\.bam$/$1.subset.fastq/;

		# Modify name if it is in reverse direction
		if ($read_direction eq "R2") {
			$new_read =~ s/_R1_/_R2_/;
		}

		# Get names of the files we want to join
		my @fastqs = @{$reduced_fastqs{$bam_file}->{$read_direction}};

		# Join the files into a single fastq
		print "cat @fastqs > '$new_read'\n";
		system("cat @fastqs > '$new_read'");
	}
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
