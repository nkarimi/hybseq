#!/usr/bin/perl
use strict;
use warnings;
use File::Path qw(remove_tree);

my @input = @ARGV;
die "You must specify fasta files for input.\n" if (!@input);

# Accession names which we need to remove from alignments
my @accessions = qw(E_S6 mad_S5 A-kilima_S12 165_S2 20_S1 A-grand2_S7 Aza_S4 A-greg2_S8 A-per2_S9 A-rub7_S11 A-suar_S3 Pseudo_S10);

# Cleanup previous runs
remove_tree("converted_nexus_out");
mkdir("converted_nexus_out");

# Loop through each file
foreach my $file (@input) {

	# Get the target name
	(my $target = $file) =~ s/\.fasta$//;	

	# Parse the alignment
	my %align = parse_fasta($file);

	# Determine how many paralogs are in the alignment
	my %paralogs = map { local $_ = $_; $_ =~ s/.*(_paralog\d+_).*/$1/; $_ } keys %align;

	# Create a separate alignment for each paralog
	my %paralog_seqs;
	foreach my $paralog (keys %paralogs) {

		# Get the sequences which describe this paralog
		my @seqs = grep { /\Q$paralog\E/ } keys %align;
		
		# Find which accession the sequence belongs to then add it to the paralog alignment
		foreach my $seq (@seqs) {
			foreach my $accession (@accessions) {
				if ($seq =~ /_\Q$accession\E$/) {
					(my $new_seq_name = $accession) =~ s/-/_/g;
					$paralog_seqs{$paralog}->{$new_seq_name} = $align{$seq};	
					last;
				}
			}
		}
	}
	
	# Create a separate file for each paralog
	foreach my $paralog (keys %paralog_seqs) {
		(my $id = $paralog) =~ s/_paralog(\d+)_/_paralog$1/;

		# So we don't have _paralog1 for targets with just a single paralog
		if (scalar(keys %paralog_seqs) == 1) {
			$id = '';
		}

		# Convert alignment to nexus
		my $out_name = "converted_nexus_out/$target$id.nex";
		write_nexus({OUT => $out_name, ALIGN => $paralog_seqs{$paralog}});
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
			#$taxon =~ s/-/_/g;
			#die "ERROR: Invalid symbol '-', present in taxon name '$taxon', file $filename, line $..\n" if ($taxon =~ /-/);
		}
		else {
			# Taxon sequence
			$align{$taxon} .= uc($line);
		}
	}
	close($alignment_file);
	
	return %align;
}

sub write_nexus {
	my $settings = shift;

	my $out_name = $settings->{'OUT'};
	my %align = %{$settings->{'ALIGN'}};

	# Determine alignment info required in NEXUS format
	my $ntaxa = scalar(values %align);
	my $nchar = length((values %align)[0]);

	open(my $out, ">", $out_name);
	print {$out} "#NEXUS\nbegin data;\n dimensions ntax=$ntaxa nchar=$nchar;\n ";
	print {$out} "format datatype=dna gap=- missing=?;\n matrix\n";
	foreach my $taxon (sort {$a cmp $b} keys %align) {
		print {$out} "  $taxon $align{$taxon}\n";
	}
	print {$out} "\n ;\nend;\n\n";
	close($out);
}
