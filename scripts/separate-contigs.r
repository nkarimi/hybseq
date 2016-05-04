#!/usr/bin/Rscript
library(ape);
library(phangorn);

max_output_clades = 2;
pat_dist_cutoff_percentile = 0.15;
terminal_branch_cutoff_percentile = 1;

# Parse command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
target_tree = args[1];
target_gene = gsub("RAxML_bestTree\\.", "", target_tree);

# Check that user gave a tree file
if (!exists("target_tree")) {
	cat("You must specify a tree file as input.\n");
	q("no", status = 1);
} else if (!file.exists(target_tree)) {
	cat(paste0("Could not locate '", target_tree, "' perhaps you made a typo?\n"));
	q("no", status = 1);
}

# Read tree in
tree = read.tree(target_tree);
#tree = midpoint(tree);
tree = unroot(tree);
plot(tree, cex = 0.3, type="phylogram", main="unmodified");

# Extract names of all species
species = gsub("_dn.*|_rep.*|_c\\d+|\\|(\\d+-\\d+)\\|$", "", tree$tip.label, perl=T);
species = unique(species);

species;

# Parse membership requirements
membership = args[2];
#membership_groups = grep("\\(.*?\\)", membership, perl=T, value=T);
membership_groups = unlist(strsplit(membership, "\\(|\\)", perl=T));
membership_groups = membership_groups[membership_groups != ""];

# Set default membership to all species if we have none
if (all(is.na(membership_groups))) {
	membership_groups = species;
}

membership_groups;

# Fetch tip labels
tips = tree$tip.label;

# Fetch terminal branch lengths
terms <- tree$edge[,2] <= Ntip(tree);
terminal_branches <- tree$edge.length[terms];

# Calculate terminal branch length cutoff
terminal_branch_cutoff_length = quantile(terminal_branches, terminal_branch_cutoff_percentile, names=F);
terminal_branch_cutoff_length;

# Calculate sd and median to help determine which branches to prune
terminal_branch_sd = sd(terminal_branches);
terminal_branch_median = median(terminal_branches);

## Recursively remove tips with terminal branch lengths greater than the threshold
#cat("removing tips with longest terminal branch lengths...\n");
#
#removed = 1;
#while (removed != 0) {
#
#	removed = 0;
#	terms <- tree$edge[,2] <= Ntip(tree);
#	terminal_branches <- tree$edge.length[terms];
#
#	for (node in length(tree$tip.label):1) {
#		contig = tree$tip.label[node];
#
#		#if (terminal_branches[node] >= terminal_branch_cutoff_length) {
#		if (terminal_branches[node] >= terminal_branch_median + 3 * terminal_branch_sd) {
#			cat("removing", contig, ", length", terminal_branches[node], "\n");
#			tree = drop.tip(tree, contig);
#			removed = removed + 1;
#		}
#	}
#}
#plot(tree, cex = 0.3, type="phylogram", main="long terminals removed");

#tree = midpoint(tree);

## Calculate patristic distances for the tree
#pat_dists = diag(vcv.phylo(tree));

## Calculate threshold patristic distance for removal
#pat_dist_cutoff_length = quantile(pat_dists, pat_dist_cutoff_percentile, names=F);

# Extract all subtrees of our tree
subtrees = subtrees(tree);

# Add in missed subtrees due to subtrees() not treating the tree as unrooted
#missed_trees;
missed_trees = vector("list", );
class(missed_trees) = "multiPhylo";
for (subtree in subtrees) {

	# Skip full tree
	if (all.equal.phylo(subtree, tree)) {
		next;
	}

#	# Determine which
#	tips_not_in_subtree = setdiff(tree$tip.label, subtree$tip.label);
#
#	if (length(tips_not_in_subtree) == 0) {
#		q();
#	}

#	print(subtree$tip.label);
#	print(tips_not_in_subtree);

#	mrca = getMRCA(tree, tips_not_in_subtree);
	#tip_tree = extract.clade(tree, mrca);
	
	# Remove this subtree's tips from the full tree
	tip_tree = drop.tip(tree, subtree$tip.label);

	# Check if this subtree exists in the current subtree list
	in_current_subtrees = 0;
	for (subtree in subtrees) {
		if (all.equal.phylo(tip_tree, subtree)) {
			in_current_subtrees = in_current_subtrees + 1;	
			break;
		}
	}

	# Add to list of missed trees if we don't have this subtree already
	if (!in_current_subtrees) {
		missed_trees[[length(missed_trees) + 1]] = tip_tree;	
	}
}
#missed_trees;
length(subtrees);
subtrees = c(subtrees, missed_trees);
length(subtrees);
#	q();

# Remove trees which don't include all species
cat(length(subtrees), "initial subtrees\n");
for (index in length(subtrees):1) {
	subtree = subtrees[[index]];
	subtree_species = gsub("_dn.*|_rep.*|_c\\d+|\\|(\\d+-\\d+)\\|$", "", subtree$tip.label, perl=T);
	subtree_species = unique(subtree_species);

#	# Subtree does not include all species, remove
#	if (length(subtree_species) != length(species)) {
#		subtrees[[index]] = NULL;
#	}
#	else {
#		print(index);
#	}

	# Check that a taxon from each membership group is present
	total_members = 0;
	for (membership_group in membership_groups) {
		taxa = unlist(strsplit(membership_group, ","));

		# Check whether a member of the group is in this subtree's species list
		has_member = 0;
		for (taxon in taxa) {
			if (taxon %in% subtree_species) {
				has_member = has_member + 1;
			}
		}
		total_members = total_members + has_member;
		
		# Remove if a member from this group is missing
		if (!has_member) {
			subtrees[[index]] = NULL;
			print("removing subtree due to lack of member(s)");
			break;
		}
	}
}
#cat(length(subtrees), "subtrees contain contigs from all species\n");
cat(length(subtrees), "subtrees contain contigs from all membership groups\n");

# Calculate smallest subtrees which don't include another subtree
subtrees_reduced = subtrees;
for (index in length(subtrees_reduced):1) {
	
	for (index2 in 1:length(subtrees_reduced)) {
		# Skip the same tree
		if (index2 == index) {
			next;
		}
	
		subtree1_tip_labels = subtrees_reduced[[index]]$tip.label;
		subtree2_tip_labels = subtrees_reduced[[index2]]$tip.label;
		
		shared_tips = intersect(subtree1_tip_labels, subtree2_tip_labels);

		if (length(shared_tips) == length(subtree1_tip_labels)) {
			#subtrees_reduced[[index]] = NULL;
			subtrees_reduced[[index2]] = NULL;
			cat("removing subtree\n");
			break;
		}
		else if (length(shared_tips) == length(subtree2_tip_labels)) {
			#subtrees_reduced[[index2]] = NULL;
			subtrees_reduced[[index]] = NULL;
			cat("removing subtree\n");
			break;
		}
#		else {
#			print(length(subtree1_tip_labels));
#			print(length(subtree2_tip_labels));
#			print(length(shared_tips));
#
#
#			cat("i think this should never happen\n");
#			q();
#		}

#		if (length(shared_tips) > 0 && length(subtree1_tip_labels) > length(subtree2_tip_labels)) {
#			subtrees_reduced[[index]] = NULL;
#			break;
#		} else if (length(shared_tips) > 0 && length(subtree1_tip_labels) < length(subtree2_tip_labels)) {
#			subtrees_reduced[[index2]] = NULL;
#			break;
#		}
#		else if (length(shared_tips) > 0) {
#			print(length(subtree1_tip_labels));
#			print(length(subtree2_tip_labels));
#			print(length(shared_tips));
#
#
#			#cat("i think this should never happen\n");
#			q();
#		}
	}
}
paste("removed", length(subtrees) - length(subtrees_reduced), "of", length(subtrees), "trees");

final_subtrees = subtrees_reduced;
final_coverages = character(length(final_subtrees));
for (index in 1:length(subtrees)) {
	subtree = subtrees[[index]];
	
	# Check that subtree doesn't include multiple reduced subtrees
	include_count = 0;
	for (reduced_subtree in subtrees_reduced) {

		# Calculate the number of shared tips between each subtree
		subtree1_tip_labels = subtree$tip.label;
		subtree2_tip_labels = reduced_subtree$tip.label;
		
		shared_tips = intersect(subtree1_tip_labels, subtree2_tip_labels);

		# Iterate counter if shared tips exist
		if (length(shared_tips) > 0) {
			include_count = include_count + 1;
		}
	}

	# Skip tree if it has multiple subtrees included
	if (include_count > 1) {
		next;
	}

	for (index in 1:length(final_subtrees)) {
		final_subtree = final_subtrees[[index]];

		# Calculate the number of shared tips between each subtree
		subtree1_tip_labels = subtree$tip.label;
		subtree2_tip_labels = final_subtree$tip.label;
		
		shared_tips = intersect(subtree1_tip_labels, subtree2_tip_labels);

		# Iterate counter if shared tips exist
		if (length(shared_tips) > 0) {
			
			# Calculate contig coverage of subtree we are checking
			subtree_full_range = numeric(0);
			subtree_coverage = gsub(".*\\|(\\d+-\\d+)\\|$", "\\1", subtree$tip.label, perl=T);
			for (range in subtree_coverage) {
				split_range = unlist(strsplit(range, "-"));
				#print(split_range);
				subtree_full_range = union(subtree_full_range, split_range[[1]]:split_range[[2]]);
			}

			# Calculate contig coverage of current final subtree
			final_subtree_full_range = numeric(0);
			final_subtree_coverage = gsub(".*\\|(\\d+-\\d+)\\|$", "\\1", final_subtree$tip.label, perl=T);
			for (range in final_subtree_coverage) {
				split_range = unlist(strsplit(range, "-"));
				final_subtree_full_range = union(final_subtree_full_range, split_range[[1]]:split_range[[2]]);
			}
			subtree_full_range = sort(subtree_full_range);
			final_subtree_full_range = sort(final_subtree_full_range);

#			cat("new:", length(subtree_full_range), "current:", length(final_subtree_full_range), "\n");

			# Subtree we are checking has better coverage, use this one instead
			if (length(subtree_full_range) > length(final_subtree_full_range)) {
				final_subtrees[[index]] = subtree;
				#final_coverages[index] = length(subtree_full_range);
				#final_coverages[index] = paste0(min(subtree_full_range), "-", max(subtree_full_range));
				
				diffs = c(1, diff(subtree_full_range));
				start_indices = c(1, which(diffs > 1));
				end_indices = c(start_indices - 1, length(subtree_full_range));
				coloned = paste(subtree_full_range[start_indices], subtree_full_range[end_indices], sep="-");
				range = paste0(coloned,collapse="_");
				final_coverages[index] = range;

				#final_coverages[index] = subtree_full_range;
			}
			else {
				#final_coverages[index] = length(final_subtree_full_range);
				#final_coverages[index] = paste0(min(final_subtree_full_range), "-", max(final_subtree_full_range));
				
				diffs = c(1, diff(final_subtree_full_range));
				start_indices = c(1, which(diffs > 1));
				end_indices = c(start_indices - 1, length(final_subtree_full_range));
				coloned = paste(final_subtree_full_range[start_indices], final_subtree_full_range[end_indices], sep="-");
				range = paste0(coloned,collapse="_");
				final_coverages[index] = range;

				#final_coverages[index] = final_subtree_full_range;
			}
#			print(range);

			break;
		}
	}
}
final_subtrees;
final_coverages;

# Output final subtrees
for (index in 1:length(final_subtrees)) {
	subtree = final_subtrees[[index]];

	# Remove contig coverage from tip labels
	subtree$tip.label = gsub("(.*)\\|\\d+-\\d+\\|$", "\\1", subtree$tip.label, perl=T);

	plot(subtree, cex = 0.3, type="phylogram", main=c("subtree", index));

	cat("subtree", index, "contains", length(subtree$tip.label), "tips\n");

	# Set file name and output tree
	filename = paste0(target_gene, "_paralog", index, ".tre");
	#filename = paste0(target_gene, "_paralog", index, "_cov", final_coverages[index], ".tre");
	print(paste0(filename, " coverage: ", final_coverages[index], "\n"));
	write.tree(subtree, file=filename);

#	# Stop output if we reach maximum
#	if (index == max_output_clades) {
#		break;
#	}
}
