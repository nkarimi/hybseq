#!/usr/bin/Rscript
library(ape);

# Parse command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
target_tree = args[1];
target_gene = gsub("RAxML_bestTree\\.", "", target_tree);

# This won't really work with the current setup of hyb-seq.pl
## Check that user gave a tree file
#if (!exists("target_tree")) {
#	cat("You must specify a tree file as input.\n");
#	q("no", status = 0);
#} else if (!file.exists(target_tree)) {
#	cat(paste0("Could not locate '", target_tree, "' perhaps you made a typo?\n"));
#	q("no", status = 0);
#}

# Read tree in
tree = read.tree(target_tree);
tree = unroot(tree);

# Determine how many paralog ids we have
paralog_ids = unique(gsub(".*(_paralog\\d+_).*", "\\1", tree$tip.label, perl=T));

# Check that paralogs are monophyletic
monophyletic_paralogs = vector("character", );
for (paralog_id in paralog_ids) {
	paralog_tips = grep(paralog_id, tree$tip.label, value=T);
	monophyletic = is.monophyletic(tree, tips=paralog_tips);

	# Add paralog to list of monophyletic paralogs
	if (monophyletic) {
		monophyletic_paralogs = c(monophyletic_paralogs, paralog_id);
	}
}

# Write outptu to csv
output = data.frame(c(target_gene), length(paralog_ids), length(monophyletic_paralogs), paste(monophyletic_paralogs, collapse=" "));

# We want to append if the file exists
if (file.exists("counts.csv")) {
	write.table(output, file="counts.csv", append=T, sep=",", col.names=F, row.names=F, quote=F);
}
else {
	write.table(output, file="counts.csv", sep=",", col.names=F, row.names=F, quote=F);
}

                                  
# Determine return status, 1 = all monophyletic, 0 = some
status = 0;
# If all paralogs are monophyletic except 1 (hence -1) it's most likely a rooting issue, count all as monophyletic
#if (length(monophyletic_paralogs) == length(paralog_ids)) {
if (length(monophyletic_paralogs) >= length(paralog_ids) - 1) {
	status = 1;
}

q(save="no", status=status);
