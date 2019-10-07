library("ggplot2")
library("phangorn")
library("phyloseq")
library("ape")
library("Rcpp")
library("readxl")
library("dplyr")

# Setting ggplot2 theme
theme_set(theme_bw())

# I want to read in the file and just figure out what it looks like after
# unpacking it from the RDS file format

seqtab <- readRDS("C:/Users/hammera/Desktop/raber_550/data/reads_290190_seqtab_final.rds")
seqtab.taxa <- readRDS("C:/Users/hammera/Desktop/raber_550/data/reads_290190_taxa_final.rds")

# check that number of ASVs and taxa rows is the same
ASV_taxa_equal <- ncol(seqtab) == nrow(seqtab.taxa)

print(ASV_taxa_equal)

# read pdf metadata in as a data_frame
meta_data <- as.data.frame(read_excel("C:/Users/hammera/Desktop/raber_550/data/550_metadata.xlsx"))
rownames(seqtab) <- meta_data$RaberSID

# set data equal to phyloseq object parameters
samples = sample_data(meta_data)
ASV = otu_table(seqtab, taxa_are_rows = FALSE)
TAX = tax_table(seqtab.taxa)

# I had an issue getting the sample names the same and into phyloseq
sample_names(samples) <- sample_names(ASV)

# create a phyloseq object
physeq = phyloseq(ASV, TAX, samples)

# if you want to prune a weird sample, look here
physeq = subset_samples(physeq, sample_names(physeq) != "93")

# start visualizing richness
plot_richness(physeq, x="Box", measures=c("Observed", "Shannon", "Simpson"))

#look at bray curtis nmds
#phyloseq.prop = transform_sample_counts(physeq, function(ASV) ASV/sum(ASV))
#ord.nmds.bray <- ordinate(phyloseq.prop, method="NMDS", distance="bray")
#plot_ordination(phyloseq.prop, ord.nmds.bray, color="HCCD68")

top15 <- names(sort(taxa_sums(physeq), decreasing = TRUE))[1:5]
physeq.top15 <- transform_sample_counts(physeq, function(ASV) ASV/sum(ASV))
physeq.top15 <- prune_taxa(top15, physeq.top15)
plot_bar(physeq, fill="Phylum")

