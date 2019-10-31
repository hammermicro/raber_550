


## R Markdown - Load your libraries here and then proceed with running functions
library("ggplot2")
library("phangorn")
library("phyloseq")
library("ape")
library("Rcpp")
library("readxl")
library("dplyr")
library("vegan")
library("DESeq2")


# I want to read in the file and just figure out what it looks like after
# unpacking it from the RDS file format

seqtab <- readRDS("C:/Users/hammera/Desktop/raber_550/data/reads_290190_seqtab_final.rds")
seqtab.taxa <- readRDS("C:/Users/hammera/Desktop/raber_550/data/reads_290190_taxa_final.rds")
seqtab <- (seqtab)
# read excel metadata in as a data_frame
meta_data <- as.data.frame(read_excel("C:/Users/hammera/Desktop/raber_550/data/550_metadata.xlsx"))
meta_data_revised <- as.data.frame(read_excel("C:/Users/hammera/Desktop/raber_550/data/550_metadata_revised.xlsx"))
meta_data_revised$Treatment <- factor(meta_data_revised$Treatment)
meta_data_envfit <- as.data.frame(read_excel("C:/Users/hammera/Desktop/raber_550/data/550_metadata_nosex.xlsx"))
seqtab <- t(seqtab)
colnames(seqtab) <- meta_data$RaberSID

# set data equal to phyloseq object parameters
# make sure to change the sample data measure based on the metadata used
samples = sample_data(meta_data_revised)
ASV = otu_table(seqtab, taxa_are_rows = TRUE)
TAX = tax_table(seqtab.taxa)

# I had an issue getting the sample names the same and into phyloseq
sample_names(samples) <- sample_names(ASV)

# create a phyloseq object rename ASVs
physeq = phyloseq(ASV, TAX, samples)
dna <- Biostrings::DNAStringSet(taxa_names(physeq))
taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))

# sample 93 was wonky in terms of # of reads, and diversity, while 59 had some metadata issues
# I ran alpha diversity measures including both of these samples and without, but there was
# no statistically significant association between radiation and alpha diversity when they were
# included vs when they weren't
physeq = subset_samples(physeq, sample_names(physeq) != "93")
physeq = subset_samples(physeq, sample_names(physeq) != "59")
phyloseq_object <- physeq



#Does rarefying change the number of significant ASVs?
set.seed(1)
phyloseq_object <- physeq
phyloseq_oject <- rarefy_even_depth(physeq)
#If you want to agglomerate to a partcular level before analysis
phyloseq_object <- tax_glom(physeq, taxrank="Family")
### DESeq2, just trying to get this to run completely
# Some covariates to consider are "Sex", "Treatment", "AMDarkAverage", "OFTotDistDay1", "OFCentDur1", "NOObjpref", "FSTTimeImmobile", "FCCuedFrzTone", "Radiated"
sigtab <- NULL
treatdds = phyloseq_to_deseq2(phyloseq_object, ~Treatment*Sex)
treatdds.deseq = DESeq(treatdds)
res = results(treatdds.deseq, cooksCutoff=FALSE)
res = data.frame(res)
alpha = 0.01
remove_nas <- is.na(res$adj)
sigtab <- subset(res, remove_nas)
res <- subset(res, !is.na(res$padj))
sigtab <- subset(res, res$padj < alpha)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq1)[rownames(sigtab), ], "matrix"))
print(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
























