
## R Markdown
library("ggplot2")
library("phangorn")
library("phyloseq")
library("ape")
library("Rcpp")
library("readxl")
library("dplyr")
library("MASS")
library("tibble")

# Setting ggplot2 theme
theme_set(theme_bw())

# I want to read in the file and just figure out what it looks like after
# unpacking it from the RDS file format

seqtab <- readRDS("C:/Users/hammera/Desktop/raber_550/data/F_280_reads_seqtab_final.rds")
seqtab.taxa <- readRDS("C:/Users/hammera/Desktop/raber_550/data/F_280_reads_taxa_final.rds")

# check that number of ASVs and taxa rows is the same
ASV_taxa_equal <- ncol(seqtab) == nrow(seqtab.taxa)

print(ASV_taxa_equal)

# read pdf metadata in as a data_frame
meta_data <- as.data.frame(read_excel("C:/Users/hammera/Desktop/raber_550/data/550_metadata.xlsx"))
meta_data_revised <- as.data.frame(read_excel("C:/Users/hammera/Desktop/raber_550/data/550_metadata_revised.xlsx"))
meta_data_envfit <- as.data.frame(read_excel("C:/Users/hammera/Desktop/raber_550/data/550_metadata_nosex.xlsx"))
rownames(seqtab) <- meta_data$RaberSID

# set data equal to phyloseq object parameters
# make sure to change the sample data measure based on the metadata used
samples = sample_data(meta_data_revised)
ASV = otu_table(seqtab, taxa_are_rows = FALSE)
TAX = tax_table(seqtab.taxa)

# I had an issue getting the sample names the same and into phyloseq
sample_names(samples) <- sample_names(ASV)

# create a phyloseq object rename ASVs
physeq = phyloseq(ASV, TAX, samples)
dna <- Biostrings::DNAStringSet(taxa_names(physeq))
physeq <- merge_phyloseq(physeq, dna)
taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))

# sample 93 was wonky in terms of # of reads, and diversity, while 59 had some metadata issues
# I ran alpha diversity measures including both of these samples and without, but there was
# no statistically significant association between radiation and alpha diversity when they were
# included vs when they weren't
physeq = subset_samples(physeq, sample_names(physeq) != "93")
physeq = subset_samples(physeq, sample_names(physeq) != "59")

# set physeq1 equal to physeq
physeq1 <- physeq

# prepare data for LEFSE analysis
# this script convert phyloseq object into lefse recognized file format

# aggregate to genus level
ps_lefse <- physeq %>% tax_glom(taxrank = 'Genus', NArm = F)

# extract taxonomic data from phyloseq object and then stored in a vector called lefse_tax
lefse_tax <- ps_lefse %>% tax_table %>% data.frame(stringsAsFactors=FALSE)
lefse_tax <- replace(lefse_tax, is.na(lefse_tax), 'Unclassified')
lefse_tax <- lefse_tax %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% mutate(id = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = "|")) %>% ungroup %>% pull(id)
  
# extract otu abundance matrix from phyloeq object and annotated with tax information
lefse_matrix <- t(otu_table(ps_lefse) %>% data.frame(stringsAsFactors = F) %>% t %>% data.frame)
colnames(lefse_matrix) <- t(lefse_tax) 

# extract sample metadata and order sample same in lefse_matrix
lefse_sample <- sample_data(ps_lefse)
  
# convert factor in lefse_sample into character in order to combine sample and abundance data
lefse_sample_isfactor <- sapply(lefse_sample, is.factor)
lefse_sample[,lefse_sample_isfactor] <- lefse_sample[,lefse_sample_isfactor] %>% lapply(as.character)
lefse_sample <- lefse_sample %>% data.frame
lefse_matrix <- data.frame(lefse_matrix)
lefse_table <- full_join(rownames_to_column(lefse_sample), rownames_to_column(lefse_matrix), by = ("rowname" = "rowname")) %>% t

lefse_table <- data.frame(lefse_table)
View(lefse_table)

write.csv2(lefse_table, file="C:/Users/hammera/Desktop/raber_550/data/lefse.xlsx", sep="\t")
write.table(lefse_table, file = "C:/Users/hammera/Desktop/raber_550/data/lefse.txt", sep = "\t",row.names = TRUE, col.names = NA)

new_lefse <- read.table(file = "C:/Users/hammera/Desktop/raber_550/data/lefse.txt", sep="\t")

write.csv2(new_lefse, file="C:/Users/hammera/Desktop/raber_550/data/lefse.xlsx", na="NA")
lefse_dataset <- (read_excel("C:/Users/hammera/Desktop/raber_550/data/lefse_prepped_final.xlsx"))
write.table(lefse_dataset, file = "C:/Users/hammera/Desktop/raber_550/data/lefse_data.txt", sep = "\t",row.names = TRUE, col.names = NA)


