#### Installing Packages

# install.packages("tidyverse")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")
# install.packages("readr")
# install.packages("seqinr")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("decontam")
# install.packages("ape")
# install.packages("vegan")
# install.packages("RColorBrewer")
# install.packages("remotes")
# remotes::install_github("microbiome/microbiome")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#   BiocManager::install("DESeq2")
  
#### Library Loaded 

library(tidyverse)
library(phyloseq)  
library(readr)  
library(seqinr)  
library(decontam)  
library(ape)  
library(vegan)  
library(RColorBrewer)  
library(remotes)  
library(DESeq2)  
library(microbiome)  

#### Importing the SRA Run Table 

SraRunTable <- read_delim("SraRunTableMod.csv", delim = ",")

#### Importing the Results from QIMME2 

count_table <- read_tsv(file="export/table/table.tsv", skip = 1)
count_table <- column_to_rownames(count_table, var = colnames(count_table)[1])
taxonomy <- read_tsv(file="export/taxonomy/taxonomy.tsv")
tree = read_tree("export/exported-tree/tree.nwk")
fasta <- read.fasta(file = "export/rep-seqs.fasta/dna-sequences.fasta")

#### Check Sequencing Depth with Rarefaction Curves

count_table_df <- as.data.frame(count_table)
rarecurve(t(count_table_df), step=100, cex=0.5, ylab="ASVs", label=T)

#### Removing Singletons 

count_table_no_singletons <- filter(count_table,rowSums(count_table)>1)

#### Modify Taxonomy Table 

head(taxonomy)
taxonomy_mod <-  taxonomy %>%
  mutate(taxonomy=str_replace_all(string=Taxon, pattern="D_\\d*\\__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species"), sep=";") %>%
  select (-Taxon, -Confidence) %>%
  column_to_rownames(var = 'Feature ID') 
head(taxonomy_mod)
head(taxa_names(TAX))

#### Put all into Phyloseq Object 

ASV =   otu_table(data.frame(count_table_no_singletons), taxa_are_rows =  TRUE)
TAX =   tax_table(as.matrix(taxonomy_mod))
META    =   sample_data(data.frame(SraRunTable, row.names = SraRunTable$"Library Name"))
head(taxa_names(TAX))
head(taxa_names(ASV))
head(taxa_names(tree))
head(sample_names(ASV))
head(sample_names(META))
ps <- phyloseq(ASV, TAX,META,tree)
rank_names(ps)
unique(tax_table(ps)[, "Domain"])
table(tax_table(ps)[, "Domain"], exclude = NULL)
ps <- subset_taxa(ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Eukaryota"))
table(tax_table(ps)[, "Domain"], exclude = NULL)
table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c(""))

pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape")
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }
my.tree <- phy_tree(ps)
out.group <- pick_new_outgroup(my.tree)

out.group

new.tree1 <- ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)
new.tree2 <- ape::multi2di(new.tree1)
phy_tree(ps) <- new.tree2
phy_tree(ps)

phylumGlommed = tax_glom(ps, "Phylum")
colourCount = length(table(tax_table(ps)[, "Phylum"], exclude = NULL))
getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
PhylaPalette = getPalette(colourCount)
plot_bar(phylumGlommed, x = "Sample", fill = "Phylum") + 
  scale_fill_manual(values = PhylaPalette)

ps_ra <- microbiome::transform(ps, transform = "compositional")
head(otu_table(ps_ra))

phylumGlommed_RA = tax_glom(ps_ra, "Phylum")
plot_bar(phylumGlommed_RA, x = "Sample", fill = "Phylum") + 
  scale_fill_manual(values = PhylaPalette)

#### The Proteobacteria 

ps_filtered <- subset_samples(ps_ra, !Library.Name %in% c("A8"))
ps_proteo_ra <- subset_taxa(ps_filtered, Phylum == "Proteobacteria")
colourCount = length(table(tax_table(ps_proteo_ra)[, "Order"], exclude = NULL))
getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
OrderPalette = getPalette(colourCount)
orderGlommed_RA = tax_glom(ps_proteo_ra, "Order")
plot_bar(orderGlommed_RA, fill = "Order") + 
  scale_fill_manual(values = OrderPalette)

#### The Actinobacteria 

ps_actino_ra <- subset_taxa(ps_filtered, Phylum == "Actinobacteria")
colourCount = length(table(tax_table(ps_actino_ra)[, "Order"], exclude = NULL))
getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
OrderPalette = getPalette(colourCount)
orderGlommed_RA = tax_glom(ps_actino_ra, "Order")
plot_bar(orderGlommed_RA,fill = "Order") + 
  scale_fill_manual(values = OrderPalette)

#### The Chlamydiae 

ps_chlamy_ra <- subset_taxa(ps_filtered, Phylum == "Chlamydiae")
colourCount = length(table(tax_table(ps_chlamy_ra)[, "Family"], exclude = NULL))
getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
FamilyPalette = getPalette(colourCount)
familyGlommed_RA = tax_glom(ps_chlamy_ra, "Family")
plot_bar(familyGlommed_RA, fill = "Family") + 
  scale_fill_manual(values = FamilyPalette)

#### Ordinations (First Transform)
ps_hellinger <- microbiome::transform(ps, transform = "hellinger")
head(otu_table(ps_hellinger))

out.pcoa <- ordinate(ps_hellinger, method = "PCoA", distance = "bray")

pcoa_plot = plot_ordination(ps_hellinger, out.pcoa, color ="Temperature", shape = "Treatment") +
  geom_point(size = 3) 

pcoa_plot

evals <- out.pcoa$values$Eigenvalues

pcoa_plot.scaled = plot_ordination(ps_hellinger, out.pcoa, color ="Temperature", shape = "Treatment") +
  geom_point(size = 3) +
  coord_fixed(sqrt(evals[2] / evals[1]))

pcoa_plot.scaled

ps_plot <- subset_samples(ps_hellinger, !Library.Name%in%c("A8")) 

out.pcoa <- ordinate(ps_plot, method = "PCoA", distance = "wunifrac")

wuf_pcoa_plot = plot_ordination(ps_plot, out.pcoa, color ="temperature", shape = "treatment") +
  geom_point(size = 3) +
  coord_fixed(sqrt(evals[2] / evals[1]))

wuf_pcoa_plot

out.nmds <- ordinate(ps_plot, method = "NMDS", distance = "bray")

nmds_plot = plot_ordination(ps_hellinger, out.nmds, color ="temperature", shape = "treatment") +
  geom_point(size = 3) 

nmds_plot

#### Test for Affects of Temperature 

ps_filtered <- subset_samples(ps, !Library.Name %in% c("A8"))

#### Differential Abundance with DeSeq2
otu_table(ps)+1
otu_table(ps) <- otu_table(ps)+1
ps_deseq <- phyloseq_to_deseq2(ps_filtered, ~temperature)
ps_deseq <- DESeq(ps_deseq)

#### Comparing 

##### First Compare Control (23°C) to (32°C)
deseq_res_temp_23_32 <- results(ps_deseq, alpha=0.01, contrast=c("temperature", "23", "32"))
summary(deseq_res_temp_23_32)

sigtab_deseq_res_temp_23_32 <- deseq_res_temp_23_32[which(deseq_res_temp_23_32$padj < 0.01), ]
sigtab_deseq_res_temp_23_32_with_tax <- cbind(as(sigtab_deseq_res_temp_23_32, "data.frame"), as(tax_table(ps_filtered)[row.names(sigtab_deseq_res_temp_23_32), ], "matrix"))
sigtab_deseq_res_temp_23_32_with_tax_ordered <- sigtab_deseq_res_temp_23_32_with_tax[order(sigtab_deseq_res_temp_23_32_with_tax$baseMean, decreasing=T), ]
sigtab_deseq_res_temp_23_32_with_tax_ordered

sig_taxa_23_32_1 <- row.names(sigtab_deseq_res_temp_23_32_with_tax_ordered)[1]
ind <- which(row.names(as(otu_table(ps_ra), "matrix")) == sig_taxa_23_32_1)
as(otu_table(ps_ra), "matrix")[ind,]

#### Significant Results 
sigtab_deseq_res_temp_23_32 <- deseq_res_temp_23_32[which(deseq_res_temp_23_32$padj < 0.01), ]
sigtab_deseq_res_temp_23_32_with_tax <- cbind(as(sigtab_deseq_res_temp_23_32, "data.frame"), as(tax_table(ps_filtered)[row.names(sigtab_deseq_res_temp_23_32), ], "matrix"))
sigtab_deseq_res_temp_23_32_with_tax[order(sigtab_deseq_res_temp_23_32_with_tax$baseMean, decreasing=T), ]

#### Plotting ASVs

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab_deseq_res_temp_23_32_with_tax$log2FoldChange, sigtab_deseq_res_temp_23_32_with_tax$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab_deseq_res_temp_23_32_with_tax$Family = factor(as.character(sigtab_deseq_res_temp_23_32_with_tax$Family), levels=names(x))
ggplot(sigtab_deseq_res_temp_23_32_with_tax, aes(x=Family, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#### Exporting 

descript_file<-writeLines(c("The phyloseq object in this file are data processed from Ramsby et.al 2018 (doi: 10.1111/mec.14544) NCBI BioProject # PRJNA407695. 
The ASV count table was produced using Qiime2 (v.2020.2), calling DADA2 for denoising, merging, and ASV inference, and a silva v132 
Naive Bayes Classifier for taxonomy calling.
Samples were collected from Little Pioneer Bay on Orpheus Island, Queensland, Australia in June 2015. Authors performed a simulated heat wave in mesocosms to test 
impacts on microbiome.
The samples were to see how the increase of temperature affected the mircrobiome community. There were 7 different types of temperature
but ran only 2 samples, a low temperature and a high temperature.
The samples corresponded from two different heat stress conditions (23˚C and 32˚C).
The metadata in this file come from the SraRunTable and have been modified to clearly indicate temperature and mortality data.
The otu_table object in the phyloseq object contains absolute abundaces of ASVs, without any tranformation or agglomeration at higher taxonomic levels.
qiime2 analysis file (done in Jupyter notebook) = PRJNA407695_QIIME_analysis.ipynb
processing of qiime2 and preliminary stats file (done in R) = seasponges.rmd"))

save(ps_plot, descript_file, file = "seasponges.RData")
