
library(phyloseq)
library(vegan)
library(proxy)
library(dendextend)
library(reshape2)
library(ggplot2)

#Using the transformed physeq object used to make the beat-diversity PCoA plot

hell.hclust = hclust(erie_PCoA_un, method="ward.D2")
plot(hell.hclust)
dend2 <- as.dendrogram(hell.hclust)
summary(dend2)

#if you have a physeq object, already trimmed
#cluster into phyla 
physeq_ACE_bac_RA = transform_sample_counts(physeq_b_rarefy, function (x) x/sum(x))
physeq_ACE_bac_phylum = tax_glom (physeq_ACE_bac_RA, taxrank = "Phylum")
#filtering out your rare taxa
physeq_dominant_3= filter_taxa(physeq_ACE_bac_phylum, function(x) mean(x) > 0.01, TRUE)
#export otu table and taxonomy
write.csv(otu_table(physeq_dominant_3), "otu_phylum_bact_pruned_revised.csv")
write.csv(tax_table(physeq_dominant_3), "tax_phylum_bact_pruned_revised.csv")
#in excel, change asv number for phyla in taxonomy file. 

#Import and format a table with phyla.
#First step, import table in which the first column has phyla names (after asvs)
otu_table_b_phylum=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
#melt it
melt_x = melt(otu_table_b_phylum)


#set up the labels in the formatted table to be identical to those in the dendogram
melt_x$variable = factor(melt_x$variable, levels = labels(dend2))

# Make the plots for the two graphs 
library(dendextend)
#plot dendogram
p3 <- ggplot(dend2, horiz = T)
#plot relative abundance table
p4 <- ggplot(melt_x, aes(y = value, x = variable, fill = Phylum )) + 
  geom_bar( stat="identity", position="fill") + coord_flip()+ scale_fill_manual(values=c("Cyanobacteria" = "green", "Acidobacteriota" = "orange", "Actinobacteriota" = "blue", "Bacteroidota" = "gold", "Chloroflexi" = "darkgreen", "Firmicutes" = "violetred1", "Gemmatimonadota" = "dark cyan", "Planctomycetota" = "yellow1", "Proteobacteria" = "red", "Myxococcota" = "aquamarine","Verrucomicrobiota" = "plum4"))

#Join the two plots together into a single image, and export as a metafile for fine tuning.
library(cowplot)
plot_grid(p3, p4, align = "h", rel_widths = c(3, 4))

#for fungi

hell.hclust_f = hclust(erie_PCoA_f_un, method="ward.D2")
plot(hell.hclust_f)
dend2_f <- as.dendrogram(hell.hclust_f)


#cluster into phyla 
physeq_ACE_bac_RA = transform_sample_counts(physeq_f_rarefy, function (x) x/sum(x))
physeq_ACE_bac_phylum = tax_glom (physeq_ACE_bac_RA, taxrank = "Phylum")
#filtering out your rare taxa
physeq_dominant_f= filter_taxa(physeq_ACE_bac_phylum, function(x) mean(x) > 0.01, TRUE)
#export otu table and taxonomy
write.csv(otu_table(physeq_dominant_f), "otu_phylum_fung_pruned_revised.csv")
write.csv(tax_table(physeq_dominant_f), "tax_phylum_fung_pruned_revised.csv")
#in excel, change asv number for phyla in taxonomy file. 

#Import and format a table with phyla.
#First step, import table in which the first column has phyla names (after asvs)
otu_table_f_phylum=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
#melt it
melt_f = melt(otu_table_f_phylum)

#Manual colors for different phya  
manual = c("Ascomycota" = "red2", "Unidentified" ="blue","Basidiomycota" = "green", "Chytridiomycota" = "gold", "Mortierellomycota" = "darkmagenta", "Rozellomycota" = "black")

#set up the labels in the formatted table to be identical to those in the dendogram
melt_f$variable = factor(melt_f$variable, levels = labels(dend2_f))

# Make the plots for the two graphs 
library(dendextend)
#plot dendogram
f3 <- ggplot(dend2_f, horiz = T)
#plot relative abundance table
f4 <- ggplot(melt_f, aes(y = value, x = variable, fill = Phylum )) + 
  geom_bar( stat="identity", position="fill") + coord_flip()+ scale_fill_manual(values=c("Ascomycota" = "red2", "Unidentified" ="blue","Basidiomycota" = "green", "Chytridiomycota" = "gold", "Mortierellomycota" = "darkmagenta", "Rozellomycota" = "black"))

#Join the two plots together into a single image, and export as a metafile for fine tuning.
library(cowplot)
plot_grid(f3, f4, align = "h", rel_widths = c(3, 4))
