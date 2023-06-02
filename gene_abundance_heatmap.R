library(phyloseq)
library(ggplot2)
library(igraph)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(pheatmap)


# First we need to import the tables
gene_table=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
taxonomy_gene = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",", quote = "" )
metadata_p = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")
metadata_p2 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")
write.csv(metadata_p, "metadata_gene_matrix.csv")

#Then we convert these tables into a phyloseq object (the tree object is not essential if we are not using unifrac)
gene_table_matrix = as.matrix(gene_table)
taxonomy_gene_matrix = as.matrix (taxonomy_gene)
ACE_gene = otu_table (gene_table_matrix, taxa_are_rows = TRUE)
ACE_tax = tax_table(taxonomy_gene_matrix)
ACE_gene_sampledata = sample_data(metadata_p)

physeq_genes = phyloseq(ACE_gene, ACE_tax, ACE_gene_sampledata)
sample_data(physeq_genes) = sample_data(metadata_p2)

manual = c("Bouvetoya" = "#004586", "Possession" = "#FA9FB5", "Bartolome" ="#FFD320","Kerguelen" = "#17BECF", "Lauff" = "#4C7022", "Maher" = "#EFEDF5", "Siple" = "#ADDD8E", "Marion" = "black", "Peter I" = "#4B1F6F", "Scott" = "#FF950E", "South Georgia" = "#C5000B")

Var1 = manual
anno_colors <- list(Var1 = Var1)

manual2 = c("Sub-Antarctic" = "red", "Continental" ="skyblue", "Maritime" = "Purple")


#Plot a PCoA ordination plot with samples colored according to geographical regions (the gene numbers are already normalized RPKM)

erie_PCoA_ordination <- ordinate(physeq = physeq_genes, method = "PCoA", distance = erie_PCoA)
PCoA = plot_ordination(physeq = physeq_genes, ordination = erie_PCoA_ordination , color = "Geography", axes = 1:2)  +scale_color_manual(values = manual2) + theme_bw()
PCoA2 = PCoA + geom_point(size = 6, alpha = 1/2) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
PCoA3  = PCoA2  + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
PCoA4  = PCoA3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +stat_ellipse()
PCoA4


#testing different groups
erie_PCoA = phyloseq::distance(physeq_genes, method = "bray")
sampledf <- data.frame(sample_data(physeq_genes))
sampledf_path <- data.frame(tax_table(physeq_genes))
adonis(erie_PCoA ~ sampledf$Geography, data = sampledf, permutations = 1000)
beta <- betadisper(erie_PCoA, sampledf$Soil_chem_group)
permutest(beta)



#Plot the heatmap from the gene abundance table

sampleDistMatrix <- as.matrix((otu_table(physeq_genes)))

pos_df = data.frame( "Island" = sampledf$Island, "Geography" = sampledf$Geography)
rownames(pos_df) = colnames(sampleDistMatrix)

color_df = data.frame("Pathways" = sampledf_path$Cycle)
rownames(color_df) = rownames(sampleDistMatrix)

pheatmap(sampleDistMatrix,clustering_distance_cols=erie_PCoA,cutree_cols = 5, cluster_rows = F, annotation = pos_df,  annotation_row = color_df, scale="row")
