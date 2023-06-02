library(ANCOMBC)
library(tibble)
library(phyloseq)



#Perform ANCOMBC at genus level
physeq_bact_genus <- tax_glom(physeq_b_rarefy, taxrank = "Genus",NArm = FALSE)
sample_data(physeq_bact_genus) = sample_data(physeq_transformed_log)

ancom_da <- ancombc(physeq_bact_genus, formula = "Group", group = "Group",
                    p_adj_method = "fdr", zero_cut = 0.90, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)



ancom_res_df <- data.frame(
  ASV = row.names(ancom_da$res_global),
  W = unlist(ancom_da$res_global$W),
  p_val = unlist(ancom_da$res_global$p_val),
  q_val = unlist(ancom_da$res_global$q_val),
  diff_abn = unlist(ancom_da$res_global$diff_abn))

Diff_rep <- ancom_res_df[ancom_res_df$q_val < 0.05,]


W = as.data.frame(ancom_da$res$W)
names(W) = 'W'
p_val = as.data.frame(ancom_da$res$p_val)
names(p_val) = 'p_val'
q_val = as.data.frame(ancom_da$res$q_val)
names(q_val) = 'q_val'
diff_abn = as.data.frame(ancom_da$res$diff_abn)
names(diff_abn) = 'diff_abn'
length(t(W))

ancom_res_df = cbind.data.frame(W, p_val, q_val, diff_abn)
ancom_res_df <- tibble::rownames_to_column(ancom_res_df, "ASV")

Diff_rep <- ancom_res_df[ancom_res_df$q_val < 0.01,]


for (n in Diff_rep$ASV) {
  Diff_rep[Diff_rep$ASV == n,c("Phylum","Class","Order","Family","Genus")] <- as.data.frame(tax_table(physeq_bact_genus))[rownames(as.data.frame(tax_table(physeq_bact_genus)))==n,c("Phylum","Class","Order","Family","Genus")]
}

Diff_rep$ASV

#Matrix of the differentially abundant ASVs
physeq_inland_genus_proportion <- transform_sample_counts(physeq_bact_genus, function (x) x/sum(x))
physeq_bact_genus_pruned <- prune_taxa(Diff_rep$ASV,physeq_inland_genus_proportion)
sampleMatrix <- as.matrix(otu_table(physeq_bact_genus_pruned))

#Cluster the ASVs by abundance
sampleMatrix_scaled <- t(scale(t(sampleMatrix)))
sampleMatrix_scaled_dist <- dist(sampleMatrix_scaled)
sampleMatrix_scaled_dist_clustered <- hclust(sampleMatrix_scaled_dist)
groups <- cutree(sampleMatrix_scaled_dist_clustered, k=3)



#Add the clusters to the differentially abundant table.
Diff_rep$Cluster <- groups
#Name the clusters
Diff_rep$Cluster <- gsub("1","Continental",Diff_rep$Cluster)
Diff_rep$Cluster <- gsub("2","Maritime",Diff_rep$Cluster)
Diff_rep$Cluster <- gsub("3","Sub-Antarctic",Diff_rep$Cluster)


#Include the RA of each genus in the different clusters
for (x in 1:nrow(Diff_rep)) {
  Diff_rep[x,"RA_Sub-Antarctic"] <- mean(sampleMatrix[row.names(sampleMatrix)==Diff_rep[x,"ASV"],colnames(sampleMatrix)%in%row.names(physeq_bact_genus_pruned@sam_data[physeq_bact_genus_pruned@sam_data$Group=="Sub-Antarctic",])])
  Diff_rep[x,"RA_Marine"] <- mean(sampleMatrix[row.names(sampleMatrix)==Diff_rep[x,"ASV"],colnames(sampleMatrix)%in%row.names(physeq_bact_genus_pruned@sam_data[physeq_bact_genus_pruned@sam_data$Group=="Maritime",])])
  Diff_rep[x,"RA_Continental"] <- mean(sampleMatrix[row.names(sampleMatrix)==Diff_rep[x,"ASV"],colnames(sampleMatrix)%in%row.names(physeq_bact_genus_pruned@sam_data[physeq_bact_genus_pruned@sam_data$Group=="Continental",])])
}


write.csv(sampleMatrix, "AncomBC_comparison_genus_bact_3_groups.csv")
write.csv(Diff_rep, "AncomBC_comparison_genus_different_rep_bact_3_groups.csv")


#Perform ANCOMBC at genus level for fungi
physeq_fungi_genus <- tax_glom(physeq_f_rarefy, taxrank = "Genus",NArm = FALSE)
sample_data(physeq_fungi_genus) = sample_data(physeq_transformed_log)

ancom_da_f <- ancombc(physeq_fungi_genus, formula = "Group", group = "Group",
                    p_adj_method = "fdr", zero_cut = 0.90, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

ancom_res_df <- data.frame(
  ASV = row.names(ancom_da_f$res_global),
  W = unlist(ancom_da_f$res_global$W),
  p_val = unlist(ancom_da_f$res_global$p_val),
  q_val = unlist(ancom_da_f$res_global$q_val),
  diff_abn = unlist(ancom_da_f$res_global$diff_abn))

Diff_rep <- ancom_res_df[ancom_res_df$q_val < 0.05,]


W = as.data.frame(ancom_da_f$res$W)
names(W) = 'W'
p_val = as.data.frame(ancom_da_f$res$p_val)
names(p_val) = 'p_val'
q_val = as.data.frame(ancom_da_f$res$q_val)
names(q_val) = 'q_val'
diff_abn = as.data.frame(ancom_da_f$res$diff_abn)
names(diff_abn) = 'diff_abn'
length(t(W))

ancom_res_df = cbind.data.frame(W, p_val, q_val, diff_abn)
ancom_res_df <- tibble::rownames_to_column(ancom_res_df, "ASV")

Diff_rep <- ancom_res_df[ancom_res_df$q_val < 0.05,]


for (n in Diff_rep$ASV) {
  Diff_rep[Diff_rep$ASV == n,c("Phylum","Class","Order","Family","Genus")] <- as.data.frame(tax_table(physeq_fungi_genus))[rownames(as.data.frame(tax_table(physeq_fungi_genus)))==n,c("Phylum","Class","Order","Family","Genus")]
}

Diff_rep$ASV

#Matrix of the differentially abundant ASVs
physeq_inland_genus_proportion <- transform_sample_counts(physeq_fungi_genus, function (x) x/sum(x))
physeq_f_genus_pruned <- prune_taxa(Diff_rep$ASV,physeq_inland_genus_proportion)
sampleMatrix <- as.matrix(otu_table(physeq_f_genus_pruned))

#Cluster the ASVs by abundance
sampleMatrix_scaled <- t(scale(t(sampleMatrix)))
sampleMatrix_scaled_dist <- dist(sampleMatrix_scaled)
sampleMatrix_scaled_dist_clustered <- hclust(sampleMatrix_scaled_dist)
plot(sampleMatrix_scaled_dist_clustered)
groups <- cutree(sampleMatrix_scaled_dist_clustered, k=3)



#Add the clusters to the differentially abundant table.
Diff_rep$Cluster <- groups
#Name the clusters
Diff_rep$Cluster <- gsub("1","Maritime",Diff_rep$Cluster)
Diff_rep$Cluster <- gsub("2","sub-Antarctic",Diff_rep$Cluster)
Diff_rep$Cluster <- gsub("3","COntinental",Diff_rep$Cluster)

#Include the RA of each genus in the different clusters
for (x in 1:nrow(Diff_rep)) {
  Diff_rep[x,"RA_Sub-Antarctic"] <- mean(sampleMatrix[row.names(sampleMatrix)==Diff_rep[x,"ASV"],colnames(sampleMatrix)%in%row.names(physeq_f_genus_pruned@sam_data[physeq_f_genus_pruned@sam_data$Group=="Sub-Antarctic",])])
  Diff_rep[x,"RA_Marine"] <- mean(sampleMatrix[row.names(sampleMatrix)==Diff_rep[x,"ASV"],colnames(sampleMatrix)%in%row.names(physeq_f_genus_pruned@sam_data[physeq_f_genus_pruned@sam_data$Group=="Maritime",])])
  Diff_rep[x,"RA_Continental"] <- mean(sampleMatrix[row.names(sampleMatrix)==Diff_rep[x,"ASV"],colnames(sampleMatrix)%in%row.names(physeq_f_genus_pruned@sam_data[physeq_f_genus_pruned@sam_data$Group=="Continental",])])
  
}

#Check table of what is in clusters
table(Diff_rep$Cluster)
table(Diff_rep[Diff_rep$Cluster=="Maritime","Phylum"])
table(Diff_rep[Diff_rep$Cluster=="sub-Antarctic","Phylum"])
table(Diff_rep[Diff_rep$Cluster=="Control","Phylum"]) #Two orders contain the vast majority of what is selected for (1 class)

write.csv(sampleMatrix, "AncomBC_comparison_genus_fungi_3_groups.csv")
write.csv(Diff_rep, "AncomBC_comparison_genus_different_rep_fungi_3_groups.csv")



table_Ancom_Inland_Marble= read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")

table_Ancom_Inland_Marble = decostand(table_Ancom_Inland_Marble, method = "log")
DF = data.frame(Groups=rep(c("sub-Antarctic","Maritime","Continental"),c(13,5,17)))
rownames(DF) = colnames(table_Ancom_Inland_Marble)

row_dend = hclust(dist(table_Ancom_Inland_Marble), method = "ward.D")

library(pheatmap)
library(dendextend)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

Rowv  <- table_Ancom_Inland_Marble %>% scale %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 3) %>% set("branches_lwd", 1.2) %>%
  ladderize

heatmap(scale(table_Ancom_Inland_Marble), Rowv = Rowv, Colv = NA,
        scale = "row")

Diff_abundance_fig <- heatmap(as.matrix(table_Ancom_Inland_Marble), Colv = NA, Rowv = row_dend)
