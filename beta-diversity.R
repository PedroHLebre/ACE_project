library(vegan)

#importing the bacterial data

otu_table_b=read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")
taxonomy_b = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

metadata_b = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

otu_table_ACE_bac = as.matrix(viral_pangenome)
taxonomy_ACE_bac = as.matrix (taxonomy_b)
ACE_tree_bac = read_tree(treefile = "./tree.nwk", errorIfNULL = FALSE)
ACE_otu_bac = otu_table (otu_table_ACE_bac, taxa_are_rows = TRUE)
ACE_tax_bac = tax_table(taxonomy_ACE_bac)
ACE_sampledata = sample_data(metadata_p)
physeq_ACE_bac = phyloseq(ACE_otu_bac, ACE_sampledata)
physeq_ACE_bac

#importing the fungal data


otu_table_f=read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")
taxonomy_f = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
metadata_f = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

otu_table_ACE_f = as.matrix(otu_table_f)
taxonomy_ACE_f = as.matrix (taxonomy_f)
ACE_tree_f = read_tree(treefile = "./fungal_tree.nwk", errorIfNULL = FALSE)
ACE_otu_f = otu_table (otu_table_ACE_f, taxa_are_rows = TRUE)
ACE_tax_f = tax_table(taxonomy_ACE_f)
taxa_names(ACE_otu_f)

physeq_ACE_f = phyloseq(ACE_otu_f, ACE_sampledata, ACE_tax_f)
physeq_ACE_f

#filtering

physeq_pruned_bac = subset_taxa(physeq_ACE_bac, tax_table(physeq_ACE_bac) != "Eukaryota")


#export tables
write.csv(otu_table(physeq_pruned_bac_20), "otu_prok_20_cuffoff.csv")
write.csv(tax_table(physeq_pruned_bac_20), "taxonomy_prok_20_cuffoff.csv")
write.csv(otu_table(physeq_pruned_fungal_20), "otu_fungal_20_cuffoff.csv")
write.csv(tax_table(physeq_pruned_fungal_20), "taxonomy_fungal_20_cuffoff.csv")

physeq_b_rarefy = rarefy_even_depth(physeq_pruned_bac, rngseed = 1, sample.size = 1*min(sample_sums(physeq_pruned_bac)), replace = F)




physeq_transformed_b = transform_sample_counts(physeq_b_rarefy, function (x) log(x+1))
physeq_transformed_b2 = transform_sample_counts(physeq_pruned_bac, function (x) x/sum(x))

sample_data(physeq_transformed_b) = sample_data(metadata_b)


sample_data(physeq_transformed_b) = sample_data(metadata_new)


erie_PCoA_b <- ordinate(
  physeq = physeq_transformed_b , 
  method = "PCoA", 
  distance = "bray"
)
summary (erie_PCoA)

manual = c("Bouvetøya" = "red2", "Possession" = "orange", "Bartolomé" ="blue","Kerguelen" = "green", "Lauff" = "skyblue4", "Maher" = "cyan", "Siple" = "palevioletred", "Marion" = "darkgreen", "Peter I" = "deeppink", "Scott" = "gold", "South Georgia" = "darkmagenta")


p1 = plot_ordination(
  physeq = physeq_transformed_b,
  ordination = erie_PCoA_b,
  color = "Island", axes = 1:2) + theme_bw() +scale_color_manual(values= manual)
p2 = p1 +  geom_point(size = 6, alpha = 1/2) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
p3 = p2 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p3 = p3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
summary(p3)
ggsave(filename = "PCOA_viral_pangenome_AMG.png", p3, width = 8, height = 6, dpi = 900)

erie_PCoA_un= phyloseq::distance(physeq_transformed_b, method = "bray")
summary(erie_PCoA_un)
sampledf <- data.frame(sample_data(physeq_transformed_b))
adonis2(erie_PCoA_un ~Geography, data = sampledf, permutations = 1000)
summary(adonis)
beta <- betadisper(erie_PCoA_un, sampledf$Geography)
permutest(beta)

manual2 = c("Sub-Antarctic" = "red", "Maritime" ="skyblue")


p1 = plot_ordination(
  physeq = physeq_transformed_b,
  ordination = erie_PCoA_b,
  color = "Geography", axes = 1:2) + theme_bw() +scale_color_manual(values= manual2)
p2 = p1 +  geom_point(size = 6, alpha = 1/2) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
p3 = p2 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p3 = p3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3




#fungal rarefy

rarecurve(t(otu_table(physeq_f_trimmed)), step = 1000, col = "blue", cex = 0.6)
physeq_f_trimmed = subset_samples ( physeq_ACE_f, sample_names(physeq_ACE_f) != "SG.17.111")
physeq_f_trimmed = subset_samples ( physeq_f_trimmed, sample_names(physeq_f_trimmed) != "CR.16.191")



physeq_f_rarefy = rarefy_even_depth(physeq_f_trimmed, rngseed = 1, sample.size = 1*min(sample_sums(physeq_f_trimmed)), replace = F)

physeq_transformed_f = transform_sample_counts(physeq_f_rarefy, function (x) log(x+1))

sample_data(physeq_transformed_f) = sample_data(metadata_new)

erie_PCoA <- ordinate(
  physeq = physeq_transformed_f, 
  method = "PCoA", 
  distance = "bray"
)

p1 = plot_ordination(
  physeq = physeq_transformed_f,
  ordination = erie_PCoA,
  color = "Island", axes = 1:2) + theme_bw() +scale_color_manual(values= manual)
p2 = p1 +  geom_point(size = 6, alpha = 1/2) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
p3 = p2 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p3 = p3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3

ggsave(filename = "PCOA_fungal_rarefy_bray.png", p3, width = 8, height = 6, dpi = 900)

erie_PCoA_un= phyloseq::distance(physeq_transformed_f, method = "bray")
sampledf <- data.frame(sample_data(physeq_transformed_f))
adonis(erie_PCoA_un ~ sampledf$Island, data = sampledf, permutations = 1000)

beta <- betadisper(erie_PCoA_un, sampledf$Island)
permutest(beta)

TukeyHSD(beta)


p1 = plot_ordination(
  physeq = physeq_transformed_f,
  ordination = erie_PCoA,
  color = "Geography", axes = 1:2) + theme_bw() +scale_color_manual(values= manual2)
p2 = p1 +  geom_point(size = 6, alpha = 1/2) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
p3 = p2 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p3 = p3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3



write.csv(table_b, "table_bacteria_transposed.csv")

table_b = as.data.frame(t(otu_table(physeq_pruned_bac_10)))
write.csv(table_b, "table_bacteria_transposed.csv")
write.csv(tax_table(physeq_pruned_bac_10), "tax_bacteria_transposed.csv")

table_b_transposed=read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")
metadata=read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")


(sim <- with(metadata, simper(table_b_transposed, Site)))
table_B_simper=summary(sim)

simper = as.data.frame(table_B_simper)

