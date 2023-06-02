library(vegan)
library(phyloseq)
library(ggplot2)

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
physeq_ACE_bac = phyloseq(ACE_otu_bac, ACE_sampledata, ACE_tax_bac)
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
tax_fungal = as.data.frame(tax_table(physeq_ACE_f))

#filtering

physeq_pruned_bac = subset_taxa(physeq_ACE_bac, tax_table(physeq_ACE_bac) != "Eukaryota")
physeq_pruned_bac = subset_taxa(physeq_pruned_bac, tax_table(physeq_pruned_bac) != "Unassigned")

physeq_pruned_bac_2 = subset_taxa(physeq_pruned_bac, Order != "Chloroplast")


#rarefy

physeq_b_rarefy = rarefy_even_depth(physeq_pruned_bac_2, rngseed = 1, sample.size = 1*min(sample_sums(physeq_pruned_bac_2)), replace = F)


#transformations
physeq_transformed_log = transform_sample_counts(physeq_b_rarefy, function (x) log(x+1))
physeq_transformed_hell = physeq_b_rarefy
otu_table(physeq_transformed_hell) <-otu_table(decostand(otu_table(physeq_b_rarefy), method = "hellinger"), taxa_are_rows=TRUE)
physeq_transformed_rel = transform_sample_counts(physeq_b_rarefy, function (x) x/sum(x))

#check which transformations give best normalization
melt_log = psmelt(physeq_transformed_log)
shapiro.test(melt_log$Abundance[10000:14990])
hist(melt_log$Abundance)

melt_hel = psmelt(physeq_transformed_hell)
shapiro.test(melt_hel$Abundance[10000:14990])
hist(melt_hel$Abundance)

melt_rel = psmelt(physeq_transformed_rel)
shapiro.test(melt_rel$Abundance[10000:14990])
hist(melt_rel$Abundance)

melt_rarefy = psmelt(physeq_b_rarefy)
shapiro.test(melt_rarefy$Abundance[10000:14990])
hist(melt_rarefy$Abundance)


#change metadata with updated climate dataset
sample_data(physeq_transformed_log) = sample_data(metadata_new)

#run ordination
erie_PCoA_b <- ordinate(
  physeq = physeq_transformed_log, 
  method = "PCoA", 
  distance = "bray"
)

#set the colours

manual = c("Bouvetoya" = "#004586", "Possession" = "#FA9FB5", "Bartolome" ="#FFD320","Kerguelen" = "#17BECF", "Lauff" = "#4C7022", "Maher" = "#EFEDF5", "Siple" = "#ADDD8E", "Marion" = "black", "Peter I" = "#4B1F6F", "Scott" = "#FF950E", "South Georgia" = "#C5000B")


p1 = plot_ordination(
  physeq = physeq_transformed_log,
  ordination = erie_PCoA_b,
  color = "Island", axes = 1:2) + theme_bw() +scale_color_manual(values= manual)
p2 = p1 +  geom_point( size = 6) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
p3 = p2 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p3 = p3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
  summary(p3)
ggsave(filename = "PCOA_viral_pangenome_AMG.png", p3, width = 8, height = 6, dpi = 900)

erie_PCoA_un= phyloseq::distance(physeq_transformed_log, method = "bray")
summary(erie_PCoA_un)
sampledf <- data.frame(sample_data(physeq_transformed_log))
adonis2(erie_PCoA_un ~Group, data = sampledf, permutations = 1000)
summary(adonis)
beta <- betadisper(erie_PCoA_un, sampledf$Geography)
permutest(beta)

anosim(erie_PCoA_un, sampledf$Group, permutations = 1000)

manual2 = c("Sub-Antarctic" = "red", "Continental" ="skyblue", "Maritime" = "Purple")


p1 = plot_ordination(
  physeq = physeq_transformed_log,
  ordination = erie_PCoA_b,
  color = "Group", axes = 1:2) + theme_bw() +scale_color_manual(values= manual2)
p2 = p1 +  geom_point(size = 6) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
p3 = p2 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p3 = p3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3




#fungal rarefy

rarecurve(t(otu_table(physeq_ACE_f)), step = 1000, col = "blue", cex = 0.6)

rainbow(11)



physeq_f_rarefy = rarefy_even_depth(physeq_ACE_f, rngseed = 1, sample.size = 1*min(sample_sums(physeq_ACE_f)), replace = F)

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
p2 = p1 +  geom_point(size = 6) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
p3 = p2 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p3 = p3 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3



erie_PCoA_f_un= phyloseq::distance(physeq_transformed_f, method = "bray")
sampledf <- data.frame(sample_data(physeq_transformed_f))
adonis2(erie_PCoA_f_un ~ Island, data = sampledf, permutations = 1000)

beta <- betadisper(erie_PCoA_un, sampledf$Island)
permutest(beta)

TukeyHSD(beta)

anosim(erie_PCoA_f_un, sampledf$Group, permutations = 1000)


p1 = plot_ordination(
  physeq = physeq_transformed_f,
  ordination = erie_PCoA,
  color = "Group", axes = 1:2) + theme_bw() +scale_color_manual(values= manual2)
p2 = p1 +  geom_point(size = 6) + theme(axis.text.x = element_text(face="bold", size=10),axis.text.y = element_text(face="bold", size=10))
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
