library(ggplot2)

b_alpha_diversity = estimate_richness(physeq_b_rarefy)
f_alpha_diversity= estimate_richness(physeq_pruned_fungal_10)

write.csv(b_alpha_diversity, "ba_alpha_rare.csv")
write.csv(f_alpha_diversity, "f_alpha.csv")

table_f_alpha = read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")

table_b_alpha = read.table (file.choose(), header = TRUE, row.names = 1 , sep = ",")
#basic boxplot
observed_compare = ggplot (table_b_alpha, aes( x= Island, y = Observed, fill = Soil_groups)) + geom_boxplot() + geom_point (aes(fill = Soil_groups), size = 3, position = position_jitterdodge())+theme_classic()
observed_compare2 = observed_compare + scale_fill_manual(values = c("Extreme" = "red", "Temperate" = "blue"))
observed_compare3 = observed_compare2 + theme(axis.text.x = element_text(size=10, angle=45), axis.text.y = element_text(face="bold",size = 14))
observed_compare4 = observed_compare3 + stat_compare_means(method = "wilcox.test", ref.group =".all.", label = "p.signif", hide.ns = TRUE)
ggsave(filename = "observed_alpha_f.png", observed_compare3, width = 8, height = 6, dpi = 900)

shapiro.test(table_f_alpha$Shannon)

kruskal.test(Shannon ~ Island, data=table_f_alpha)
pairwise.wilcox.test(table_f_alpha$Shannon, table_f_alpha$Island, p.adjust.method="fdr")

aov.observed.b = aov(Shannon ~ Island, data=table_f_alpha)
summary(aov.observed.b)

TukeyHSD(aov.observed.b)

manual = (values = c("Bouvetøya" = "red", "Possession" = "yellow", "Bartolomé" ="blue","Kerguelen" = "green", "Lauff" = "darkorange1", "Maher" = "skyblue1", "Siple" = "palevioletred", "Marion" = "darkgreen", "Peter I" = "mediumaquamarine", "Scott" = "gold", "South Georgia" = "darkmagenta"))
c("Extreme" = "red", "Temperate" = "blue")