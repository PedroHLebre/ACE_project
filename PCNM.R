library(vegan)

#export table from transformed phyloseq
asv_table_bact = as.data.frame(otu_table(physeq_transformed_log))

#import coordinates and create PCNM matrix
coordinates_xy = read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
coor_pcnm = pcnm(dist(coordinates_xy))
table_pcnm = as.data.frame(coor_pcnm$vectors)

set.seed(5)
upr = capscale(t(asv_table_bact) ~ ., data=table_pcnm, dist="bray", add=TRUE)
lwr = capscale(t(asv_table_bact) ~ 1, data=table_pcnm, dist="bray", add=TRUE)

mods = ordistep(lwr, scope = formula(upr), trace = TRUE, perm.max = 999, direction = c("forward"), R2scope=TRUE)
anova(mods, permutations = 999)

vif.cca(mods)

table_pcnm_sig <- table_pcnm[,c("PCNM1", "PCNM2","PCNM3","PCNM4", "PCNM5", "PCNM6", "PCNM7")]

dbRDA = capscale(t(otu_table_gene_rda) ~ ., data = table_pcnm_sig, dist="bray", add=TRUE) 
plot(dbRDA)
summary(dbRDA)
vif.cca(dbRDA)
dbRDA_anova = anova(dbRDA, by="terms", permu=999)
summary(dbRDA_anova)
R2adj <- RsquareAdj(dbRDA)$adj.r.squared
summary(R2adj, by= "terms")

site <- factor(sample_data(physeq_transformed_log)$Island)
col.gr <- c("blue", "red2", "green", "skyblue4", "cyan", "darkgreen" ,
            "deeppink", "orange", "gold", "palevioletred", "darkmagenta")

gr.use <- factor(site)



plot(dbRDA , scaling=1, cex.lab =1.3, cex.axis=1.3, 
     font.lab = 1, display = "sites", 
     xlim=c(-0.75,0.75), ylim=c(-0.75,0.75),
     type="n", xlab=c("CAP1 (13.6%)"), ylab=c("CAP2 (8.3%)"))
points(scores(dbRDA , display="sites", choices=c(1,2), scaling=1),
       pch = 21, col = "black", bg = col.gr[gr.use], cex=2.5, lwd = 1.3)
legend(1.75, 2, legend = levels(gr.use), col= "black", pt.bg = col.gr, pch=c(21), cex=0.75, pt.cex = 1.5)

text(-1.0,-0.5, "p-value < 0.001, R2 = 0.26", font = 2)

text(dbRDA, display = "bp", 
     col="black", cex=1,)

write.csv(table_pcnm_sig, "PCNM_sig_bacteria.csv")

#For fungi

#export table from transformed phyloseq
asv_table_fung = as.data.frame(otu_table(physeq_transformed_f))

upr = capscale(t(asv_table_fung) ~ ., data=table_pcnm, dist="bray", add=TRUE)
lwr = capscale(t(asv_table_fung) ~ 1, data=table_pcnm, dist="bray", add=TRUE)

mods_f = ordistep(lwr, scope = formula(upr), trace = TRUE, perm.max = 999, direction = c("forward"), R2scope=TRUE)
anova(mods_f, permutations = 999)

vif.cca(mods_f)


table_pcnm_sig <- table_pcnm[,c("PCNM1", "PCNM2","PCNM3","PCNM4", "PCNM5", "PCNM6", "PCNM7")]

dbRDA = capscale(t(asv_table_fung) ~ ., data = table_pcnm_sig, dist="bray", add=TRUE) 
plot(dbRDA)
summary(dbRDA)
vif.cca(dbRDA)
dbRDA_anova = anova(dbRDA, by="terms", permu=999)
summary(dbRDA_anova)
R2adj <- RsquareAdj(dbRDA)$adj.r.squared


site <- factor(sample_data(physeq_transformed_f)$Island)
col.gr <- c("blue", "red2", "green", "skyblue4", "cyan", "darkgreen" ,
            "deeppink", "orange", "gold", "palevioletred", "darkmagenta")

gr.use <- factor(site)



plot(dbRDA , scaling=1, cex.lab =1.3, cex.axis=1.3, 
     font.lab = 1, display = "sites", 
     xlim=c(-0.75,0.75), ylim=c(-0.75,0.75),
     type="n", xlab=c("CAP1 (13.6%)"), ylab=c("CAP2 (8.3%)"))
points(scores(dbRDA , display="sites", choices=c(1,2), scaling=1),
       pch = 21, col = "black", bg = col.gr[gr.use], cex=2.5, lwd = 1.3)
legend(1.75, 2, legend = levels(gr.use), col= "black", pt.bg = col.gr, pch=c(21), cex=0.75, pt.cex = 1.5)

text(-1.0,-0.5, "p-value < 0.001, R2 = 0.26", font = 2)

text(dbRDA, display = "bp", 
     col="black", cex=1,)
