library(phyloseq)
library(ggrepel)

#import log transformed otu table

asv_table_bact = as.data.frame(otu_table(physeq_transformed_log))
#change metadata to include the lat and long 

#import the env table 
env_table=as.data.frame(sample_data(physeq_transformed_log))
#standardize the env_table
env.bio <- env_table[,c(1:21)]
env.z.L5 <- decostand(env.bio, method="standardize")

set.seed(5)
upr = capscale(t(asv_table_bact) ~ ., data=env.z.L5, dist="bray", add=TRUE)
lwr = capscale(t(asv_table_bact) ~ 1, data=env.z.L5, dist="bray", add=TRUE)


mods = ordistep(lwr, scope = formula(upr), trace = TRUE, perm.max = 999, direction = c("forward"), R2scope=TRUE)
anova(mods, permutations = 999) # provides statistical signficance of the model
vif.cca(mods)


dbRDA = capscale(t(asv_table_bact) ~ pH..KCl. + K + Mean.Annual.Temp + Temp.Seasonality + Mg + CEC + OM + NH4.N, data = env.z.L5, dist="bray", add=TRUE)
vif.cca(dbRDA)

dbRDA = capscale(t(asv_table_bact) ~ pH..KCl.+ Mean.Annual.Temp + K + Mg + CEC + OM + NH4.N, data = env.z.L5, dist="bray", add=TRUE)
vif.cca(dbRDA)

dbRDA = capscale(t(asv_table_bact) ~ pH..KCl.+ Mean.Annual.Temp + Mg + CEC + OM + NH4.N, data = env.z.L5, dist="bray", add=TRUE)
vif.cca(dbRDA)
dbRDA_anova = anova(dbRDA, by="terms", permu=999)
summary(dbRDA)
R2adj <- RsquareAdj(dbRDA)$adj.r.squared

summary(dbRDA)

site <- factor(env_table$Island)
col.gr <- c("#CCFF00", "blue", "#D2B48C", "#7FD2FF", "#FF5500", "#9248A7" ,
            "#1C8859", "#333333", "#8B4513", "#49F149", "#FF0000")

gr.use <- factor(site)



plot(dbRDA , scaling=1, cex.lab =1.3, cex.axis=1.3, 
     font.lab = 1, display = "sites", 
     xlim=c(-0.75,0.75), ylim=c(-0.75,0.75),
     type="n", xlab=c("RDA1 (13.9%)"), ylab=c("RDA2 (10.2%)"))
points(scores(dbRDA , display="sites", choices=c(1,2), scaling=1),
       pch = 21, col = "black", bg = col.gr[gr.use], cex=2.5, lwd = 1.3)
legend(0.5, 0.8, legend = levels(gr.use), col= "black", pt.bg = col.gr, pch=c(21), cex=0.75, pt.cex = 1.5)

text(-0.5,-0.5, "p-value < 0.001, R2 = 0.26", font = 2)

text(dbRDA, display = "bp", 
     col="black", cex=1,)

env.z.L5_sig <- env.z.L5[,c("pH..KCl.", "Mean.Annual.Temp", "Mg", "CEC", "OM", "NH4.N")]

pcnm_soilchem = merge(env.z.L5_sig, table_pcnm_sig, by = 'row.names' )
pcnm_soilchem_2 = pcnm_soilchem[,-1]
rownames(pcnm_soilchem_2) <- pcnm_soilchem[,1]

#partition variation analysis using the same data, after calculating the explanatory variables with the ordistep modelling.Variables should be clustered into 2 to 4 groups, for instance, climatic (temp) vs chemical (Ca) 
mod=varpart(t(asv_table_bact), ~  Mg + CEC + OM + NH4.N + pH..KCl., ~ Mean.Annual.Temp, ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7, data= pcnm_soilchem_2 ,scale = FALSE)
#plotting as venn diagram
plot(mod, cutoff = -Inf, cex = 0.7, bg=2:5)


#fungal

otu_table_gene_rda_f=otu_table(physeq_transformed_f)


#change metadata to include the lat and long 
metadata_new= read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")
sample_data(physeq_transformed_f) = sample_data(metadata_new)

#import the env table 
env_table_f=as.data.frame(sample_data(physeq_transformed_f))
#standardize the env_table
env.bio_f <- env_table_f[,c(1:21)]
env.z.L5_f <- decostand(env.bio_f, method="standardize")

DCA = decorana(t(otu_table_gene_rda_f))



upr = capscale(t(otu_table_gene_rda_f) ~ ., data=env.z.L5_f, dist="bray", add=TRUE)
lwr = capscale(t(otu_table_gene_rda_f) ~ 1, data=env.z.L5_f, dist="bray", add=TRUE)
set.seed(5)

mods = ordistep(lwr, scope = formula(upr), trace = TRUE, perm.max = 999, direction = c("forward"), R2scope=TRUE)
anova(mods, permutations = 999) # provides statistical signficance of the model
vif.cca(mods)

dbRDA = capscale(t(otu_table_gene_rda_f) ~ K + Temp.Seasonality + Mean.Annual.Temp + Density + Ca  + P, data = env.z.L5_f, dist="bray", add=TRUE) # add a constant (add = TRUE) to prevent negative eigenvalues
vif.cca(dbRDA)
dbRDA = capscale(t(otu_table_gene_rda_f) ~ K +  Mean.Annual.Temp + Density + Ca  + P, data = env.z.L5_f, dist="bray", add=TRUE) # add a constant (add = TRUE) to prevent negative eigenvalues
vif.cca(dbRDA)

summary(dbRDA)

dbRDA_anova = anova(dbRDA, by="terms", permu=999)
R2adj <- RsquareAdj(dbRDA)$adj.r.squared

#choose significant factors 
data.env.subs.L5_f <- env.z.L5_f[,c("K", "Mean.Annual.Temp",  "Digtheid", "Ca", "PBray1")]

pcnm_soilchem = merge(data.env.subs.L5_f, table_pcnm_sig, by = 'row.names' )
pcnm_soilchem_2 = pcnm_soilchem[,-1]
rownames(pcnm_soilchem_2) <- pcnm_soilchem[,1]



#partition variation analysis using the same data, after calculating the explanatory variables with the ordistep modelling.Variables should be clustered into 2 to 4 groups, for instance, climatic (temp) vs chemical (Ca) 
mod=varpart(t(otu_table_gene_rda_f), ~ Ca + K + Digtheid + PBray1, ~ Mean.Annual.Temp, ~ PCNM1+ PCNM2+ PCNM3+ PCNM4 + PCNM5 + PCNM6 + PCNM7, data= pcnm_soilchem_2 ,scale = FALSE)
mod
#plotting as venn diagram
plot(mod, cutoff = -Inf, cex = 0.7, bg=2:5)

site <- factor(env_table$Island)
col.gr <- c("#CCFF00", "blue", "#D2B48C", "#7FD2FF", "#FF5500", "#9248A7" ,
            "#1C8859", "#333333", "#8B4513", "#49F149", "#FF0000")

gr.use <- factor(site)



plot(dbRDA , scaling=1, cex.lab =1.3, cex.axis=1.3, 
     font.lab = 1, display = "sites", 
     xlim=c(-0.75,0.75), ylim=c(-0.75,0.75),
     type="n", xlab=c("RDA1 (8.5%)"), ylab=c("RDA2 (7.5%)"))
points(scores(dbRDA , display="sites", choices=c(1,2), scaling=1),
       pch = 21, col = "black", bg = col.gr[gr.use], cex=2.5, lwd = 1.3)
legend(1.75, 2, legend = levels(gr.use), col= "black", pt.bg = col.gr, pch=c(21), cex=0.75, pt.cex = 1.5)

text(-0.5,-0.5, "p-value < 0.001, R2 = 0.15", font = 2)

text(dbRDA, display = "bp", 
     col="black", cex=1,)

#GLMs to confirm the RDA results - using the rarefied counts table
library(mvabund)
install.packages('Rcpp')
library(Rcpp)
library(dplyr)

#for bacteria

physeq_bact_genus = tax_glom(physeq_b_rarefy, taxrank = "Genus", NArm = F)
otu_table_gene_rda=otu_table(physeq_bact_genus)

abund_table = t(otu_table_gene_rda)

abund_table_m = mvabund(abund_table)


env.z.L5_sig <- env.z.L5[,c("pH..KCl.", "Mean.Annual.Temp", "Mg", "CEC", "OM", "NH4.N")]

physeq_f_genus = tax_glom(physeq_f_rarefy, taxrank = "Genus", NArm = F)
env_table_f = sample_data(physeq_f_genus)

otu_table_gene_rda_f=otu_table(physeq_f_genus)
abund_table_f = t(otu_table_gene_rda_f)
abund_table_m_f = mvabund(abund_table_f)

ft_pH_b = manyglm(abund_table_m ~ env_table$pH..KCl. + env_table$Mean.Annual.Temp + env_table$Mg + env_table$CEC + env_table$OM + env_table$NH4.N  , family = "negative.binomial")
plot(ft_pH_b)

bact_pH_GLM = anova.manyglm(ft_pH_b)

ft_Temp = manyglm(abund_table_m ~ env_table$Mean.Annual.Temp, family = "negative.binomial")
bact_Temp_GLM = anova.manyglm(ft_Temp)

ft_Mg = manyglm(abund_table_m ~ env_table$Mg, family = "negative.binomial")
bact_Mg_GLM = anova.manyglm(ft_Mg)

ft_CEC = manyglm(abund_table_m ~ env_table$CEC, family = "negative.binomial")
bact_CEC_GLM = anova.manyglm(ft_CEC)

ft_OM = manyglm(abund_table_m ~ env_table$OM, family = "negative.binomial")
bact_OM_GLM = anova.manyglm(ft_OM)

ft_NH4_N = manyglm(abund_table_m ~ env_table$NH4.N, family = "negative.binomial")
bact_NH4_GLM = anova.manyglm(ft_NH4_N)


#for fungi


ft_Lat_f = manyglm(abund_table_m_f ~ env.bio_f$Lat, family = "negative.binomial")
anova.manyglm(ft_Lat_f)

ft_Long_f = manyglm(abund_table_m_f ~ env.bio_f$Long, family = "negative.binomial")
anova.manyglm(ft_Long_f)
summary(ft_Long_f)

ft_Ca_f = manyglm(abund_table_m_f ~ env.bio_f$Ca, family = "negative.binomial")
anova.manyglm(ft_Ca_f)

ft_Mg_f = manyglm(abund_table_m_f ~ env.bio_f$Mg, family = "negative.binomial")
anova.manyglm(ft_Mg_f)

ft_D_f = manyglm(abund_table_m_f ~ env.bio_f$Digtheid, family = "negative.binomial")
anova.manyglm(ft_D_f)

ft_OM_f = manyglm(abund_table_m_f ~ env.bio_f$OM, family = "negative.binomial")
anova.manyglm(ft_OM_f)

ft_MAT_f = manyglm(abund_table_m_f ~ env.bio_f$Mean.Annual.Temp, family = "negative.binomial")
anova.manyglm(ft_MAT_f)

ft_AP_f = manyglm(abund_table_m_f ~ env.bio_f$Annual.Precipitation, family = "negative.binomial")
anova.manyglm(ft_AP_f)

summary(ft_Long_f, test = "LR")
