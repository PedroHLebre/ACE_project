library(Hmisc)
library(corrplot)
library(vegan)
library(ggbiplot)
library(FactoMineR)
library(factoextra)
library(grid)
library(phyloseq)

#upload the soil chemistry and climatic data 

table=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

#select only the numerical columns

env.bio <- table[,c(1:19)]
env <- table[,20:21]
hist(chemistry_norm)

#Standardize the data and calculate the pearson correlation between variables.

chemistry_norm<-decostand(env.bio, "stand") 
res_pearson <- rcorr(as.matrix(chemistry_norm),type = c("pearson"))


#Plot correlation bubble plot

flattenCorrMatrix(res2$r, res2$P)
corrplot(res_pearson$r, type="upper", order="hclust", 
         p.mat = res_pearson$P, sig.level = 0.01, insig = "n")

#Calculate PCA distance between samples and plot as an ordination plot

prC <- prcomp(chemistry_norm, center=TRUE, scale=TRUE) 
prc = PCA(chemistry_norm)
eig.val <- get_eigenvalue(prc)
summary(prC)

manual = c("Bouvetoya" = "#004586", "Possession" = "#FA9FB5", "Bartolome" ="#FFD320","Kerguelen" = "#17BECF", "Lauff" = "#4C7022", "Maher" = "#EFEDF5", "Siple" = "#ADDD8E", "Marion" = "black", "Peter I" = "#4B1F6F", "Scott" = "#FF950E", "South Georgia" = "#C5000B")

head(prc$coord)
plot(prcomp_soilchem, type = "l") 
my_text <- "R2 = 0.6, p-value > 0.01"
my_grob = grid.text(my_text, x=0.3,  y=0.8, gp=gpar(col="firebrick", fontsize=12, fontface="bold"))

print(prcomp_soilchem)
prC <- prcomp(chemistry_norm, center=TRUE, scale=TRUE)

g <- ggbiplot(prc, groups = table$Island,ellipse=T, var.scale = 1, var.axes = F,circle = F) + theme_bw()
g2 <- g 
g3 <- g2 + scale_color_manual(values = manual)
g5 = g3 + theme_bw()
p6 = g5 + theme(axis.title.x = element_text(size=12, face="bold", margin = margin(t = 20,r = 20, b = 20,l = 0)),axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 20)))
p7 = p6 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + annotation_custom(my_grob)
p7

#Test the significance of the sample distances according to groups

dtp <- data.frame('Site' = as.character(table$Island),
                  prC$x[,1:2])

PCAcoords <- dtp[,c("PC1","PC2")]


PCAdist <- dist(PCAcoords)


adonis(PCAdist  ~ Site, data = dtp, permutations = 1000)

beta <- betadisper(PCAdist, dtp$Site)
permutest(beta)
