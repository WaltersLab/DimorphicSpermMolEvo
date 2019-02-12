#Read in FPKM data 
setwd("~/Documents/ExpressionWork/Danaus_Tissue_Specificity/Input_Data/")
monarch.fpkm.raw <- read.table("danaus.fpkm.tmm.matrix.txt", header=TRUE, as.is = T)
str(monarch.fpkm.raw)

setwd("~/Documents/ExpressionWork/Danaus_Tissue_Specificity/")
#Identify genes with no expression (i.e rows with all zeros)
monarch.fpkm.zero <- monarch.fpkm.raw[apply(monarch.fpkm.raw[,-c(1:3)], 1, function(x) all(x==0)),]

#List of genes which have been excluded due to no expression 
write(unlist(lapply(monarch.fpkm.zero[,1], paste, collapse="")), file="Unexpressed.Genes", ncolumns = 1, sep="/t")

#Remove rows of all zeros (i.e. no expression in any tissue)
monarch.fpkm <- monarch.fpkm.raw[apply(monarch.fpkm.raw[,-c(1:3)], 1, function(x) !all(x==0)),]

#Add "1" to all values 
monarch.fpkm[4:30] <- monarch.fpkm[4:30] + 1

#Log transform all values 
monarch.fpkm[4:30] <- log(monarch.fpkm[4:30])

#Separate data into tissue type and sex 
fpkm.head.female <- monarch.fpkm[4:6]
fpkm.head.male <- monarch.fpkm[7:9]
fpkm.mg.female <- monarch.fpkm[10:12]                          
fpkm.mg.male <- monarch.fpkm[13:15]
fpkm.thorax.male <- monarch.fpkm[16:18]
fpkm.thorax.female <- monarch.fpkm[19:21]
fpkm.ovary <- monarch.fpkm[22:24]
fpkm.testes <- monarch.fpkm[25:27]
fpkm.ag <- monarch.fpkm[28:30]



#Average expression across biological replicates and for separate sexes
head.male <- rowMeans(fpkm.head.male)
mg.male <- rowMeans(fpkm.mg.male)
thorax.male <- rowMeans(fpkm.thorax.male)
fpkm.testes.avg <- rowMeans(fpkm.testes)
fpkm.ag.avg <- rowMeans(fpkm.ag)


#Form a dataframe with averages for tissue and sex 
fpkm.avg.sex <- cbind(monarch.fpkm[,1:3], head.male, mg.male,  thorax.male,  fpkm.testes.avg, fpkm.ag.avg)


#SPM shows, for each tissue separately, how specific a gene is to that tissue
##Calculating SPM
###Square every value 
fpkm.avg.square <- (fpkm.avg.sex[,4:8])^2
###Divide by squared sum of an entire row (gene)
SPM <- fpkm.avg.square/rowSums(fpkm.avg.square)
###Switch out of scientific notation
options(scipen=8)
###Recombine into full dataframe with gene information
SpmGenes <- cbind(fpkm.avg.sex[,1:3], SPM)
head(SpmGenes)  ##SpmGenes is the dataset for SPM values 
###Write a .csv table of SpmGenes
write.csv(SpmGenes, file="SpmGenesAJM.csv")



#SPM
##SPM multi-panel histograms
###Separate into SPM by tissue 
Spm.head.male <- as.data.frame((SpmGenes[,4]))
Spm.head.female <- as.data.frame(SpmGenes[,5])
Spm.mg.male <- as.data.frame(SpmGenes[,6])
Spm.mg.female <- as.data.frame(SpmGenes[,7])
Spm.thorax.male <- as.data.frame(SpmGenes[,8])
Spm.thorax.female <- as.data.frame(SpmGenes[,9])
Spm.ovary <- as.data.frame(SpmGenes[,10])
Spm.testes <- as.data.frame(SpmGenes[,11])
Spm.ag <- as.data.frame(SpmGenes[,12])

###Create histograms for each dataset 
####Male head 
jpeg(file="Male_Head_Spm_Hist.jpeg")
Male.Head.Spm.Hist <- ggplot(Spm.head.male, aes(x=Spm.head.male, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="darkolivegreen", fill="darkolivegreen3") + 
  geom_density(alpha=0.2, fill="darkolivegreen", col="darkolivegreen") +
  labs(title="Head (Male)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Female head
jpeg(file="Female_Head_Spm_Hist.jpeg")
Female.Head.Spm.Hist <- ggplot(Spm.head.female, aes(x=Spm.head.female, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="darkolivegreen", fill="darkolivegreen1") + 
  geom_density(alpha=0.2, fill="darkolivegreen", col="darkolivegreen") +
  labs(title="Head (Female) ", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Male midgut 
jpeg(file="Male_Mg_Spm_Hist.jpeg")
Male.Mg.Spm.Hist <- ggplot(Spm.mg.male, aes(x=Spm.mg.male, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="burlywood4", fill="burlywood3") + 
  geom_density(alpha=0.2, fill="burlywood4", col="burlywood4") +
  labs(title="Midgut (Male)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Female Midgut
jpeg(file="Female_Mg_Spm_Hist.jpeg")
Female.Mg.Spm.Hist <- ggplot(Spm.mg.male, aes(x=Spm.mg.male, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="burlywood3", fill="burlywood1") + 
  geom_density(alpha=0.2, fill="burlywood3", col="burlywood3") +
  labs(title="Midgut (Female)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Male Thorax
jpeg(file="Male_Thorax_Spm_Hist.jpeg")
Male.Thorax.Spm.Hist <- ggplot(Spm.thorax.male, aes(x=Spm.thorax.male, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="gold4", fill="gold2") + 
  geom_density(alpha=0.2, fill="gold4", col="gold4") +
  labs(title="Thorax (Male) ", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Female Thorax
jpeg(file="Female_Thorax_Spm_Hist.jpeg")
Female.Thorax.Spm.Hist <- ggplot(Spm.thorax.female, aes(x=Spm.thorax.female, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="gold3", fill="lightgoldenrod") + 
  geom_density(alpha=0.2, fill="gold3", col="gold3") +
  labs(title="Thorax (Female)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Testes 
jpeg(file="Testes_Spm_Hist.jpeg")
Testes.Spm.Hist <- ggplot(Spm.testes, aes(x=Spm.testes, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="skyblue4", fill="skyblue3") + 
  geom_density(alpha=0.2, fill="skyblue4", col="skyblue4") +
  labs(title="Testes", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Ovaries
jpeg(file="Ovaries_Spm_Hist.jpeg")
Ovary.Spm.Hist <- ggplot(Spm.ovary, aes(x=Spm.ovary, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="lightsteelblue4", fill="lightsteelblue1") + 
  geom_density(alpha=0.2, fill="lightsteelblue4", col="lightsteelblue4") +
  labs(title="Ovaries", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Accessory Glands
jpeg(file="Ag_Head_Spm_Hist.jpeg")
Ag.Spm.Hist <- ggplot(Spm.ag, aes(x=Spm.ag, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="indianred4", fill="indianred1") + 
  geom_density(alpha=0.2, fill="indianred4", col="indianred4") +
  labs(title="Accessory Glands", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
dev.off()

####Combine density histograms
#####Install ggpubr to easily format multipanel plots 
install.packages("ggpubr")
library("cowplot")
pdf(file="Combined_Spm_Hist.pdf")
Combined.Spm.Hist <- plot_grid(Male.Head.Spm.Hist, Female.Head.Spm.Hist,
         Male.Mg.Spm.Hist, Female.Mg.Spm.Hist,
         Male.Thorax.Spm.Hist, Female.Thorax.Spm.Hist,
         Testes.Spm.Hist, Ovary.Spm.Hist,
         Ag.Spm.Hist, 
         ncol=2, nrow=5,
         labels="auto")
dev.off()

##SPM violin plots 
###"Melt" data into a "variable" column and a "value" column 
library(reshape2)
SpmMelt <- melt(SpmGenes[,4:12])
str(SpmMelt)

###Initiate plot 
p <- ggplot(SpmMelt, aes(x=variable, y=value, fill=variable))

###Customize plot 
pdf(file="VioPlotSPM.pdf")
VioPlotLabs <- c("Head (Male)", "Head (Female)", "Midgut (Male)",
                 "Midgut (Female)", "Thorax (Male)", "Thorax (Female)", 
                 "Ovary", "Testes", "Accessory Glands")
VioPlotSPM <- p + geom_violin() + 
  labs(x="Tissue", y="SPM Tissue Specificity") + 
  theme(axis.text.x = element_text(angle=60, hjust=1)) + 
  scale_x_discrete(labels=VioPlotLabs) + 
  theme(legend.position = "none") 
VioPlotSPM + scale_fill_brewer(palette="greens")
dev.off()

##SPM density curve 
###Rename variables to correspond to how you want them labeled in the plot
SpmMelt$variable <- factor(SpmMelt$variable, levels=c("head.male", "head.female",
                                                      "mg.male", "mg.female", 
                                                      "thorax.male", "thorax.female",
                                                      "fpkm.ovary.avg", "fpkm.testes.avg",
                                                      "fpkm.ag.avg"), 
                                              labels=c("Head (Male)", "Head (Female)", "Midgut (Male)",
                                                       "Midgut (Female)", "Thorax (Male)", "Thorax (Female)", 
                                                       "Ovary", "Testes", "Accessory Glands"))
###Create density plot of all tissues 
pdf("DensityPlotSPM.pdf")
SPMDensity <- ggplot(SpmMelt, aes(x=value, color=variable, ..density..)) + 
  geom_density() + 
  labs(x="SPM Tissue Specificity Values", y="Density", color="Tissue") + 
  theme(legend.key.size = unit(0.4, "cm"), legend.position = c(0.8, 0.75)) + 
  scale_color_brewer(palette="Set1")
dev.off()

#ks.test comparing tissue by tissue - applied Bonferroni correction by dividing p value by 81 
f1 <- function(x,y,...,
              alternative="two.sided", exact=TRUE) {
  p <- ks.test(x,y,..., alternative=alternative, exact=exact)$p.value
  p/81
}
ks.results.monarch <- apply(SPM, 2, function(x) apply(SPM, 2, function(y) f(x,y)))

#repeat for Wilcox test (also w/Bonferroni correction) 
f2 <- function(x,y,...,
              alternative="two.sided", exact=TRUE) {
  p <- wilcox.test(x,y,..., alternative=alternative, exact=exact)$p.value
  p/81
}
wilcox.results.monarch <- apply(SPM, 2, function(x) apply(SPM, 2, function(y) f(x,y)))

#######

