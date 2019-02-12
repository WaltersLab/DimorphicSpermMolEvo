setwd("~/Documents/ExpressionWork/Danaus_Tissue_Specificity")

#Read in manduca FPKM data
manduca.fpkm.raw <- read.csv("./Input_Data/manduca_fpkm.csv", header=TRUE)

#Take out gene.description 
manduca.fpkm <- manduca.fpkm.raw[,-2]

#Identify genes with no expression
manduca.fpkm.zero <- manduca.fpkm[apply(manduca.fpkm[,-1], 1, function(x) all(x==0)),]

#List of genes which have been excluded due to no expression 
write(unlist(lapply(manduca.fpkm.zero[,1], paste, collapse="")), file="Unexpressed.Manduca.Genes", ncolumns = 1, sep="/t")

#Remove rows of all zeros (i.e. no expression in any tissue)
manduca.fpkm <- manduca.fpkm[apply(manduca.fpkm[,-1], 1, function(x) !all(x==0)),]

#Add "1" to all values 
manduca.fpkm[,c(2:68)] <- manduca.fpkm[,c(2:68)] + 1

#Log transform all values 
manduca.fpkm[,c(2:68)] <- log(manduca.fpkm[,c(2:68)])

#Group into tissue and developmental stage and average across biological replicates for tissues w/more than one  
headLarvae <- rowMeans(manduca.fpkm[,c(2:8)])
headPupae <- manduca.fpkm[,9]
headAdult <- rowMeans(manduca.fpkm[,c(10:12)])
fatLarvae <- rowMeans(manduca.fpkm[,c(13:16)])
fatPupae <- rowMeans(manduca.fpkm[,c(17:18)])
fatAdult <- rowMeans(manduca.fpkm[,c(19:20)])
gutLarvae <- rowMeans(manduca.fpkm[,c(26:35)])
gutPupae <- rowMeans(manduca.fpkm[,c(36:37)])
gutAdult <- manduca.fpkm[,38]
mtLarvae <- manduca.fpkm[,39]
mtAdult <- rowMeans(manduca.fpkm[,c(40:41)])
muscleLarvae <- rowMeans(manduca.fpkm[,c(42:48)])
testesPupae <- rowMeans(manduca.fpkm[,c(49:50)])
testesAdult <- manduca.fpkm[,51]
ovary.Pupae <- manduca.fpkm[,52]
ovary.Adult <- manduca.fpkm[,53]
headAdultFemale <- rowMeans(manduca.fpkm[,c(54:57)])
headAdultMale <- rowMeans(manduca.fpkm[,c(58:61)])
antennaeLarvae <- rowMeans(manduca.fpkm[,c(62:64)])
antennaeAdultFemale <- rowMeans(manduca.fpkm[,c(65:67)])
antennaeAdultMale <- manduca.fpkm[,68]

#Create dataframe with all tissues averaged across replciates  
geneID <- manduca.fpkm[,1]
manduca.fpkm.tissue <- cbind.data.frame(as.character(geneID), headLarvae, headPupae, headAdult, fatLarvae, fatPupae, fatAdult, 
                                        gutLarvae, gutPupae, gutAdult, mtLarvae, mtAdult, muscleLarvae, testesPupae,
                                        testesAdult, ovary.Pupae, ovary.Adult, headAdultFemale, headAdultMale, antennaeLarvae,
                                        antennaeAdultFemale, antennaeAdultMale)



##We want to make an apples-to-apples comparison with our expression data. SPM is going to be very sensitive to the number of tissues sampled, so because we only have a a few tissues in monarchs (and all adult male), we want to down-sample the Manduca RNA data. 

manduca.fpkm.tissue<-cbind.data.frame(as.character(geneID),manduca.fpkm.tissue$headAdultMale,manduca.fpkm.tissue$gutAdult, manduca.fpkm.tissue$mtAdult, manduca.fpkm.tissue$testesAdult)
colnames(manduca.fpkm.tissue)<-c("Gene","headMale","gutAdult","mtAdult","testesAdult")

#Calculating SPM
##Square every value 
manduca.fpkm.tissue.square <- (manduca.fpkm.tissue[,2:5])^2
##Divide by squared sum of an entire row (gene)
SPMManduca <- manduca.fpkm.tissue.square/rowSums(manduca.fpkm.tissue.square)
##Recombine into full dataframe with gene information
SpmGenesManduca <- cbind.data.frame(as.character(geneID), SPMManduca)
###Write a .csv table of SpmGenes
write.csv(SpmGenesManduca, file="SpmGenesManduca.csv")

#Tau Histogram 
library(ggplot2)


#SPM Histograms 
Spm.head.larvae <- as.data.frame(SPMManduca[,1])
Spm.head.pupae <- as.data.frame(SPMManduca[,2])
Spm.head.adult <- as.data.frame(SPMManduca[,3])
Spm.fat.larvae <- as.data.frame(SPMManduca[,4])
Spm.fat.pupae <- as.data.frame(SPMManduca[,5])
Spm.fat.adult <- as.data.frame(SPMManduca[,6])
Spm.gut.larvae <- as.data.frame(SPMManduca[,7])
Spm.gut.pupae <- as.data.frame(SPMManduca[,8])
Spm.gut.adult <- as.data.frame(SPMManduca[,9])
Spm.mt.larvae <- as.data.frame(SPMManduca[,10])
Spm.mt.adult <- as.data.frame(SPMManduca[,11])
Spm.muscle.larvae <- as.data.frame(SPMManduca[,12])
Spm.testes.pupae <- as.data.frame(SPMManduca[,13])
Spm.testes.adult <- as.data.frame(SPMManduca[,14])
Spm.ovary.pupae <- as.data.frame(SPMManduca[,15])
Spm.ovary.adult <- as.data.frame(SPMManduca[,16])
Spm.head.adult.female <- as.data.frame(SPMManduca[,17])
Spm.head.adult.male <- as.data.frame(SPMManduca[,18])
Spm.antennae.larvae <- as.data.frame(SPMManduca[,19])
Spm.antennae.adult.female <- as.data.frame(SPMManduca[,20])
Spm.antennae.adult.male <- as.data.frame(SPMManduca[,21])

Spm.head.larvae.plot <- ggplot(Spm.head.larvae, aes(x=Spm.head.larvae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="darkolivegreen", fill="darkolivegreen3") + 
  geom_density(alpha=0.2, fill="darkolivegreen", col="darkolivegreen") +
  labs(title="Head (Larvae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.head.pupae.plot <- ggplot(Spm.head.pupae, aes(x=Spm.head.pupae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="olivedrab3", fill="olivedrab1") + 
  geom_density(alpha=0.2, fill="olivedrab3", col="olivedrab1") +
  labs(title="Head (Pupae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.head.adult.plot <- ggplot(Spm.head.adult, aes(x=Spm.head.adult, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="seagreen4", fill="seagreen3") + 
  geom_density(alpha=0.2, fill="seagreen4", col="seagreen3") +
  labs(title="Head (Adult)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.fat.larvae.plot <- ggplot(Spm.fat.larvae, aes(x=Spm.fat.larvae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="peachpuff4", fill="peachpuff3") + 
  geom_density(alpha=0.2, fill="peachpuff4", col="peachpuff3") +
  labs(title="Fat (Larvae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.fat.pupae.plot <- ggplot(Spm.fat.pupae, aes(x=Spm.fat.pupae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="peachpuff2", fill="peachpuff") + 
  geom_density(alpha=0.2, fill="peachpuff2", col="peachpuff") +
  labs(title="Fat (Pupae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.fat.adult.plot <- ggplot(Spm.fat.adult, aes(x=Spm.fat.adult, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="wheat3", fill="wheat1") + 
  geom_density(alpha=0.2, fill="wheat3", col="wheat1") +
  labs(title="Fat (Adult)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.gut.larvae.plot <- ggplot(Spm.gut.larvae, aes(x=Spm.gut.larvae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="palevioletred4", fill="palevioletred") + 
  geom_density(alpha=0.2, fill="palevioletred4", col="palevioletred") +
  labs(title="Gut (Larvae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.gut.pupae.plot <- ggplot(Spm.gut.pupae, aes(x=Spm.gut.pupae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="lightpink4", fill="lightpink2") + 
  geom_density(alpha=0.2, fill="lightpink4", col="lightpink2") +
  labs(title="Gut (Pupae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.gut.adult.plot <- ggplot(Spm.gut.adult, aes(x=Spm.gut.adult, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="pink3", fill="pink1") + 
  geom_density(alpha=0.2, fill="pink3", col="pink1") +
  labs(title="Gut (Adult)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.mt.larvae.plot <- ggplot(Spm.mt.larvae, aes(x=Spm.mt.larvae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="mediumpurple3", fill="mediumpurple1") + 
  geom_density(alpha=0.2, fill="mediumpurple3", col="mediumpurple1") +
  labs(title="Malphigian Tubules (Larvae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.mt.adult.plot <- ggplot(Spm.mt.adult, aes(x=Spm.mt.adult, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="darkorchid4", fill="darkorchid") + 
  geom_density(alpha=0.2, fill="darkorchid4", col="darkorchid") +
  labs(title="Malphigian Tubules (Adult)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.muscle.larvae.plot <- ggplot(Spm.muscle.larvae, aes(x=Spm.muscle.larvae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="gold4", fill="gold") + 
  geom_density(alpha=0.2, fill="gold4", col="gold") +
  labs(title="Muscle (Larvae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.testes.pupae.plot <- ggplot(Spm.testes.pupae, aes(x=Spm.testes.pupae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="dodgerblue4", fill="dodgerblue3") + 
  geom_density(alpha=0.2, fill="dodgerblue4", col="dodgerblue3") +
  labs(title="Testes (Pupae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.testes.adult.plot <- ggplot(Spm.testes.adult, aes(x=Spm.testes.adult, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="skyblue4", fill="skyblue1") + 
  geom_density(alpha=0.2, fill="skyblue4", col="skyblue1") +
  labs(title="Testes (Adult)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.ovary.pupae.plot <-  ggplot(Spm.ovary.pupae, aes(x=Spm.ovary.pupae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="plum4", fill="plum1") + 
  geom_density(alpha=0.2, fill="plum4", col="plum1") +
  labs(title="Ovary (Pupae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.ovary.adult.plot <- ggplot(Spm.ovary.adult, aes(x=Spm.ovary.adult, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="pink3", fill="pink") + 
  geom_density(alpha=0.2, fill="pink3", col="pink") +
  labs(title="Ovary (Adult)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.head.adult.female.plot <- ggplot(Spm.head.adult.female, aes(x=Spm.head.adult.female, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="slategray4", fill="slategray2") + 
  geom_density(alpha=0.2, fill="slategray4", col="slategray2") +
  labs(title="Head (Adult Female)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.head.adult.male.plot <- ggplot(Spm.head.adult.male, aes(x=Spm.head.adult.male, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="lightcyan4", fill="lightcyan") + 
  geom_density(alpha=0.2, fill="lightcyan4", col="lightcyan") +
  labs(title="Head (Adult Male)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.antennae.larvae.plot <- ggplot(Spm.antennae.larvae, aes(x=Spm.antennae.larvae, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="lightgoldenrod4", fill="lightgoldenrod") + 
  geom_density(alpha=0.2, fill="lightgoldenrod4", col="lightgoldenrod") +
  labs(title="Antennae (Larvae)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.antennae.adult.female.plot <- ggplot(Spm.antennae.adult.female, aes(x=Spm.antennae.adult.female, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="khaki3", fill="khaki1") + 
  geom_density(alpha=0.2, fill="khaki3", col="khaki1") +
  labs(title="Antennae (Adult Female)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
Spm.antennae.adult.male.plot <- ggplot(Spm.antennae.adult.male, aes(x=Spm.antennae.adult.male, y=..density..)) + 
  geom_histogram(binwidth=0.02, col="darkgoldenrod4", fill="darkgoldenrod3") + 
  geom_density(alpha=0.2, fill="darkgoldenrod4", col="darkgoldenrod3") +
  labs(title="Antennae (Adult Male)", x="SPM Tissue Specificity", y="Density") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_minimal() + 
  scale_x_continuous()
#Combine to multipanel histogram 
install.packages("ggpubr")
library("cowplot")
Combined.Spm.Manduca <- plot_grid(Spm.head.larvae.plot, Spm.head.pupae.plot, Spm.head.adult.plot,
                                  Spm.head.adult.female.plot, Spm.head.adult.male.plot, 
                                Spm.fat.larvae.plot, Spm.fat.pupae.plot, Spm.fat.adult.plot,
                               Spm.gut.larvae.plot, Spm.gut.pupae.plot, Spm.gut.adult.plot,
                              Spm.mt.larvae.plot, Spm.mt.adult.plot,
                                Spm.muscle.larvae.plot,
                               Spm.testes.pupae.plot, Spm.testes.adult.plot,
                               Spm.ovary.pupae.plot, Spm.ovary.adult.plot, 
                               Spm.antennae.larvae.plot, Spm.antennae.adult.female.plot, 
                               Spm.antennae.adult.male.plot,
                               ncol=3, nrow=8,
                               labels="auto")

#ks.test comparing tissue by tissue - applied Bonferroni correction by dividing p value by 81 
f1 <- function(x,y,...,
              alternative="two.sided", exact=TRUE) {
  p <- ks.test(x,y,..., alternative=alternative, exact=exact)$p.value
  p/81
}
ks.results.manduca <- apply(SPMManduca, 2, function(x) apply(SPMManduca, 2, function(y) f1(x,y)))


#Wilcox test
f2 <- function(x,y,...,
              alternative="two.sided", exact=TRUE) {
  p <- wilcox.test(x,y,..., alternative=alternative, exact=exact)$p.value
  p/81
}
wilcox.results.manduca <- apply(SPMManduca, 2, function(x) apply(SPMManduca, 2, function(y) f2(x,y)))

