#let's play with the 0x, 4x stuff 
#first let's validate jamie's work with the snpEff annotations
#here's Jamie's file
load("/Users/Andrew/Downloads/Dplex_0x4x.Rdata")
#it's called my.coords.table
jstart<-as.data.frame(my.coords.table)
colnames(jstart)<-c("nth.gene","Scaffold","Strand","POS","coord.cds","transcript","codon","codon.base",
"AminoAcid","Codon.POS","Degeneracy")
#here's mine
start<-read.table(file="/Users/Andrew/Downloads/SNPs/filteredDnDs/newCounts/wholegenomerelevantpolysnps.txt",header=T, stringsAsFactors=F)
#examen<-merge(jstart,start,by=c("Scaffold","POS"))
##indeed, the 1x sites are all changes (missense, stop gained, etc)
#the 2 and 3xs are a mix of synon and non synon
#the 4xs are all synonymous 
#also of note, there are 20922015 CDS sites in the Dplex v3 genome
#and with 15130 genes, that's an average of 1382.817 CDS bases per gene, 
#and a length of 460.9389 amino acids per gene
#
#
#now we can proceed to new downstream analyses
#goals:
#split out autosomes vs Z
#split out sperm proteome and non-sperm
#to do so we need to change the naming scheme (i.e. remove the "-TA")
jstart$transcript<-substr(as.character(jstart$transcript),1,nchar(as.character(jstart$transcript))-3)
#now we can merge to get sperm proteome status and genomic location info

#get 0x and 4x
fulldegens<-jstart[which(jstart$Degeneracy == "4"),]
nodegens<-jstart[which(jstart$Degeneracy == "0"),]
#in the whole geneome there are 13423553 0x sites and 3040820 4x sites
#meaning there are 4.4x as many 0xs as 4xs



