##A file for figures for the Sperm Comp manuscript
##a to-do list from Jamie's comments:
#1. remove Manduca commixed but not localized proteins

#first some file loading and data slicing
library("data.table")
manmaster<-read.csv(file="/Users/Andrew/Documents/ManducaSextaFullGenomeWithZorA.csv", header=T)
man_append<-read.csv(file="/Users/Andrew/Documents/manducaorthosfordanaus.csv",header=T)
colnames(man_append)<-c("Gene","Chromosome","Danaus.Sperm.Ortho")
manmaster<-as.data.frame(manmaster, stringsAsFactors=F)
danmaster<-read.csv(file="/Users/Andrew/Downloads/evoconfiles/NewKanzen.csv",header=T,stringsAsFactors=F)
dan_append<-read.csv(file="/Users/Andrew/Documents/danausorthologsformanduca.csv",header=T)
colnames(dan_append)<-c("Gene","ZorA","Chromosome","Manduca.Sperm.Ortho")
danmaster<-as.data.frame(danmaster,stringsAsFactors=F)
mantranscend<-merge(manmaster,man_append,by=c("Gene","Chromosome"),all.x=T, all.y=T)
dantranscend<-merge(danmaster,dan_append,by=c("Gene","Chromosome","ZorA"),all.x=T, all.y=T)
#Okay, now we have master files with ortholog annotation
#let's grab the definitely non-sperm
mtbg<-mantranscend[which(mantranscend$Class=="non-sperm"),]
dtbg<-dantranscend[which(dantranscend$gene.class=="non-sperm"),]
#and just the sperm proteome
mtsperm<-as.data.frame(mantranscend[which(mantranscend$Class!="non-sperm"),])
dtsperm<-as.data.frame(dantranscend[which(dantranscend$Apyrene==1 | dantranscend$Eupyrene==1),])
#some have NAs still...because some genes have no identifiable ortho
#we will count these as not shared
mtsperm[is.na(mtsperm)]<-0
dtsperm[is.na(dtsperm)]<-0
#we put aside the Z for now
mtsauto<-mtsperm[which(mtsperm$Chromosome!=1),]
dtsauto<-dtsperm[which(dtsperm$ZorA=="A"),]
mtbgauto<-mtbg[which(mtbg$Chromosome!=1),]
dtbgauto<-dtbg[which(dtbg$ZorA=="A"),]

#let's look at dimorphic subsets
dtsae<-dtsauto[which(dtsauto$gene.class=="eupyrene"),]
dtsaa<-dtsauto[which(dtsauto$gene.class=="apyrene"),]
dtsas<-dtsauto[which(dtsauto$gene.class=="0"),]

mtsae<-mtsauto[which(mtsauto$Class=="eupyrene"),]
mtsaa<-mtsauto[which(mtsauto$Class=="apyrene"),]
mtsas<-mtsauto[which(mtsauto$Class=="shared"),]

par(mfrow=c(1,2))
boxplot(dtbgauto$pN/dtbgauto$Non.sites,outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3),main="Polymorphism vs Dimorphism")
boxplot(dtsae$pN/dtsae$Non.sites, add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(dtsas$pN/dtsas$Non.sites, add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(dtsaa$pN/dtsaa$Non.sites, add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)
text(1.5,0.025, "*",cex=2)
text(2,0.025, "*",cex=2)

boxplot(mtbgauto$pN/mtbgauto$Non.sites,outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3))
boxplot(mtsae$pN/mtsae$Non.sites, add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(mtsas$pN/mtsas$Non.sites, add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(mtsaa$pN/mtsaa$Non.sites, add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)

wilcox.test((dtbgauto$pN/dtbgauto$Non.sites),(dtsae$pN/dtsae$Non.sites))
#danaus Pn, eu vs bg, significant
#W = 1291700, p-value = 0.0002815
wilcox.test((dtbgauto$pN/dtbgauto$Non.sites),(dtsas$pN/dtsas$Non.sites))
#danaus Pn, shared vs bg, significant
#W = 1486100, p-value = 1.167e-08
wilcox.test((dtbgauto$pN/dtbgauto$Non.sites),(dtsaa$pN/dtsaa$Non.sites))
#danaus Pn, ap vs bg, not-significant!
#W = 284620, p-value = 0.1684

#for funsies
wilcox.test((dtsae$pN/dtsae$Non.sites),(dtsaa$pN/dtsaa$Non.sites))
#not significant

wilcox.test((mtbgauto$pN/mtbgauto$Non.sites),(mtsae$pN/mtsae$Non.sites))
#manduca Pn, eu vs bg, marginal, but NS
#W = 995190, p-value = 0.08571

wilcox.test((mtbgauto$pN/mtbgauto$Non.sites),(mtsas$pN/mtsas$Non.sites))
#manduca Pn, shared vs bg, v not significant
#W = 1470300, p-value = 0.3277
wilcox.test((mtbgauto$pN/mtbgauto$Non.sites),(mtsaa$pN/mtsaa$Non.sites))
#manduca Pn, shared vs bg, not-significant!
#W = 548570, p-value = 0.4974

#Dn time
wilcox.test((dtbgauto$dN/dtbgauto$Non.sites),(dtsae$dN/dtsae$Non.sites))
#danaus Dn, eu vs bg, it's significant
#W = 1021800, p-value = 0.01351
wilcox.test((dtbgauto$dN/dtbgauto$Non.sites),(dtsas$dN/dtsas$Non.sites))
#danaus Dn, shared vs bg, super not significant
#W = 1218000, p-value = 0.9185
wilcox.test((dtbgauto$dN/dtbgauto$Non.sites),(dtsaa$dN/dtsaa$Non.sites))
#danaus Dn, shared vs bg, not-significant!
#W = 266580, p-value = 0.6042

wilcox.test((mtbgauto$dN/mtbgauto$Non.sites),(mtsae$dN/mtsae$Non.sites))
#manduca Dn, eu vs bg, it's significant
#W = 966630, p-value = 0.3183
wilcox.test((mtbgauto$dN/mtbgauto$Non.sites),(mtsas$dN/mtsas$Non.sites))
#manduca Dn, shared vs bg, super not significant
#W = 1540700, p-value = 0.6883
wilcox.test((mtbgauto$dN/mtbgauto$Non.sites),(mtsaa$dN/mtsaa$Non.sites))
#manduca Dn, shared vs bg, not-significant!
#W = 560910, p-value = 0.2741

par(mfrow=c(1,2))
boxplot(dtbgauto$dN/dtbgauto$Non.sites,outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3))
boxplot(dtsae$dN/dtsae$Non.sites, add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(dtsas$dN/dtsas$Non.sites, add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(dtsaa$dN/dtsaa$Non.sites, add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)

boxplot(mtbgauto$dN/mtbgauto$Non.sites,outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3))
boxplot(mtsae$dN/mtsae$Non.sites, add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(mtsas$dN/mtsas$Non.sites, add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(mtsaa$dN/mtsaa$Non.sites, add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)

#Ds
par(mfrow=c(1,2))
boxplot(dtbgauto$dS/dtbgauto$Syn.sites,outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3))
boxplot(dtsae$dS/dtsae$Syn.sites, add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(dtsas$dS/dtsas$Syn.sites, add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(dtsaa$dS/dtsaa$Syn.sites, add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)

boxplot(mtbgauto$dS/mtbgauto$Syn.sites,outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3))
boxplot(mtsae$dS/mtsae$Syn.sites, add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(mtsas$dS/mtsas$Syn.sites, add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(mtsaa$dS/mtsaa$Syn.sites, add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)

#Ds time
wilcox.test((dtbgauto$dS/dtbgauto$Syn.sites),(dtsae$dS/dtsae$Syn.sites))
#danaus Dn, eu vs bg, it's marginally not significant!
#W = 1056000, p-value = 0.09275
wilcox.test((dtbgauto$dS/dtbgauto$Syn.sites),(dtsas$dS/dtsas$Syn.sites))
#danaus Dn, shared vs bg, super not significant
#W = 1209000, p-value = 0.7665
wilcox.test((dtbgauto$dS/dtbgauto$Syn.sites),(dtsaa$dS/dtsaa$Syn.sites))
#danaus Dn, shared vs bg, not-significant!
#W = 279420, p-value = 0.2594

wilcox.test((mtbgauto$dS/mtbgauto$Syn.sites),(mtsae$dS/mtsae$Syn.sites))
#manduca Ds, eu vs bg, it's significant
#W = 963410, p-value = 0.3596
wilcox.test((mtbgauto$dS/mtbgauto$Syn.sites),(mtsas$dS/mtsas$Syn.sites))
#manduca Ds, shared vs bg, marginally not significant
#W = 1437100, p-value = 0.103
wilcox.test((mtbgauto$dS/mtbgauto$Syn.sites),(mtsaa$dS/mtsaa$Syn.sites))
#manduca Ds, shared vs bg, not-significant!
#W = 495180, p-value = 0.2653


#what if we do Dn/Ds all together?
#Dn time
wilcox.test(((dtbgauto$dN/dtbgauto$Non.sites)/(dtbgauto$dS/dtbgauto$Syn.sites)),((dtsae$dN/dtsae$Non.sites)/(dtsae$dS/dtsae$Syn.sites)))
#danaus Dn/Ds, eu vs bg, it's juuuuust shy of significant
#W = 1047500, p-value = 0.05934
wilcox.test(((dtbgauto$dN/dtbgauto$Non.sites)/(dtbgauto$dS/dtbgauto$Syn.sites)),((dtsas$dN/dtsas$Non.sites)/(dtsas$dS/dtsas$Syn.sites)))
#danaus Dn, shared vs bg, super not significant
#W = 1240000, p-value = 0.7065
wilcox.test(((dtbgauto$dN/dtbgauto$Non.sites)/(dtbgauto$dS/dtbgauto$Syn.sites)),((dtsaa$dN/dtsaa$Non.sites)/(dtsaa$dS/dtsaa$Syn.sites)))
#danaus Dn/Ds, shared vs bg, not-significant!
#W = 273710, p-value = 0.389

wilcox.test(((mtbgauto$dN/mtbgauto$Non.sites)/(mtbgauto$dS/mtbgauto$Syn.sites)),((mtsae$dN/mtsae$Non.sites)/(mtsae$dS/mtsae$Syn.sites)))
#manduca, eu vs bg, is also marginal?
#W = 996490, p-value = 0.07971
wilcox.test(((mtbgauto$dN/mtbgauto$Non.sites)/(mtbgauto$dS/mtbgauto$Syn.sites)),((mtsas$dN/mtsas$Non.sites)/(mtsas$dS/mtsas$Syn.sites)))
#manduca, shared vs bg, marginal again?
#W = 1613300, p-value = 0.06797
wilcox.test(((mtbgauto$dN/mtbgauto$Non.sites)/(mtbgauto$dS/mtbgauto$Syn.sites)),((mtsaa$dN/mtsaa$Non.sites)/(mtsaa$dS/mtsaa$Syn.sites)))
#manduca, shared vs bg, not-significant
#W = 573340, p-value = 0.1307

par(mfrow=c(2,2))
# bottom, left, top, right margins
par(mai=c(0.12,0.48,0.12,0.1))
#Pn/Ps time
boxplot(((dtbgauto$pN/dtbgauto$Non.sites)/(dtbgauto$pS/dtbgauto$Syn.sites)),outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3),main ="Monarch Pn/Ps",ylim=c(0,0.40))
boxplot(((dtsae$pN/dtsae$Non.sites)/(dtsae$pS/dtsae$Syn.sites)), add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(((dtsas$pN/dtsas$Non.sites)/(dtsas$pS/dtsas$Syn.sites)), add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(((dtsaa$pN/dtsaa$Non.sites)/(dtsaa$pS/dtsaa$Syn.sites)), add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)
text(1.5, 0.29,"*",cex=1.8)
text(1.5, 0.275,"*",cex=1.8)
text(2, 0.295,"*",cex=1.8)
text(2, 0.275,"*",cex=1.8)
text(2, 0.255,"*",cex=1.8)
text(2, 0.315,"*",cex=1.8)

boxplot(((mtbgauto$pN/mtbgauto$Non.sites)/(mtbgauto$pS/mtbgauto$Syn.sites)),outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3), main="Carolina sphinx Pn/Ps")
boxplot(((mtsae$pN/mtsae$Non.sites)/(mtsae$pS/mtsae$Syn.sites)), add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(((mtsas$pN/mtsas$Non.sites)/(mtsas$pS/mtsas$Syn.sites)), add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(((mtsaa$pN/mtsaa$Non.sites)/(mtsaa$pS/mtsaa$Syn.sites)), add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)
text(2, 1.05,"*",cex=1.8)


#Dn/Ds time?
boxplot(((dtbgauto$dN/dtbgauto$Non.sites)/(dtbgauto$dS/dtbgauto$Syn.sites)),outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3))
boxplot(((dtsae$dN/dtsae$Non.sites)/(dtsae$dS/dtsae$Syn.sites)), add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(((dtsas$dN/dtsas$Non.sites)/(dtsas$dS/dtsas$Syn.sites)), add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(((dtsaa$dN/dtsaa$Non.sites)/(dtsaa$dS/dtsaa$Syn.sites)), add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)
text(1.3)

boxplot(((mtbgauto$dN/mtbgauto$Non.sites)/(mtbgauto$dS/mtbgauto$Syn.sites)),outline=F, col=, notch=T, las=1,cex.axis=1.3,xlim=c(0.8,3))
boxplot(((mtsae$dN/mtsae$Non.sites)/(mtsae$dS/mtsae$Syn.sites)), add=T, at = 1.5,outline=F, col="blue", notch=T, las=1,cex.axis=1.3)
boxplot(((mtsas$dN/mtsas$Non.sites)/(mtsas$dS/mtsas$Syn.sites)), add=T, at = 2,outline=F, col="purple", notch=T, las=1,cex.axis=1.3)
boxplot(((mtsaa$dN/mtsaa$Non.sites)/(mtsaa$dS/mtsaa$Syn.sites)), add=T, at = 2.5,outline=F, col="red", notch=T, las=1,cex.axis=1.3)



#what if we do pn/ps all together?
wilcox.test(((dtbgauto$pN/dtbgauto$Non.sites)/(dtbgauto$pS/dtbgauto$Syn.sites)),((dtsae$pN/dtsae$Non.sites)/(dtsae$pS/dtsae$Syn.sites)))
#danaus Pn/Ps, eu vs bg, it's real significant
#W = 1269000, p-value = 0.001815
wilcox.test(((dtbgauto$pN/dtbgauto$Non.sites)/(dtbgauto$pS/dtbgauto$Syn.sites)),((dtsas$pN/dtsas$Non.sites)/(dtsas$pS/dtsas$Syn.sites)))
#danaus pn/ps, shared vs bg, super the significantest
#W = 1477100, p-value = 3.566e-08
wilcox.test(((dtbgauto$pN/dtbgauto$Non.sites)/(dtbgauto$pS/dtbgauto$Syn.sites)),((dtsaa$pN/dtsaa$Non.sites)/(dtsaa$pS/dtsaa$Syn.sites)))
#danaus pn/ps, shared vs bg, not-significant!
#W = 273590, p-value = 0.3951

wilcox.test(((mtbgauto$pN/mtbgauto$Non.sites)/(mtbgauto$pS/mtbgauto$Syn.sites)),((mtsae$pN/mtsae$Non.sites)/(mtsae$pS/mtsae$Syn.sites)))
#manduca pn/ps, eu vs bg, is also marginal?
#W = 992280, p-value = 0.09998
wilcox.test(((mtbgauto$pN/mtbgauto$Non.sites)/(mtbgauto$pS/mtbgauto$Syn.sites)),((mtsas$pN/mtsas$Non.sites)/(mtsas$pS/mtsas$Syn.sites)))
#manduca shared vs bg, is significant?
#W = 1623000, p-value = 0.04403
wilcox.test(((mtbgauto$pN/mtbgauto$Non.sites)/(mtbgauto$pS/mtbgauto$Syn.sites)),((mtsaa$pN/mtsaa$Non.sites)/(mtsaa$pS/mtsaa$Syn.sites)))
#manduca Dn, shared vs bg, not-significant but marginal
#W = 576470, p-value = 0.1062


#Ds time
wilcox.test((dtbgauto$dS/dtbgauto$Syn.sites),(dtsae$dS/dtsae$Syn.sites))
#danaus Dn, eu vs bg, it's marginally not significant!
#W = 1056000, p-value = 0.09275
wilcox.test((dtbgauto$dS/dtbgauto$Syn.sites),(dtsas$dS/dtsas$Syn.sites))
#danaus Dn, shared vs bg, super not significant
#W = 1209000, p-value = 0.7665
wilcox.test((dtbgauto$dS/dtbgauto$Syn.sites),(dtsaa$dS/dtsaa$Syn.sites))
#danaus Dn, shared vs bg, not-significant!
#W = 279420, p-value = 0.2594

wilcox.test((mtbgauto$dS/mtbgauto$Syn.sites),(mtsae$dS/mtsae$Syn.sites))
#manduca Ds, eu vs bg, it's significant
#W = 963410, p-value = 0.3596
wilcox.test((mtbgauto$dS/mtbgauto$Syn.sites),(mtsas$dS/mtsas$Syn.sites))
#manduca Ds, shared vs bg, marginally not significant
#W = 1437100, p-value = 0.103
wilcox.test((mtbgauto$dS/mtbgauto$Syn.sites),(mtsaa$dS/mtsaa$Syn.sites))
#manduca Ds, shared vs bg, not-significant!
#W = 495180, p-value = 0.2653

#we also need to sort to get rid of genes that lack variation in certain places that explode the calculations
#specifically if pS or dS are zero things get wonky
mtbgauto<-mtbgauto[which(mtbgauto$pS!=0 & mtbgauto$dS!=0),]
dtbgauto<-dtbgauto[which(dtbgauto$pS!=0 & dtbgauto$dS!=0),]
mtsauto<-mtsauto[which(mtsauto$pS!=0 & mtsauto$dS!=0),]
#here we subset out those 40 orphan commixed proteins
mtsauto<-mtsauto[which(!(mtsauto$Commixed=="1"&mtsauto$Eupyrene=="0"&mtsauto$Apyrene=="0")),]
dtsauto<-dtsauto[which(dtsauto$pS!=0 & dtsauto$dS!=0),]
#this gives us:
#Species: Unique  Ortho
#Manduca: 401     236
#Danaus:  294     216
mtsaUn<-mtsauto[which(mtsauto$Danaus.Sperm.Ortho==0),]
mtsaOr<-mtsauto[which(mtsauto$Danaus.Sperm.Ortho==1),]
dtsaUn<-dtsauto[which(dtsauto$Manduca.Sperm.Ortho==0),]
dtsaOr<-dtsauto[which(dtsauto$Manduca.Sperm.Ortho==1),]

#what if we naively partition alpha by these distinctions?
NItgCalc<-function(dn,ds,pn,ps)
{
	unbiasedNI<-sum((ds*pn)/(ps+ds))/sum((ps*dn)/(ps+ds))
	return(unbiasedNI)
}

#here's some background vs sperm proteome stuff
mbg<-1-NItgCalc(mtbgauto$dN,mtbgauto$dS,mtbgauto$pN,mtbgauto$pS)
dbg<-1-NItgCalc(dtbgauto$dN,dtbgauto$dS,dtbgauto$pN,dtbgauto$pS)
mbs<-1-NItgCalc(mtsauto$dN,mtsauto$dS,mtsauto$pN,mtsauto$pS)
dbs<-1-NItgCalc(dtsauto$dN,dtsauto$dS,dtsauto$pN,dtsauto$pS)
#Point est: BG            Sperm
#Manduca:   -0.0767369   -0.05070336
#Danaus:    -0.1383413    0.06568262


#Fig 1 
#first panel whole sperm vs genome background both species

#we'll do autosomal
par(mfrow(1,1))
plot(-9,9,ylim=c(-0.21,0.2),xlim=c(0.05,1.9),ylab="",
xaxt="n",xlab="",main="Molecular evolution of the sperm proteome", cex.lab=2, cex.axis=1.3, las=1, cex.main=1.5)
mtext(expression(alpha),side=2,las=1,line=3,cex=2)

abline(h=0)
#here is the genome background
segments(0.45, -0.05847823, 0.65, -0.05847823, lwd=3, col="dimgray")
segments(0.55, -0.05847823, 0.55, -0.09532876, lwd=3, col="dimgray")
segments(0.45, -0.09532876, 0.65, -0.09532876, lwd=3, col="dimgray")
points(0.55,-0.07673690,pch=0, cex=1.9, col="dimgray")
points(0.55,-0.07673690,pch=0, cex=2, col="dimgray")
points(0.55,-0.07673690,pch=0, cex=2.1, col="dimgray")
points(0.55,-0.07673690,pch=0, cex=2.2, col="dimgray")

#here is all autosomal sperm genes
segments(0.15, -0.12922974, 0.35, -0.12922974, lwd=3, col="dimgray")
segments(0.25, -0.12922974, 0.25, 0.02966341, lwd=3, col="dimgray")
segments(0.15, 0.02966341 , 0.35, 0.02966341, lwd=3, col="dimgray")
points(0.25,-0.04407427,pch=16, cex=2.5, col="dimgray")

segments(0.25,-0.19,0.55,-0.19)
segments(0.25,-0.17,0.25,-0.19)
segments(0.55,-0.17,0.55,-0.19)
text(0.4,-0.206,"p = 0.409",cex=1.4)


#now for Monarchs
#here is the genome background
# auto :   1.043279   1.070868    1.098379
segments(1.45, 1-1.118138 , 1.65, 1-1.118138, lwd=3, col="orange")
segments(1.55, 1-1.118138 , 1.55, 1-1.158768, lwd=3, col="orange")
segments(1.45, 1-1.158768, 1.65, 1-1.158768, lwd=3, col="orange")
points(1.55,1-1.138341,pch=0, cex=1.9, col="orange")
points(1.55,1-1.138341,pch=0, cex=2, col="orange")
points(1.55,1-1.138341,pch=0, cex=2.1, col="orange")
points(1.55,1-1.138341,pch=0, cex=2.2, col="orange")

#here is all autosomal sperm genes

#sperm: 0.8349017 0.9342584 1.0439800
segments(1.15, 1-0.8172010, 1.35, 1-0.8172010, lwd=3, col="orange")
segments(1.25, 1-0.8172010, 1.25, 1-1.0025506, lwd=3, col="orange")
segments(1.15, 1-1.0025506, 1.35, 1-1.0025506, lwd=3, col="orange")
points(1.25,1-0.9342584,pch=16, cex=2.5, col="orange")

segments(1.25,-0.19,1.55,-0.19)
segments(1.25,-0.17,1.25,-0.19)
segments(1.55,-0.17,1.55,-0.19)
text(1.4,-0.206,"p < 0.00001", cex=1.4,font=2)
text(1.25,-0.02,"*",cex=4)
text(1.25,-0.05,"*",cex=4)
text(1.25,-0.08,"*",cex=4)
text(1.25,-0.11,"*",cex=4)
legend(0.05,0.2,pch=c(16,0),c("Sperm proteome","Background genome"),col=c("black","black"),cex=1.2)
mtext(side=1, adj=0.1, expression(italic("Manduca sexta")), cex=1.4, line=1)
mtext(side=1, adj=0.8, expression(italic("Danaus plexippus")), cex = 1.4, line=1)
mtext(side=3, adj=-0.1, "A", cex = 2, line=1)


#second panel, decomposing alpha into pn,ps,dn,ds
#we can do signifcance testing right here
#Polymorphism
wilcox.test((mtbgauto$pN/mtbgauto$Non.sites),(mtsauto$pN/mtsauto$Non.sites))
#manduca Pn
#W = 3014100, p-value = 0.5964
wilcox.test((dtbgauto$pN/dtbgauto$Non.sites),(dtsauto$pN/dtsauto$Non.sites))
#danaus Pn, this one's the significant one
#W = 3062400, p-value = 3.224e-11
wilcox.test((mtbgauto$pS/mtbgauto$Syn.sites),(mtsauto$pS/mtsauto$Syn.sites))
#manduca PS
#W = 2879300, p-value = 0.183
wilcox.test((dtbgauto$pS/dtbgauto$Syn.sites),(dtsauto$pS/dtsauto$Syn.sites))
#danaus PS
#W = 2684200, p-value = 0.272

#Divergence
wilcox.test((mtbgauto$dN/mtbgauto$Non.sites),(mtsauto$dN/mtsauto$Non.sites))
#manduca Dn
#W = W = 3068300, p-value = 0.2009
wilcox.test((dtbgauto$dN/dtbgauto$Non.sites),(dtsauto$dN/dtsauto$Non.sites))
#danaus Dn, this is marginal
#W = 2506400, p-value = 0.13
wilcox.test((mtbgauto$dS/mtbgauto$Syn.sites),(mtsauto$dS/mtsauto$Syn.sites))
#manduca dS
#W = 2895700, p-value = 0.2686
wilcox.test((dtbgauto$dS/dtbgauto$Syn.sites),(dtsauto$dS/dtsauto$Syn.sites))
#danaus dS
#W = 2544400, p-value = 0.3437


#bottom left, Pn/Ps for both
par(mfrow=c(2,2))
# bottom, left, top, right margins
par(mai=c(0.22,0.94,0.44,0.18))
#pN
boxplot(mtsauto$pN/mtsauto$Non.sites,xlim=c(0.5,4),ylim=c(0,0.038),outline=F, col="dimgray", notch=T, main = "Non-syn. Polymorphism", ylab="",cex.lab=1.5, las=1, cex.main=1.5,cex.axis=1.3)
segments(1,0.035,1.5,0.035)
segments(1.5,0.033,1.5,0.035)
segments(1,0.033,1,0.035)
text(1.2,0.0365, "p = 0.596", cex=1.4)
boxplot(mtbgauto$pN/mtbgauto$Non.sites, add=T, at = 1.5, outline=F, notch=T, las=1,cex.axis=1.3)
boxplot(dtsauto$pN/dtsauto$Non.sites,add=T, at = 3,outline=F, col="orange", notch=T, las=1,cex.axis=1.3)
boxplot(dtbgauto$pN/dtbgauto$Non.sites, add=T, at = 3.5, outline=F, notch=T, las=1,cex.axis=1.3)
segments(3,0.035,3.5,0.035)
segments(3.5,0.033,3.5,0.035)
segments(3,0.033,3,0.035)
text(3.2,0.0365, "p < 0.00001", cex=1.4)
text(3,0.025, "*", cex=4)
text(3,0.029, "*", cex=4)
mtext(side=3, adj=-0.2, "B", cex = 2, line=0)
mtext(side=2, "Variants per Site", cex = 1.5, line=4)
#pS
par(mai=c(0.22,0.42,0.42,0.12))
boxplot(mtsauto$pS/mtsauto$Syn.sites,xlim=c(0.5,4),ylim=c(0,0.29),outline=F, col="dimgray", notch=T, main = "Syn. Polymorphism",las=1, cex.main=1.5,cex.axis=1.3)
boxplot(mtbgauto$pS/mtbgauto$Syn.sites, add=T, at = 1.5, outline=F, notch=T, las=1,cex.axis=1.3)
segments(1,0.15,1.5,0.15)
segments(1.5,0.13,1.5,0.15)
segments(1,0.13,1,0.15)
text(1.2,0.165, "p = 0.183", cex=1.4)
boxplot(dtsauto$pS/dtsauto$Syn.sites,add=T, at = 3,outline=F, col="orange", notch=T, las=1,cex.axis=1.3)
boxplot(dtbgauto$pS/dtbgauto$Syn.sites, add=T, at = 3.5, outline=F, notch=T, las=1,cex.axis=1.3)
segments(3,0.26,3.5,0.26)
segments(3.5,0.24,3.5,0.26)
segments(3,0.24,3,0.26)
text(3.2,0.275, "p = 0.272", cex=1.4)

#dN
# bottom, left, top, right margins
par(mai=c(0.52,0.94,0.42,0.18))
boxplot(mtsauto$dN/mtsauto$Non.sites,xlim=c(0.5,4),ylim=c(0,0.024),outline=F, col="dimgray", notch=T, main = "Non-syn. Divergence",ylab="",cex.lab=1.5,las=1,cex.main=1.5, ytics=seq(0,1,0.05),cex.axis=1.3)
boxplot(mtbgauto$dN/mtbgauto$Non.sites, add=T, at = 1.5, outline=F, notch=T, las=1,cex.axis=1.3)
segments(1,0.022,1.5,0.022)
segments(1.5,0.02,1.5,0.022)
segments(1,0.02,1,0.022)
text(1.2,0.0235, "p = 0.201", cex=1.4)
boxplot(dtsauto$dN/dtsauto$Non.sites,add=T, at = 3,outline=F, col="orange", notch=T, las=1,cex.axis=1.3)
boxplot(dtbgauto$dN/dtbgauto$Non.sites, add=T, at = 3.5, outline=F, notch=T, las=1,cex.axis=1.3)
segments(3,0.019,3.5,0.019)
segments(3.5,0.017,3.5,0.019)
segments(3,0.017,3,0.019)
text(3.2,0.0205, "p = 0.130", cex=1.4)
mtext(side=1, adj=0, expression(italic("M. sexta")), cex=1.2, line=1)
mtext(side=1, adj=1, expression(italic("D. plexippus")), cex = 1.2, line=1)
mtext(side=2, "Variants per Site", cex = 1.5, line=4)

#dS
par(mai=c(0.52,0.42,0.42,0.12))
boxplot(mtsauto$dS/mtsauto$Syn.sites,xlim=c(0.5,4),ylim=c(0,0.22),outline=F, col="dimgray", notch=T, main = "Syn. Divergence",las=1, cex.main=1.5,cex.axis=1.3)
boxplot(mtbgauto$dS/mtbgauto$Syn.sites, add=T, at = 1.5, outline=F, notch=T, las=1,cex.axis=1.3)
segments(1,0.11,1.5,0.11)
segments(1.5,0.09,1.5,0.11)
segments(1,0.09,1,0.11)
text(1.2,0.12, "p = 0.269", cex=1.4)
boxplot(dtsauto$dS/dtsauto$Syn.sites,add=T, at = 3,outline=F, col="orange", notch=T, las=1,cex.axis=1.3)
boxplot(dtbgauto$dS/dtbgauto$Syn.sites, add=T, at = 3.5, outline=F, notch=T, las=1,cex.axis=1.3)
segments(3,0.19,3.5,0.19)
segments(3.5,0.17,3.5,0.19)
segments(3,0.17,3,0.19)
text(3.2,0.2, "p = 0.344", cex=1.4)
mtext(side=1, adj=0, expression(italic("M. sexta")), cex=1.2, line=1)
mtext(side=1, adj=1, expression(italic("D. plexippus")), cex = 1.2, line=1)


#these don't add much visually
#pN/pS
#boxplot((mtsauto$pN/mtsauto$Non.sites)/(mtsauto$pS/mtsauto$Syn.sites),xlim=c(0.7,4),ylim=c(0,1.54),outline=F, col="dimgray", notch=T, main = "Scaled Polymorphism")
#boxplot((mtbgauto$dN/mtbgauto$Non.sites)/(mtbgauto$pS/mtbgauto$Syn.sites), add=T, at = 1.4, outline=F, notch=T)
#boxplot((dtsauto$dN/dtsauto$Non.sites)/(dtsauto$pS/dtsauto$Syn.sites),add=T, at = 3,outline=F, col="orange", notch=T)
#boxplot((dtbgauto$dN/dtbgauto$Non.sites)/(dtbgauto$pS/dtbgauto$Syn.sites), add=T, at = 3.5, outline=F, notch=T)
# Dn/Ds for both
#dN/dS
#boxplot((mtsauto$dN/mtsauto$Non.sites)/(mtsauto$dS/mtsauto$Syn.sites),xlim=c(0.7,4),ylim=c(0,1),outline=F, col="dimgray", notch=T, main = "Scaled Divergence")
#boxplot((mtbgauto$dN/mtbgauto$Non.sites)/(mtbgauto$dS/mtbgauto$Syn.sites), add=T, at = 1.4, outline=F, notch=T)
#boxplot((dtsauto$dN/dtsauto$Non.sites)/(dtsauto$dS/dtsauto$Syn.sites),add=T, at = 3,outline=F, col="orange", notch=T)
#boxplot((dtbgauto$dN/dtbgauto$Non.sites)/(dtbgauto$dS/dtbgauto$Syn.sites), add=T, at = 3.5, outline=F, notch=T)



#Fig 2
#ortho vs not comparison
#bootstrap p vals can be found in spermsharedornot script

plot(-9,9,ylim=c(-0.2,.41),xlim=c(.1,.7),ylab="",
xaxt="n",xlab="",main="Adaptive evolution accounting for orthology in the sperm proteins",cex.lab=2, las=1)
mtext(expression(alpha),side=2,las=1,line=3,cex=2)
abline(h=0)
#Manduca
#unique
points(0.2, -0.02897364, pch=1, col="dimgray",cex=1.8 )
points(0.2, -0.02897364, pch=1, col="dimgray",cex=2 )
points(0.2, -0.02897364, pch=1, col="dimgray",cex=2.2 )
segments(0.2, 0.06224726, 0.2, -0.14513430, col="dimgray",lwd=3)
segments(0.15, 0.06224726, 0.25, 0.06224726, col="dimgray",lwd=3)
segments(0.15, -0.14513430, 0.25, -0.14513430, col="dimgray",lwd=3)
#shared
points(0.3,  -0.07190798, pch=19, col="dimgray", cex=2 )
segments(0.3, 0.03045990, 0.3, -0.18798665, col="dimgray",lwd=3)
segments(0.25, -0.18798665, 0.35, -0.18798665, col="dimgray",lwd=3)
segments(0.25, 0.03045990, 0.35, 0.03045990, col="dimgray",lwd=3)

#the p-val stuff
segments(0.2, 0.08, 0.2, 0.1, col="black",lwd=2)
segments(0.3, 0.08, 0.3, 0.1, col="black",lwd=2)
segments(0.2, 0.1, 0.3, 0.1, col="black",lwd=2)
text(0.25, 0.12, "p = 0.617")

#Danaus
points(0.6,  0.01140614, pch=1, col="orange",cex=1.8)
points(0.6,  0.01140614, pch=1, col="orange",cex=2)
points(0.6,  0.01140614, pch=1, col="orange",cex=2.2)
segments(0.6, 0.14263047, 0.6, -0.13503302, col="orange",lwd=3)
segments(0.55, 0.14263047, 0.65, 0.14263047, col="orange",lwd=3)
segments(0.55, -0.13503302, 0.65,-0.13503302, col="orange",lwd=3)
points(0.5, 0.2395017, pch=19, col="orange", cex=2)
segments(0.5, 0.3386478, 0.5, 0.1213998, col="orange",lwd=3)
segments(0.45, 0.1213998, 0.55, 0.1213998, col="orange",lwd=3)
segments(0.45, 0.3386478, 0.55, 0.3386478, col="orange",lwd=3)
#the p-val stuff
segments(0.5, 0.36, 0.5, 0.39, col="black",lwd=2)
segments(0.6, 0.36, 0.6, 0.39, col="black",lwd=2)
segments(0.5, 0.39, 0.6, 0.39, col="black",lwd=2)
text(0.55, 0.41, "p = 0.0372")
text(0.6, 0.2395, "*", cex=4)
mtext(side=1, adj=0, expression(italic("      Manduca sexta")), cex=1.5, line=1)
mtext(side=1, adj=1, expression(italic("Danaus plexippus")), cex = 1.5, line=1)
legend(0.1,0.4, c("Unique proteins","Sperm homologs"),pch=c(1,19),col=c("black","black"),cex=1.4)
mtext(side=3, adj=-0.1, "B", cex = 2, line=1)

#Fig 3 DFEs
#(see SfsVisualizationAndComparisons.R)


#Fig 4 Dimorphism
#p-vals come from bootstrapping in sharedornot file
par(mfrow=c(1,2))
# bottom, left, top, right margins
par(mai=c(0.12,0.82,0.32,0.14))
#per the top of the figure, Manduca should go first on the left
plot(-9,9,ylim=c(-0.35,0.5),xlim=c(0.08,0.8),ylab="",
xaxt="n",xlab="",main="Evolution of dimorphic sperm in Carolina sphinx",cex.lab=1.5,las=1,
panel.first = rect(-9,-0.09532876,9,-0.05847823, col='lightgray', border=NA))
abline(h=0)
mtext(expression(alpha),side=2,las=1,line=3,cex=2)

#here is the genome background
#this is correct
segments(-1, -0.09532876, 2, -0.09532876, lty=4, col="black")
#segments(-1, -0.07673690, 2, -0.07673690, lty=3, col="black")
segments(-1, -0.05847823, 2, -0.05847823, lty=4, col="black")

#we  lead with apyrene
points(0.75,1-1.0806311,pch=16, cex=2,col="firebrick3")
segments(0.75,1-0.9029121, 0.75, 1-1.3256465,lwd=2, col="firebrick3")
segments(0.7,1-0.9029121, 0.8, 1-0.9029121,lwd=2, col="firebrick3")
segments(0.7,1-1.3256465, 0.8, 1-1.3256465,lwd=2, col="firebrick3")

#then shared
points(0.5,1-1.03721, pch=17,cex=2,col="maroon4")
segments(0.5, 1-0.9413469, 0.5, 1-1.1568803,lwd=2, col="maroon4")
segments(0.45,1-1.1568803 , 0.55, 1-1.1568803,lwd=2, col="maroon4")
segments(0.45,1-0.9413469,0.55, 1-0.9413469,lwd=2, col="maroon4")

#then eupyrene
points(0.25,1-1.0369739,pch=15, col="deepskyblue",cex=2)
segments(0.25,1-0.9100191,0.25, 1- 1.1901942,lwd=2, col="deepskyblue")
segments(0.2, 1-1.1901942, 0.3, 1-1.1901942,lwd=2, col="deepskyblue")
segments(0.2,1-0.9100191 ,0.3, 1-0.9100191,lwd=2, col="deepskyblue")

legend(.075,0.5,pch=c(15,17,16),c("Eupyrene-specific","Shared","Apyrene-specific"), col = c("deepskyblue","maroon4","firebrick3"))



#followed by monarchs on the right
plot(-9,9,ylim=c(-0.35,0.5),xlim=c(0.08,0.8),ylab="",
xaxt="n",xlab="",main="Evolution of dimorphic sperm in monarch",cex.lab=1.5,las=1,
panel.first = rect(-9,-0.1181787,9,-0.1587623, col='lightgray', border=NA))
mtext(expression(alpha),side=2,las=1,line=3,cex=2)
#here are the limits for the genome background
abline(h=0)
segments(-1, -0.1181787, 2, -0.1181787, lty=4, col="black")
#segments(-1, 1-1.138341, 2, 1-1.138341, lty=3, col="black")
segments(-1, -0.1587623, 2, -0.1587623, lty=4, col="black")

segments(1.45, 1-1.118138 , 1.65, 1-1.118138, lwd=3, col="black")
segments(1.55, 1-1.118138 , 1.55, 1-1.158768, lwd=3, col="black")
segments(1.45, 1-1.158768, 1.65, 1-1.158768, lwd=3, col="black")
points(1.55,1-1.138341,pch=15, cex=2, col="black")

#here are the eupyrene + CIs
points(0.25,1-0.8747993 ,pch=15,cex=2,col="deepskyblue")
segments(0.25,1-0.7927828,0.25,1-0.9625260,lwd=2,col="deepskyblue")
segments(0.2,1-0.9625260,0.3,1-0.9625260,lwd=2,col="deepskyblue")
segments(0.2,1-0.7927828,0.3,1-0.7927828,lwd=2,col="deepskyblue")
#the shared set
points(0.5,1-0.9195244, pch=17, cex=2, col="maroon4")
segments(0.5, 1-0.7708373, 0.5, 1-1.0758177,lwd=2, col="maroon4")
segments(0.45,1-1.0758177 , 0.55, 1-1.0758177,lwd=2, col="maroon4")
segments(0.45, 1-0.7708373,0.55, 1-0.7708373,lwd=2, col="maroon4")
#apyrene
points(0.75,1-1.0250200,pch=16,cex=2, col="firebrick3")
segments(0.75,1-0.9123317, 0.75, 1-1.1545854,lwd=2, col="firebrick3")
segments(0.7,1-0.9123317, 0.8, 1-0.9123317,lwd=2, col="firebrick3")
segments(0.7,1-1.1545854, 0.8, 1-1.1545854,lwd=2, col="firebrick3")

legend(.075,0.5,pch=c(15,17,16),c("Eupyrene-specific","Shared","Apyrene-specific"), col = c("deepskyblue","maroon4","firebrick3"))

#within-sperm comparisons
#segments(0.25,0.26,0.48,0.26)
#segments(0.52,0.26,0.75,0.26)
segments(0.25,0.31,0.75,0.31)
#segments(0.25,0.23,0.25,0.26)
#segments(0.48,0.23,0.48,0.26)
#segments(0.52,0.23,0.52,0.26)
#segments(0.75,0.23,0.75,0.26)
segments(0.25,0.28,0.25,0.31)
segments(0.75,0.28,0.75,0.31)
#text(0.62,0.27,"p = 0.633")
#text(0.37,0.27,"p = 0.997")
text(0.48,0.32,"p = 0.098")

#sperm type vs the world
#ap
#segments(0.69,1-1.0250200,0.69,1-1.138341)
#segments(0.69,1-1.0250200,0.71,1-1.0250200)
#segments(0.69,1-1.138341,0.71,1-1.138341)
#text(0.6,-0.09,"p = 0.559")

#shared
segments(0.44,1-0.9195244,0.44,1-1.138341)
segments(0.44,1-0.9195244,0.46,1-0.9195244)
segments(0.44,1-1.138341,0.46,1-1.138341)
text(0.35,-0.04, "p = 0.010")
text(0.5,-0.1,"*",cex=4)
#eu
segments(0.19,1-0.9195244,0.19,1-1.138341)
segments(0.19,1-0.9195244,0.21,1-0.9195244)
segments(0.19,1-1.138341,0.21,1-1.138341)
text(0.12,0.02, "p < 0.001")
text(0.25,-0.06,"*",cex=4)
text(0.25,-0.02,"*",cex=4)

