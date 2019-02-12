###Expression Data Analysis
NItgCalc<-function(dn,ds,pn,ps)
{
	unbiasedNI<-sum((ds*pn)/(ps+ds))/sum((ps*dn)/(ps+ds))
	return(unbiasedNI)
}

#Here we will tackle the question of tissue specificity of gene expression
#with two ultimate goals:
#1. to subset the sperm proteome based on genes that appear to be truly sperm (testes) specific
#2. to test the Dapper & Wade prediction that sex-limited expression has a noticeable impact on 
#gene evolution
#we start with the gold standard Kanzen files and will append to them as necessary
library("data.table")
manmaster<-read.csv(file="/Users/Andrew/Documents/ManducaSextaFullGenomeWithZorA.csv", header=T)
man_append<-read.csv(file="/Users/Andrew/Documents/manducaorthosfordanaus.csv",header=T)
colnames(man_append)<-c("Gene","Chromosome","Danaus.Sperm.Ortho")
manmaster<-as.data.frame(manmaster, stringsAsFactors=F)
danmaster<-read.csv(file="/Users/Andrew/Downloads/evoconfiles/NewKanzen.csv",header=T)
dan_append<-read.csv(file="/Users/Andrew/Documents/danausorthologsformanduca.csv",header=T)
colnames(dan_append)<-c("Gene","ZorA","Chromosome","Manduca.Sperm.Ortho")
danmaster<-as.data.frame(danmaster,stringsAsFactors=F)
mantranscend<-merge(manmaster,man_append,by=c("Gene","Chromosome"),all.x=T, all.y=T)
dantranscend<-merge(danmaster,dan_append,by=c("Gene","Chromosome","ZorA"),all.x=T, all.y=T)
#first let's tackle the tau data...it's not ultimately that useful as tau will be large
#if a gene has specific expression in ANY tissue...
#but we can use it to test the prediction that genes with high specificity evolve more adaptively
man_tau<-read.csv(file="/Users/Andrew/Documents/ExpressionWork/Danaus_Tissue_Specificity/Output/TauGenesManduca.csv",header=T)
dan_tau<-read.csv(file="/Users/Andrew/Documents/ExpressionWork/Danaus_Tissue_Specificity/Output/TauGenes.csv",header=T)
mantran2<-merge(mantranscend,man_tau,by="Gene",all.x=T, all.y=T)
dantran2<-merge(dantranscend,dan_tau,by="Gene",all.x=T, all.y=T)
#we also need a simple measure of alpha
dan_alpha_simple<-1-((dantran2$dS*dantran2$pN)/(dantran2$dN*dantran2$pS))
man_alpha_simple<-1-((mantran2$dS*mantran2$pN)/(mantran2$dN*mantran2$pS))
par(mfrow=c(1,2))
plot(mantran2$Tau,man_alpha_simple,ylim=c(-10,1))
plot(dantran2$Tau,dan_alpha_simple,ylim=c(-10,1))
#Yeah, so that wasn't very helpful after all...what about Dn/Ds?
dan_dnds<-(dantran2$dN/dantran2$Non.sites)/(dantran2$dS/dantran2$Syn.sites)
man_dnds<-(mantran2$dN/mantran2$Non.sites)/(mantran2$dS/mantran2$Syn.sites)
par(mfrow=c(1,2))
plot(mantran2$Tau,man_dnds,ylim=c(0,10.5))
plot(dantran2$Tau,dan_dnds,ylim=c(0,10.5))
##okay enough silly buggers, let's move onto the more advanced metrics of tissue specificity
#Namely, we want the SPM values, which ares specific to tissue
dan_spm<-read.csv(file="/Users/Andrew/Documents/ExpressionWork/Danaus_Tissue_Specificity/SpmGenesAJM.csv",header=T)
#we need to get the names back onto the genes in Manduca, so I've hacked it back into Megan's file
man_spm<-read.csv(file="/Users/Andrew/Documents/ExpressionWork/Danaus_Tissue_Specificity/SpmGenesManduca.csv",header=T)
mantran3<-merge(mantran2,man_spm,by="Gene",all.x=T, all.y=T)
mantran4<-cbind(mantran3,man_dnds)
write.csv(mantran4, file="/Users/Andrew/Documents/ExpressionWork/Andrews_analysis/MsextaSHINKANZEN2.csv", row.names=F,quote=F)
#gods forgive me i'm going to edit the CV in excel for names and column order
dantran3<-merge(dantran2,dan_spm,by="Gene",all.x=T, all.y=T)
dantran4<-cbind(dantran3,dan_dnds)
write.csv(dantran4, file="/Users/Andrew/Documents/ExpressionWork/Andrews_analysis/DplexSHINKANZEN2.csv", row.names=F,quote=F)
###okay there was some editing in excel to rename and move around columns
##the real dark souls starts here:
#E
#	x
#		p
#			r
#				e
#					s
#						s
#							
#							Y
#								o
#									s
#										e
#											l	
#												f
#
manfull<-read.csv(file="/Users/Andrew/Documents/ExpressionWork/Andrews_analysis/MsextaSHINKANZEN2.csv",stringsAsFactors=F)
danfull<-read.csv(file="/Users/Andrew/Documents/ExpressionWork/Andrews_analysis/DplexSHINKANZEN.csv",stringsAsFactors=F)



############################################working line#########################################################
#okay, so Aloy doesn't want to give up his full RNA-seq dataset, but would be willing to spare males
#what if we just do male in both species, ignoring female-limited expression and saving it for another day

tstset<-danfull[which(danfull$dS !=0 & danfull$pS !=0 ),]
spmind<-seq(0,1,0.05)
slide_alphas<-rep(0,length(spmind))
for (i in 1:length(spmind))
{
passers<-tstset[which(tstset$head.male >spmind[i] |
tstset$mg.male >spmind[i] | tstset$thorax.male>spmind[i]  | tstset$fpkm.testes.avg >spmind[i] | tstset$fpkm.ag.avg>spmind[i]),]
slide_alphas[i]<-1-NItgCalc(passers$dN,passers$dS,passers$pN,passers$pS)
}
plot(spmind,slide_alphas,ylab="",xlab="SPM Threshold",cex.lab=2,
main="Effect of increasing tissue specificity on monarch results",ylim=c(-0.15,0.46), las=1, xlim=c(0,0.95),cex=2, cex.axis=1.5, cex.main=1.5)
mtext(expression(alpha),side=2,las=1,line=3,cex=2)
abline(h=0)


###what about if we do the same, but for sperm?
tstset<-danfull[which(danfull$dS !=0 & danfull$pS !=0 ),]
tstset<-tstset[which(tstset$gene.class=="eupyrene" | tstset$gene.class=="apyrene"),]
spmind<-seq(0,1,0.05)
slide_alphas<-rep(0,length(spmind))
for (i in 1:length(spmind))
{
passers<-tstset[which(tstset$fpkm.testes.avg >spmind[i] | tstset$fpkm.ag.avg>spmind[i]),]
slide_alphas[i]<-1-NItgCalc(passers$dN,passers$dS,passers$pN,passers$pS)
}
points(spmind,slide_alphas,pch=16,col="blue",cex=2)
#legend(0,0.25,pch=c(1,16),col=c("black","blue"),c("all genes","sperm genes"))


tstset<-danfull[which(danfull$dS !=0 & danfull$pS !=0 ),]
tstset<-tstset[which(tstset$gene.class=="non-sperm"),]
spmind<-seq(0,1,0.05)
slide_alphas<-rep(0,length(spmind))

for (i in 1:length(spmind))
{
passers<-tstset[which(tstset$fpkm.testes.avg >spmind[i] | tstset$fpkm.ag.avg>spmind[i]),]
slide_alphas[i]<-1-NItgCalc(passers$dN,passers$dS,passers$pN,passers$pS)
}
points(spmind,slide_alphas,pch=18,col="red",cex=2)

legend(-0.01,0.455,pch=c(1,16,18),col=c("black","blue","red"),c("all genes","sperm genes","male-limited (non-sperm)"),pt.cex=2,cex=1.5)
mtext(side=3, adj=-0.05, "B", cex = 2, line=1)


#
#now Manduca
#SPM Cut-off graph
#
tstset<-manfull[which(manfull$dS !=0 & manfull$pS !=0 ),]
spmind<-seq(0,1,0.05)
slide_alphas<-rep(0,length(spmind))
for (i in 1:length(spmind))
{
passers<-tstset[which(tstset$mtAdult>spmind[i] |
tstset$gutAdult>spmind[i] | tstset$headMale | tstset$testesAdult>spmind[i]),]
slide_alphas[i]<-1-NItgCalc(passers$dN,passers$dS,passers$pN,passers$pS)
}

plot(spmind,slide_alphas,ylab="",xlab="SPM Threshold",cex.lab=2,
main="Effect of increasing tissue specificity on Carolina sphinx moth results", 
ylim=c(-0.13, 0.062),las=1, xlim=c(0,0.95),cex=2, cex.axis=1.5, cex.main=1.5)
mtext(expression(alpha),side=2,las=1,line=3,cex=2)
abline(h=0)

tstset<-manfull[which(manfull$dS !=0 & manfull$pS !=0 ),]
tstset<-tstset[which(tstset$Class=="eupyrene" | tstset$Class=="apyrene" | tstset$Class=="shared"),]
spmind<-seq(0,1,0.05)
slide_alphas<-rep(0,length(spmind))
for (i in 1:length(spmind))
{
passers<-tstset[which(tstset$testesAdult>spmind[i]),]
slide_alphas[i]<-1-NItgCalc(passers$dN,passers$dS,passers$pN,passers$pS)
}
points(spmind,slide_alphas,pch=16,col="blue", cex=2)

tstset<-manfull[which(manfull$dS !=0 & manfull$pS !=0 ),]
tstset<-tstset[which(tstset$Class=="non-sperm"),]
spmind<-seq(0,1,0.05)
slide_alphas<-rep(0,length(spmind))
for (i in 1:length(spmind))
{
passers<-tstset[which(tstset$testesAdult>spmind[i]),]
slide_alphas[i]<-1-NItgCalc(passers$dN,passers$dS,passers$pN,passers$pS)
}
points(spmind,slide_alphas,pch=18,col="red", cex=2)

legend(-0.01,0.065,pch=c(1,16,18),col=c("black","blue","red"),c("all genes","sperm genes","male-limited (non-sperm)"),pt.cex=2, cex=1.5)
mtext(side=3, adj=-0.05, "A", cex = 2, line=1)



###Response to reviewer comments


#histogram of sperm proteome SPM


dstset<-danfull[which(danfull$dS !=0 & danfull$pS !=0 ),]
dsperm<-dstset[which(dstset$gene.class != "non-sperm"),]

mtstset<-manfull[which(manfull$dS !=0 & manfull$pS !=0 ),]
msperm<-mtstset[which(mtstset$Class !="non-sperm"),]

#Figure S3 
par(mfrow=c(2,2))
#note that monarch SPMs are incorrectly labeled fpkms because i didn't rename my columns, but the values are still correct
hist(dstset$fpkm.testes.avg,xlab="SPM",main="Specificity of all genes expressed in testes of monarchs")
hist(mtstset$testesAdult,xlab="SPM",main="Specificity of all genes expressed in testes of sphinx moth")

hist(dsperm$fpkm.testes.avg,xlab="SPM",main="Specificity of sperm proteins in testes of monarch butterfly")
hist(msperm$testesAdult,xlab="SPM",main="Specificity of sperm proteins in testes of sphinx moth")
#this also helps us explain why SPM has a hump in alpha in the middle?
#the hump is most in Manduca, and occurs in the range (0.4,0.6) for sperm genes, but monarchs also see a sharp increase around the (0.1,0.2) range
#these regions line up with switch points between distributions of broadly expressed and specifically expressed genes, so it shouldn't suprise that there's variablity in there



#Signal peptides?






#Alpha for head/thorax/abdomen/ vs background?

par(mfrow=c(2,3))
hist(dstset$head.male,xlab="SPM",main="head monarchs",ylim=c(0,100))
hist(dstset$thorax.male,xlab="SPM",main="thorax monarchs",ylim=c(0,100))
hist(dstset$mg.male,xlab="SPM",main="gut monarchs",ylim=c(0,100))

hist(mtstset$headMale,xlab="SPM",main="head sphinx moth",ylim=c(0,600))
hist(mtstset$mtAdult,xlab="SPM",main="thorax sphinx moth",ylim=c(0,100))
hist(mtstset$gutAdult,xlab="SPM",main="gut sphinx moth",ylim=c(0,100))

#okay, per the bimodal distribution, let's bin genes with SPM < 0.5 and those with SPM > 0.5 
#remove genes that explode alpha calc
tstset<-danfull[which(danfull$dS !=0 & danfull$pS !=0 ),]
tstset<-tstset[which(tstset$ZorA=="A"),]
spnot<-tstset[which(tstset$gene.class =="non-sperm"),]
spsp<-tstset[which(tstset$gene.class !="non-sperm"),]
spmind<-0.5
#monarch non-sperm
#parse the low bin from the hi
passerslo<-spnot[which(spnot$thorax.male<spmind |
spnot$mg.male<spmind | spnot$head.male<spmind ),]
passershi<-spnot[which(spnot$thorax.male>spmind |
spnot$mg.male>spmind | spnot$head.male>spmind ),]
alo<-1-NItgCalc(passerslo$dN,passerslo$dS,passerslo$pN,passerslo$pS)
ahi<-1-NItgCalc(passershi$dN,passershi$dS,passershi$pN,passershi$pS)


size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
for(i in 1:size)
	{
		bsaph<-passershi[sample(nrow(passershi),size=nrow(passershi),replace=T),]
		bsapl<-passerslo[sample(nrow(passerslo),size=nrow(passerslo),replace=T),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		
	}
hist(bshi)

#hist(bsstata)
hdist<-sort(bshi)
acish<-c(hdist[250],ahi,hdist[9750])
acish
#for monarch non-sperm hi specificity
#2.5%         pt. est.     97.5%
#-0.25352949 -0.15942813 -0.07326127

ldist<-sort(bslo)
acisl<-c(ldist[250],alo,ldist[9750])
acisl
#for monarch non-sperm lo specificity
#2.5%         pt. est.     97.5%
#-0.1567134 -0.1359472 -0.1160281

#significance test
testor<-ahi-alo
size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
bsdif<-rep(0,size)
bswhole<-rbind(passershi,passerslo)
for(i in 1:size)
	{
		bsap<-bswhole[sample(nrow(bswhole),replace=F),]
		bsaph<-bsap[1:nrow(passershi),]
		bsapl<-bsap[(nrow(passershi)+1):nrow(bsap),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		bsdif[i]<-bshi[i]-bslo[i]
		
	}
hist(bsdif)
abline(v=testor)
bsdist<-sort(bsdif)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
# 0.6831

#sperm time
#parse the low bin from the hi
passerslo<-spsp[which(spsp$thorax.male<spmind |
spsp$mg.male<spmind | spsp$head.male<spmind | spsp$fpkm.testes.avg<spmind | spsp$fpkm.ag.avg<spmind),]
passershi<-spsp[which(spsp$thorax.male>spmind |
spsp$mg.male>spmind | spsp$head.male>spmind |spsp$fpkm.testes.avg>spmind |spsp$fpkm.ag.avg>spmind),]
alo<-1-NItgCalc(passerslo$dN,passerslo$dS,passerslo$pN,passerslo$pS)
ahi<-1-NItgCalc(passershi$dN,passershi$dS,passershi$pN,passershi$pS)


size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
for(i in 1:size)
	{
		bsaph<-passershi[sample(nrow(passershi),size=nrow(passershi),replace=T),]
		bsapl<-passerslo[sample(nrow(passerslo),size=nrow(passerslo),replace=T),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		
	}
hist(bshi)

#hist(bsstata)
hdist<-sort(bshi)
acish<-c(hdist[250],ahi,hdist[9750])
acish
#for monarch sperm hi specificity
#2.5%         pt. est.     97.5%
#0.1278896    0.2099777    0.2822741

ldist<-sort(bslo)
acisl<-c(ldist[250],alo,ldist[9750])
acisl
#for monarch sperm lo specificity
#2.5%         pt. est.     97.5%
#-0.04628542  0.06574155   0.16336396

#significance test
testor<-ahi-alo
size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
bsdif<-rep(0,size)
bswhole<-rbind(passershi,passerslo)
for(i in 1:size)
	{
		bsap<-bswhole[sample(nrow(bswhole),replace=F),]
		bsaph<-bsap[1:nrow(passershi),]
		bsapl<-bsap[(nrow(passershi)+1):nrow(bsap),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		bsdif[i]<-bshi[i]-bslo[i]
		
	}
hist(bsdif)
abline(v=testor)
bsdist<-sort(bsdif)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
# 0.0242


##let's try for our male-limited proxy##############
#aka non-sperm from the testes perspective##

passerslo<-spnot[which(spnot$thorax.male<spmind |
spnot$mg.male<spmind | spnot$head.male<spmind | spnot$fpkm.testes.avg<spmind | spnot$fpkm.ag.avg<spmind),]
passershi<-spnot[which(spnot$thorax.male>spmind |
spnot$mg.male>spmind | spnot$head.male>spmind |spnot$fpkm.testes.avg>spmind |spnot$fpkm.ag.avg>spmind),]
alo<-1-NItgCalc(passerslo$dN,passerslo$dS,passerslo$pN,passerslo$pS)
ahi<-1-NItgCalc(passershi$dN,passershi$dS,passershi$pN,passershi$pS)


size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
for(i in 1:size)
	{
		bsaph<-passershi[sample(nrow(passershi),size=nrow(passershi),replace=T),]
		bsapl<-passerslo[sample(nrow(passerslo),size=nrow(passerslo),replace=T),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		
	}
hist(bshi)

#hist(bsstata)
hdist<-sort(bshi)
acish<-c(hdist[250],ahi,hdist[9750])
acish
#for monarch testes hi specificity
#2.5%         pt. est.     97.5%
#-0.12138052  -0.07480633  -0.03030386

ldist<-sort(bslo)
acisl<-c(ldist[250],alo,ldist[9750])
acisl
#for monarch testes lo specificity
#2.5%         pt. est.     97.5%
#-0.1569366   -0.1359472   -0.1157395

##this is close....we can significance test it


testor<-ahi-alo
size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
bsdif<-rep(0,size)
bswhole<-rbind(passershi,passerslo)
for(i in 1:size)
	{
		bsap<-bswhole[sample(nrow(bswhole),replace=F),]
		bsaph<-bsap[1:nrow(passershi),]
		bsapl<-bsap[(nrow(passershi)+1):nrow(bsap),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		bsdif[i]<-bshi[i]-bslo[i]
		
	}
hist(bsdif)
abline(v=testor)
bsdist<-sort(bsdif)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.0137

###Manduca#######
#remove genes that explode alpha calc
tstset<-manfull[which(manfull$dS !=0 & manfull$pS !=0 ),]
tstset<-tstset[which(tstset$ZorA=="A"),]
mtstset<-tstset

#separate by class
spnot<-tstset[which(tstset$Class =="non-sperm"),]
spsp<-tstset[which(tstset$Class !="non-sperm"),]
#set our cut-off
spmind<-0.5
#manduca non-sperm
#parse the low bin from the hi
passerslo<-spnot[which(spnot$mtAdult<spmind |
spnot$gutAdult<spmind | spnot$headMale<spmind),]
passershi<-spnot[which(spnot$mtAdult>spmind |
spnot$gutAdult>spmind | spnot$headMale>spmind),]
alo<-1-NItgCalc(passerslo$dN,passerslo$dS,passerslo$pN,passerslo$pS)
ahi<-1-NItgCalc(passershi$dN,passershi$dS,passershi$pN,passershi$pS)


size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
for(i in 1:size)
	{
		bsaph<-passershi[sample(nrow(passershi),size=nrow(passershi),replace=T),]
		bsapl<-passerslo[sample(nrow(passerslo),size=nrow(passerslo),replace=T),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		
	}
hist(bshi)

#hist(bsstata)
hdist<-sort(bshi)
acish<-c(hdist[250],ahi,hdist[9750])
acish
#for manduca non-sperm hi specificity
#2.5%         pt. est.     97.5%
#-0.10656197 -0.07180327 -0.03868588

ldist<-sort(bslo)
acisl<-c(ldist[250],alo,ldist[9750])
acisl
#for manduca non-sperm lo specificity
#2.5%         pt. est.     97.5%
#-0.09546457 -0.07697234 -0.05869831

testor<-ahi-alo
size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
bsdif<-rep(0,size)
bswhole<-rbind(passershi,passerslo)
for(i in 1:size)
	{
		bsap<-bswhole[sample(nrow(bswhole),replace=F),]
		bsaph<-bsap[1:nrow(passershi),]
		bsapl<-bsap[(nrow(passershi)+1):nrow(bsap),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		bsdif[i]<-bshi[i]-bslo[i]
		
	}
hist(bsdif)
abline(v=testor)
bsdist<-sort(bsdif)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.3868

#now for sperm
passerslo<-spsp[which(spsp$mtAdult<spmind |
spsp$gutAdult<spmind | spsp$headMale<spmind | spsp$testesAdult<spmind),]
passershi<-spsp[which(spsp$mtAdult>spmind |
spsp$gutAdult>spmind | spsp$headMale>spmind |spsp$testesAdult>spmind),]
alo<-1-NItgCalc(passerslo$dN,passerslo$dS,passerslo$pN,passerslo$pS)
ahi<-1-NItgCalc(passershi$dN,passershi$dS,passershi$pN,passershi$pS)

size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
for(i in 1:size)
	{
		bsaph<-passershi[sample(nrow(passershi),size=nrow(passershi),replace=T),]
		bsapl<-passerslo[sample(nrow(passerslo),size=nrow(passerslo),replace=T),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		
	}
hist(bshi)

#hist(bsstata)
hdist<-sort(bshi)
acish<-c(hdist[250],ahi,hdist[9750])
acish
#for manduca sperm hi specificity
#2.5%         pt. est.     97.5%
#-0.11509348  -0.02440375  0.05130695

ldist<-sort(bslo)
acisl<-c(ldist[250],alo,ldist[9750])
acisl
#for manduca sperm lo specificity
#2.5%         pt. est.     97.5%
#-0.13026454  -0.05070336  0.02033531

testor<-ahi-alo
size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
bsdif<-rep(0,size)
bswhole<-rbind(passershi,passerslo)
for(i in 1:size)
	{
		bsap<-bswhole[sample(nrow(bswhole),replace=F),]
		bsaph<-bsap[1:nrow(passershi),]
		bsapl<-bsap[(nrow(passershi)+1):nrow(bsap),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		bsdif[i]<-bshi[i]-bslo[i]
		
	}
hist(bsdif)
abline(v=testor)
bsdist<-sort(bsdif)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.3248

####Manduca non-sperm but testes
#now for sperm
passerslo<-spnot[which(spnot$mtAdult<spmind |
spnot$gutAdult<spmind | spnot$headMale<spmind | spnot$testesAdult<spmind),]
passershi<-spnot[which(spnot$mtAdult>spmind |
spnot$gutAdult>spmind | spnot$headMale>spmind |spnot$testesAdult>spmind),]
alo<-1-NItgCalc(passerslo$dN,passerslo$dS,passerslo$pN,passerslo$pS)
ahi<-1-NItgCalc(passershi$dN,passershi$dS,passershi$pN,passershi$pS)

size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
for(i in 1:size)
	{
		bsaph<-passershi[sample(nrow(passershi),size=nrow(passershi),replace=T),]
		bsapl<-passerslo[sample(nrow(passerslo),size=nrow(passerslo),replace=T),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		
	}
hist(bshi)

#hist(bsstata)
hdist<-sort(bshi)
acish<-c(hdist[250],ahi,hdist[9750])
acish
#for manduca sperm hi specificity
#2.5%          pt. est.      97.5%
#-0.10638291   -0.07957497   -0.05300785

ldist<-sort(bslo)
acisl<-c(ldist[250],alo,ldist[9750])
acisl
#for manduca sperm lo specificity
#2.5%         pt. est.     97.5%
#-0.09541171  -0.07697234  -0.05881828

testor<-ahi-alo
size<-10000
bslo<-rep(0,size)
bshi<-rep(0,size)
bsdif<-rep(0,size)
bswhole<-rbind(passershi,passerslo)
for(i in 1:size)
	{
		bsap<-bswhole[sample(nrow(bswhole),replace=F),]
		bsaph<-bsap[1:nrow(passershi),]
		bsapl<-bsap[(nrow(passershi)+1):nrow(bsap),]
		bshi[i]<-1-NItgCalc(bsaph$dN,bsaph$dS,bsaph$pN,bsaph$pS)
		bslo[i]<-1-NItgCalc(bsapl$dN,bsapl$dS,bsapl$pN,bsapl$pS)
		bsdif[i]<-bshi[i]-bslo[i]
		
	}
hist(bsdif)
abline(v=testor)
bsdist<-sort(bsdif)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val


#we can make a dummy set of max SPMs
dspmdummy<-cbind(dstset[,22],dstset[,24],dstset[,26],dstset[,29:30])
dspmdummy$max<-apply(dspmdummy, 1, max)
mspmdummy<-cbind(mtstset[,18:21])
mspmdummy$max<-apply(mspmdummy, 1, max)

#########proposed new  figure ################################
#Figure S3 
par(mfrow=c(3,2))
# bottom, left, top, right margins
par(mai=c(0.26,0.64,0.46,0.06))
#note that monarch SPMs are incorrectly labeled fpkms because i didn't rename my columns, but the values are still correct
hist(dspmdummy$max,xlab="SPM",main="Specificity of all genes in current dataset",las=1)
mtext("Monarch butterfies", side = 3, adj = 0.5, line = 2)
hist(mspmdummy$max,xlab="SPM",main="Specificity of all genes expressed in current dataset",las=1)
mtext("Sphinx moths", side = 3, adj = 0.5, line = 2)
par(mai=c(0.26,0.64,0.26,0.06))
hist(dsperm$fpkm.testes.avg,xlab="SPM",main="Specificity of sperm proteins in testes",las=1)
abline(v=0.5,lty=2)
hist(msperm$testesAdult,xlab="SPM",main="Specificity of sperm proteins in testes")
abline(v=0.5,lty=2)
#this also helps us explain why SPM has a hump in alpha in the middle?
#the hump is most in Manduca, and occurs in the range (0.4,0.6) for sperm genes, but monarchs also see a sharp increase around the (0.1,0.2) range
#these regions line up with switch points between distributions of broadly expressed and specifically expressed genes, so it shouldn't suprise that there's variablity in there


# bottom, left, top, right margins
par(mai=c(0.3,0.64,0.06,0.06))
plot(-9,-9,xlim=c(0,1),ylim=c(-0.28,0.4),ylab="",xlab="",xaxt='n',las=1)
mtext(expression(alpha),side=2,las=1,line=3,cex=2)
points(0.22,-0.1359472 , pch=0, cex=1.3, col="orange")
points(0.22,-0.1359472 , pch=0, cex=1.1, col="orange")
points(0.22,-0.1359472 , pch=0, cex=0.9, col="orange")
segments(0.22,-0.1569141,0.22,-0.1155204,lwd=2,col="orange")
segments(0.21,-0.1569141,0.23,-0.1569141,lwd=2,col="orange")
segments(0.21,-0.1155204,0.23,-0.1155204,lwd=2,col="orange")
#sperm proteins lo
points(0.27,0.06574155, pch=16, cex=1.3, col="orange")
segments(0.27,-0.04628542,0.27,0.16336396,lwd=2,col="orange")
segments(0.26,-0.04628542,0.28,-0.04628542,lwd=2,col="orange")
segments(0.26,0.16336396,0.28,0.16336396,lwd=2,col="orange")
#testes lo
points(0.32, -0.1359472 , pch=17, cex=1.3, col="orange")
segments(0.32,-0.1569366,0.32,-0.1157395,lwd=2,col="orange")
segments(0.31,-0.1569366,0.33,-0.1569366,lwd=2,col="orange")
segments(0.31,-0.1157395,0.33,-0.1157395,lwd=2,col="orange")

abline(v=0.5,lty=2)

#non-sperm hi
points(0.72,-0.15942813, pch=0, cex=1.3, col="orange")
points(0.72,-0.15942813, pch=0, cex=1.1, col="orange")
points(0.72,-0.15942813, pch=0, cex=0.9, col="orange")
segments(0.72,-0.25352949,0.72, -0.07326127,lwd=2,col="orange")
segments(0.71,-0.25352949,0.73,-0.25352949,lwd=2,col="orange")
segments(0.71,  -0.07326127,0.73,-0.07326127,lwd=2,col="orange")
#sperm proteins hi
points(0.77,0.2099777 , pch=16, cex=1.3, col="orange")
segments(0.77, 0.1278896, 0.77,  0.2822741,lwd=2,col="orange")
segments(0.76, 0.1278896, 0.78, 0.1278896,lwd=2,col="orange")
segments(0.76,  0.2822741,0.78, 0.2822741,lwd=2,col="orange")
#testes
points(0.82,-0.07480633  , pch=17, cex=1.3, col="orange")
segments(0.82, -0.12138052, 0.82,  -0.03030386,lwd=2,col="orange")
segments(0.81, -0.12138052, 0.83, -0.12138052,lwd=2,col="orange")
segments(0.81,  -0.03030386,0.83, -0.03030386,lwd=2,col="orange")
mtext("Low-specificity", side = 1, adj = 0.1, line = 1)
mtext("High-specificity", side = 1, adj = 0.9, line = 1)
legend(0,0.4,pch=c(0,16,17),col="orange",c("background","sperm","male-limited"))
segments(0.32, 0.01, 0.82, 0.01)
segments(0.82, 0.01, 0.82, -0.01)
segments(0.32, 0.01, 0.32, -0.1)
text(0.62,0.04,"p = 0.0137")
segments(0.27, 0.31, 0.77, 0.31)
segments(0.77, 0.31, 0.77, 0.295)
segments(0.27, 0.31, 0.27, 0.21)
text(0.62,0.34,"p = 0.0242")


# bottom, left, top, right margins
par(mai=c(0.3,0.64,0.06,0.06))
plot(-9,-9,xlim=c(0,1),ylim=c(-0.15,0.1),ylab="",xlab="",xaxt='n',las=1)
mtext(expression(alpha),side=2,las=1,line=3,cex=2)
points(0.22,-0.07697234, pch=0, cex=1.3, col="dimgray")
points(0.22,-0.07697234, pch=0, cex=1.1, col="dimgray")
points(0.22,-0.07697234, pch=0, cex=0.9, col="dimgray")
segments(0.22,-0.09529088,0.22,-0.05868484,lwd=2,col="dimgray")
segments(0.21,-0.09529088,0.23,-0.09529088,lwd=2,col="dimgray")
segments(0.21,-0.05868484,0.23,-0.05868484,lwd=2,col="dimgray")
#sperm proteins lo
points(0.27,-0.05070336, pch=16, cex=1.3, col="dimgray")
segments(0.27,-0.13026454,0.27,0.02033531,lwd=2,col="dimgray")
segments(0.26,-0.13026454,0.28,-0.13026454,lwd=2,col="dimgray")
segments(0.26, 0.02033531,0.28,0.02033531,lwd=2,col="dimgray")
#non-sperm testes
points(0.32,-0.07697234, pch=17, cex=1.3, col="dimgray")
segments(0.32,-0.09541171,0.32,-0.05881828,lwd=2,col="dimgray")
segments(0.31,-0.09541171,0.33,-0.09541171,lwd=2,col="dimgray")
segments(0.31, -0.05881828,0.33,-0.05881828,lwd=2,col="dimgray")


abline(v=0.5,lty=2)

#non-sperm hi
points(0.72,-0.07180327, pch=0, cex=1.3, col="dimgray")
points(0.72,-0.07180327, pch=0, cex=1.1, col="dimgray")
points(0.72,-0.07180327, pch=0, cex=0.9, col="dimgray")
segments(0.72,-0.10656197,0.72,-0.03868588,lwd=2,col="dimgray")
segments(0.71,-0.10656197,0.73,-0.10656197,lwd=2,col="dimgray")
segments(0.71, -0.03868588,0.73,-0.03868588,lwd=2,col="dimgray")
#sperm proteins hi
points(0.77,-0.02440375, pch=16, cex=1.3, col="dimgray")
segments(0.77, -0.11509348, 0.77,  0.05130695,lwd=2,col="dimgray")
segments(0.76, -0.11509348, 0.78, -0.11509348,lwd=2,col="dimgray")
segments(0.76,  0.05130695,0.78, 0.05130695,lwd=2,col="dimgray")
#non-sperm testes
points(0.82,-0.07957497, pch=17, cex=1.3, col="dimgray")
segments(0.82, -0.10638291, 0.82,  -0.05300785,lwd=2,col="dimgray")
segments(0.81, -0.10638291, 0.83, -0.10638291,lwd=2,col="dimgray")
segments(0.81,  -0.05300785,0.83, -0.05300785,lwd=2,col="dimgray")
mtext("Low-specificity", side = 1, adj = 0.1, line = 1)
mtext("High-specificity", side = 1, adj = 0.9, line = 1)
legend(0,0.05,pch=c(0,16,17),col="dimgray",c("background","sperm", "male-limited"))
