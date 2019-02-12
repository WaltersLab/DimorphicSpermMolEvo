
###this is Manduca
ancuf<-c(2364810.285072, 49911.892188, 23946.818661, 15429.569841, 11245.885870, 8916.528336, 7367.203951, 6535.537473, 5893.685237, 5176.135592, 4750.860467, 4494.885802, 5333.897160, 3709.265661, 3471.305452, 3288.573917, 3445.394245, 3271.845586, 3470.677227, 3557.369325, 3786.061613, 4153.048505, 5375.609119, 8691.121864, 6106.541834)
 
angsdupdateuf<-c(106713.424293, 2036.725628, 987.382202, 595.918337, 400.249217, 356.019683, 303.950451, 250.439892, 248.645846, 202.742447, 201.998323, 181.655934, 274.948221, 174.464247, 145.351390, 170.611281, 165.563983, 153.681046, 179.642826, 163.829448, 190.452388, 236.934044, 329.030935, 600.345759, 11747.992177)
 
 par(mfrow=c(2,1))
 barplot(ancuf,ylim=c(0,59911))
 barplot(angsdupdateuf, ylim=c(0,2500))
 
 
anc0x<-c(11024505.834231, 47313.978325, 19925.634164, 11591.726433, 7804.615239, 5654.296256, 4536.397722, 3686.457362, 3253.512785, 2852.943596, 2571.934696, 2439.119335, 3543.368427, 1662.943771, 1485.031909, 1369.058073, 1325.487508, 1331.971052, 1354.332952, 1290.529149, 1379.923371, 1600.646383, 1910.507283, 3275.299825, 2918.450152) 

 par(mfrow=c(2,1))
 barplot(ancuf,ylim=c(0,18800),main="4x SFS")
 barplot(anc0x, ylim=c(0,18800),main="0x SFS")
 
 msws4<-c(74967.556671,1414.576730,743.726048,458.534529,322.710610,243.460475,249.107901,178.452967,189.141892,157.117390,182.655963,110.875740,195.817078,137.872758,103.785837,108.669637,144.915788,89.310980,121.806868,110.606810,95.446545,160.635618,161.772834,251.942257,154.500077)
msws0<-c(358600.499608,1401.249383,588.837838,366.927651,188.763856,140.018237,121.773763,96.485492,98.007318,69.708963,65.009793,71.899383,93.009940,153.943236,52.124584,37.923004,39.425422,55.470704,52.077660,36.267697,49.317008,85.280098,80.740385,88.787692,62.451288) 
 
 
 msor4<-c(54278.541851,1161.304685,618.156786,388.109574,378.549938,261.064111,220.991836,189.649967,176.041798,108.124203,107.110063,117.885731,126.156372,98.264338,96.417167,85.516140,120.020847,93.798766,126.301880,113.325213,155.375091,126.598149,143.638695,227.526802,171.529994)
 msor0<-c(267913.665852,825.368412,336.205784,202.117454,171.934703,98.414104,60.467343,68.733167,51.791134,34.372997,42.601834,30.642066,37.754274,25.631152,24.321999,10.284862,45.380622,26.446030,26.219529,19.368183,42.645118,34.837718,28.938425,59.818248,28.038989)
 
 
 ###this is monarchs
 non4auto<-c(2291724.807804,122782.506115,42881.776205,20862.556514,13629.857569,9912.499805,7942.488050,6562.470018,5448.170195,4857.374588,4500.819163,4273.589985,4442.533759,3231.523896,2961.740312,2954.356253,2822.115448,2660.536977,2831.065478,2890.261595,3095.761358,3594.956500,4990.830043,8796.560808,4136.841562)
 non0auto<-c(11161470.126366,98803.655118,26246.705715,11092.728792,7274.423362,5093.870138,4056.982927,3291.253772,3104.740602,2663.690263,2481.075028,2414.553587,3416.612710,1475.519022,1330.188635,1301.264125,1170.595637,1123.101881,1195.338379,1147.301089,1281.688260,1386.461293,1776.288605,3145.445529,2952.389163)
 eu4auto<-c(56951.274437,2880.252168,961.285373,525.634410,321.605738,249.351843,164.409406,143.356734,112.477320,105.186811,122.305607,89.582024,117.775204,75.326046,59.216540,59.392357,71.488200,67.847153,54.627099,67.209530,66.321716,88.371059,97.405364,197.026721,75.271140)
 eu0auto<-c(289782.460325,1820.025595,495.241705,166.755978,122.823620,122.889806,68.184197,63.517308,52.733162,43.887176,41.994751,36.438407,48.976581,30.225990,23.622510,26.607080,15.633988,23.713915,27.495679,30.506042,23.702628,25.987868,36.571814,69.807588,45.196288)
 
 ap4auto<-c(30734.788425,1585.203287,570.517987,278.662559,164.378512,137.849859,90.512959,96.740853,65.296427,70.699096,59.772687,43.825447,55.845290,43.702094,49.563308,34.251322,35.204424,35.138151,31.171577,43.399891,59.726729,51.054129,77.270997,117.467206,42.956782)
 
 ap0auto<-c(155189.521284,957.087661,238.034476,89.393631,66.304412,43.521491,36.990691,27.341399,30.995023,22.940171,13.308236,25.748697,17.633388,12.698053,19.327227,22.360706,10.651689,10.823704,16.392034,11.437906,14.590831,14.332996,22.564839,38.599877,23.399579)
 
par(mfrow=c(2,1))
barplot(non4auto,ylim=c(0,48800),main="4x SFS")
barplot(non0auto, ylim=c(0,48800),main="0x SFS")
 
par(mfrow=c(2,1))
barplot(eu4auto,ylim=c(0,580),main="4x SFS")
barplot(eu0auto, ylim=c(0,580),main="0x SFS")


par(mfrow=c(2,1))
barplot(non4auto-non0auto,ylim=c(0,1800),main="BG diff")
barplot(eu4auto-eu0auto,ylim=c(0,350),main="eu diff")

#monarch whole sperm
mws0<-c(765247.147894,4234.416765,1039.390045,374.253612,244.325669,199.476121,136.271205,128.655907,113.858912,87.765367,88.257793,88.542176,87.889342,58.410380,68.058492,53.378964,50.823623,46.519601,54.544558,62.066428,65.286813,77.629993,96.700048,180.413582,91.916710)
mws4<-c(148380.508891,7230.657692,2518.151156,1238.635366,845.343202,617.992628,440.628017,357.334831,279.035814,270.188700,280.872484,227.162587,241.706245,198.538093,177.091073,148.628430,173.087569,179.859928,158.950857,173.257595,216.001001,236.892039,301.644556,564.875516,210.955728)

#the ortho sfses
dpor0<-c(265889.777291,1023.381089,227.793753,81.778009,55.092408,45.591866,23.034220,36.267946,25.539556,14.994370,23.006628,19.272103,22.639831,19.230783,14.658388,6.963023,17.448562,7.136567,10.133673,13.287020,14.816144,18.709867,19.728662,42.599557,23.118682)
dpor4<-c(51458.782736,2546.414005,910.157682,422.817018,312.422129,230.933287,165.412621,123.924214,84.288690,99.608123,105.103622,77.268662,81.447178,76.966031,54.791798,46.704007,73.040245,63.135975,62.497499,59.845653,70.274005,84.478762,91.015888,205.768273,76.901897)

par(mfrow=c(2,1))
barplot(dpor4,ylim=c(0,1500),main="4x SFS")
barplot(dpor0, ylim=c(0,200),main="0x SFS")



par(mfrow=c(2,1))
barplot(euAFRdfe)

barplot(eudfe)

par(mfrow=c(2,1))
barplot(bgCBFGdfe,ylim=c(0,0.85),main="Distribution of fitness effects of new mutations in genome BG",ylab="Freq",xlab="s")
barplot(euCBFGdfe,ylim=c(0,0.85),main="Distribution of fitness effects of new mutations in eupyrene sperm",ylab="Freq",xlab="s")

###Actual manuscript code starts here

#how does bootstrapping work?
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/bootstraps/Dplex_eu_data_file", rep = 100)
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/bootstraps/Msexta_bg_file", rep = 100)
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/bootstraps/Dplex_ortho_file", rep = 130)
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/bootstraps/Msexta_wsperm_file", rep = 130)
bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/bootstraps/Msexta_ortho_file", rep = 130)


bootstrapData(inputfile="/Users/Andrew/Documents/polyDFE/bootstraps/Dplex_whole_sperm", rep = 130)

#can we iterate loading files into a list for bootstrapping
#we can, but only when there are no failed runs, "max iterations reached" or whatever
setwd("/Users/Andrew/Documents/polyDFE/bsData")
#eudataFiles <- lapply(Sys.glob("Dplex_eu_out_strap*"), parseOutput)
msbg <- lapply(Sys.glob("Msexta_bg_out_strap*"), parseOutput)
dpbg<- lapply(Sys.glob("Dplex_bg_out_strap*"), parseOutput)
dpor<- lapply(Sys.glob("dplex_ortho_out_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Ms_ws/")
msws<- lapply(Sys.glob("Msexta_wsperm_out_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Dp_ws/")
dpws<- lapply(Sys.glob("Dplex_wsperm_out_strap*"), parseOutput)
setwd("/Users/Andrew/Documents/polyDFE/bsData/Ms_or/")
msor<- lapply(Sys.glob("Msexta_ortho_out_strap*"), parseOutput)

#edf<-eudataFiles



#modify_depth can be used to subset nested lists, but requires the stupidly named "purrr" package
library("purrr")



#alphs<-unlist(modify_depth(edf,2,"alpha"))

dpbga<-unlist(modify_depth(dpbg,2,"alpha"))
msbga<-unlist(modify_depth(msbg,2,"alpha"))
dpora<-unlist(modify_depth(dpor,2,"alpha"))
mswsa<-unlist(modify_depth(msws,2,"alpha"))
dpwsa<-unlist(modify_depth(dpws,2,"alpha"))
msora<-unlist(modify_depth(msor,2,"alpha"))

#let's bracket these bois
#sem<-sd(alphs)/sqrt(length(alphs))
#cis<-c(mean(alphs)-2*sem,mean(alphs)+2*sem)


dpbgasem<-sd(dpbga)/sqrt(length(dpbga))
dpbgacis<-c(mean(dpbga)-2*dpbgasem,mean(dpbga)+2*dpbgasem)

msbgasem<-sd(msbga)/sqrt(length(msbga))
msbgacis<-c(mean(msbga)-2*msbgasem,mean(msbga)+2*msbgasem)

mswsasem<-sd(mswsa)/sqrt(length(mswsa))
mswsacis<-c(mean(mswsa)-2*mswsasem,mean(mswsa)+2*mswsasem)

dpwsasem<-sd(dpwsa)/sqrt(length(dpwsa))
dpwsacis<-c(mean(dpwsa)-2*dpwsasem,mean(dpwsa)+2*dpwsasem)

dporasem<-sd(dpora)/sqrt(length(dpora))
dporacis<-c(mean(dpora)-2*dporasem,mean(dpora)+2*dporasem)

msorasem<-sd(msora)/sqrt(length(msora))
msoracis<-c(mean(msora)-2*msorasem,mean(msora)+2*msorasem)

#let's do an alpha graphic
dumxm<-c(0.1,0.2,0.3)
dumxd<-c(0.5,0.6,0.7)
plot(-9,-9,xlim=c(0,1),ylim=c(0,1))
points(dumxm,c(mean(msbga),mean(mswsa),mean(msora)))
segments(dumxm,c(msbgacis[1],mswsacis[1],msoracis[1]),dumxm,c(msbgacis[2],mswsacis[2],msoracis[2]))
points(dumxd,c(mean(dpbga),mean(dpwsa),mean(dpora)))
segments(dumxd,c(dpbgacis[1],dpwsacis[1],dporacis[1]),dumxd,c(dpbgacis[2],dpwsacis[2],dporacis[2]))


#can we get the bootstrapped DFE?
dpbgtops<-length(dpbg)
dpbgdfemat<-data.frame(matrix(NA, nrow = dpbgtops, ncol = 7))
for(i in 1:dpbgtops)
{
dpbgdfemat[i,]<-unlist(lapply(dpbg[[i]],getDiscretizedDFE))
}

dpbgbsest<-colMeans(dpbgdfemat)

names(dpbgbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")
#names(bsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpbgsdsdfe<-apply(dpbgdfemat, 2, sd)
dpbgsemdfe<-dpbgsdsdfe/sqrt(dpbgtops)

dumx<-seq(0.5,7.5,1)
dpbgciHi<-dpbgbsest+2*dpbgsemdfe
dpbgciLo<-dpbgbsest-2*dpbgsemdfe




msbgtops<-length(msbg)
msbgdfemat<-data.frame(matrix(NA, nrow = msbgtops, ncol = 7))
for(i in 1:msbgtops)
{
msbgdfemat[i,]<-unlist(lapply(msbg[[i]],getDiscretizedDFE))
}

msbgbsest<-colMeans(msbgdfemat)

names(msbgbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

msbgsdsdfe<-apply(msbgdfemat, 2, sd)
msbgsemdfe<-msbgsdsdfe/sqrt(msbgtops)

dumx<-seq(0.5,7.5,1)
msbgciHi<-msbgbsest+2*msbgsemdfe
msbgciLo<-msbgbsest-2*msbgsemdfe

dpwstops<-length(dpws)
dpwsdfemat<-data.frame(matrix(NA, nrow = dpwstops, ncol = 7))
for(i in 1:dpwstops)
{
dpwsdfemat[i,]<-unlist(lapply(dpws[[i]],getDiscretizedDFE))
}

dpwsbsest<-colMeans(dpwsdfemat)

names(dpwsbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")

dpwssdsdfe<-apply(dpwsdfemat, 2, sd)
dpwssemdfe<-dpwssdsdfe/sqrt(dpwstops)

dumx<-seq(0.5,7.5,1)
dpwsciHi<-dpwsbsest+2*dpwssemdfe
dpwsciLo<-dpwsbsest-2*dpwssemdfe


dportops<-length(dpor)
dpordfemat<-data.frame(matrix(NA, nrow = dportops, ncol = 7))
for(i in 1:dportops)
{
dpordfemat[i,]<-unlist(lapply(dpor[[i]],getDiscretizedDFE))
}

dporbsest<-colMeans(dpordfemat)

names(dporbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")


dporsdsdfe<-apply(dpordfemat, 2, sd)
dporsemdfe<-dporsdsdfe/sqrt(dportops)

dumx<-seq(0.5,7.5,1)
dporciHi<-dporbsest+2*dporsemdfe
dporciLo<-dporbsest-2*dporsemdfe


###now the Manduca side
mswstops<-length(msws)
mswsdfemat<-data.frame(matrix(NA, nrow = mswstops, ncol = 7))
for(i in 1:mswstops)
{
mswsdfemat[i,]<-unlist(lapply(msws[[i]],getDiscretizedDFE))
}

mswsbsest<-colMeans(mswsdfemat)

names(mswsbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")


mswssdsdfe<-apply(mswsdfemat, 2, sd)
mswssemdfe<-mswssdsdfe/sqrt(mswstops)

dumx<-seq(0.5,7.5,1)
mswsciHi<-mswsbsest+2*mswssemdfe
mswsciLo<-mswsbsest-2*mswssemdfe

msortops<-length(msor)
msordfemat<-data.frame(matrix(NA, nrow = msortops, ncol = 7))
for(i in 1:msortops)
{
msordfemat[i,]<-unlist(lapply(msor[[i]],getDiscretizedDFE))
}

msorbsest<-colMeans(msordfemat)

names(msorbsest)<-c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<")


msorsdsdfe<-apply(msordfemat, 2, sd)
msorsemdfe<-msorsdsdfe/sqrt(msortops)

dumx<-seq(0.5,7.5,1)
msorciHi<-msorbsest+2*msorsemdfe
msorciLo<-msorbsest-2*msorsemdfe

ebx1<-seq(0.3,7.3,1)
ebx2<-seq(0.7,7.7,1)

#the DFE comparison plots

library("plotrix")
#we have to define our gap bounds
from<-0.175
to<-0.675

effcols<-c("brown4","brown3","brown1","gray66","darkseagreen2","darkseagreen3","darkseagreen4")

par(mfrow=c(3,2))
# bottom, left, top, right margins
par(mai=c(0.26,0.14,0.42,0.16))
gap.barplot(msbgbsest,xlim=c(0,7.2),main="Carolina sphinx moth
Genome backgrond",gap=c(from,to),col=effcols,ytics=seq(0,1,0.05),
xaxlab=names(msbgbsest),xlab="Effect class",ylim=c(0,0.44),ylab="")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2, at=from)
segments(dumx+0.5,msbgciHi,dumx+0.5,msbgciLo,lwd=2)
segments(ebx1+0.5,msbgciHi,ebx2+0.5,msbgciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,msbgciHi[1]-(to-from),ebx2[1]+0.5,msbgciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,msbgciHi[1]-(to-from),dumx[1]+0.5,msbgciLo[1]-(to-from),lwd=2)

gap.barplot(dpbgbsest,xlim=c(0,7.2),main="Monarch butterfly
Genome background",gap=c(from,to),col=effcols,ytics=seq(0,1,0.05),
xaxlab=names(msbgbsest),xlab="Effect class",ylim=c(0,0.44),ylab="")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2, at=from)
segments(dumx+0.5,dpbgciHi,dumx+0.5,dpbgciLo,lwd=2)
segments(ebx1+0.5,dpbgciHi,ebx2+0.5,dpbgciHi,lwd=2)
segments(ebx1[1]+0.5,dpbgciHi[1]-(to-from),ebx2[1]+0.5,dpbgciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpbgciHi[1]-(to-from),dumx[1]+0.5,dpbgciLo[1]-(to-from),lwd=2)

# bottom, left, top, right margins
par(mai=c(0.26,0.14,0.12,0.16))

gap.barplot(mswsbsest,xlim=c(0,7.2),main="Whole sperm proteome",gap=c(from,to),col=effcols,ytics=seq(0,1,0.05),
xaxlab=names(mswsbsest),xlab="Effect class",ylim=c(0,0.44),ylab="")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2, at=from)
segments(dumx+0.5,mswsciHi,dumx+0.5,mswsciLo,lwd=2)
segments(ebx1+0.5,mswsciHi,ebx2+0.5,mswsciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,mswsciHi[1]-(to-from),ebx2[1]+0.5,mswsciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,mswsciHi[1]-(to-from),dumx[1]+0.5,mswsciLo[1]-(to-from),lwd=2)

gap.barplot(dpwsbsest,xlim=c(0,7.2),main="Whole sperm proteome",gap=c(from,to),col=effcols,ytics=seq(0,1,0.05),
xaxlab=names(mswsbsest),xlab="Effect class",ylim=c(0,0.44),ylab="")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2, at=from)
segments(dumx+0.5,dpwsciHi,dumx+0.5,dpwsciLo,lwd=2)
segments(ebx1+0.5,dpwsciHi,ebx2+0.5,dpwsciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dpwsciHi[1]-(to-from),ebx2[1]+0.5,dpwsciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dpwsciHi[1]-(to-from),dumx[1]+0.5,dpwsciLo[1]-(to-from),lwd=2)

gap.barplot(msorbsest,xlim=c(0,7.2),main="Shared sperm proteins",gap=c(from,to),col=effcols,ytics=seq(0,1,0.05),
xaxlab=names(msorbsest),xlab="Effect class",ylim=c(0,0.44),ylab="")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2, at=from)
segments(dumx+0.5,msorciHi,dumx+0.5,msorciLo,lwd=2)
segments(ebx1+0.5,msorciHi,ebx2+0.5,msorciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,msorciHi[1]-(to-from),ebx2[1]+0.5,msorciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,msorciHi[1]-(to-from),dumx[1]+0.5,msorciLo[1]-(to-from),lwd=2)

gap.barplot(dporbsest,xlim=c(0,7.2),main="Shared sperm proteins",gap=c(from,to),col=effcols,ytics=seq(0,1,0.05),
xaxlab=names(dporbsest),xlab="Effect class",ylim=c(0,0.44),ylab="")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2, at=from)
segments(dumx+0.5,dporciHi,dumx+0.5,dporciLo,lwd=2)
segments(ebx1+0.5,dporciHi,ebx2+0.5,dporciHi,lwd=2)
#we have to fudge that first error bar due to the nature of the plot-break
segments(ebx1[1]+0.5,dporciHi[1]-(to-from),ebx2[1]+0.5,dporciHi[1]-(to-from),lwd=2)
segments(dumx[1]+0.5,dporciHi[1]-(to-from),dumx[1]+0.5,dporciLo[1]-(to-from),lwd=2)






par(mfrow=c(3,2))
barplot(msbgbsest,xlim=c(0,7.2),main="Msexta Genome BG",width=1,space=0,ylim=c(0,0.15),col=effcols)
segments(dumx,msbgciHi,dumx,msbgciLo,lwd=2)
segments(ebx1,msbgciHi,ebx2,msbgciHi,lwd=2)

barplot(dpbgbsest,xlim=c(0,7.2),main="Dplex Genome BG",width=1,space=0,ylim=c(0,0.15),col=effcols)
segments(dumx,dpbgciHi,dumx,dpbgciLo,lwd=2)
segments(ebx1,dpbgciHi,ebx2,dpbgciHi,lwd=2)

barplot(mswsbsest,xlim=c(0,7.2),main="Msexta Sperm Proteome",width=1,space=0,ylim=c(0,0.15),col=effcols)
segments(dumx,mswsciHi,dumx,mswsciLo,lwd=2)
segments(ebx1,mswsciHi,ebx2,mswsciHi,lwd=2)

barplot(dpwsbsest,xlim=c(0,7.2),main="Dplex Sperm Proteome",width=1,space=0,ylim=c(0,0.15),col=effcols)
segments(dumx,dpwsciHi,dumx,dpwsciLo,lwd=2)
segments(ebx1,dpwsciHi,ebx2,dpwsciHi,lwd=2)

barplot(msorbsest, space = 0, main="Msex Sperm Orthologs",names=c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<"), ylim=c(0,0.15),col=effcols)
segments(dumx,msorciHi,dumx,msorciLo,lwd=2)
segments(ebx1,msorciHi,ebx2,msorciHi,lwd=2)

barplot(dporbsest, space = 0, main="Dplex Sperm Orthologs",names=c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<"), ylim=c(0,0.15),col=effcols)
segments(dumx,dporciHi,dumx,dporciLo,lwd=2)
segments(ebx1,dporciHi,ebx2,dporciHi,lwd=2)





#barplot(bsest,xlim=c(0,7.2),width=1,space=0,ylim=c(0,0.09),main="Dplex Eupyrene")
#segments(dumx,ciHi,dumx,ciLo,lwd=4)

#Manduca




test<-c(bsest,dpbgbsest)
barplot(test, beside=T)

tops<-length(edf)
dfemat<-data.frame(matrix(NA, nrow = tops, ncol = 7))
for(i in 1:tops)
{
dfemat[i,]<-unlist(lapply(edf[[i]],getDiscretizedDFE))
}

bsest<-colMeans(dfemat)

sdsdfe<-apply(dfemat, 2, sd)
semdfe<-sdsdfe/sqrt(tops)

dumx<-seq(0.5,7.5,1)
ciHi<-bsest+2*semdfe
ciLo<-bsest-2*semdfe



#####scraps####
#now we fire up the polydfePost.R script to get the functions we need to analyze our data


dplexbgpoly<-parseOutput("/Users/Andrew/Documents/polyDFE/Dplex_bg_test_out_pr_B_no_div", init = FALSE)
dplexeupoly<-parseOutput("/Users/Andrew/Documents/polyDFE/Dplex_eu_test_out_pr_B_no_div", init = FALSE)
euAFR<-parseOutput("/Users/Andrew/Documents/polyDFE/Dplex_eu_test_out_pr_A_no_div", init = FALSE)
euCBFG<-parseOutput("/Users/Andrew/Documents/polyDFE/Dplex_eu_test_outCbfgs", init = FALSE)
bgCBFG<-parseOutput("/Users/Andrew/Documents/polyDFE/Dplex_bg_test_outCbfgs", init = FALSE)

euCBFGdfe<-getDiscretizedDFE(euCBFG[[1]], sRanges = c(-100, -10, -1, 0, 1, 10))
bgCBFGdfe<-getDiscretizedDFE(bgCBFG[[1]], sRanges = c(-100, -10, -1, 0, 1, 10))
euAFRdfe<-getDiscretizedDFE(euAFR[[1]], sRanges = c(-100, -10, -1, 0, 1, 10))
bgdfe<-getDiscretizedDFE(dplexbgpoly[[1]], sRanges = c(-100, -10, -1, 0, 1, 10))
eudfe<-getDiscretizedDFE(dplexeupoly[[1]], sRanges = c(-100, -10, -1, 0, 1, 10))

par(mfrow=c(2,1))
barplot(bsest,xlim=c(0,7.2),names=c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<"),width=1,space=0,ylim=c(0,0.89))
segments(dumx,ciHi,dumx,ciLo,lwd=4)

barplot(dpbgbsest,xlim=c(0,7.2),names=c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<"),width=1,space=0,ylim=c(0,0.89))
segments(dumx,dpbgciHi,dumx,dpbgciLo,lwd=4)


par(mfrow=c(2,1))
barplot(msbgbsest,xlim=c(0,7.2),names=c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<"),width=1,space=0,ylim=c(0,0.89))
segments(dumx,msbgciHi,dumx,msbgciLo,lwd=4)

barplot(dpbgbsest,xlim=c(0,7.2),names=c("<-100","(-100,-10)","(-10,-1)","(-1,0)","(0,1)","(1,10)","10<"),width=1,space=0,ylim=c(0,0.89))
segments(dumx,dpbgciHi,dumx,dpbgciLo,lwd=4)


