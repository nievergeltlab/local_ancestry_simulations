echo ""afeffectsize" "eureffectsize" "afmaf" "eurmaf" "admixturegroup" "N" "prevalence" "power_newanc" "power_oldanc"" > header.txt 

cat sims_jan2/power500sims*.txt | grep -v afeffectsize | cat header.txt  - |  awk '{if (NR==1 || $1 != ""afeffectsize"") print}' > jan10_powersims.txt


R

library(plyr)
powergivenn1 <- read.table('jan10_powersims.txt',header=T,stringsAsFactors=F)
names(powergivenn)[8:9] 

powergivenn <- aggregate(cbind(power_oldanc,power_newanc) ~afeffectsize + eureffectsize +afmaf + eurmaf+ admixturegroup + N +prevalence ,data=powergivenn1,FUN=mean)

powergivenn_alt <-ddply(powergivenn1, ~ afeffectsize + eureffectsize +afmaf + eurmaf+ admixturegroup + N +prevalence  ,colwise(mean,c("power_oldanc","power_newanc"),na.rm=T ))


#wesanderson::wes_palette
fineness=.005


#Subset data to one comparison

#I dont think prevalence matters.
#compare aggregate and ddply functions as a test
#effect sizes:  euronly matching     null stronger   weaker
#eur maf options: larger smaller same (60%?)
#afr maf options 0.1 0.2 0.3 0.4
#admixture group : 0.8 0.5



d1_alt <- subset(powergivenn_alt,eureffectsize == "matching"&afmaf==0.1&eurmaf=="smaller"&admixturegroup==0.5&prevalence==0.1)
d1_alt <- d1[order(d1_alt$afeffectsize),]



d1 <- subset(powergivenn,eureffectsize == "null"&afmaf==0.2&eurmaf=="same"&admixturegroup==0.5&prevalence==0.2)
d1 <- d1[order(d1$afeffectsize),]



d1 <- subset(powergivenn,eureffectsize == "null"&afmaf==0.4&eurmaf=="same"&admixturegroup==0.8&prevalence==0.1)
d1 <- subset(powergivenn,eureffectsize == "weaker"&afmaf==0.1&eurmaf=="larger"&admixturegroup==0.8&prevalence==0.2)



##Just vary the maf
efparm='matching'
afmafparm='varying'
eurmafparm='same'
admixtureparm=080
prevalenceparm=010

pdf(paste(efparm,'_',afmafparm,'_',eurmafparm,'_',admixtureparm,'_',prevalenceparm, '.pdf',sep=''),8,4)
for (maf in c(0.1,0.2,0.3,0.4))
{
 d1 <- subset(powergivenn,eureffectsize == "matching"&afmaf==maf&eurmaf=="same"&admixturegroup==0.8&prevalence==0.1)
 d1 <- d1[order(d1$afeffectsize),]

 plot(d1$afeffectsize,100*d1$power_newanc/100,type='l',lwd=2,cex.axis=1.25,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,100),main=paste("Overall MAF set to", maf))
 lines(d1$afeffectsize,100*d1$power_old/100,type='l',lwd=2,col='darkgrey',lty=1)
 legend("topleft",legend=c("12000 Cases, LANC","12000 Cases, GLOB" ), col=c('black','darkgrey'),lty=c(2,1),bty='n',cex=.5, pt.cex = .4)

}
dev.off()

#vary the maf but have null euro effect
efparm='null'
afmafparm='varying'
eurmafparm='same'
admixtureparm=080
prevalenceparm=010

pdf(paste(efparm,'_',afmafparm,'_',eurmafparm,'_',admixtureparm,'_',prevalenceparm, '.pdf',sep=''),8,4)
for (maf in c(0.1,0.2,0.3,0.4))
{
 d1 <- subset(powergivenn,eureffectsize == "null"&afmaf==maf&eurmaf=="same"&admixturegroup==0.8&prevalence==0.1)
 d1 <- d1[order(d1$afeffectsize),]

 plot(d1$afeffectsize,100*d1$power_newanc/100,type='l',lwd=2,cex.axis=1.25,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,100),main=paste("Overall MAF set to", maf))
 lines(d1$afeffectsize,100*d1$power_old/100,type='l',lwd=2,col='darkgrey',lty=1)
 legend("topleft",legend=c("12000 Cases, LANC","12000 Cases, GLOB" ), col=c('black','darkgrey'),lty=c(2,1),bty='n',cex=.5, pt.cex = .4)

}
dev.off()

##Conclusion : MAF works like in a regular GWAS

##vary the effect size
#vary the maf but have null euro effect
efparm='stronger'
afmafparm='020'
eurmafparm='same'
admixtureparm=080
prevalenceparm=010

pdf(paste('varying','_',afmafparm,'_',eurmafparm,'_',admixtureparm,'_',prevalenceparm, '.pdf',sep=''),8,4)
for (efparm in c("euronly", "matching"  ,  "null", "stronger" ,  "weaker"))
{
 d1 <- subset(powergivenn,eureffectsize == efparm &afmaf==0.2&eurmaf=="same"&admixturegroup==0.8&prevalence==0.1)
 d1 <- d1[order(d1$afeffectsize),]

 plot(d1$afeffectsize,100*d1$power_newanc/100,type='l',lwd=2,cex.axis=1.25,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,100),main=paste("Effect size in europeans:", efparm))
 lines(d1$afeffectsize,100*d1$power_old/100,type='l',lwd=2,col='darkgrey',lty=1)
 legend("topleft",legend=c("12000 Cases, LANC","12000 Cases, GLOB" ), col=c('black','darkgrey'),lty=c(2,1),bty='n',cex=.5, pt.cex = .4)

}
dev.off()

##vary  admixture and effect
efparm='stronger'
afmafparm='020'
eurmafparm='same'
admixtureparm=080
prevalenceparm=010

pdf(paste('varyingef','_',afmafparm,'_',eurmafparm,'_',"varyingadmix",'_',prevalenceparm, '.pdf',sep=''),8,4)
for (efparm in c("euronly", "matching"  ,  "null", "stronger" ,  "weaker"))
{
 for (admixparm in c(0.8, 0.5))
 {
 d1 <- subset(powergivenn,eureffectsize == efparm &afmaf==0.2&eurmaf=="same"&admixturegroup==admixparm&prevalence==0.1)
 d1 <- d1[order(d1$afeffectsize),]

 plot(d1$afeffectsize,100*d1$power_newanc/100,type='l',lwd=2,cex.axis=1.25,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,100),main=paste("Effect size in europeans:", efparm, ".Admixture rate:", admixparm))
 lines(d1$afeffectsize,100*d1$power_old/100,type='l',lwd=2,col='darkgrey',lty=1)
 legend("topleft",legend=c("12000 Cases, LANC","12000 Cases, GLOB" ), col=c('black','darkgrey'),lty=c(2,1),bty='n',cex=.5, pt.cex = .4)
 } 
}
dev.off()



##vary  prevalence and effect
efparm='varyingef'
afmafparm=0.2
eurmafparm='same'
admixparm=0.8
prevalenceparm=varying

pdf(paste(efparm,'_',afmafparm,'_',eurmafparm,'_',admixtureparm,'_',prevalenceparm, '.pdf',sep=''),8,4)
for (efparm in c("euronly", "matching"  ,  "null", "stronger" ,  "weaker"))
{
 for (prevalenceparm in c(0.1,0.2))
 {
 d1 <- subset(powergivenn,eureffectsize == efparm &afmaf==afmafparm&eurmaf==eurmafparm&admixturegroup==admixparm&prevalence==prevalenceparm)
 d1 <- d1[order(d1$afeffectsize),]

 plot(d1$afeffectsize,100*d1$power_newanc/100,type='l',lwd=2,cex.axis=1.25,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,100),main=paste("Effect size in europeans:", efparm, ". Prevalence:", prevalenceparm))
 lines(d1$afeffectsize,100*d1$power_old/100,type='l',lwd=2,col='darkgrey',lty=1)
 legend("topleft",legend=c("12000 Cases, LANC","12000 Cases, GLOB" ), col=c('black','darkgrey'),lty=c(2,1),bty='n',cex=.5, pt.cex = .4)
 } 
}
dev.off()



 ##vary  maf
efparm='varyingef'
afmafparm=0.2
eurmafparm='same'
admixparm=0.8
prevalenceparm=0.2

pdf(paste(efparm,'_',afmafparm,'_',eurmafparm,'_',admixtureparm,'_',prevalenceparm, '.pdf',sep=''),8,4)
for (efparm in c("euronly", "matching"  ,  "null", "stronger" ,  "weaker"))
{
 for (eurmafparm in c("larger", "smaller","same"))
 {
 d1 <- subset(powergivenn,eureffectsize == efparm &afmaf==afmafparm&eurmaf==eurmafparm&admixturegroup==admixparm&prevalence==prevalenceparm)
 d1 <- d1[order(d1$afeffectsize),]

 plot(d1$afeffectsize,100*d1$power_newanc/100,type='l',lwd=2,cex.axis=1.25,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,100),main=paste("Effect size in europeans:", efparm, ". Maf difference:", eurmafparm))
 lines(d1$afeffectsize,100*d1$power_old/100,type='l',lwd=2,col='darkgrey',lty=1)
 legend("topleft",legend=c("12000 Cases, LANC","12000 Cases, GLOB" ), col=c('black','darkgrey'),lty=c(2,1),bty='n',cex=.5, pt.cex = .4)
 } 
}
dev.off()



 
 

 
 
 
 
 
 
 


#Graph w/o green lne
pdf('powersimg2_scen1_v2.pdf',5.5,4.5)

plot(efseq,powergivenn4000[,1]/10,type='l',lwd=2,cex.axis=1.15,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,110),yaxt='n',bty='l')

abline(a=80,b=0,lty=3)
axis(2,c(0,20,40,60,80,100),cex.axis=1.15) #jesus christ you fuck jkust tra
lines(efseq,powergivenn4000[,3]/10,type='l',lwd=2,col='black',lty=1)

lines(efseq,powergivenn12000[,1]/10,type='l',lwd=2,col='blue',lty=2)
lines(efseq,powergivenn12000[,3]/10,type='l',lwd=2,col='blue',lty=1)

lines(efseq,powergivenn13[,1]/10,type='l',lwd=2,col='red',lty=2)
lines(efseq,powergivenn13[,3]/10,type='l',lwd=2,col='red',lty=1)

legend("topleft",legend=c("12000 Cases, LANC, MAF20%","12000 Cases, GLOB" ,"12000 Cases, LANC, MAF10% AFR, 30% EUR" ,"12000 Cases, GLOB" ,"4000 Cases, LANC, MAF20%" ,"4000 Cases, GLOB" ), col=rep(c('blue','blue','red','red','black','black'),3),lty=rep(c(2,1),3),bty='n',cex=.5, pt.cex = .4)

dev.off()



#Graph with black and blue (grant). Put green color first in legend
pdf('powersimg2_scen2.pdf',5.5,4.5)
 par(mar=c(5, 4, 4, 2) + 0.5)
plot(efseq,powergivenn4000[,1]/10,type='l',lwd=2,cex.axis=1.4,cex.lab=1.8,xlab="", ylab="",lty=2,ylim=c(0,110),yaxt='n',bty='l')

mtext(side=1,"Odds Ratio",cex=1.8,line=2.6)
mtext(side=2,"Power",cex=1.8,line=2.6)
 
 
abline(a=80,b=0,lty=3)
axis(2,c(0,20,60,100),cex.axis=1.4) 
axis(2,c(40,80),cex.axis=1.4) 

lines(efseq,powergivenn4000[,3]/10,type='l',lwd=2,col='black',lty=1)

lines(efseq,powergivenn12000[,1]/10,type='l',lwd=2,col='blue',lty=2)
lines(efseq,powergivenn12000[,3]/10,type='l',lwd=2,col='blue',lty=1)

lines(efseq,powergivenn13[,1]/10,type='l',lwd=2,col='green',lty=2)
lines(efseq,powergivenn13[,3]/10,type='l',lwd=2,col='green',lty=1)
par(xpd=TRUE)
legend(x=1.035,y=130,legend=c("LANC, MAF10% AFR, 30% EUR" ,"GLOB, MAF10% AFR, 30% EUR"  ,"LANC, MAF20%","GLOB, MAF20%" ,"LANC, MAF20%" ,"GLOB, MAF 20%" ), col=rep(c('green','green','blue','blue','black','black'),3),lty=rep(c(2,1),3),bty='n',cex=.85, pt.cex = .85)
#legend("topleft",legend=c("12000 Cases, LANC, MAF10% AFR, 30% EUR" ,"12000 Cases, GLOB"  ,"12000 Cases, LANC, MAF20%","12000 Cases, GLOB" ,"4000 Cases, LANC, MAF20%" ,"4000 Cases, GLOB" ), col=rep(c('green','green','blue','blue','black','black'),3),lty=rep(c(2,1),3),bty='n',cex=.5, pt.cex = .4)

dev.off()

#Graph with red and green

pdf('powersimg2_scen3_v2.pdf',5.5,4.5)

plot(efseq,powergivenn4000[,1]/10,type='l',lwd=2,cex.axis=1.15,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,110),yaxt='n',bty='l',col='white')
axis(2,c(0,20,40,60,80,100),cex.axis=1.15) 

abline(a=80,b=0,lty=3)
lines(efseq,powergivenn24[,1]/10,type='l',lwd=2,col='green',lty=2)
lines(efseq,powergivenn24[,3]/10,type='l',lwd=2,col='green',lty=1)

lines(efseq,powergivenn13[,1]/10,type='l',lwd=2,col='red',lty=2)
lines(efseq,powergivenn13[,3]/10,type='l',lwd=2,col='red',lty=1)

legend("topleft",legend=c("12000 Cases, LANC, MAF20% AFR, 40% EUR","12000 Cases, GLOB" ,"12000 Cases, LANC, MAF10% AFR, 30% EUR" ,"12000 Cases, GLOB" ), col=rep(c('green','green','red','red'),3),lty=rep(c(2,1),2),bty='n',cex=.5, pt.cex = .4)

dev.off()




#graph all 
pdf('powersimg2_202020401030.pdf',5.5,4.5)

plot(efseq,powergivenn4000[,1]/10,type='l',lwd=2,cex.axis=1.15,cex.lab=1.45,xlab="Odds Ratio", ylab="Power",lty=2,ylim=c(0,110),yaxt='n',bty='l')
axis(2,c(0,20,40,60,80,100),cex.axis=1.15) #jesus christ you fuck jkust tra
lines(efseq,powergivenn4000[,3]/10,type='l',lwd=2,col='black',lty=1)
abline(h=0)

lines(efseq,powergivenn12000[,1]/10,type='l',lwd=2,col='blue',lty=2)
lines(efseq,powergivenn12000[,3]/10,type='l',lwd=2,col='blue',lty=1)

lines(efseq,powergivenn24[,1]/10,type='l',lwd=2,col='red',lty=2)
lines(efseq,powergivenn24[,3]/10,type='l',lwd=2,col='red',lty=1)

lines(efseq,powergivenn13[,1]/10,type='l',lwd=2,col='green',lty=2)
lines(efseq,powergivenn13[,3]/10,type='l',lwd=2,col='green',lty=1)
legend("topleft",legend=c("12000 Cases, asaMap","12000 Cases, standard" ), col=rep(c('red','pink'),3),lty=c(2,1),bty='n')


dev.off()

