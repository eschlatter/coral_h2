setwd("/Users/Ted/Desktop/RNAseq_Larvae/Orpheus/Orpheus/protein")

st=read.csv("Protein_Orpheus.csv")
attach(st)
summary(st)
st$rep=as.factor(st$rep)
st$loss=((st$wgg-st$day3)/st$wgg)
summary(st)
hist(st$loss)

summary(st)
library(ggplot2)
source("summarySE.R")
st$culture=paste(st$mom,st$dad,st$rep, sep="")
all.marg=summarySE(st,measurevar="loss",groupvars=c("culture"), na.rm=T)
all.marg
write.csv(all.marg, file = "protein_av.csv")


library(lme4) 
m0=lmer(loss~(1|mom), st)
m1=lmer(loss~(1|mom)+(1|dad),st)
m2=lmer(loss~(1|mom)*(1|dad),st)
anova(m0,m1,m2)
summary(m1)
#Linear mixed model fit by REML 
#Formula: loss ~ (1 | mom) + (1 | dad) 
#   Data: st 
#    AIC    BIC logLik deviance REMLdev
# -261.7 -252.6  134.9   -275.7  -269.7
#Random effects:
# Groups   Name        Variance   Std.Dev.
# dad      (Intercept) 0.00029229 0.017097
# mom      (Intercept) 0.00070482 0.026548
# Residual             0.00103201 0.032125
#Number of obs: 72, groups: dad, 9; mom, 2

lm.va=VarCorr(m1)$dad[1]
vres= as.numeric ( attributes ( summary(m1))$REmat [ attributes (summary(m1)) $REmat [,1] == "Residual" , 3 ] ) 
> h2lmer=2*lm.va/(lm.va+vres)  
> h2lmer 
# dad explains 44% of the variation!! [1] 0.4414314

par(mfrow=c(2,2))
plot(aov1)


aov1=aov(loss~dad, data=st)
summary(aov1)
TukeyHSD(aov1)
plot(loss~dad, data=st)

# Analysing only Mom A across all dads
a=subset(st, mom=="A")
summary(a)
A1=aov(loss~dad, data=a)
summary(A1)
TukeyHSD(A1)
plot(loss~dad, data=a)


# Analysing only Mom C across all dads
c=subset(st, mom=="C")
summary(c)

amil2=aov(loss~dad, data=c)
summary(amil2)
TukeyHSD(amil2)
plot(loss~dad, data=c)


library(ggplot2)
source("summarySE.R")
all.marg=summarySE(st,measurevar="loss",groupvars=c("dad", "mom"), na.rm=T)
all.marg
pd=position_dodge(.3)
ggplot(all.marg,aes(x=dad,y=loss,colour=mom,group=mom))+
 	geom_errorbar(aes(ymin=loss+se,ymax=loss-se),lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=mom,pch=mom),position=pd,size=2.5)+
	theme_bw()

qplot(dad, loss, data = st, geom = "boxplot", xlab = "", ylab = "Protein Loss") +
theme_bw() +
opts(axis.title.x =theme_text(size = 14), axis.title.y = theme_text(size= 14, angle = 90, vjust = 0.4), axis.text.x = theme_text(size = 12), axis.text.y = theme_text(size = 12))


qplot(dad, loss, data = st, geom = "boxplot", xlab = "", ylab = "Protein Loss") +
theme_bw() +
opts(axis.title.x =theme_text(size = 14), axis.title.y = theme_text(size= 14, angle = 90, vjust = 0.4), axis.text.x = theme_text(size = 12), axis.text.y = theme_text(size = 12)) +
ylim(0,0.15)
par(mfrow=c(2,2))
qplot(dad, loss, data = a, geom = "boxplot", xlab = "", ylab = "Protein Loss") +
theme_bw() +
opts(axis.title.x =theme_text(size = 14), axis.title.y = theme_text(size= 14, angle = 90, vjust = 0.4), axis.text.x = theme_text(size = 12), axis.text.y = theme_text(size = 12)) +
ylim(0,0.15)
qplot(dad, loss, data = c, geom = "boxplot", xlab = "", ylab = "Protein Loss") +
theme_bw() +
opts(axis.title.x =theme_text(size = 14), axis.title.y = theme_text(size= 14, angle = 90, vjust = 0.4), axis.text.x = theme_text(size = 12), axis.text.y = theme_text(size = 12)) +
ylim(0,0.15)

# REML modeling with proportions, sire effects for Mom A at 48hrs
library(lme4) 
m0=lmer(asin(sqrt(prop))~(1|pl),a48)
m1=lmer(asin(sqrt(prop))~(1|pl)+(1|dad),a48)
anova(m0,m1)
summary(m1)
#Data: a48
#Models:
#m0: asin(sqrt(prop)) ~ (1 | pl)
#m1: asin(sqrt(prop)) ~ (1 | pl) + (1 | dad)
#   Df    AIC    BIC  logLik  Chisq Chi Df Pr(>Chisq)   
#m0  3 22.576 28.859 -8.2878                            
#m1  4 17.663 26.041 -4.8316 6.9124      1    0.00856 **
#Significant effect of dad
#Now plot the dad to see
library(ggplot2)
qplot(dad, asin(sqrt(prop)), data=a)
qplot(dad, asin(sqrt(prop)), data = a48, geom = "boxplot", xlab = "Sire", ylab = "Proportion of Settlement", main = "Dame A") +
opts(axis.text.x = theme_text(angle = 90, hjust = 1, size = 8)) +
theme_bw()

summary(a)
amil1=aov(asin(sqrt(prop))~time+dad, data=a)
summary(amil1)
TukeyHSD(amil1)

# REML modeling with proportions, sire effects for Mom C at 48hrs
m2=lmer(asin(sqrt(prop))~(1|pl),c48)
m3=lmer(asin(sqrt(prop))~(1|pl)+(1|dad),c48)
anova(m2,m3)
summary(m3)
#Data: c48
#Models:
#m2: asin(sqrt(prop)) ~ (1 | pl)
#m3: asin(sqrt(prop)) ~ (1 | pl) + (1 | dad)
#   Df      AIC    BIC logLik Chisq Chi Df Pr(>Chisq)  
#m2  3 -0.34682 5.9362 3.1734                          
#m3  4 -1.88878 6.4886 4.9444 3.542      1    0.05983 .

qplot(dad, asin(sqrt(prop)), data=c48)
qplot(dad, asin(sqrt(prop)), data = c48, geom = "boxplot", xlab = "Sire", ylab = "Proportion of Settlement", main = "Dame C") +
opts(axis.text.x = theme_text(angle = 90, hjust = 1, size = 8)) +
theme_bw()

#include repeated measures
m10=lmer(asin(sqrt(prop))~(1|pl),a)
m11=lmer(asin(sqrt(prop))~(1|pl)+(1|time), a)
m12=lmer(asin(sqrt(prop))~(1|pl)+(1|time)+(1|dad),a)
anova(m10,m11,m12)
summary(m12)
##not working
library(multcomp)
aTukey<- glht(m12, linfct=mcp(varfert="Tukey")) 
summary(aTukey)

#------------------------------
# MCMC modeling with proportions, sire effects

library(MCMCglmm) 
mpr=MCMCglmm(asin(sqrt(prop))~mom,random=~pl+dad,data=st,verbose=F,family="gaussian",nitt=100000,thin=20,burnin=70000)
summary(mpr)
mc.va=posterior.mode(mpr$VCV[,'dad'])
h2p=2*mpr$VCV[,'dad']/(mpr$VCV[,'dad']+mpr$VCV[,'pl']+mpr$VCV[,'units']) # heritability
plot(h2p)
posterior.mode(h2p)  # 0.45 - 0.5
HPDinterval(h2p) # 0.07-0.96
evol.mc=mc.va/(mset^2)
evol.mc # evolvability = 0.028

#------------------------------
# REML modeling with proportions, sire effects
library(lme4) 
m0=lmer(asin(sqrt(prop))~(1|pl),st)
m1=lmer(asin(sqrt(prop))~mom+(1|pl),st)
m2=lmer(asin(sqrt(prop))~mom+(1|pl)+(1|dad),st)
anova(m0,m1,m2)
summary(m2)
lm.va=VarCorr(m2)$dad[1]
vpl=VarCorr(m2)$pl[1]
vres= as.numeric ( attributes ( summary(m2))$REmat [ attributes (summary(m2)) $REmat [,1] == "Residual" , 3 ] ) 
h2lmer=2*lm.va/(lm.va+vpl+vres)  
h2lmer # 0.45
evol.lm=lm.va/(mset^2) # evolvability = 0.036


# Analysing only dad V and W 
vw=read.csv("Protein.test.csv")
summary(vw)
#st$Rep=as.factor(st$Rep)
vw$loss=vw$wgg-vw$day3
summary(vw)
hist(vw$loss)
par(mfrow=c(2,2))
plot(aov1)

aov1=aov(loss~dad, data=vw)
summary(aov1)
TukeyHSD(aov1)
plot(loss~dad, data=st)

# Analysing only Mom A across all dads
a=subset(vw, mom=="A")
summary(a)
A1=aov(loss~dad, data=a)
summary(A1)
TukeyHSD(A1)
plot(loss~dad, data=a)


# Analysing only Mom C across all dads
c=subset(vw, mom=="C")
summary(c)

amil2=aov(loss~dad, data=c)
summary(amil2)
TukeyHSD(amil2)
plot(loss~dad, data=c)

#regressing settlement to prop protein loss
set=read.csv("Settlement_Orpheus48.csv")
set$rep=as.factor(set$rep)
summary(set)
set$prop=(set$rec/(set$rec+set$lar))
library(ggplot2)
source("summarySE.R")
all.marg=summarySE(st,measurevar="loss",groupvars=c("dad", "mom", "rep"), na.rm=T)
all.marg
all.marg1=summarySE(set,measurevar="prop",groupvars=c("dad", "mom", "rep"), na.rm=T)
all.marg1
plot1=lm(all.marg$loss~all.marg1$prop)
summary(plot1)
plot(all.marg$loss~all.marg1$prop, xlab="Protein Loss", ylab="Settlement", cex=1.5, col="blue", cex.axis=1, cex.lab=1.2)
abline(plot1 ,col="red")

pd=position_dodge(.3)
ggplot(all.marg,aes(x=dad,y=loss,colour=mom,group=mom))+
 	geom_errorbar(aes(ymin=loss+se,ymax=loss-se),lwd=0.4,width=0.3,position=pd)+
	geom_line(aes(group=mom,linecol=mom),position=pd)+
	geom_point(aes(group=mom,pch=mom),position=pd,size=2.5)+
	theme_bw()
	
	
#############New in July
st=read.csv("Protein_Orpheus.csv")
attach(st)
summary(st)
st$rep=as.factor(st$rep)
st$loss=((st$wgg-st$day3)/st$wgg)
summary(st)
str(st)
all.margrg=summarySE(st,measurevar="loss",groupvars=c("mom", "dad"), na.rm=T)
all.margrg
pd=position_dodge(.3)
ggplot(all.margrg,aes(x=dad,y=loss,colour=mom,group=mom))+
 	geom_errorbar(aes(ymin=loss+se,ymax=loss-se), lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=mom,pch=mom),position=pd,size=2.5)+
	xlab("Sire")+
	ylab("Protein Loss (mg)")+
	theme_bw()	
	
# REML modeling with proportions, sire effects for Mom A at 48hrs
library(lme4) 
st$culture=paste(st$mom,st$dad,st$rep,sep="")

m2=lmer(loss~(1|dad),st)
m3=lmer(loss~(1|mom),st)
m4=lmer(loss~(1|dad)+(1|mom),st)
m5=lmer(loss~(1|dad)*(1|mom), st)
m6=lmer(loss~(1|dad)*(1|mom)+culture, st)
anova(m2,m3,m4,m5,m6)

lm.va=VarCorr(m6)$dad[1]
vres= as.numeric ( attributes ( summary(m6))$REmat [ attributes (summary(m6)) $REmat [,1] == "Residual" , 3 ] ) 
h2lmer=2*lm.va/(lm.va+vres)  
h2lmer # 0.49999

#------------------------------
# MCMC modeling
summary(st)
library(MCMCglmm) 
hist(log(st$loss+1))

mc=MCMCglmm(loss~1,random=~dad+mom+mom:dad+culture,family="gaussian",data=st)

summary(mc)
mmm=mc
names(colMeans(mmm$VCV))

bigH=c()
 for(i in 1:length(mmm$VCV[,1])) {
 	bigH=append(bigH,sum(mmm$VCV[i,1:3])/sum(mmm$VCV[i,1:5]))
  }
 quantile(bigH,0.025) #1.595138e-09       
 quantile(bigH,0.975) #0.7732985 
   
mean(bigH) # 0.1383126


mom=c()
for(i in 1:length(mmm$VCV[,1])) {
	mom=append(mom,sum(mmm$VCV[i,2])/sum(mmm$VCV[i,1:5]))
}
quantile(mom,0.025) # super small
quantile(mom,0.975) # 0.74
mean(mom) # 0.05

dad=c()
for(i in 1:length(mmm$VCV[,1])) {
	dad=append(dad,(mmm$VCV[i,1])/sum(mmm$VCV[i,1:5]))
}
quantile(dad,0.025) # really small
quantile(dad,0.975) # 0.09
mean(dad) # 0.007

int=c()
for(i in 1:length(mmm$VCV[,1])) {
	int=append(int,(sum(mmm$VCV[i,3]))/sum(mmm$VCV[i,1:5]))
}
quantile(int,0.025) # really small
quantile(int,0.975) # .025
mean(int) # 0.002

h.means=c(mean(mom),mean(dad),mean(int))
h.lower=c(quantile(mom,0.025),quantile(dad,0.025),quantile(int,0.025))
h.upper=c(quantile(mom,0.975),quantile(dad, 0.975),quantile(int, 0.975))
genetics=data.frame("variance.explained"=h.means,h.lower,h.upper)
genetics=cbind("component"=c("Dam","Sire","Interaction"),genetics)

qplot(component,variance.explained, data=genetics,geom="bar",fill=component)+geom_errorbar(aes(ymin=h.lower,ymax=h.upper),lwd=0.4,width=0.3)+theme_bw()+scale_x_discrete(breaks=NULL)

#getting average for WGCNA
all.margp=summarySE(st,measurevar="loss",groupvars=c("culture"), na.rm=T)
all.margp
write.csv(all.margp, file = "proteinloss_av_new.csv")



