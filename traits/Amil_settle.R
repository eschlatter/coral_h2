setwd("/Users/Ted/Desktop/RNAseq_Larvae/Orpheus/Orpheus/settlement")
library(lme4)
st=read.csv("A.millepora_Settlement_2011.csv")
attach(st)
st$rep=as.factor(st$rep)
st$time=as.factor(st$time)
summary(st)
head(st)

# proportions
prop=st$rec/(st$rec+st$lar)
prop[prop==0]=0.0001
prop[prop==1]=0.9999
st$prop=prop

source("summarySE.R")
st$culture=paste(st$mom,st$dad,st$rep, sep="")
all.marg=summarySE(st,measurevar="prop",groupvars=c("culture"), na.rm=T)
all.marg
library(ggplot2)
pd=position_dodge(.3)
ggplot(all.marg,aes(x=culture,y=prop))+
 	geom_errorbar(aes(ymin=prop+se,ymax=prop-se),lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=culture),position=pd,size=2.5)+
	theme_bw()+
	scale_colour_grey(start = 0, end = .7)

write.csv(all.marg, file = "settle_av.csv")



aov1=lmer(asin(sqrt(prop))~mom+(1|dad), data=st)

# Analysing only Mom A across all dads
a=subset(st, mom=="A")
summary(a)
amil1=aov(asin(sqrt(prop))~time+dad, data=a)
summary(amil1)
TukeyHSD(amil1)

# Analysing only Mom A across all dads at time point 48hrs
a48=subset(a, time=="48")
summary(a48)

la48=lm(asin(sqrt(prop))~pl+dad,a48)
par(mfrow=c(2,2))
plot(la48) #settlement_bysire_diagnostics.pdf
summary(la48)  # dadW is the only significant one...

# Analysing only Mom A across all dads at time point 24hrs
a24=subset(a, time=="24")
summary(a24)

la24=lm(asin(sqrt(prop))~pl+dad,a24)
par(mfrow=c(2,2))
plot(la24) #settlement_bysire_diagnostics.pdf
summary(la24)  # still only dadW is the only significant one...


# Analysing only Mom C across all dads
c=subset(st, mom=="C")
summary(c)

amil2=aov(asin(sqrt(prop))~time+dad, data=c)
summary(amil2)
TukeyHSD(amil2)


# Analysing only Mom C across all dads at time point 48hrs
c48=subset(c, time=="48")
summary(c48)

lc48=lm(asin(sqrt(prop))~pl+dad,c48)
par(mfrow=c(2,2))
plot(lc48) #settlement_bysire_diagnostics.pdf
summary(lc48)  # dadW is the only significant one...

# Analysing only Mom C across all dads at time point 24hrs
c24=subset(c, time=="24")
summary(c24)

lc24=lm(asin(sqrt(prop))~pl+dad,c24)
par(mfrow=c(2,2))
plot(lc24) #settlement_bysire_diagnostics.pdf
summary(lc24)  # still only dadW is the only significant one...

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
qplot(dad, prop, data = a48, geom = "boxplot", xlab = "Sire", ylab = "Proportion of Settlement", main = "Dame A") +
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



#############New in July
summary(st)
str(st)
all.margrg=summarySE(st,measurevar="prop",groupvars=c("mom", "dad"), na.rm=T)
all.margrg
pd=position_dodge(.3)
ggplot(all.margrg,aes(x=dad,y=prop,colour=mom,group=mom))+
 	geom_errorbar(aes(ymin=prop+se,ymax=prop-se),lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=mom,pch=mom),position=pd,size=2.5)+
	xlab("Sire")+
	ylab("Proportion of Settlement")+
	theme_bw()	
	
# REML modeling with proportions, sire effects for Mom A at 48hrs
library(lme4) 
hist(asin(sqrt(st$prop)))
st$culture=paste(st$mom,st$dad,st$rep,sep="")

m1=lmer(asin(sqrt(prop))~(1|pl)+time,st)
m2=lmer(asin(sqrt(prop))~(1|pl)+time+(1|dad),st)
m3=lmer(asin(sqrt(prop))~(1|pl)+time+(1|mom),st)
m4=lmer(asin(sqrt(prop))~(1|pl)+time+(1|dad)+(1|mom),st)
m5=lmer(asin(sqrt(prop))~(1|pl)+time+(1|dad)*(1|mom), st)
m6=lmer(asin(sqrt(prop))~(1|pl)+time+(1|dad)*(1|mom)+culture, st)
anova(m1,m2,m3,m4,m5,m6)

lm.va=VarCorr(m6)$dad[1]
vpl=VarCorr(m6)$pl[1]
vres= as.numeric ( attributes ( summary(m6))$REmat [ attributes (summary(m6)) $REmat [,1] == "Residual" , 3 ] ) 
h2lmer=2*lm.va/(lm.va+vpl+vres)  
h2lmer # 0.158

#------------------------------
# MCMC modeling
summary(st)
library(MCMCglmm) 

mc=MCMCglmm(asin(sqrt(prop))~time,random=~idh(time):pl+idh(time):dad+idh(time):mom+idh(time):mom:dad,rcov=~idh(time):units,family="gaussian", data=st)

summary(mc)
mmm=mc
names(colMeans(mmm$VCV))

mom=c()
for(i in 1:length(mmm$VCV[,1])) {
	mom=append(mom,sum(mmm$VCV[i,c(5,6)])/sum(mmm$VCV[i,c(1:10)]))
}
quantile(mom,0.025) # small
quantile(mom,0.975) # 0.98
mean(mom) # 0.23

dad=c()
for(i in 1:length(mmm$VCV[,1])) {
	dad=append(dad,sum(mmm$VCV[i,c(3,4)])/sum(mmm$VCV[i,c(1:10)]))
}
quantile(dad,0.025) # really small
quantile(dad,0.975) # 0.307
mean(dad) # 0.045

int=c()
for(i in 1:length(mmm$VCV[,1])) {
	int=append(int,sum(mmm$VCV[i,8])/sum(mmm$VCV[i,c(4,6,8,10, 12)]))
}
quantile(int,0.025) # small
quantile(int,0.975) # .45
mean(int) # 0.11

h.means=c(mean(mom),mean(dad),mean(int))
h.lower=c(quantile(mom,0.025),quantile(dad,0.025),quantile(int,0.025))
h.upper=c(quantile(mom,0.975),quantile(dad, 0.975),quantile(int, 0.975))
genetics=data.frame("variance.explained"=h.means,h.lower,h.upper)
genetics=cbind("component"=c("Dam","Sire","Interaction"),genetics)

qplot(component,variance.explained, data=genetics,geom="bar",fill=component)+geom_errorbar(aes(ymin=h.lower,ymax=h.upper),lwd=0.4,width=0.3)+theme_bw()+scale_x_discrete(breaks=NULL)


# MCMC modeling without plate effect
summary(st)
library(MCMCglmm) 
prop=st$rec/(st$rec+st$lar)
prop[prop==0]=0.0001
prop[prop==1]=0.9999
st$prop=prop

mc=MCMCglmm(asin(sqrt(prop))~time,random=~idh(time):dad+idh(time):mom+idh(time):mom:dad+idh(time):pl,data=st)

summary(mc)
mmm=mc
names(colMeans(mmm$VCV))

bigH=c()
for(i in 1:length(mmm$VCV[,1])) {
	bigH=append(bigH,sum(mmm$VCV[i,1:6])/sum(mmm$VCV[i,1:9]))
}
quantile(bigH,0.025) # 0.1835714
quantile(bigH,0.975) # 0.9949671 
mean(bigH) # 0.5570522

smallH=c()
for(i in 1:length(mmm$VCV[,1])) {
	smallH=append(smallH,2*sum(mmm$VCV[i,1:2])/sum(mmm$VCV[i,1:9]))
}
quantile(smallH,0.025) #0.0003272761
quantile(smallH,0.975) #0.9100735 
mean(smallH) # 0.3268541

moms=c()
for(i in 1:length(mmm$VCV[,1])) {
	moms=append(moms,(sum(mmm$VCV[i,2:4])-sum(mmm$VCV[i,1:2]))/sum(mmm$VCV[i,1:9]))
}
quantile(moms,0.025) # -0.03923223
quantile(moms,0.975) # 0.9928374
mean(moms) #  0.3650897

h.means=c(mean(bigH),mean(smallH),mean(moms),mean(inter))
h.lower=c(quantile(bigH,0.025),quantile(smallH,0.025),quantile(moms,0.025),quantile(inter,0.025))
h.upper=c(quantile(bigH,0.975),quantile(smallH, 0.975),quantile(moms, 0.975),quantile(inter, 0.975))
genetics=data.frame("variance.explained"=h.means,h.lower,h.upper)
genetics=cbind("component"=c("H","h2","maternal","non-additive"),genetics)

qplot(component,variance.explained, data=genetics,geom="bar",fill=component)+geom_errorbar(aes(ymin=h.lower,ymax=h.upper),lwd=0.4,width=0.3)+theme_bw()+scale_x_discrete(breaks=NULL)

############# FINAL MCMC modeling
summary(st)
library(MCMCglmm) 
st$pl=as.factor(st$pl)
#with plate and both times
mc=MCMCglmm(asin(sqrt(prop))~time,random=~idh(time):pl+idh(time):dad+idh(time):mom+idh(time):mom:dad,rcov=~idh(time):units,family="gaussian", data=st)
summary(mc)
mmm=mc
names(colMeans(mmm$VCV))

mom=c()
for(i in 1:length(mmm$VCV[,1])) {
	mom=append(mom,sum(mmm$VCV[i,c(5,6)])/sum(mmm$VCV[i,c(1:10)]))
}
quantile(mom,0.025) # small 4.371004e-13
quantile(mom,0.975) # 0.9795231
mean(mom) # 0.2528941

dad=c()
for(i in 1:length(mmm$VCV[,1])) {
	dad=append(dad,sum(mmm$VCV[i,c(3,4)])/sum(mmm$VCV[i,c(1:10)]))
}
quantile(dad,0.025) # 4.118331e-06
quantile(dad,0.975) # 0.3433702 
mean(dad) # 0.1037364

int=c()
for(i in 1:length(mmm$VCV[,1])) {
	int=append(int,sum(mmm$VCV[i,c(7,8)])/sum(mmm$VCV[i,c(1:10)]))
}
quantile(int,0.025) # 2.049027e-10 
quantile(int,0.975) # 0.2525508 
mean(int) # 0.04480941

bigH=c()
for(i in 1:length(mmm$VCV[,1])) {
	bigH=append(bigH,sum(mmm$VCV[i,2:7])/sum(mmm$VCV[i,1:10]))
}
quantile(bigH,0.025) #0.01318835
quantile(bigH,0.975) # 0.9832271 
mean(bigH) # 0.369348

#####Now with only 48 hours with plate USE THIS ONE
summary(st)
st$time = as.factor(st$time)
t48=subset(st , time=="48")
summary(t48)

mc2=MCMCglmm(asin(sqrt(prop))~1,random=~pl+dad+mom+mom:dad,family="gaussian", data=t48)
summary(mc2)
mmm=mc2
names(colMeans(mmm$VCV))

mom=c()
for(i in 1:length(mmm$VCV[,1])) {
	mom=append(mom,sum(mmm$VCV[i, 3])/sum(mmm$VCV[i,c(1:5)]))}
quantile(mom,0.025) # 1.010283e-13
quantile(mom,0.975) #0.9948508
mean(mom) #0.2818381

dad=c()
for(i in 1:length(mmm$VCV[,1])) {
	dad=append(dad,sum(mmm$VCV[i,2])/sum(mmm$VCV[i,c(1:5)]))}
quantile(dad,0.025) #1.401965e-14
quantile(dad,0.975) # 0.502977 
mean(dad) #0.1597487

int=c()
for(i in 1:length(mmm$VCV[,1])) {
	int=append(int,sum(mmm$VCV[i,4])/sum(mmm$VCV[i,c(1:5)]))}
quantile(int,0.025) # 2.182642e-15 
quantile(int,0.975) #0.395365  
mean(int) #  0.05186796

bigH=c()
for(i in 1:length(mmm$VCV[,1])) {
	bigH=append(bigH,sum(mmm$VCV[i,2:4])/sum(mmm$VCV[i,1:5]))}
quantile(bigH,0.025) #0.1419714
quantile(bigH,0.975) # 0.9958634 
mean(bigH) # 0.4934548
