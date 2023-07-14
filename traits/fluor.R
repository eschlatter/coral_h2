setwd("/Users/Ted/Desktop/RNAseq_Larvae/Orpheus/Orpheus/fluor")
st=read.csv("fluor_final_nosize.csv")
summary(st)
library(ggplot2)
source("summarySE.R")
st$dam=as.factor(st$dam)
st$sire=as.factor(st$sire)


#looking only at relative red
st$relred=st$red/(st$geren+st$red+st$blue)
hist(asin(sqrt(st$relred)))
red=aov(asin(sqrt(st$relred))~sire, data=st)
summary(red)
TukeyHSD(red)
plot(red~sire, data=st)

all.margr=summarySE(st,measurevar="relred",groupvars=c("dam", "sire"), na.rm=T)
all.margr
pd=position_dodge(.3)
ggplot(all.margr,aes(x=sire,y=relred,colour=dam,group=dam))+
 	geom_errorbar(aes(ymin=relred+se,ymax=relred-se), lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=dam,pch=dam),position=pd,size=2.5)+
	theme_bw()

#looking only at relative green
st$relgreen=st$geren/(st$geren+st$red+st$blue)
hist(asin(sqrt(st$relgreen)))
green=aov(asin(sqrt(st$relgreen))~sire, data=st)
summary(green)
TukeyHSD(green)
plot(relgreen~sire, data=st)

all.margg=summarySE(st,measurevar="relgreen",groupvars=c("dam", "sire"), na.rm=T)
all.margg
pd=position_dodge(.3)
ggplot(all.margg,aes(x=sire,y=relgreen,colour=dam,group=dam))+
 	geom_errorbar(aes(ymin=relgreen+se,ymax=relgreen-se), lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=dam,pch=dam),position=pd,size=2.5)+
	theme_bw()
	
#looking only at relative blue
st$relblue=st$blue/(st$geren+st$red+st$blue)
hist(asin(sqrt(st$relblue)))
blue=aov(asin(sqrt(st$relblue))~sire, data=st)
summary(blue)
TukeyHSD(blue)
plot(relblue~sire, data=st)

all.margb=summarySE(st,measurevar="relblue",groupvars=c("dam", "sire"), na.rm=T)
all.margb
pd=position_dodge(.3)
ggplot(all.margb,aes(x=sire,y=relblue,colour=dam,group=dam))+
 	geom_errorbar(aes(ymin=relblue+se,ymax=relblue-se), lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=dam,pch=dam),position=pd,size=2.5)+
	theme_bw()
	
#What about total fluorescence?
st$all=st$blue+st$red+st$geren
head(st)
hist(st$all)
allf=aov(st$all~sire, data=st)
summary(allf)
TukeyHSD(allf)
plot(all~sire, data=st)

all.marga=summarySE(st,measurevar="all",groupvars=c("dam", "sire"), na.rm=T)
all.marga
pd=position_dodge(.3)
ggplot(all.marga,aes(x=sire,y=all,colour=dam,group=dam))+
 	geom_errorbar(aes(ymin=all+se,ymax=all-se), lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=dam,pch=dam),position=pd,size=2.5)+
	theme_bw()
	
####Final plot for poster
all.margr=summarySE(st,measurevar="relred",groupvars=c("dam", "sire"), na.rm=T)
all.margr
pd=position_dodge(.3)
ggplot(all.margr,aes(x=sire,y=relred,colour=dam,group=dam))+
 	geom_errorbar(aes(ymin=relred+se,ymax=relred-se), lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=dam,pch=dam),position=pd,size=2.5)+
	xlab("Sire")+
	ylab("Relative Red Fluorescence")+
	theme_bw()

# REML modeling with proportions, sire effects for Mom A at 48hrs
library(lme4) 
st$culture=as.factor(st$culture)
#m0=lmer(log(rg)~1,st)
m1=lmer(asin(sqrt(relred))~1+(1|sire),st)
m3=lmer(asin(sqrt(relred))~1+(1|dam),st)
m4=lmer(asin(sqrt(relred))~1+(1|sire)+(1|dam),st)
m5=lmer(asin(sqrt(relred))~1+(1|sire)*(1|dam),st)
m6=lmer(asin(sqrt(relred))~1+(1|sire)*(1|dam)+(1|culture), st)
anova(m1,m3,m4,m5,m6)
summary(m6)
# Data: st
# Models:
# m1: log(rg) ~ 1 + (1 | sire)
# m3: log(rg) ~ 1 + (1 | dam)
# m4: log(rg) ~ 1 + (1 | sire) + (1 | dam)
# m5: log(rg) ~ 1 + (1 | sire) * (1 | dam)
# m6: log(rg) ~ 1 + (1 | sire) * (1 | dam) + culture
   # Df     AIC     BIC logLik  Chisq Chi Df Pr(>Chisq)    
# m1  3 -1276.9 -1263.9 641.43                             
# m3  3 -1423.1 -1410.2 714.57 146.28      0     <2e-16 ***
# m4  4 -1590.1 -1572.9 799.05 168.96      1     <2e-16 ***
# m5  4 -1590.1 -1572.9 799.05   0.00      0          1    
# m6 40 -1691.5 -1519.3 885.74 173.38     36     <2e-16 ***    
#Significant effect of dad
#Now plot the dad to see
qplot(sire, asin(sqrt(relred)), data=st)
qplot(sire, asin(sqrt(relred)), data =st, geom = "boxplot", xlab = "Sire", ylab = "Red/Green+Red") +
opts(axis.text.x = theme_text(angle = 90, hjust = 1, size = 8)) +
theme_bw()

lm.va=VarCorr(m6)$sire[1]
vres= as.numeric ( attributes ( summary(m6))$REmat [ attributes (summary(m6)) $REmat [,1] == "Residual" , 3 ] ) 
h2lmer=2*lm.va/(lm.va+vres)  
h2lmer # 0.282

#------------------------------
# MCMC modeling with red/green, sire effects
summary(st)
names(st)
head(st)
library(MCMCglmm) 
mpr=MCMCglmm(relred~1, random=~dam+sire+dam:sire+culture,data=st, family="gaussian")
summary(mpr)
mmm=mpr
names(colMeans(mmm$VCV))

bigH=c()
for(i in 1:length(mmm$VCV[,1])) {
	bigH=append(bigH,sum(mmm$VCV[i,1:3])/sum(mmm$VCV[i,1:5]))
}
quantile(bigH,0.025) # 0.371
quantile(bigH,0.975) # 0.999
mean(bigH) # 0.760

# smallH=c()
# for(i in 1:length(mmm$VCV[,1])) {
	# smallH=append(smallH,2*sum(mmm$VCV[i,2])/sum(mmm$VCV[i,1:5]))
# }
 # quantile(smallH,0.025) # really small
 # quantile(smallH,0.975) # 0.635
 # mean(smallH) # 0.162

moms=c()
for(i in 1:length(mmm$VCV[,1])) {
	moms=append(moms,(mmm$VCV[i,1])/sum(mmm$VCV[i,1:5]))}
quantile(moms,0.025) # 0.12
quantile(moms,0.975) # .999
mean(moms) # 0.649

dads=c()
for(i in 1:length(mmm$VCV[,1])) {
	dads=append(dads,(mmm$VCV[i,2])/sum(mmm$VCV[i,1:5]))}
quantile(dads,0.025) # very small
quantile(dads,0.975) # 0.39
mean(dads) # 0.09

int=c()
for(i in 1:length(mmm$VCV[,1])) {
	int=append(int,(mmm$VCV[i,3])/sum(mmm$VCV[i,1:5]))}
quantile(int,0.025) # very small
quantile(int,0.975) # 0.265
mean(int) # 0.027

h.means=c(mean(bigH), mean(moms),mean(dads),mean(int))
h.lower=c(quantile(bigH,0.025), quantile(moms,0.025),quantile(dads,0.025),quantile(int,0.025))
h.upper=c(quantile(bigH,0.975), quantile(moms,0.975),quantile(dads, 0.975),quantile(int, 0.975))
genetics=data.frame("variance.explained"=h.means,h.lower,h.upper)
genetics=cbind("component"=c("BigH","Dam","Sire","Interaction"),genetics)

qplot(component,variance.explained, data=genetics,geom="bar",fill=component)+geom_errorbar(aes(ymin=h.lower,ymax=h.upper),lwd=0.4,width=0.3)+theme_bw()+scale_x_discrete(breaks=NULL)



#regressing green to red
summary(st)
all.marg2=summarySE(st,measurevar="geren",groupvars=c("culture"), na.rm=T)
all.marg2
all.marg3=summarySE(st,measurevar="red",groupvars=c("culture"), na.rm=T)
all.marg3
plot1=lm(all.marg2$geren~all.marg3$red)
summary(plot1)
plot(all.marg2$geren~all.marg3$red, xlab="rel green", ylab="rel red", cex=1.5, col="blue", cex.axis=1, cex.lab=1.2)
abline(plot1 ,col="red")

#getting average for WGCNA
all.marg=summarySE(st,measurevar="relred",groupvars=c("culture"), na.rm=T)
all.marg
write.csv(all.marg, file = "relred_av_new.csv")
all.marg1=summarySE(st,measurevar="relgreen",groupvars=c("culture"), na.rm=T)
all.marg1
write.csv(all.marg1, file = "relgreen_av_new.csv")
all.marg2=summarySE(st,measurevar="relblue",groupvars=c("culture"), na.rm=T)
all.marg2
write.csv(all.marg2, file = "relblue_av_new.csv")


