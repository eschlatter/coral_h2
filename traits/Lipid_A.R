setwd("/Users/Ted/Desktop/RNAseq_Larvae/Orpheus/Orpheus/lipid")
lipid=read.csv("Lipid_2011_AC.csv")
attach(lipid)
summary(lipid)
lipid$time=as.factor(lipid$time)
summary(lipid)
lipid$loss=((lipid$lipid-lipid$egg)/(lipid$egg))
hist(lipid$loss)

aov1=aov(abs(loss)~dad, data=lipid)
summary(aov1)

library(ggplot2)

#Averaged change over time
qplot(dad, loss, data = lipid, geom = "boxplot", xlab = "", ylab = "Lipid Loss") +
theme_bw() +
opts(axis.title.x =theme_text(size = 14), axis.title.y = theme_text(size= 14, angle = 90, vjust = 0.4), axis.text.x = theme_text(size = 12), axis.text.y = theme_text(size = 12))

qplot(mom, loss, data = lipid, geom = "boxplot", xlab = "", ylab = "Lipid Loss") +
theme_bw() +
opts(axis.title.x =theme_text(size = 14), axis.title.y = theme_text(size= 14, angle = 90, vjust = 0.4), axis.text.x = theme_text(size = 12), axis.text.y = theme_text(size = 12))

# Testing for heritibility
summary(lipid)
library(lme4) 
m0=lmer(loss))~(1|mom), lipid)
m1=lmer(loss))~(1|mom)+(1|dad),lipid)
m2=lmer(loss))~(1|mom)*(1|dad),lipid)
anova(m0,m1,m2)
summary(m1)

lm.va=VarCorr(m1)$dad[1]
vres= as.numeric ( attributes ( summary(m1))$REmat [ attributes (summary(m1)) $REmat [,1] == "Residual" , 3 ] ) 
h2lmer=2*lm.va/(lm.va+vres)  
> h2lmer 
[1] 0.4253725
par(mfrow=c(2,2))
plot(m1)












#change through time
g=ggplot(st,aes(time,(log(lipid.loss+100),group=1,colour=dad))
g+geom_smooth(aes(group=dad))+theme_bw()
g+stat_summary(aes(group=dad),fun.data="median_hilow",conf.int=0.5,geom="smooth")

#Lipid change at three days
s3=subset(st, st$time=="3")
summary(s3)
aov2=aov(abs(lipid.loss)~dad, data=s3)
summary(aov2)


#Lipid change at five days
s5=subset(lipid, lipid$time=="5")
summary(s5)
aov3=aov(abs(lipid.loss)~dad, data=s5)
summary(aov3)
library(ggplot2)
qplot(dad, abs(lipid.loss/100), data = s5, geom = "boxplot", xlab = "", ylab = "5-day Lipid Loss") +
theme_bw() +
opts(axis.title.x =theme_text(size = 14), axis.title.y = theme_text(size= 14, angle = 90, vjust = 0.4), axis.text.x = theme_text(size = 12), axis.text.y = theme_text(size = 12))

#regressing settlement and lipid
set=read.csv("Settlement_Orpheus48.csv")
set$rep=as.factor(set$rep)
summary(set)
set$prop=(set$rec/(set$rec+set$lar))
library(ggplot2)
source("summarySE.R")
all.marg=summarySE(lipid,measurevar="loss",groupvars=c("dad"), na.rm=T)
all.marg
all.marg1=summarySE(set,measurevar="prop",groupvars=c("dad"), na.rm=T)
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



#regressing protein and lipid
st=read.csv("Protein_Orpheus.csv")
attach(st)
summary(st)
st$rep=as.factor(st$rep)
st$loss=st$wgg-st$day3
summary(st)

source("summarySE.R")
all.marg=summarySE(lipid,measurevar="loss",groupvars=c("dad"), na.rm=T)
all.marg
all.marg1=summarySE(st,measurevar="loss",groupvars=c("dad"), na.rm=T)
all.marg1
plot1=lm(all.marg$loss~all.marg1$loss)
summary(plot1)
plot(log(abs(all.marg$loss))~log(abs(all.marg1$loss)), xlab="Protein Loss", ylab="Settlement", cex=1.5, col="blue", cex.axis=1, cex.lab=1.2)
abline(plot1 ,col="red")

pd=position_dodge(.3)
ggplot(all.marg,aes(x=dad,y=loss,colour=mom,group=mom))+
 	geom_errorbar(aes(ymin=loss+se,ymax=loss-se),lwd=0.4,width=0.3,position=pd)+
	geom_line(aes(group=mom,linecol=mom),position=pd)+
	geom_point(aes(group=mom,pch=mom),position=pd,size=2.5)+
	theme_bw()



	
#############New in July
setwd("/Users/Ted/Desktop/RNAseq_Larvae/Orpheus/Orpheus/lipid")
lipid=read.csv("Lipid_2011_AC.csv")
attach(lipid)
summary(lipid)
lipid$time=as.factor(lipid$time)
summary(lipid)
lipid$loss=((lipid$egg-lipid$lipid)/(lipid$egg))
hist(lipid$loss)
str(lipid)
all.margrg=summarySE(lipid,measurevar="loss",groupvars=c("mom", "dad"), na.rm=T)
all.margrg
pd=position_dodge(.3)
ggplot(all.margrg,aes(x=dad,y=loss,colour=mom,group=mom))+
 	geom_errorbar(aes(ymin=loss+se,ymax=loss-se), lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=mom,pch=mom),position=pd,size=2.5)+
	xlab("Sire")+
	ylab("Lipid Loss (mg)")+
	theme_bw()	
	
# REML modeling with proportions, sire effects for Mom A at 48hrs
library(lme4) 
lipid$sample=as.factor(lipid$sample)
m2=lmer(loss~(1|dad),lipid)
m3=lmer(loss~(1|mom),lipid)
m4=lmer(loss~(1|dad)+(1|mom),lipid)
m5=lmer(loss~(1|dad)*(1|mom), lipid)
anova(m2,m3,m4,m5)

lm.va=VarCorr(m5)$dad[1]
vres= as.numeric ( attributes ( summary(m6))$REmat [ attributes (summary(m5)) $REmat [,1] == "Residual" , 3 ] ) 
h2lmer=2*lm.va/(lm.va+vres)  
h2lmer # 0.597

#------------------------------
# MCMC modeling
summary(st)
library(MCMCglmm) 
hist(log(st$loss+1))

mc=MCMCglmm(loss~1,random=~dad+mom+mom:dad,family="gaussian",data=lipid)

summary(mc)
mmm=mc
names(colMeans(mmm$VCV))

bigH=c()
 for(i in 1:length(mmm$VCV[,1])) {
	bigH=append(bigH,sum(mmm$VCV[i,1:3])/sum(mmm$VCV[i,1:4]))
 }
quantile(bigH,0.025) #     2.5% 
# 5.620227e-11
    
quantile(bigH,0.975) #     97.5% 
# 0.4476038   

mean(bigH) #0.0443763




dads=c()
for(i in 1:length(mmm$VCV[,1])) {
	dads=append(dads,sum(mmm$VCV[i,1])/sum(mmm$VCV[i,1:4]))
}
quantile(dads,0.025) # super small
quantile(dads,0.975) # 0.46
mean(dads) # 0.08

moms=c()
for(i in 1:length(mmm$VCV[,1])) {
	moms=append(moms,(sum(mmm$VCV[i,2]))/sum(mmm$VCV[i,1:4]))
}
quantile(moms,0.025) # very small
quantile(moms,0.975) # 0.44
mean(moms) # 0.03

int=c()
for(i in 1:length(mmm$VCV[,1])) {
	int=append(int,sum(mmm$VCV[i,3])/sum(mmm$VCV[i,1:4]))
}
quantile(int,0.025) # really small
quantile(int,0.975) # 0.44
mean(int) # 0.04

h.means=c(mean(moms),mean(dads),mean(int))
h.lower=c(quantile(moms,0.025),quantile(dads,0.025),quantile(int,0.025))
h.upper=c(quantile(moms,0.975),quantile(dads, 0.975),quantile(int, 0.975))
genetics=data.frame("variance.explained"=h.means,h.lower,h.upper)
genetics=cbind("component"=c("Dam","sire","interaction"),genetics)

qplot(component,variance.explained, data=genetics,geom="bar",fill=component)+geom_errorbar(aes(ymin=h.lower,ymax=h.upper),lwd=0.4,width=0.3)+theme_bw()+scale_x_discrete(breaks=NULL)



##########getting data out for WGCNA Traits
setwd("/Users/Ted/Desktop/Orpheus/Orpheus/lipid")
lipid=read.csv("Lipid_2011_AC.csv")
attach(lipid)
summary(lipid)
lipid$time=as.factor(lipid$time)
summary(lipid)
lipid$loss=((lipid$egg-lipid$lipid)/(lipid$egg))
hist(lipid$loss)
str(lipid)
all.margrg=summarySE(lipid,measurevar="loss",groupvars=c("mom", "dad"), na.rm=T)
all.margrg
pd=position_dodge(.3)
ggplot(all.margrg,aes(x=dad,y=loss,colour=mom,group=mom))+
 	geom_errorbar(aes(ymin=loss+se,ymax=loss-se), lwd=0.4,width=0.3,position=pd)+
	geom_point(aes(group=mom,pch=mom),position=pd,size=2.5)+
	xlab("Sire")+
	ylab("Lipid Loss (mg)")+
	theme_bw()	
	
#getting average for WGCNA
all.margl=summarySE(lipid,measurevar="loss",groupvars=c("sample"), na.rm=T)
all.margl
write.csv(all.margl, file = "lipidloss_av_new.csv")

