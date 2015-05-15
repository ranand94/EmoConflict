#Emo Conflict Analysis
#Read in data
library(sas7bdat)
library(ez) 
library (ggplot2)
library(lme4)
library(plotrix)
library(plyr)
library(Hmisc)
library(multcomp)
library(reshape)
library(foreign)
source(lmerCellMeans.R)
#Read in data
#conflict<-read.table("C:/Users/ranan_000/Dropbox/Conflict Analysis/econflict_merged/data/all.csv", header=T, sep=",")
conflict<-read.table("C:/Users/anandr2/dropbox/Conflict Analysis/econflict_merged/data/all.csv", header=T, sep=",")
#conflict<-read.table("/Users/michael/Dropbox/Conflict Analysis/econflict_merged/data/all.csv", header=T, sep=",")

#head(conflict)
conflict<-within(conflict, {Subject<-factor(Subject)})
#str(conflict)
conflict<-within(conflict, {Cong<-factor(Cong,
                                         levels=c(0:1),
                                         labels=c("incong",
                                                  "cong"))})

conflict<-within(conflict, {Diag<-factor(BPD,
                                         levels=c(0:1),
                                         labels=c("control",
                                                  "BPD"))})
#Need a method to make emotion pairs. if face emo=1 and word emo=1-> hap-hap pair
conflict$EmoPair <- paste(conflict$FaceEmo,conflict$WordEmo, sep = '-') 
#convert char string to factor in EmoPair
conflict<-within(conflict, {EmoPair<-factor(EmoPair)})
#Next step is not necessary but makes plots better organized
#recategorize vars
conflict$EP2[conflict$EmoPair%in%c(' angry- angry')]<-1
conflict$EP2[conflict$EmoPair%in%c(' angry- fearful')]<-2
conflict$EP2[conflict$EmoPair%in%c(' fearful- fearful')]<-3
conflict$EP2[conflict$EmoPair%in%c(' fearful- happy')]<-4
conflict$EP2[conflict$EmoPair%in%c(' happy- happy')]<-5
conflict$EP2[conflict$EmoPair%in%c(' happy- fearful')]<-6
#Rename levels
conflict<-within(conflict,{EP2<-factor(EP2, levels=c(1:6), labels=c("Ang-Ang", "Ang-Fear", "Fear-Fear", "Fear-Hap", "Hap-Hap", "Hap-Fear"))})

#Vars of Interest= Subject, Trial, FaceEmo, WordEmo, Cong, ACC, RT, Diag, EmoPair or EP2

#Replace outlier values of subject trial means with NA in new col- RT_Trim
conflict <- ddply(conflict, .(Subject), function(subdf) {
  pct75 <- quantile(subdf$RT, probs=0.75, na.rm=T)
  iqr <- IQR(subdf$RT, na.rm=T)
  upperCut <- pct75 + 3*iqr
  subdf$RT_Trim <- subdf$RT
  subdf$RT_Trim[which(subdf$RT_Trim > upperCut)] <- NA
  return(subdf)
})
#NO NA's returned for outliers because the game has a very short RT cut-off
#a<-boxplot(conflict$RT_Trim): there is one 0 RT. need to drop
#subset data where RT_Trim does not =0
conflict2<-subset(conflict, RT_Trim != 0)
conflict<-subset(conflict2, ACC != 0)
rtcollapse <- ddply(conflict, c("Subject" , "Cong" , "EP2", "Diag" ), summarise, rtavg=mean(RT_Trim, na.rm=TRUE), rtsd=sd(RT_Trim, na.rm=TRUE))

# rtcollapse <- ddply(conflict, c("Subject", "EP2"), function(subdf) {
#   retdf <- data.frame(subdf$Subject[1], subdf$EmoPair[1], rtavg=mean(subdf$RT_Trim, na.rm=TRUE))
#   return(retdf)
# })



#b<-boxplot(conflict$RT_Trim)

#Group Descriptives Table
getDescriptivesByGroup <- function(df, outcomes, group, zscale=FALSE, tscale=FALSE) {
  require(plotrix)
  oneVar <- function(df, outcome) {
    stats <- c("mean", "sd", "std.error", "median", "min", "max")
    bigMatrix <- c()
    
    for (stat in stats) {
      if (zscale) {
        bigMatrix <- cbind(bigMatrix, tapply(as.vector(scale(df[[outcome]])), df[[group]], stat, na.rm=TRUE))  
      } else if (tscale) {
        scaleVec <- as.vector(scale(df[[outcome]]))
        scaleVec <- 50 + 10*scaleVec
        bigMatrix <- cbind(bigMatrix, tapply(scaleVec, df[[group]], stat, na.rm=TRUE))
      } else {
        bigMatrix <- cbind(bigMatrix, tapply(df[[outcome]], df[[group]], stat, na.rm=TRUE))
      }
      
      colnames(bigMatrix)[ncol(bigMatrix)] <- stat
    }
    
    bigMatrix <- as.data.frame(bigMatrix)
    bigMatrix$var <- outcome
    bigMatrix$group <- factor(rownames(bigMatrix))
    rownames(bigMatrix) <- NULL
    
    return(bigMatrix)    
  }
  
  allVars <- c()
  for (outcome in outcomes) {
    allVars <- rbind(allVars, oneVar(df, outcome))
  }
  allVars <- plyr::rename(allVars, c(group=group))
  return(allVars)
}
#end descriptives function
#type in vars of interest for descriptives table
desc<-getDescriptivesByGroup(conflict, c("ACC", "RT"), "Subject")
accdesc<-subset(desc, var=="ACC")
attach(accdesc)
boxplot(mean)
print(boxplot(mean))
#across subs, acc below 75 is considered outlier so drop scores below 75 (only 1 at 72)

accdesc<-subset(desc, var=="RT")
attach(accdesc)
boxplot(mean)
print(boxplot(mean))


#NOW plot observed values
#PLot observed values

#Rt Observed
g<-ggplot(conflict, aes(x=EmoPair, y=RT_Trim, color=FaceEmo, shape=Diag), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)
g<-ggplot(conflict, aes(x=FaceEmo, y=RT_Trim, color=EP2, shape=Diag), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)

h<-ggplot(conflict, aes(x=EmoPair, y=ACC, color=FaceEmo, shape=Diag), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)
h<-ggplot(conflict, aes(x=FaceEmo, y=ACC, color=EP2, shape=Diag), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)

#Plot Observed RT in pretty graph
j<-ggplot(rtcollapse, aes(x=EP2, y=rtavg, color=Diag, shape=Cong))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emo Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict RT Per Person")+xlab("Subject\n")+ylab("Mean RTs")
j+ggtitle("Emo Conflict RT Per Person")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("Cong")


#Or plot by TargetEmotion
i<-ggplot(conflict, aes(x=FaceEmo, y=RT_Trim, color=EP2, shape=Diag))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emo Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict RT")+xlab("Target Emo\n")+ylab("Mean RTs")
i+ggtitle("Emo Conflict RT")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("Diag")

#Plot Observed ACC in pretty graph
dd<-ggplot(conflict, aes(x=EmoPair, y=ACC, color=Diag))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Diag", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict Acc")+xlab("Emotion Pair\n")+ylab("Mean ACC (%)")
dd+ggtitle("Emo Conflict Acc")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))

#Or plot by target emo
dd<-ggplot(conflict, aes(x=FaceEmo, y=ACC, color=EmoPair, shape=Diag))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emo Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict Acc")+xlab("Target Emo\n")+ylab("Mean ACC (%)")
dd+ggtitle("Emo Conflict Acc")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("Diag")

#We can think of different ways to set up these plots but interesting that on incongruent and congruent negative
#words, BPD's perform better compared to happy target trials. effect should be stronger with more subjects hopefully

#Inidividual graphs by diag
rtcollapseBPD<-subset (rtcollapse, Diag=="BPD")
str(rtcollapseBPD)

pdf("subject_emopair_rts.bpd.pdf", width=10, height=8)

bysubj <- split(rtcollapseBPD, rtcollapseBPD$Subject)
for (subject in bysubj) {
  if(nrow(subject)==0L) next;
  subjgraph <- ggplot(subject, aes(x=EP2, y=rtavg, color=Cong, ymin=rtavg-rtsd, ymax=rtavg+rtsd)) + geom_pointrange(size=3) + ggtitle(paste("Subject", subject$Subject[1]))
  plot(subjgraph)
}

dev.off()

rtcollapseCon<-subset (rtcollapse, Diag=="control")

pdf("subject_emopair_rts.con.pdf", width=10, height=8)

bysubj <- split(rtcollapseCon, rtcollapseCon$Subject)
for (subject in bysubj) {
  if(nrow(subject)==0L) next;
  subjgraph <- ggplot(subject, aes(x=EP2, y=rtavg, color=Cong, ymin=rtavg-rtsd, ymax=rtavg+rtsd)) + geom_pointrange(size=3) + ggtitle(paste("Subject", subject$Subject[1]))
  plot(subjgraph)
}

dev.off()

#Average RT for each Diag by EP2
rtagg <- ddply(rtcollapse, c("EP2", "Diag"), summarise, rtavg=mean(rtavg, na.rm=TRUE), rtsd=sd(rtavg, na.rm=TRUE)) #function(subdf) { browser () })
               
#Individual Differences between Cong and Incong

rtindiv<-ddply(rtcollapse, c("Cong", "Diag", "Subject"), summarise, rtavg=mean(rtavg, na.rm=TRUE), rtsd=sd(rtavg, na.rm=TRUE))

Irt<-subset (rtindiv, Cong=="incong")
rtindivfinal<- subset (rtindiv, Cong=="cong")
names(Irt)[4] <- "IRT"

rtindivfinal$IRT<-Irt$IRT
names(rtindivfinal)[4] <- "CRT" 

rtindivfinal$RTdiff<- (rtindivfinal$IRT - rtindivfinal$CRT)
 #done with individual difference values


#congruency effect
incong<-as.vector(conflict$Cong)

priorincong <- Hmisc::Lag(incong, shift = 1)
conflict$priorincong<- priorincong
conflict<-within(conflict, {priorincong<-factor(priorincong,
                                         levels=c("cong", "incong"),
                                         labels=c(0:1))})
conflict[["priorincong"]][is.na(conflict[["priorincong"]])] <- 2
conflict$priorincong[is.NA(conflict$priorincong)]<-2
which(is.na(conflict$priorincong))
#basic linear modeling

rt.noint <- lmer(RT_Trim ~ FaceEmo + Cong + EP2 + (1 | Subject), conflict)
summary(rt.noint)
library(afex)
#mixed(RT_Trim ~ FaceEmo + Cong + EP2 + (1 | Subject), conflict)

#effect of run?
rt.noint.runfixed <- lmer(RT_Trim ~ FaceEmo + Cong + EP2 + Trials + (1|Subject), conflict)
summary(rt.noint.runfixed)
anova(rt.noint, rt.noint.runfixed)

rt.noint.runrand <- lmer(RT_Trim ~ FaceEmo + Cong + EP2 + (1 + Trials | Subject), conflict)
nlme_rt.noint.runrand <- nlme::lme(conflict.RT_Trim ~ FaceEmo + Cong+ EP2, random=~1|Subject, subData.all, na.action=na.exclude )
anova(nlme_rt.noint.runrand)
summary(rt.noint.runrand)

anova(rt.noint, rt.noint.runrand)

rt.noint.runboth<-lmer(RT_Trim~FaceEmo*aceEmo/EP2*Cong*priorincong+Diag+(1|Subject)+Trials+(1+Trials|Subject), data=conflict)
summary(rt.noint.runboth)

rt.int <- lmer(RT_Trim ~ FaceEmo*Cong + (1|Subject), conflict)
summary(rt.int)
print((pv <- pvals.fnc(rt.int))$fixed)

anova(rt.noint, rt.int)
anova(rt.noint, rt.noint.runfixed, rt.noint.runrand, rt.noint.runboth) 

# Conclusion: Including both fixed and random effects of trial fits best.

# Data: conflict
# Models:
# rt.noint: RT_Trim ~ FaceEmo + Cong + EP2 + (1 | Subject)
# rt.noint.runfixed: RT_Trim ~ FaceEmo + Cong + EP2 + Trials + (1 | Subject)
# rt.noint.runrand: RT_Trim ~ FaceEmo + Cong + EP2 + (1 + Trials | Subject)
# rt.noint.runboth: RT_Trim ~ FaceEmo + Cong + EP2 + Trials + (1 + Trials | Subject)
# Df   AIC   BIC logLik deviance   Chisq Chi Df Pr(>Chisq)    
# rt.noint           8 70433 70486 -35209    70417                              
# rt.noint.runfixed  9 70422 70482 -35202    70404 12.6847      1  0.0003687 ***
# rt.noint.runrand  10 70382 70449 -35181    70362 41.9896      1  9.176e-11 ***
# rt.noint.runboth  11 70380 70453 -35179    70358  4.0509      1  0.0441466 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#lmer cell means RT
rtcm <- lmerCellMeans(rt.noint.runboth)

# pdf("rtMeans.pdf", width=8, height=6) #units="in", res=300)
# ggplot(rtcm, aes(x=FaceEmo, y=RT_Trim, color=Cong)) + #geom_bar(stat="identity", position="dodge") +
#   geom_pointrange(aes(ymin=RT_Trim - se, ymax=RT_Trim + se), size=1.8, position=position_dodge()) +
#   theme_bw(base_size=22) + xlab("\nFacial Emotion") + ylab("Reaction Time (ms)\n") + ylim(400, 1000) + 
#   scale_color_brewer("Face-Word\nCongruence", palette="Set1") + theme(axis.title.x=element_text(vjust=0),legend.key.size=unit(2, "lines")) +
#   ggtitle("Fig 2. Reaction Times Across Conditions") #+ facet_wrap (~EP2)
# dev.off()

###this is where we were plotting rtcm
ggplot(rtcm, aes(x=EP2, y=RT_Trim, color=Cong, shape=Diag,
                 ymin=RT_Trim-se, ymax=RT_Trim+se))+
  geom_pointrange(size=1.0, position=position_dodge(width=.5))+
  scale_color_brewer("Face-Word\nCongruence", palette="Set2")+
  xlab("\nFacial Emotion") + ylab("Reaction Time (ms)\n") +
  theme(axis.title.x=element_text(vjust=0),legend.key.size=unit(2, "lines")) +
  ggtitle("Fig 2. Reaction Times Across Conditions")

a<-ggplot(rtcm, aes(x=EP2, y=RT_Trim, color=Diag, shape=priorincong))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emo Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict RT Model")+xlab("Subject\n")+ylab("Mean RTs")

##there is an effect of prior congruency
h<-ggplot(conflict, aes(x=priorincong, y=RT_Trim, color=Cong), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)

#ACC modeling
rt.nointa <- glmer(ACC ~ FaceEmo + Cong + EP2 + (1 | Subject), data=conflict2, family="binomial")
summary(rt.nointa)
library(afex)
#mixed(RT_Trim ~ FaceEmo + Cong + EP2 + (1 | Subject), conflict)

#effect of run?
rt.noint.runfixeda <- glmer(ACC ~ FaceEmo + Cong + EP2 + Trials + (1|Subject), data=conflict2, family="binomial")
summary(rt.noint.runfixeda)
anova (rt.nointa, rt.noint.runfixeda)

rt.noint.runranda <- glmer(ACC ~ FaceEmo + Cong + EP2 + (1 + Trials | Subject), data=conflict2, family="binomial")
nlme_rt.noint.runrand <- nlme::lme(conflict.RT_Trim ~ FaceEmo + Cong+ EP2, random=~1|Subject, subData.all, na.action=na.exclude )
summary(rt.noint.runranda)



rt.noint.runbotha <- glmer(ACC ~ FaceEmo + Cong + EP2 + Trials + (1 + Trials |Subject), data=conflict2, family="binomial")
summary(rt.noint.runbotha)

rt.inta <- glmer(ACC ~ FaceEmo*Cong + (1|Subject), data=conflict2, family="binomial")
summary(rt.inta)

anova (rt.noint.runfixeda, rt.noint.runbotha)
anova(rt.nointa, rt.noint.runfixeda, rt.noint.runranda, rt.noint.runbotha) 
 
#conclusion - runfixeda is best trials as just fixed  
# Models:
#   rt.nointa: ACC ~ FaceEmo + Cong + EP2 + (1 | Subject)
# rt.noint.runfixeda: ACC ~ FaceEmo + Cong + EP2 + Trials + (1 | Subject)
# rt.noint.runranda: ACC ~ FaceEmo + Cong + EP2 + (1 + Trials | Subject)
# rt.noint.runbotha: ACC ~ FaceEmo + Cong + EP2 + Trials + (1 + Trials | Subject)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# rt.nointa           7 3965.7 4013.0 -1975.8   3951.7                         
# rt.noint.runfixeda  8 3927.6 3981.7 -1955.8   3911.6 40.059      1  2.464e-10
# rt.noint.runranda   9 3955.8 4016.7 -1968.9   3937.8  0.000      1          1
# rt.noint.runbotha  10 3928.9 3996.5 -1954.5   3908.9 28.891      1  7.655e-08
# 
# rt.nointa             
# rt.noint.runfixeda ***
#   rt.noint.runranda     
# rt.noint.runbotha  ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1