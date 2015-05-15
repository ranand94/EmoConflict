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

i<-ggplot(conflict, aes(x=FaceEmo, y=RT_Trim, color=EP2, shape=Diag))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emo Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict RT")+xlab("Target Emo\n")+ylab("Mean RTs")+ facet_wrap(~priorincong)
i+ggtitle("Emo Conflict RT")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("Diag") 
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

#congruency effect
incong<-as.vector(conflict$Cong)

priorincong <- Hmisc::Lag(incong, shift = 1)
conflict$priorincong<- priorincong
conflict<-within(conflict, {priorincong<-factor(priorincong,
                                                levels=c("cong", "incong"),
                                                labels=c(0:1))})

#basic linear modeling RT
rt.noint.runboth<-lmer(RT_Trim~FaceEmo*FaceEmo/EP2*Cong/priorincong+Diag+(1|Subject)+Trials+(1+Trials|Subject), data=conflict)
summary(rt.noint.runboth)

cm <- lmerCellMeans(rt.noint.runboth)
#check about trials**********

#TO SUBSET OUT NONSENSE VARS, merge emo and emo pair and remove nonsense

cm$EmoEmoPair <- paste(cm[,1], cm[,5], sep=" ")
cm<-subset(cm, !(EmoEmoPair %in% c(" fearful Ang-Ang", " happy Ang-Ang", " fearful Ang-Ang",
                                   " happy Fear-Fear" , " fearful Ang-Fear", " happy Ang-Fear",
                                   " angry Fear-Fear" , " angry Hap-Hap", " fearful Hap-Hap",
                                   " angry Hap-Hap", " angry Hap-Fear", " angry Fear-Hap")))
cm$CongEmoPair <- paste(cm[,2], cm[,5], sep=" ")
cm<-subset(cm, !(CongEmoPair %in% c("incong Ang-Ang", "cong Ang-Fear",  "incong Fear-Fear",
                                    "cong Fear-Hap" ,  "incong Hap-Hap" , "cong Hap-Fear")))


ggplot(cm, aes(x=FaceEmo, y=RT_Trim, color=EP2, shape=Diag))+
  geom_pointrange(aes(ymin=RT_Trim - se, ymax=RT_Trim + se))+
  scale_color_brewer("Emo Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict RT Model")+xlab("Subject\n")+ylab("Mean RTs")

d <- ggplot(cm, aes(x=FaceEmo, y=RT_Trim,
                    color=EP2, shape=Diag,
                    ymin=RT_Trim-se, ymax=RT_Trim+se)) +
  geom_point( size=1,position=position_dodge(width=.5)) +
  scale_color_brewer("Emotion Pair", palette="Set2")+
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())

ggplot(cm, aes(x=FaceEmo, y=RT_Trim, color=EP2, shape=Diag, ymin=RT_Trim-se, ymax=RT_Trim+se))+
geom_pointrange(size = 1, position=position_dodge(width=.5))+
  scale_color_brewer("Emo Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict RT")+xlab("Target Emo\n")+ylab("Mean RTs")+ facet_wrap(~priorincong)+
ggtitle("Emo Conflict RT")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("Diag") 





#linear model ACC
rt.noint.runbotha <- glmer(ACC ~ FaceEmo*FaceEmo/EP2*Cong/priorincong+Diag+(1|Subject)+Trials+(1+Trials|Subject), data=conflict, family="binomial")
summary(rt.noint.runbotha)

rtcmacc <- lmerCellMeans(rt.noint.runbotha)

b<-ggplot(rtcmacc, aes(x=EP2, y=RT_Trim, color=Diag, shape=priorincong))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emo Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Conflict RT Model")+xlab("Subject\n")+ylab("Mean RTs")
