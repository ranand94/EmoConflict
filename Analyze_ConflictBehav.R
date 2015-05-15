setwd(file.path(getMainDir(), "K01_Study", "conflictTask"))
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
load(file="emotConflictBehav6Jun2012.RData")
library(foreign)
library(lme4)
library(ggplot2)
library(languageR)
library(grid)

library(plyr)
#pull together a long dataset (all subjects, all runs)
subData.all <- c()

for (s in 1:length(subData)) {
  sub <- subData[[s]]
  subData.all <- c(subData.all, sub)
}

subData.all <- do.call(rbind.fill, subData.all)

subData.all$subid <- factor(subData.all$subid)
subData.all$runnum <- factor(subData.all$runnum)
subData.all$congruent <- factor(subData.all$congruent, levels=c(0,1), labels=c("Incongruent", "Congruent"))
subData.all$faceemot <- factor(subData.all$faceemot, levels=c("angry", "fearful", "happy"), labels=c("Angry", "Fearful", "Happy"))

#drop unused cols
subData.all <- subset(subData.all, select=-c(Running.Trial., Procedure.Trial., test_img, img, iti, emoSlide.DurationError, 
        emoSlide.OnsetTime, emoSlide.OnsetTime.Abs, emoSlide.RESP, emoSlide.CRESP, emoSlide.RTTime, emoSlide.RTTime.Abs))

#data cleaning
#dev.new()
#ggplot(subData.all, aes(x=emoSlide.RT)) + geom_histogram() + facet_grid(subid ~ runnum) + ggtitle("Histograms for reaction times by subject and run")

#recode missed responses (0 RT) as NA
subData.all$emoSlide.RT <- sapply(subData.all$emoSlide.RT, function(x) ifelse(x > 0, x, NA))

#recode RTs < 150ms as missing (NA) -- quantile < 5%
quantile(subData.all$emoSlide.RT, c(.00, .05, .10, .15), na.rm=T)
subData.all[which(subData.all$emoSlide.RT < 150),c("emoSlide.ACC", "emoSlide.RT")] <- NA

#for now, recode RTs > 2000ms as inaccurate (NA) -- quantile > 99.5 
#quantile(subData.all$emoSlide.RT, c(.90, .95, .99, .995, .999, 1), na.rm=T)
#subData.all[which(subData.all$emoSlide.RT > 2000),c("emoSlide.ACC", "emoSlide.RT")] <- NA

#remove RTs for trials where person's response was incorrect (ANT_ACC = 0)
subData.all[which(subData.all$emoSlide.ACC == 0), "emoSlide.RT"] <- NA

#descriptive accuracy statistics per subject and run
with(subData.all, table(runnum, emoSlide.ACC, subid)) #counts of correct trials
with(subData.all, prop.table(table(runnum, emoSlide.ACC, subid), margin=c(1,3))) #same as proportion correct per run

#atrocious performance by 005_sb
#also, head movement was quite bad
#drop
#, , subid = 005_sb
#
#emoSlide.ACC
#runnum          
#           0          1
#1 0.31944444 0.68055556
#2 0.43055556 0.56944444
#3 0.48611111 0.51388889
#4 0.59722222 0.40277778

subData.all <- gdata::drop.levels(subset(subData.all, subid!="005_sb"))

#basic lmer
options(width=130)

rt.noint <- lmer(emoSlide.RT ~ faceemot + congruent + (1|subid), subData.all)
summary(rt.noint)
print((pv <- pvals.fnc(rt.noint))$fixed)

#effect of run?
rt.noint.runfixed <- lmer(emoSlide.RT ~ faceemot + congruent + runnum + (1|subid), subData.all)
summary(rt.noint.runfixed)

rt.noint.runrand <- lmer(emoSlide.RT ~ faceemot + congruent + (1 + runnum |subid), subData.all)
nlme_rt.noint.runrand <- nlme::lme(emoSlide.RT ~ faceemot * congruent, random=~1|subid, subData.all, na.action=na.exclude )
anova(nlme_rt.noint.runrand)
summary(rt.noint.runrand)

rt.noint.runboth <- lmer(emoSlide.RT ~ faceemot + congruent + runnum + (1 + runnum |subid), subData.all)
summary(rt.noint.runboth)

rt.int <- lmer(emoSlide.RT ~ faceemot*congruent + (1|subid), subData.all)
summary(rt.int)
print((pv <- pvals.fnc(rt.int))$fixed)

anova(rt.noint, rt.int)
anova(rt.noint, rt.noint.runfixed, rt.noint.runrand, rt.noint.runboth) 

#conclusion: accept the model with unique random intercepts per run
# looks like performance was much more variable on run 1
# and in general, people were faster in runs 2-4.
cm <- lmerCellMeans(rt.noint.runrand)

pdf("rtMeans.pdf", width=8, height=6) #units="in", res=300)
ggplot(cm, aes(x=faceemot, y=emoSlide.RT, color=congruent)) + #geom_bar(stat="identity", position="dodge") +
    geom_pointrange(aes(ymin=emoSlide.RT - se, ymax=emoSlide.RT + se), size=1.8, position=position_dodge()) +
    theme_bw(base_size=22) + xlab("\nFacial Emotion") + ylab("Reaction Time (ms)\n") + ylim(700, 1000) + 
    scale_color_brewer("Face-Word\nCongruence", palette="Set1") + theme(axis.title.x=element_text(vjust=0),legend.key.size=unit(2, "lines")) +
    ggtitle("Fig 2. Reaction Times Across Conditions")
dev.off()

library(multcomp)

#Pairwise contrasts for emotion 
summary(glht(rt.noint.runrand, linfct = mcp(faceemot = "Tukey")))
summary(glht(rt.noint.runrand, linfct = mcp(faceemot = "Tukey")), test=Chisqtest())
summary(glht(rt.noint.runrand, linfct = mcp(congruent = "Tukey")), test=Chisqtest())

#conclusion: no congruence x emotion interaction
# ME of emotion: Angry > Fearful > Happy
# ME of congruence: incongruent slowed RT in general


#accuracy analysis
acc.noint <- glmer(emoSlide.ACC ~ faceemot + congruent + (1|subid), subData.all, family=binomial)
summary(acc.noint)

acc.int <- glmer(emoSlide.ACC ~ faceemot * congruent + (1|subid), subData.all, family=binomial)

#no evidence of interaction
anova(acc.noint, acc.int)

#effect of run?
acc.noint.runfixed <- glmer(emoSlide.ACC ~ faceemot + congruent + runnum + (1|subid), subData.all, family=binomial)
acc.noint.runrand <- glmer(emoSlide.ACC ~ faceemot + congruent + (1 + runnum |subid), subData.all, family=binomial)
acc.noint.runboth <- glmer(emoSlide.ACC ~ faceemot + congruent + runnum + (1 + runnum |subid), subData.all, family=binomial)

#no evidence that modeling run improves fit
anova(acc.noint, acc.noint.runfixed, acc.noint.runrand, acc.noint.runboth)

cm <- lmerCellMeans(acc.noint)

cm$prob.emoSlide.ACC <- gtools::inv.logit(cm$emoSlide.ACC)
cm$sehigh <- gtools::inv.logit(cm$emoSlide.ACC+cm$se)
cm$selow <- gtools::inv.logit(cm$emoSlide.ACC-cm$se)

pdf("accMeans.pdf", width=8, height=6)
ggplot(cm, aes(x=faceemot, y=prob.emoSlide.ACC, color=congruent, ymin=selow, ymax=sehigh)) + 
    geom_pointrange(size=1.7, position=position_dodge(width=0.2)) + theme_bw(base_size=22) +
    ylim(0.5, 1) + scale_color_brewer("Face-Word\nCongruence", palette="Set1") + theme(axis.title.x=element_text(vjust=0),legend.key.size=unit(2, "lines")) +
    ylab("Probability of Correct Response\n") + xlab("\nFacial Emotion") + 
    ggtitle("Fig 4. Accuracy Across Conditions")#geom_errorbar(position=position_dodge(width=1), width=0.4) + 
dev.off()


summary(glht(acc.noint, linfct = mcp(faceemot = "Tukey")))
summary(glht(acc.noint, linfct = mcp(faceemot = "Tukey")), test=Chisqtest())

#now look at how behavior is influenced by psychological variables

#psych data import
gpsbpd <- read.spss("PGS fMRI Data.sav", to.data.frame=TRUE)
#dump weird padding
gpsbpd$PTNUM <- factor(as.numeric(as.character(gpsbpd$PTNUM)))
myids <- data.frame(PTNUM=factor(c(
        50199, 50413, 50498, 50589, 50696, 60456, 60718, 60743, 70778, 70812)),
  subid=factor(c("009_cm", "003_jd", "005_sb", "004_dm", "010_nf", "001_yc", "002_zz", "006_jn", "007_at", "008_ss")))

gpsbpd <- merge(gpsbpd, myids, by="PTNUM")

gpsbpd <- gpsbpd[order(gpsbpd$subid),]
gpsbpd$bpddiagScore <- apply(gpsbpd[,paste("BORDL", 1:9, sep="")], 1, function(row) {
      numcode <- sapply(row, function(x) {
            if (x=="Not present") 0
            else if (x=="Subthreshold") 1
            else if (x=="Present") 2
            else if (x=="Strongly Present") 3
            else stop("cannot match:", x)
          })

      sum(numcode)
    })

#depressed mood, and anhedonia... leave of irritability/anger for now 
gpsbpd$depanh <- apply(gpsbpd[,c("DEPR1", "DEPR3")], 1, function(row) {
      numcode <- sapply(row, function(x) {
            if (x=="Not present") 0
            else if (x=="Subthreshold") 1
            else if (x=="Threshold") 2
            else if (x=="Strongly Present") 3
            else stop("cannot match:", x)
          })
      sum(numcode)
    })

gpsbpd$anh <- sapply(gpsbpd$DEPR3, function(x) {
      if (x=="Not present") 0
      else 1
    })

#it ain't pretty, but 001_yc is missing the PAI scores from the initial visit
#use the PAI AI score from baseline (screen)
gpsbpd[which(gpsbpd$subid == "001_yc"), "paiaffe"] <- gpsbpd[which(gpsbpd$subid == "001_yc"), "PAIscoreBS"] 

#even uglier, just take my best guess at the self-harm score based on other clinical information
#score of 9 is comparable to others at this level of overall BPD symptomatology and given that 001_yc endorsed self-harm
gpsbpd[which(gpsbpd$subid == "001_yc"), "paiharm"] <- 9

gpsbpd[,c("PTNUM", "subid", "Age", "PAIscoreBS", "PAIscoreW1", "paiaffe", "paiharm", "bpddiagScore", "depanh")]  #, paste("BORDL", 1:9, sep=""))]

cor(gpsbpd[,c("Age", "paiaffe", "paiharm", "bpddiagScore", "depanh")])  #, paste("BORDL", 1:9, sep=""))]

subData.all <- merge(subData.all, gpsbpd[,c("subid", "Age", "paiaffe", "paiharm", "bpddiagScore", "depanh")], by="subid")
subData.all$bpddiagScore <- as.vector(scale(subData.all$bpddiagScore, center=TRUE, scale=FALSE)) 

#add bpd diagnostic score to RT analysis
#pulling out the extra random effects per run here to reduce parameters (stabilize convergence)
rt.noint.cov <- lmer(emoSlide.RT ~ faceemot*congruent*bpddiagScore + (1|subid), subData.all) #faceemot + congruent + faceemot*bpddiagScore + congruent*bpddiagScore + (1 |subid), subData.all)
summary(rt.noint.cov)
print((pv <- pvals.fnc(rt.noint.cov))$fixed)

cm <- lmerCellMeans(rt.noint.cov)

#evidence of a BPD x Emotion interaction
anova(rt.noint.cov)

ggplot(cm, aes(x=bpddiagScore, y=emoSlide.RT)) + geom_line() + facet_grid(congruent~faceemot)

#since there is no emotion x congruence interaction, just look at emotion x bpd for a minute
rt.cov_noint <- lmer(emoSlide.RT ~ congruent + faceemot + bpddiagScore + (1|subid), subData.all) #faceemot + congruent + faceemot*bpddiagScore + congruent*bpddiagScore + (1 |subid), subData.all)
summary(rt.cov_noint)

#add faceemot:bpddiagScore interaction to test for improvement
rt.cov <- lmer(emoSlide.RT ~ congruent + faceemot + bpddiagScore + faceemot:bpddiagScore + (1|subid), subData.all) #faceemot + congruent + faceemot*bpddiagScore + congruent*bpddiagScore + (1 |subid), subData.all)
summary(rt.cov)

print((pv <- pvals.fnc(rt.cov))$fixed)

anova(rt.cov_noint, rt.cov)


rt.cov_nocong <- lmer(emoSlide.RT ~ faceemot + bpddiagScore + faceemot:bpddiagScore + (1|subid), subData.all) #faceemot + congruent + faceemot*bpddiagScore + congruent*bpddiagScore + (1 |subid), subData.all)
summary(rt.cov_nocong)

rt.cov_nocong_run <- lmer(emoSlide.RT ~ faceemot + bpddiagScore + faceemot:bpddiagScore + (1+runnum|subid), subData.all) #faceemot + congruent + faceemot*bpddiagScore + congruent*bpddiagScore + (1 |subid), subData.all)
summary(rt.cov_nocong_run)

anova(rt.cov_nocong, rt.cov_nocong_run)

cm <- lmerCellMeans(rt.cov_nocong_run)
ggplot(cm, aes(x=bpddiagScore, y=emoSlide.RT)) + geom_line() + facet_wrap(~faceemot)

cm <- lmerCellMeans(rt.cov_nocong_run, divide="bpddiagScore")
cm$bpddiagScore <- factor(cm$bpddiagScore, levels=c("bpddiagScore: -1 SD", "bpddiagScore: M", "bpddiagScore: +1 SD"), labels=c("Low\n(-1 SD)", "Medium\n(M)", "High\n(+1 SD)"))
ggplot(cm, aes(x=bpddiagScore, y=emoSlide.RT, fill=faceemot)) + geom_bar(stat="identity", position="dodge")

pdf("rtModBPDHapAng.pdf", width=8, height=6) #units="in", res=300)
ggplot(cm, aes(x=bpddiagScore, y=emoSlide.RT, color=faceemot)) + #geom_bar(stat="identity", position="dodge") +
    geom_pointrange(aes(ymin=emoSlide.RT - se, ymax=emoSlide.RT + se), size=1.8, position=position_dodge(width=0.5)) +
    theme_bw(base_size=22) + xlab("\nBPD Symptom Severity") + ylab("Reaction Time (ms)\n") + ylim(700, 1050) + 
    scale_color_brewer("Facial Emotion", palette="Dark2") + theme(axis.title.x=element_text(vjust=0), legend.key.size=unit(2, "lines")) +
    ggtitle("Fig 3. Modulation of RTs by BPD Symptoms") + geom_line(aes(x=as.numeric(bpddiagScore)), position=position_dodge(width=0.5))
dev.off()


#try testing slope for each emotion
#ref is angry
conmat <- rbind("bpdDiag effect for Fearful"=c(
        0, #Intercept
        0, #faceemot fearful
        0, #faceemot happy
        1, #bpd diag score (angry)
        1, #fearful x bpddiag
        0  #happy x bpddiag
    ),
    "bpdDiag effect for Angry"=c(
        0, #Intercept
        0, #faceemot fearful
        0, #faceemot happy
        1, #bpd diag score (angry)
        0, #fearful x bpddiag
        0  #happy x bpddiag
    ),
    "bpdDiag effect for Happy"=c(
        0, #Intercept
        0, #faceemot fearful
        0, #faceemot happy
        1, #bpd diag score (angry)
        0, #fearful x bpddiag
        1  #happy x bpddiag
    ),
    "bpdDiag effect for Fear - Angry"=c(
        0, #Intercept
        0, #faceemot fearful
        0, #faceemot happy
        -1, #bpd diag score (angry)
        0, #fearful x bpddiag
        1  #happy x bpddiag
    ),
    "bpdDiag effect for M(Angry+Happy) - Fear"=c(
        0, #Intercept
        0, #faceemot fearful
        0, #faceemot happy
        0, #bpd diag score (angry)
        -1, #fearful x bpddiag
        0.5  #happy x bpddiag
    ))

summary(glht(rt.cov_nocong_run, linfct=conmat))



rt.noint.cov <- lmer(emoSlide.RT ~ faceemot + congruent + faceemot*paiharm + faceemot*paiaffe + (1|subid), subData.all) #faceemot + congruent + faceemot*bpddiagScore + congruent*bpddiagScore + (1 |subid), subData.all)
summary(rt.noint.cov)


#accuracy model including bpd
acc.noint <- glmer(emoSlide.ACC ~ faceemot + congruent + bpddiagScore + (1|subid), subData.all, family=binomial)

#no evidence of 3-way interaction or face x congruence interaction
acc.int <- glmer(emoSlide.ACC ~ faceemot * congruent * bpddiagScore + (1|subid), subData.all, family=binomial)

#pare down bpd model
acc.int <- glmer(emoSlide.ACC ~ faceemot + congruent + bpddiagScore*faceemot + (1|subid), subData.all, family=binomial, nAGQ=20)
summary(acc.int)

anova(acc.noint, acc.int)

#drop congruence for a minute
acc.int <- glmer(emoSlide.ACC ~ faceemot + bpddiagScore*faceemot + (1|subid), subData.all, family=binomial, nAGQ=20) #+ depanh*faceemot 
summary(acc.int)

#cm <- lmerCellMeans(acc.int)
cm <- lmerCellMeans(acc.int, divide="bpddiagScore")
cm$prob.emoSlide.ACC <- gtools::inv.logit(cm$emoSlide.ACC)
cm$sehigh <- gtools::inv.logit(cm$emoSlide.ACC+cm$se)
cm$selow <- gtools::inv.logit(cm$emoSlide.ACC-cm$se)
cm$bpddiagScore <- factor(cm$bpddiagScore, levels=c("bpddiagScore: -1 SD", "bpddiagScore: M", "bpddiagScore: +1 SD"), labels=c("Low\n(-1 SD)", "Medium\n(M)", "High\n(+1 SD)"))

pdf("accModBPDHapAng.pdf", width=8, height=6) #units="in", res=300)
ggplot(cm, aes(x=bpddiagScore, y=prob.emoSlide.ACC, color=faceemot)) + #geom_bar(stat="identity", position="dodge") +
    geom_pointrange(aes(ymin=selow, ymax=sehigh), size=1.7, position=position_dodge(width=0.5)) +
    theme_bw(base_size=22) + xlab("\nBPD Symptom Severity") + ylab("Probability of Correct Response\n") + ylim(0.5, 1.0) + 
    scale_color_brewer("Facial Emotion", palette="Dark2") + theme(axis.title.x=element_text(vjust=0), legend.key.size=unit(2, "lines")) +
    ggtitle("Fig 5. Modulation of RTs by BPD Symptoms") + geom_line(aes(x=as.numeric(bpddiagScore)), position=position_dodge(width=0.5))
dev.off()

#pdf("accBPD.pdf", width=8, height=6)
#ggplot(cm, aes(x=bpddiagScore, y=prob.emoSlide.ACC, color=bpddiagScore, ymin=selow, ymax=sehigh)) + 
#    geom_pointrange(size=1.7, position=position_dodge(width=0.2)) + theme_bw(base_size=22) +
#    ylim(0.5, 1) + scale_color_brewer("Face-Word\nCongruence", palette="Set1") + theme(axis.title.x=element_text(vjust=0),legend.key.size=unit(2, "lines")) +
#    ylab("Probability of Correct Response\n") + xlab("\nFacial Emotion") + 
#    ggtitle("Fig 4. Accuracy Across Conditions")#geom_errorbar(position=position_dodge(width=1), width=0.4) + 
#dev.off()
