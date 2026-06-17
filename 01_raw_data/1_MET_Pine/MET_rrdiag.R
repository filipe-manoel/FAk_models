##################################
## Multi-Environment Trials - MET
##  Conifer Dataset (4 trials)
##################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/FilipeF/MET_Pine")
library(asreml)
library(ASRtriala)
library(ggplot2)

# Reading data
datam<-read.table("TRIALS4.txt",h=T)
head(datam)
datam$Rep<-as.factor(datam$Rep)
datam$Iblock<-as.factor(datam$Iblock)
datam$Test<-as.factor(datam$Test)
datam$Genotype<-as.factor(datam$Genotype)
str(datam)

# Performing some EDA; Exploratory data analysis
boxplot(HT~Test,data=datam)
aggregate(HT~Test,var,data=datam)
mean(aggregate(HT~Test,mean,data=datam)[,2])
table(datam$Test,datam$Genotype)
meanst<-aggregate(HT~Genotype+Test,mean,data=datam)
Int.Plot<-ggplot(data=meanst, aes(x=Test, y=HT, group=Genotype))+
  geom_line(linewidth=1, aes(color=Genotype))+
  ylab("HT")+xlab("Test")
Int.Plot

##########################
# Model for a single test

# BV = 2*GCA
# Var(BV) = Var(a) = Var(2*GCA) = 4*Var(GCA) 

# y = mu + Rep + Rep:Iblock + Genotype + e
levels(datam$Test)
model1<-asreml(fixed=HT~Rep,
	             random=~Rep:Iblock+Genotype,
               residual=~units, # ar1v(row):ar1(col)
	             subset=Test==3, #test only for site 1
               data=datam)
summary(model1)$varcomp
plot(model1)
(h2b<-vpredict(model1, h2b~4*V1/(V1+V2+V3)))

#############################
### Fitting MET models 

# Connectivity
? audit.met
ck <- audit.met(data=datam, trial='Test', gen='Genotype', resp='HT')
ck$met.stats
ck$gen.inc
levels(datam$Test)

############
# US  -78061.03

modelu<-asreml(fixed=HT~Test+Test:Rep,
                random=~at(Test):Rep:Iblock+corgh(Test):Genotype,
                residual=~dsum(~units|Test),
                data=datam)
summary(modelu)$varcomp
BLUP <- summary(modelu, coef=TRUE)$coef.random
head(BLUP)
tail(BLUP)
pred_site <- predict.asreml(modelu, classify="Test:Genotype")$pvals
pred_all <- predict.asreml(modelu, classify="Genotype")$pvals
head(pred_site)
head(pred_all)

############
# Factor Analytic - fa1

initg<-c(450,450,450,450,0.5,0.5,0.5,0.5)
model2f<-asreml(fixed=HT~Test+Test:Rep,
                random=~at(Test):Rep:Iblock+
                  fa(Test,1,init=initg):Genotype,
                residual=~dsum(~units|Test),data=datam)
summary(model2f)$varcomp

(R<-summary(model2f)$varcomp$component[1:4])
(L<-summary(model2f)$varcomp$component[5:8])
(V<-L%*%t(L)+diag(R))
(CORR<-cov2cor(V))

BLUP2f<-summary(model2f,coef=TRUE)$coef.random
View(BLUP2f)

############
# Factor Analytic - fa2
initg<-c(197,74,129,185,18,22,16,13,0.01,0.01,0.01,0.01)
model2f2<-asreml(fixed=HT~Test+Test:Rep,
                random=~at(Test):Rep:Iblock+
                  fa(Test,2,init=initg):Genotype,
                residual=~dsum(~units|Test),data=datam)
model2f2<-update.asreml(model2f2)
summary(model2f2)$varcomp

(R<-summary(model2f2)$varcomp$component[1:4])
(L1<-summary(model2f2)$varcomp$component[5:8])
(L2<-summary(model2f2)$varcomp$component[9:12])
(L<-cbind(L1,L2))
V<-L%*%t(L)+diag(R)
R <- cov2cor(V)
colnames(V) <- levels(datam$Test)
rownames(V) <- levels(datam$Test)
colnames(R) <- levels(datam$Test)
rownames(R) <- levels(datam$Test)
V
R

###################################
# Evaluating rr2+diag

summary(model2f2)

# Factor Analytic - fa2 - as rr2 + daig
#initg<-c(197,74,129,185,18,22,16,13,0.01,0.01,0.01,0.01)
model2rrd<-asreml(fixed=HT~Test+Test:Rep,
                 random=~at(Test):Rep:Iblock+
                   rr(Test,2):Genotype + diag(Test):Genotype,
                 residual=~dsum(~units|Test),data=datam)
model2rrd<-update.asreml(model2rrd)
summary(model2rrd)

(R<-summary(model2rrd)$varcomp$component[1:4])
(L1<-summary(model2rrd)$varcomp$component[9:12])
(L2<-summary(model2rrd)$varcomp$component[13:16])
(L<-cbind(L1,L2))
V<-L%*%t(L)+diag(R)
R <- cov2cor(V)
colnames(V) <- levels(datam$Test)
rownames(V) <- levels(datam$Test)
colnames(R) <- levels(datam$Test)
rownames(R) <- levels(datam$Test)
V
R

#######################
#######################
# Comparing some definitions

# corgh -78061.03 
model2f<-asreml(fixed=HT~Test+Test:Rep,
                random=~at(Test):Rep:Iblock+
                  corgh(Test):Genotype,
                residual=~dsum(~units|Test),data=datam)
summary(model2f)$varcomp

# corh -78066.73 
model3<-asreml(fixed=HT~Test+Test:Rep,
                random=~at(Test):Rep:Iblock+
                  corh(Test):Genotype,
                residual=~dsum(~units|Test),data=datam)
summary(model3)$varcomp
preds3 <- predict(model3, classify="Test:Genotype")$pvals

lrt.asreml(model2f, model3, boundary=FALSE)

# corv -78068.76 
model4<-asreml(fixed=HT~Test+Test:Rep,
               random=~at(Test):Rep:Iblock+
                 corv(Test):Genotype,
               residual=~dsum(~units|Test),data=datam)
summary(model4)$varcomp

# cs -78068.76  
model5<-asreml(fixed=HT~Test+Test:Rep,
               random=~at(Test):Rep:Iblock+
                 Genotype + Test:Genotype,
               residual=~dsum(~units|Test),data=datam)
summary(model5)$varcomp
vpredict(model5, r2 ~ V1/(V1+V2))

# cs + diag  -78068.59 
model6<-asreml(fixed=HT~Test+Test:Rep,
               random=~at(Test):Rep:Iblock+
                 Genotype + diag(Test):Genotype,
               residual=~dsum(~units|Test),data=datam)
summary(model6)$varcomp
vpredict(model6, r2 ~ V1/(V1+(V2+V3+V4+V5)/4))
preds6 <- predict(model6, classify="Test:Genotype")$pvals

lrt.asreml(model5, model6, boundary=FALSE)

# corgv (failing)
asreml.options(ai.sing=TRUE)
initg<-c(0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 100)
model7<-asreml(fixed=HT~Test+Test:Rep,
               random=~at(Test):Rep:Iblock+
                 corgv(Test, init=initg):Genotype,
               residual=~dsum(~units|Test),data=datam)
summary(model7)$varcomp
