##################################
## Multi-Environment Trials - MET
##  Conifer Dataset (4 trials)
##################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Workshops ASReml-R/ASReml_2025/Material/Session10")
library(asreml)
library(ASRtriala)
library(ggplot2)

# Reading data
datam<- read.table("01_raw_data/1_MET_Pine/TRIALS4.txt",h=T)

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

# y = mu + Test + Test:Rep + at(Test):Rep:Iblock + Genotype + Genotype:Test + dsum(e|Test)

### GxE as Interaction
modelF<-asreml(fixed=HT~Test+Rep:Test,
               random=~at(Test,c(1,2,3,4)):Rep:Iblock +
                 Genotype+Genotype:Test, # Interaction
               residual=~dsum(~units|Test),
               data=datam)
summary(modelF)$varcomp

# Simple Model for all sites #
model2<-asreml(fixed=HT~Test+Test:Rep,
               random=~at(Test):Rep:Iblock #if all sites have incomplete blocks
                      +Genotype+Test:Genotype,
               residual=~dsum(~units|Test),
               #residual=~dsum(~ar1v(row):ar1(col)|Test,level='1') +
               #          dsum(~ar1v(row):id(col)|Test,level='2') +
               #          dsum(~idv(units)|Test,level=c('3','4'),
               data=datam)
summary(model2)$varcomp
plot(model2)

vpredict(model2,r2B~V1/(V1+V2)) #Type-B correlation
vpredict(model2,h2~4*V1/(V1+V2+(V3+V4+V5+V6)/4+(V7+V8+V9+V10)/4)) #heritability across environments

# BLUPs (we can reconstruct responses)
BLUP<-summary(model2,coef=TRUE)$coef.random
View(BLUP)

# Some predictions
ppG2<-predict(model2,classify="Genotype")$pvals
head(ppG2,20)
ppGE2<-predict(model2,classify="Test:Genotype")$pvals
head(ppGE2)
View(ppGE2)

######################
### GxE as Nested Eff.

# Simple corv similar to GxE-Int
initg<-c(0.65,458)
model2b<-asreml(fixed=HT~Test+Test:Rep,
               random=~at(Test):Rep:Iblock+
                      corv(Test,init=initg):Genotype,
                residual=~dsum(~units|Test),data=datam)
summary(model2b)$varcomp
plot(model2b)

BLUP2b<-summary(model2b,coef=TRUE)$coef.random
head(BLUP2b)
tail(BLUP2b)

ppGE2b<-predict(model2b,classify="Test:Genotype")$pvals
head(ppGE2b)

# Complex corgh
initg<-c(0.65,0.65,0.65,0.65,0.65,0.65,450,450,450,450)
model2c<-asreml(fixed=HT~Test+Test:Rep,
                random=~at(Test):Rep:Iblock+
                  corgh(Test,init=initg):Genotype,
                residual=~dsum(~units|Test),data=datam)
summary(model2c)$varcomp
BLUP2c<-summary(model2c,coef=TRUE)$coef.random
head(BLUP2c)

# Comparing models 2b & 2c
lrt.asreml(model2b,model2c,boundary=FALSE)

ppGE2c<-predict(model2c,classify="Test:Genotype")$pvals
head(ppGE2c)
View(ppGE2c)

preds <-predict(model2c,classify="Genotype")$pvals
head(preds)

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

# Comparing Models (only nested ones)
lrt(model2f2,model2f,boundary=FALSE)

# AIC and BIC (the more negative the better)
c(summary(model2f)$bic, summary(model2f)$aic)  # BIC, AIC
c(summary(model2f2)$bic, summary(model2f2)$aic)

preds <- predict.asreml(model2f2, classify='Genotype')$pvals
head(preds,10)
ppGE2c<-predict(model2f2,classify="Test:Genotype")$pvals
head(ppGE2c)

BLUP2 <- summary(model2f2,coef=TRUE)$coef.random
View(BLUP2)

# Some stability tools
? ASRtriala::stability
stab.index <- stability(data = ppGE2c,
                        trial = "Test", gen = "Genotype", resp = "predicted.value",
                        method = "static", best = "max", plot = TRUE, top = TRUE,
                        bottom = FALSE, percentage = 5)
stab.index$stability.plot

# Biplot
? ASRtriala::gbiplot
gbiplot(data = ppGE2c, vcov.g = V,
        vector = "Test", unit = "Genotype", resp = "predicted.value",
        unit.label = TRUE)

