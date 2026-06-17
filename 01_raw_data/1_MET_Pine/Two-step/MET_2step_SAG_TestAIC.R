##########################################
## Using ASRtriala for 2-stage analyses
##  Wheat Dataset
##########################################

rm(list=ls())
setwd("C:/Users/sgeza/OneDrive/Desktop/Workshops ASReml-R/ASReml_2025/Material/Session10")
library(asreml)   
library(ASRtriala)

# Reading dataset
? pheno.wheat        # raw dataset
? pheno.wheat.preds  # predicted from step 1

# Pre-fitted data from pheno.wheat.
pheno.wheat.preds <- ASRtriala::pheno.wheat.preds
str(pheno.wheat.preds)
head(pheno.wheat.preds,20)
tail(pheno.wheat.preds,20)
levels(pheno.wheat.preds$location)
pheno.wheat.preds[pheno.wheat.preds$gen == 'NE16401',]

# Audit MET
audit.all <- audit.met(
  data = pheno.wheat.preds, gen = "gen",
  trial = "location", resp = "predicted.value")
audit.all$met.stats 
audit.all$gen.inc 
audit.all$gen.inc.perc

# NEED to make a decision on what to do with Checks! (drop or not drop?)

# Selecting best MET model from a set of models
? select.met

# Get statistics for several GxE variance-covariance structures.
met.stats <- select.met(
  data = pheno.wheat.preds, trial = "location", gen = "gen",
  resp = "predicted.value", weight = "weight", type.trial = 'fixed',
  vc.models = c('corv', 'corh', 'fa1', 'fa2', 'fa3'),
  criteria = "AIC")
met.stats$gof.stats
met.stats$best.call

# Now fitting the selected MET model
? fit.met

# Fit factor analytic of 3rd order.
met.model.fa <- fit.met(
  data = pheno.wheat.preds, gen = "gen", trial = "location",
  resp = "predicted.value", weight = "weight",
  type.gen = "random", type.trial = "fixed", vc.model = "fa1")
met.model.fa$call # ASReml-R call.
met.model.fa$gof.stats # g-o-f statistics.
#head(met.model.fa$predictions) # Predictions of gen.
#met.model.fa$vcov.g # GxE variance-covariance matrix.
#met.model.fa$corr.g # GxE correlation matrix.
#met.model.fa$fa.loadings # The loading of each factor (in FA).

# Looking at the call
str(pheno.wheat.preds)
is.na(pheno.wheat.preds$weight)
# Data needs to be filtered by missing values
pheno2 <- pheno.wheat.preds[!is.na(pheno.wheat.preds$weight),]
mymod <- asreml(fixed = predicted.value ~ 1 + location, 
               random = ~fa(location, 1):id(gen), 
               na.action = list(x = "include", y = "include"), 
               weights = weight,
               family = asreml::asr_gaussian(dispersion = 1), 
               data = pheno2, 
               workspace = 1.28e+08)
mymod <- update.asreml(mymod)

# Comparing varcomp
summary(mymod)$varcomp
summary(met.model.fa$mod)$varcomp

# Comparing logL
mymod$loglik
met.model.fa$mod$loglik

summary(mymod)$bic
summary(met.model.fa$mod)$bic


