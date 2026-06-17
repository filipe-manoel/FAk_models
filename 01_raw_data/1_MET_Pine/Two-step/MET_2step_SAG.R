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
  type.gen = "random", type.trial = "fixed", vc.model = "fa3")
met.model.fa$call # ASReml-R call.
met.model.fa$gof.stats # g-o-f statistics.
head(met.model.fa$predictions) # Predictions of gen.
met.model.fa$vcov.g # GxE variance-covariance matrix.
met.model.fa$corr.g # GxE correlation matrix.
met.model.fa$fa.loadings # The loading of each factor (in FA).

BLUP <- summary(met.model.fa$mod,coef=TRUE)$coef.random
head(BLUP)
met.model.fa$blups

library(ASRgenomics)
kinship.heatmap(K=met.model.fa$corr.g)

# Expand FA for additional output
? fa.summary
sx <- fa.summary(met.model.fa, gen.id = "Freeman", rotation='svd')
sx <- fa.summary(met.model.fa, gen.id = "NE16516", rotation='svd', type.resp = "blup")
sx$comp.gen.rot
sx$fa.plot
sx$fa.loadings.rot
sx$cum.var
sx$cum.var.trial.rot

sx <- fa.summary(met.model.fa, gen.id = "Freeman", rotation='svd', type.resp = "blup")
sx$fa.plot
sx <- fa.summary(met.model.fa, gen.id = "NE16402", rotation='svd', type.resp = "blup")
sx$fa.plot

# Biplot
gbiplot(data = met.model.fa$predictions, vcov.g = met.model.fa$vcov.g,
        vector = "location", unit = "gen", resp = "predicted.value")

# Stability Index
stab.index <- stability(
  data = met.model.fa$predictions, trial = "location", gen = "gen",
  resp = "predicted.value", method = "static", top = TRUE)
stab.index$stability.plot

