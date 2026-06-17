# FA and Variants ---------------------------------------------------------
rm(list=ls()) 

library(asreml)
library(ASRtriala)
library(janitor)
library(tidyverse)

# Phenotypic data --------------------------------------------------------------------
# Load the raw data (rdata)
rdata = read.table("01_raw_data/1_MET_Pine/TRIALS4.txt",h=T)
str(rdata)

# Tidy the data (dat)
dat = rdata %>% 
  janitor::clean_names() %>% 
  rename(env = test,
         gen = genotype,
         col = column) %>% 
  mutate_at(vars(idd:col), as.factor)%>% 
  mutate_at(vars(surv:ht), as.numeric) 

write_csv(dat, "02_data/tidy_data.csv")


# Fit the US model --------------------------------------------------------
asreml.options(ai.sing = TRUE)

m_us = asreml(
  fixed   =  ht ~ env + env:rep,
  random  = ~ gen + us(env):id(gen) + at(env):rep:iblock,
  residual = ~ dsum(~units|env),
  workspace = "4gb",
  data = dat)

vc_us = summary(m_us)$varcomp

# Fit the corgh model --------------------------------------------------------
asreml.options(ai.sing = TRUE)

m_corgh = asreml(
  fixed   =  ht ~ env + env:rep,
  random  = ~ gen + corgh(env):id(gen) + at(env):rep:iblock,
  residual = ~ dsum(~units|env),
  workspace = "4gb",
  data = dat)

vc_corgh = summary(m_corgh)$varcomp

# Fit the rr + diag 1 model --------------------------------------------------------
asreml.options(ai.sing = TRUE)

m_rrdiag1 = asreml(
  fixed   =  ht ~ env + env:rep,
  random  = ~ gen + rr(env,1):gen + diag(env):gen + at(env):rep:iblock,
  residual = ~ dsum(~units|env),
  workspace = "4gb",
  data = dat)

m_rrdiag1 = update(m_rrdiag1)
m_rrdiag1 = update(m_rrdiag1)
m_rrdiag1 = update(m_rrdiag1)
m_rrdiag1 = update(m_rrdiag1)
m_rrdiag1 = update(m_rrdiag1)

vc_rrdiag1 = summary(m_rrdiag1)$varcomp

# Fit the rr + diag 2 model --------------------------------------------------------
asreml.options(ai.sing = TRUE)

m_rrdiag2 = asreml(
  fixed   =  ht ~ env + env:rep,
  random  = ~ gen + rr(env,2):gen + diag(env):gen + at(env):rep:iblock,
  residual = ~ dsum(~units|env),
  workspace = "4gb",
  data = dat)

m_rrdiag2 = update(m_rrdiag2)
summary(m_rrdiag2)$varcomp

# Fit the FA1 model --------------------------------------------------------
asreml.options(ai.sing = TRUE)

m_fa1 = asreml(
  fixed   =  ht ~ env + env:rep,
  random  = ~ gen + fa(env,1):gen + at(env):rep:iblock,
  residual = ~ dsum(~units|env),
  workspace = "4gb",
  data = dat)

m_fa1 = update(m_fa1)
m_fa1 = update(m_fa1)
m_fa1 = update(m_fa1)
m_fa1 = update(m_fa1)
m_fa1 = update(m_fa1)

vc_fa1 = summary(m_fa1)$varcomp

# Fit the FA2 model --------------------------------------------------------
asreml.options(ai.sing = TRUE)

m_fa2 = asreml(
  fixed   =  ht ~ env + env:rep,
  random  = ~ gen + fa(env,2):gen + at(env):rep:iblock,
  residual = ~ dsum(~units|env),
  workspace = "4gb",
  data = dat)

m_fa2 = update(m_fa2)

vc_fa2 = summary(m_fa2)$varcomp


# Fit the facv 1 model --------------------------------------------------------

asreml.options(ai.sing = TRUE)

m_facv1 = asreml(
  fixed   =  ht ~ env + env:rep,
  random  = ~ gen + facv(env,1):gen + at(env):rep:iblock,
  residual = ~ dsum(~units|env),
  workspace = "4gb",
  data = dat)

m_facv1 = update(m_facv1)
m_facv1 = update(m_facv1)
m_facv1 = update(m_facv1)

vc_facv1 = summary(m_facv1)$varcomp


# Fit the sfa 1 model --------------------------------------------------------

asreml.options(ai.sing = TRUE)

m_sfa1 = asreml(
  fixed   =  ht ~ env + env:rep,
  random  = ~ gen + sfa(env,1):gen + at(env):rep:iblock,
  residual = ~ dsum(~units|env),
  workspace = "4gb",
  data = dat)

m_sfa1 = update(m_sfa1)
m_sfa1 = update(m_sfa1)
m_sfa1 = update(m_sfa1)
m_sfa1 = update(m_sfa1)
m_sfa1 = update(m_sfa1)

vc_sfa1 = summary(m_sfa1)$varcomp

# Componentes de Variâncias -----------------------------------------------

write.csv(vc_us, "04_output/varcomp_us.csv")
write.csv(vc_corgh, "04_output/varcomp_us.csv")
write.csv(vc_rrdiag1, "04_output/varcomp_us.csv")
write.csv(vc_rrdiag2, "04_output/varcomp_us.csv")
write.csv(vc_fa1, "04_output/varcomp_us.csv")
write.csv(vc_fa2, "04_output/varcomp_us.csv")
write.csv(vc_facv1, "04_output/varcomp_us.csv")
write.csv(vc_sfa1, "04_output/varcomp_us.csv")


























# Fit the US model --------------------------------------------------------

names(dat)
traits = names(dat)[8:10]

results_us = list()

for (i in 1:length(traits)) {
  
  dat = dat %>% arrange(column, row)

    # Fórmula fixa dinâmica
  fixed_formula = as.formula(
    paste0(traits[i], 
           " ~ test")
  )
  
  #Modelo
  #asreml.options(gammaPar = FALSE, ai.sing = TRUE)
  m_us = asreml(
    #fixed   = fixed_formula,
    fixed   = fixed_formula,
    random  = ~ rep:iblock + us(test):id(genotype),
    residual = ~ dsum(~units|test),
    na.action = na.method(x = "include", y = "include"),
    workspace = "4gb",
    data = dat,
    sing.ai = TRUE
  )
  
  vc  = summary(m_us)$varcomp
  AIC = summary(m_us)$aic
  BIC = summary(m_us)$bic
  
  H2_us = vpredict(m_us, H2 ~ V2 / (V1 + V2 + V3 + (V5+V7+V10+V14)/4))
  sum.us = summary(m_us)$varcomp
  
  mvd.us = mean(((predict.asreml(m_us, classify = "at(tipo, 'Teste'):gen", sed = TRUE)$sed)^2)[
    upper.tri((predict.asreml(m_us, classify = "at(tipo, 'Teste'):gen", sed = TRUE)$sed)^2, diag = FALSE)
  ])
  PEV.us = mean(diag((predict.asreml(m_us, classify = "at(tipo, 'Teste'):gen", vcov = TRUE)$vcov)))
  acc.us = sqrt(1 - (PEV.us / sum.us["at(tipo, 'Teste'):gen","component"]))
  rel.us = 1 - (PEV.us / sum.us["at(tipo, 'Teste'):gen","component"])
  
  H2_Cullis_us = 1 - (mvd.us / (2 * sum.us["at(tipo, 'Teste'):gen","component"]))
  
  # BLUEs
  BLUE = data.frame(summary(m_us, coef = TRUE)$coef.fixed)
  BLUE = BLUE[grep("gen_T", rownames(BLUE)), ]
  
  # BLUPs
  BLUP = data.frame(summary(m_us, coef = TRUE)$coef.random)
  BLUP = BLUP[grep("gen_", rownames(BLUP)), ]
  
  pop_mean = mean(dat[[traits[i]]], na.rm = TRUE)
  blup.m0 = predict.asreml(m_us, classify = "at(tipo, 'Teste'):gen")$pvals
  sel_mean = mean(blup.m0[order(blup.m0$predicted.value, decreasing = TRUE),]$predicted.value[1:10])
  
  gain.us = (sel_mean - pop_mean) / pop_mean * 100
  
  # Guardar tudo num único objeto
  resultados_us[[traits[i]]] = list(
    vc      = list(varcomp = vc, AIC = AIC, BIC = BIC),
    H2      = list(
      H2_asreml = H2_us,
      H2_Cullis = H2_Cullis_us,
      PEV       = PEV.us,
      acc       = acc.us,
      rel       = rel.us
    ),
    gain    = list(
      pop_mean = pop_mean,
      sel_mean = sel_mean,
      gain     = gain.us
    ),
    effects = list(
      BLUE = BLUE,
      BLUP = BLUP
    )
  )
}






head(datam)
datam$Rep=as.factor(datam$Rep)
datam$Iblock=as.factor(datam$Iblock)
datam$Test=as.factor(datam$Test)
datam$Genotype=as.factor(datam$Genotype)
str(datam)

# Performing some EDA; Exploratory data analysis
boxplot(HT~Test,data=datam)
aggregate(HT~Test,var,data=datam)
mean(aggregate(HT~Test,mean,data=datam)[,2])
table(datam$Test,datam$Genotype)
meanst=aggregate(HT~Genotype+Test,mean,data=datam)
Int.Plot=ggplot(data=meanst, aes(x=Test, y=HT, group=Genotype))+
  geom_line(linewidth=1, aes(color=Genotype))+
  ylab("HT")+xlab("Test")
Int.Plot