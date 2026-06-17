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



