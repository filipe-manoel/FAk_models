fit_asreml_model <- function(structure, k = NULL, dat, max_update = 5) {
  
  # Garante que singularidades não parem o ajuste
  asreml.options(ai.sing = TRUE)  
  
  # Define a estrutura aleatória
  rand_str <- switch(structure,
                     "rr"    = paste0("~ gen + rr(env,", k, "):gen + diag(env):gen + at(env):rep:iblock"),
                     "fa"    = paste0("~ gen + fa(env,", k, "):gen + at(env):rep:iblock"),
                     "sfa"   = paste0("~ gen + sfa(env,", k, "):gen + at(env):rep:iblock"),
                     "facv"  = paste0("~ gen + facv(env,", k, "):gen + at(env):rep:iblock"),
                     "us"    = "~ gen + us(env):id(gen) + at(env):rep:iblock",
                     "corgh" = "~ gen + corgh(env):id(gen) + at(env):rep:iblock",
                     stop("Unknown structure")
  )
  
  rand_term <- as.formula(rand_str)
  
  # Ajusta modelo inicial
  m <- asreml(
    fixed   = ht ~ env + env:rep,
    random  = rand_term,
    residual = ~ dsum(~units|env),
    workspace = "4gb",
    data = dat
  )
  
  # Atualiza apenas se necessário
  if (!m$converge) {
    count <- 0
    while (!m$converge && count < max_update) {
      m <- update(m)
      count <- count + 1
      print(paste("Update number", count))      
    }
    if (!m$converge) {
      warning(paste("Modelo", structure, "não convergiu após", max_update, "updates"))
    }
  }
  
  return(m)
}


# Modelo US
modelo_us <- fit_asreml_model("us", dat = dat)
summary(modelo_us)$varcomp

# Modelo FA(2)
m_fa2 <- fit_asreml_model("fa", k = 2, dat = dat)
summary(m_fa2)$varcomp
