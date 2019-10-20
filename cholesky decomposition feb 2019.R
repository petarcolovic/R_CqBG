# cholesky decomposition - preliminary syntax

# arguments:

# data: dataset (data frame); the function should convert the data into data frame

# first, second: variables for the first / second twin. Variable names and order should be stated exactly as they are in the dataset
  # the "first" and "second" arguments should be stated as e.g. first = c("V1_tw1", "V2_tw1", "V3_tw1"), second = c("V1_tw2", "V2_tw2", "V3_tw2")

# zygosity: zygosity variable; preferably, the cases with the missing zygosity data (if any) should be excluded prior to the analysis

# Currently, only Cholesky ACE (aditive genetic factor, shared environment (C), non-shared environment (E)) and AE (additive genetic influences, non-shared environment) models, as the most frequent, are available. Other models, if there is need for them, can be added.

cholfun <- function(data, first, second, zygo){
  
  library(lavaan)
  library(semPlot)
  library(stringr)
  library(stringi)
  library(car)
  library(dplyr)
  library(psych)
  
  data <- as.data.frame(data)
  
  zygosity <- zygo
  
  zygosity <- subset(data, select = c(zygosity))
  
  colnames(zygosity) <- "zygosity"
  
  str.sub.first <- stri_extract_all_boundaries(first, simplify = TRUE)
  
  str.sub.second <- stri_extract_all_boundaries(second, simplify = TRUE)
  
  nvar <- length(first)
  
  first <- function(nvar){for(i in 1:nvar)
    (cat("V",i,".twin1", sep = "", fill = TRUE))
  }
  
  second <- function(nvar){for(i in 1:nvar)
    (cat("V",i,".twin2", sep = "", fill = TRUE))
  }
  
  firstlabels <- capture.output(first(nvar))
  
  secondlabels <- capture.output(second(nvar))
  
  for(i in 1:nvar){
    labels.first <- (str_sub(str.sub.first))
    labels.second <- (str_sub(str.sub.second))
    twinvar <- subset(data, select = c(labels.first, labels.second))}
  twinvard <- as.data.frame(twinvar)
  twinvardren <- as.data.frame(twinvard)
  colnames(twinvardren) <- c(firstlabels,  secondlabels)
  complete <- cbind(twinvardren, zygosity)  

  cholesky.loadings.model.a<- function(nvar){for(i in (1:nvar)){for (j in (1:2)){
    cat(sprintf("A%s.%s=~c(%s, %s)*V%s.twin%s", i, j, paste("a",i,seq(from=i, to=nvar, by=1), sep=""), paste("a",i,seq(from=i, to=nvar, by=1), sep = ""), seq(from=i, to=nvar, by=1), j), sep = "\n", fill = FALSE)}}}
  
  cholesky.loadings.model.c<- function(nvar){for(i in (1:nvar)){for (j in (1:2)){
    cat(sprintf("C%s.%s=~c(%s, %s)*V%s.twin%s", i, j, paste("c",i,seq(from=i, to=nvar, by=1), sep=""), paste("c",i,seq(from=i, to=nvar, by=1), sep = ""), seq(from=i, to=nvar, by=1), j), sep = "\n", fill = FALSE)}}}
  
  cholesky.loadings.model.e<- function(nvar){for(i in (1:nvar)){for (j in (1:2)){
    cat(sprintf("E%s.%s=~c(%s, %s)*V%s.twin%s", i, j, paste("e",i,seq(from=i, to=nvar, by=1), sep=""), paste("e",i,seq(from=i, to=nvar, by=1), sep = ""), seq(from=i, to=nvar, by=1), j), sep = "\n", fill = FALSE)}}}
  
  cholesky.a <- as.data.frame(capture.output(cholesky.loadings.model.a(nvar)), ncol = 1, rownames = FALSE, optional = TRUE)
  cholesky.c <- as.data.frame(capture.output(cholesky.loadings.model.c(nvar)), ncol = 1, rownames = FALSE, optional = TRUE)
  cholesky.e <- as.data.frame(capture.output(cholesky.loadings.model.e(nvar)), ncol = 1, rownames = FALSE, optional = TRUE)
  
  
  
  chol.cor.ace <- function(nvar){for(i in 1:nvar)(cat(
    "A",i,".1","~~c(1,.5)*A",i,".2", "\n",
    "C",i,".1","~~c(1,1)*C",i,".2", "\n",
    "E",i,".1","~~c(0,0)*E",i,".2", "\n",
    "V",i,".twin1~~c(0,0)*V",i,".twin1","\n",
    "V",i,".twin2~~c(0,0)*V",i,".twin2", "\n",
    sep="", fill=FALSE))}
  
  chol.cor.ae <- function(nvar){for(i in 1:nvar)(cat(
    "A",i,".1","~~c(1,.5)*A",i,".2", "\n",
    "E",i,".1","~~c(0,0)*E",i,".2", "\n",
    "V",i,".twin1~~c(0,0)*V",i,".twin1","\n",
    "V",i,".twin2~~c(0,0)*V",i,".twin2", "\n",
    sep="", fill=FALSE))}
  
  
  chol.cor.ace.out <- as.data.frame(capture.output(chol.cor.ace(nvar)), ncol = 1, optional = TRUE)
  
  chol.cor.ae.out <- as.data.frame(capture.output(chol.cor.ae(nvar)), ncol = 1, optional = TRUE)
  
  colnames(cholesky.a) <- "col"
  colnames(cholesky.c) <- "col"
  colnames(cholesky.e) <- "col"
  colnames(chol.cor.ace.out) <- "col"
  colnames(chol.cor.ae.out) <- "col"
  
  cholesky.ace.mod <- rbind(cholesky.a, cholesky.c, cholesky.e, chol.cor.ace.out)
  
  cholesky.ae.mod <- rbind(cholesky.a, cholesky.e, chol.cor.ae.out)
  
  colnames(cholesky.ace.mod) <- NULL
  
  colnames(cholesky.ae.mod) <- NULL
  
  cholesky.ace.mod.def <- print(cholesky.ace.mod)
  
  cholesky.ae.mod.def <- print(cholesky.ae.mod)
  
  
  
  cholesky.ace.fit <- cfa(as.matrix(cholesky.ace.mod), data = complete, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "default")   
  
  cholesky.ae.fit <- cfa(as.matrix(cholesky.ae.mod), data = complete, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "default")
  
  # parameters Cholesky
  
  ekst.par.chol.stand <- function(modelfit){
    
    ord <- lavInspect(modelfit, what = "std", add.labels = TRUE)
    
    parameters <- ord$`1`$lambda
    
    parameters.data <- as.data.frame(parameters)
    
    columns <- select(parameters.data, ends_with(".1"))
    
    columns.t.data <- as.data.frame(t(columns))
    
    rows <- select(columns.t.data, ends_with(".twin1"))
    
    rows.squared <- rows^2
    
    rows.squared.data <- as.data.frame(rows.squared)
    
    rows.squared.data.coln <- rows.squared.data
    
    colnames(rows.squared.data.coln) <- (str.sub.first)
    
    return(rows.squared.data.coln)}
  
  parameters.stand.cholesky.ace <- ekst.par.chol.stand(cholesky.ace.fit)
  
  parameters.stand.cholesky.ae <- ekst.par.chol.stand(cholesky.ae.fit)
  
  # unstandardized parameters
  
  ekst.par.chol.unstand <- function(modelfit){
    
    ord <- lavInspect(modelfit, what = "free", add.labels = TRUE)
    
    parameters <- ord$`1`$lambda
    
    parameters.data <- as.data.frame(parameters)
    
    columns <- select(parameters.data, ends_with(".1"))
    
    columns.t.data <- as.data.frame(t(columns))
    
    rows <- select(columns.t.data, ends_with(".twin1"))
    
    rows.squared <- rows^2
    
    rows.squared.data <- as.data.frame(rows.squared)
    
    rows.squared.data.coln <- rows.squared.data
    
    colnames(rows.squared.data.coln) <- (str.sub.first)
    
    return(rows.squared.data.coln)}
  
  parameters.unstand.cholesky.ace <- ekst.par.chol.unstand(cholesky.ace.fit)
  
  parameters.unstand.cholesky.ae <- ekst.par.chol.unstand(cholesky.ae.fit)
  
  
  
  # calculation of genetic and environmental correlations
  
  # ace model
  
  # cholesky ace genetic and environmental correlations
  
  ord.ace <- lavInspect(cholesky.ace.fit, what = "est", add.labels = TRUE)
  
  parameters.ace <- ord.ace$`1`$lambda
  
  parameters.ace.sub <- parameters.ace[1:nvar,]
  
  parameters.ace.sub.data <- as.data.frame(parameters.ace.sub)
  
  parameters.sub.ace.a <- parameters.ace.sub.data
  
  parameters.sub.ace.c <- parameters.ace.sub.data
  
  parameters.sub.ace.e <- parameters.ace.sub.data
  
  parameters.sub.ace.a <- select(parameters.sub.ace.a, starts_with("A")) 
  
  parameters.sub.ace.a <- select(parameters.sub.ace.a, contains(".1"))
  
  parameters.sub.ace.c <- select(parameters.sub.ace.c, starts_with("C")) 
  
  parameters.sub.ace.c <- select(parameters.sub.ace.c, contains(".1"))
  
  parameters.sub.ace.e <- select(parameters.sub.ace.e, starts_with("E")) 
  
  parameters.sub.ace.e <- select(parameters.sub.ace.e, contains(".1"))
  
  parameters.sub.ace.ce <- cbind(parameters.sub.ace.c, parameters.sub.ace.e)
  
  gen.cov.ace.a <- as.matrix(parameters.sub.ace.a)%*%t(as.matrix(parameters.sub.ace.a))
  
  env.cov.ace.ce <- as.matrix(parameters.sub.ace.ce)%*%t(as.matrix(parameters.sub.ace.ce))
  
  gen.cor.ace.a <- cov2cor(as.matrix(parameters.sub.ace.a)%*%t(as.matrix(parameters.sub.ace.a)))
  
  env.cor.ace.ce <- cov2cor(as.matrix(parameters.sub.ace.ce)%*%t(as.matrix(parameters.sub.ace.ce)))
  
  colnames(gen.cor.ace.a) <- str.sub.first
  
  rownames(gen.cor.ace.a) <- str.sub.first
  
  colnames(env.cor.ace.ce) <- str.sub.first
  
  rownames(env.cor.ace.ce) <- str.sub.first
  
  # ae model
  
  ord.ae <- lavInspect(cholesky.ae.fit, what = "est", add.labels = TRUE)
  
  parameters.ae <- ord.ae$`1`$lambda
  
  parameters.ae.sub <- parameters.ae[1:nvar,]
  
  parameters.ae.sub.data <- as.data.frame(parameters.ae.sub)
  
  parameters.sub.ae.a <- parameters.ae.sub.data
  
  parameters.sub.ae.e <- parameters.ae.sub.data
  
  parameters.sub.ae.a <- select(parameters.sub.ae.a, starts_with("A")) 
  
  parameters.sub.ae.a <- select(parameters.sub.ae.a, contains(".1"))
  
  parameters.sub.ae.e <- select(parameters.sub.ae.e, starts_with("E")) 
  
  parameters.sub.ae.e <- select(parameters.sub.ae.e, contains(".1"))
  
  gen.cov.ae.a <- as.matrix(parameters.sub.ae.a)%*%t(as.matrix(parameters.sub.ae.a))
  
  env.cov.ae.e <- as.matrix(parameters.sub.ae.e)%*%t(as.matrix(parameters.sub.ae.e))
  
  gen.cor.ae.a <- cov2cor(as.matrix(parameters.sub.ae.a)%*%t(as.matrix(parameters.sub.ae.a)))
  
  env.cor.ae.e <- cov2cor(as.matrix(parameters.sub.ae.e)%*%t(as.matrix(parameters.sub.ae.e)))
  
  colnames(gen.cor.ae.a) <- str.sub.first
  
  rownames(gen.cor.ae.a) <- str.sub.first
  
  colnames(env.cor.ae.e) <- str.sub.first
  
  rownames(env.cor.ae.e) <- str.sub.first
  
  # comparison of cholesky models
  
  comparison <- as.matrix(anova(cholesky.ace.fit, cholesky.ae.fit))
  
  comp.sub <- as.matrix(subset(comparison, select=c("BIC")))
  
  comp.sub.sort <- as.matrix(sort(comp.sub, decreasing = FALSE))
  
  return(list(cholesky.ace = cholesky.ace.fit, cholesky.ae = cholesky.ae.fit, par.ace = parameters.stand.cholesky.ace, par.ae = parameters.stand.cholesky.ae, genetic.corr.ace = gen.cor.ace.a, envir.corr.ace = env.cor.ace.ce, genetic.cov.ace.a = gen.cov.ace.a, envir.cov.ace.ce = env.cov.ace.ce, genetic.corr.ae = gen.cor.ae.a, envir.corr.ae = env.cor.ae.e, genetic.cov.ae = gen.cov.ae.a, envir.cov.ae = env.cov.ae.e, comparison.anova = comparison, comparison.bic = comp.sub.sort))}




# specifying the analysis - example

# cholesky.trial <- cholfun(data = twindata, first = c("bis_1", "bas_1", "fight_1", "flight_1", "freeze_1"), second = c("bis_2", "bas_2", "fight_2", "flight_2", "freeze_2"), zygo = "zygo")

# results

# cholesky.trial$cholesky.ace # basic info on ACE model fit (for more detailed information, the lavaan function is e.g. summary(cholesky.trial$cholesky.ace, standardized = TRUE, fit.measures = TRUE)
# cholesky.trial$cholesky.ae # the same for the AE model
# cholesky.trial$par.ace # squared parameters from the ACE model; the parameters for each variable, or for the entire table, can be summed up to check if they sum to 1
# cholesky.trial$par.ae # the same for the AE model
# cholesky.trial$genetic.corr.ace # genetic correlations - ACE
# cholesky.trial$genetic.corr.ae # genetic correlations - AE
# cholesky.trial$envir.corr.ace # environmental correlations - ACE
# cholesky.trial$envir.corr.ae #environmental correlations - AE
# cholesky.trial$comparison.anova # model comparison (ACE/AE) - ANOVA
# cholesky.trial$comparison.bic # model comparison (ACE/AE) - BIC only (extracted from the ANOVA object)
