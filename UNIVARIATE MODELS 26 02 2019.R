
# univariate models - preliminary syntax

# for the common pathways models, the function permits to set multiple common factors, but not to precisely determine their loadings 

# in other words, all common factors will load on all variables

# this may be obsolete, but it seems that it could be helpful in determination of the number of non-trivial common factors

# the aforementioned assumption is, of course, to be further explored by simulation studies

# arguments:

# data: dataset (data frame); the function should convert the data into data frame

# first, second: variables for the first / second twin - exact names and order as in the dataset
# the "first" and "second" arguments should be stated as e.g. first = c("V1_tw1", "V2_tw1", "V3_tw1"), second = c("V1_tw2", "V2_tw2", "V3_tw2")

# important: the function provides the results for each variable separately, in form of a list; since 


# zygosity: zygosity variable; preferably, the cases with the missing zygosity data (if any) should be excluded prior to the analysis



# univariate models - preliminary syntax

# for the common pathways models, the function permits to set multiple common factors, but not to precisely determine their loadings 

# in other words, all common factors will load on all variables

# this may be obsolete, but it seems that it could be helpful in determination of the number of non-trivial common factors

# the aforementioned assumption is, of course, to be further explored by simulation studies

# arguments:

# data: dataset (data frame); the function should convert the data into data frame

# first, second: variables for the first / second twin - exact names and order as in the dataset
# the "first" and "second" arguments should be stated as e.g. first = c("V1_tw1", "V2_tw1", "V3_tw1"), second = c("V1_tw2", "V2_tw2", "V3_tw2")

# important: the function provides the results for each variable separately, in form of a list; since 


# zygosity: zygosity variable; preferably, the cases with the missing zygosity data (if any) should be excluded prior to the analysis



# univariate models - preliminary syntax

# for the common pathways models, the function permits to set multiple common factors, but not to precisely determine their loadings 

# in other words, all common factors will load on all variables

# this may be obsolete, but it seems that it could be helpful in determination of the number of non-trivial common factors

# the aforementioned assumption is, of course, to be further explored by simulation studies

# arguments:

# data: dataset (data frame); the function should convert the data into data frame

# first, second: variables for the first / second twin - exact names and order as in the dataset
# the "first" and "second" arguments should be stated as e.g. first = c("V1_tw1", "V2_tw1", "V3_tw1"), second = c("V1_tw2", "V2_tw2", "V3_tw2")

# important: the function provides the results for each variable separately, in form of a list; since 


# zygosity: zygosity variable; preferably, the cases with the missing zygosity data (if any) should be excluded prior to the analysis



univariate <- function(data, first, second, zygo, estimator){
  
  library(lavaan)
  library(semPlot)
  library(stringr)
  library(stringi)
  library(car)
  library(dplyr)
  library(psych)
  library(colr)
  
  
  if (length(first)!=length(second)){
    stop ('Number of variables is different for twin1 and twin2')}
  else
    
    data <- as.data.frame(data)
  
  zygosity <- zygo
  
  
  
  if (missing(estimator)){estimator <- "ML"} else {estimator <- estimator}
  
  zygosity <- subset(data, select = c(zygosity))
  
  colnames(zygosity) <- "zygosity"
  
  str.sub.first <- stri_extract_all_words(first)
  
  str.sub.second <- stri_extract_all_words(second)
  
  nvar <- length(first)
  
  first <- function(nvar){for(i in 1:nvar)
    (cat("V",i,".twin1", sep = "", fill = TRUE))
  }
  
  second <- function(nvar){for(i in 1:nvar)
    (cat("V",i,".twin2", sep = "", fill = TRUE))
  }
  
  firstlabel <- capture.output(first(nvar))
  
  secondlabel <- capture.output(second(nvar))
  
  for(i in 1:nvar){
    labels.first <- (str_sub(str.sub.first))
    labels.second <- (str_sub(str.sub.second))
    twinvar <- subset(data, select = c(labels.first, labels.second))}
  twinvard <- as.data.frame(twinvar)
  twinvardren <- as.data.frame(twinvard)
  colnames(twinvardren) <- c(firstlabel, secondlabel)
  completed <- cbind(twinvardren, zygosity)
  
  
  
  listanaziva <- vector("list")  
  
  for(i in 1:nvar){listanaziva[[i]] <- c(firstlabel[i], secondlabel[i], "zygosity")}
  
  listavarijabli <- vector("list")
  
  for(i in 1:nvar){listavarijabli[[i]] <- subset(completed, select= c(listanaziva[[i]]))}
  
  for(i in 1:nvar){colnames(listavarijabli[[i]]) <- c("V1.twin1", "V1.twin2", "zygosity")}
  
  
  
  
  
  univfunc <- function(data){ 
    univariate.ace <- function(nvar){
      for(i in 1:1)(cat("A.1","=~ c(a",i,",a",i,")*V",i,".twin1", "\n", 
                        "A.2","=~ c(a",i,",a",i,")*V",i,".twin2", "\n",
                        "C.1","=~ c(c",i,",c",i,")*V",i,".twin1", "\n", 
                        "C.2","=~ c(c",i,",c",i,")*V",i,".twin2", "\n",
                        "E.1","=~ c(e",i,",e",i,")*V",i,".twin1", "\n", 
                        "E.2","=~ c(e",i,",e",i,")*V",i,".twin2", "\n",
                        "A.1~~c(1,.5)*A.2", "\n",
                        "C.1~~c(1,1)*C.2", "\n",
                        "E.1~~c(0,0)*E.2", "\n",
                        "V",i,".twin1~~c(0,0)*V",i,".twin1","\n",
                        "V",i,".twin2~~c(0,0)*V",i,".twin2", "\n", 
                        "V",i,".twin1~~c(0,0)*V",i,".twin2", "\n",
                        sep="", fill=FALSE))}
    
    univariate.ace.out <- capture.output(univariate.ace(1))
    
    univariate.ace.mod <- unique(univariate.ace.out)
    
    univariate.ade.mod <- univariate.ace.mod
    
    univariate.ade.mod <- gsub("C", "D", x = univariate.ade.mod)
    
    univariate.ade.mod <- gsub("C", "D", x = univariate.ade.mod)
    
    univariate.ade.mod <- gsub("1,1", "1,.25", x = univariate.ade.mod)
    
    univariate.ae.mod <- univariate.ace.mod[grep("C", x= univariate.ace.mod, invert = TRUE)]
    
    univariate.de.mod <- univariate.ade.mod[grep("A", x= univariate.ade.mod, invert = TRUE)]
    
    univariate.ce.mod <- univariate.ace.mod[grep("A", x= univariate.ace.mod, invert = TRUE)]
    
    univariate.e.mod <- univariate.ce.mod[grep("C", x= univariate.ce.mod, invert = TRUE)]
    
    
    univariate.ace.fit <- cfa(univariate.ace.mod, data = data, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "direct")   
    
    univariate.ade.fit <- cfa(univariate.ade.mod, data = data, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "direct")   
    
    univariate.ae.fit <- cfa(univariate.ae.mod, data = data, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "direct")
    
    univariate.de.fit <- cfa(univariate.de.mod, data = data, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "direct")
    
    univariate.ce.fit <- cfa(univariate.ce.mod, data = data, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "direct")   
    
    univariate.e.fit <- cfa(univariate.e.mod, data = data, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "direct")   
    
    
    # funkcija za ekstrakciju parametara univariate
    
    ext.par.ind <- function(modelfit){
      
      ord <- lavInspect(modelfit, what = "std", add.labels = TRUE)
      
      parameters <- ord$`1`$lambda
      
      parameters.data <- as.data.frame(parameters)
      
      columns <- select(parameters.data, ends_with(".1"))
      
      columns.t.data <- as.data.frame(t(columns))
      
      rows <- select(columns.t.data, ends_with(".twin1"))
      
      rows.squared <- rows^2
      
      rows.squared.data <- as.data.frame(rows.squared)
      
      rows.squared.data.coln <- rows.squared.data
      
      return(rows.squared.data.coln)}
    
    
    par.univariate.ace.fit <- ext.par.ind(univariate.ace.fit)
    
    par.univariate.ae.fit <- ext.par.ind(univariate.ae.fit)
    
    par.univariate.ade.fit <- ext.par.ind(univariate.ade.fit)
    
    par.univariate.de.fit <- ext.par.ind(univariate.de.fit)
    
    par.univariate.ce.fit <- ext.par.ind(univariate.ce.fit)
    
    par.univariate.e.fit <- ext.par.ind(univariate.e.fit)
    
    
    fitted.models.list <- list(univariate.ace.fit, univariate.ae.fit, univariate.ade.fit, univariate.de.fit, univariate.ce.fit, univariate.e.fit)
    
    comparison.constrained <- lapply(FUN = anova, X = fitted.models.list)
    
    names(comparison.constrained) <- c("saturated vs. univariate ACE", "saturated vs. univariate AE", "saturated vs. univariate ADE", "saturated vs. univariate DE", "saturated vs. univariate CE", "saturated vs. univariate E")
    
    comparison <- as.matrix(anova(univariate.ace.fit, univariate.ae.fit, univariate.ade.fit, univariate.de.fit, univariate.ce.fit, univariate.e.fit, method = "satorra.2000"))
    
    comp.sub.bic <- as.matrix(subset(comparison, select=c("BIC")))
    
    colnames(comp.sub.bic) <- c("bic")
    
    comp.sub.mat <- comp.sub.bic
    
    row.names(comp.sub.mat) <- row.names(comparison)
    
    colnames(comp.sub.mat) <- c("bic")
    
    #csdat.names <- as.data.frame(rownames(comp.sub.dat))
    
    #colnames(csdat.names) <- c("models")
    
    #csdatsve <- cbind(csdat.names, comp.sub.dat)
    
    #rownames(csdatsve) <- NULL
    
    comp.sub.sort <- comp.sub.mat[order(comp.sub.mat[,1]),]
    
    comp.sub.sort <- as.data.frame(comp.sub.sort)
    
    colnames(comp.sub.sort) <- c("bic")
    
    bic2 <- as.matrix(c(fitmeasures(univariate.ace.fit, fit.measures = "bic2"), fitmeasures(univariate.ae.fit, fit.measures = "bic2"), fitmeasures(univariate.ade.fit, fit.measures = "bic2"), fitmeasures(univariate.de.fit, fit.measures = "bic2"), fitmeasures(univariate.ce.fit, fit.measures = "bic2"), fitmeasures(univariate.e.fit, fit.measures = "bic2")))
    
    bic2mat <- as.matrix(bic2)
    
    bic2dat <- as.data.frame(bic2mat, row.names = c("univariate.ACE", "univariate.AE", "univariate.ADE", "univariate.DE", "univariate.CE", "univariate.E"), colnames = "bic2")
    
    bic2datnames <- as.data.frame(rownames(bic2dat))
    
    colnames(bic2datnames) <- c("models")
    
    bic2datsve <- cbind(bic2datnames, bic2dat)
    
    colnames(bic2datsve) <- c("models", "bic2")
    
    rownames(bic2datsve) <- NULL
    
    
    rownames(bic2mat) <- c("univariate.ACE", "univariate.AE", "univariate.ADE", "univariate.DE", "univariate.CE", "univariate.E")
    
    colnames(bic2mat) <- c("bic2")
    
    
    
    comp.sub.sort.2 <- as.data.frame(bic2datsve[order(bic2datsve$bic2),])
    
    
    return(list(names = colnames(twinvardren), orig = as.data.frame(twinvar), renamed = twinvardren, complete.matrix = completed, uni.ace = univariate.ace.fit, uni.ae = univariate.ae.fit, uni.ade = univariate.ade.fit, uni.de = univariate.de.fit,
                uni.ce = univariate.ce.fit, uni.e = univariate.e.fit, constrained = comparison.constrained, comparison.anova = comparison, comparison.bic = comp.sub.sort, comparison.bic2 = comp.sub.sort.2, par.uni.ace = par.univariate.ace.fit, par.uni.ae = par.univariate.ae.fit, par.uni.ade = par.univariate.ade.fit, par.uni.de = par.univariate.de.fit, par.uni.ce = par.univariate.ce.fit, par.uni.e = par.univariate.e.fit))}
  
  rezultat <- lapply(FUN = univfunc, X = listavarijabli)
  
  names(rezultat) <- str.sub.first
  
  return(rezultat)}


# example - specifying analysis parameters:



#univariate.fit <- univariate(data = uop29_data, first = c("bis_parc", "bas_parc", "bor_parc", "bez_parc", "blo_parc"), second = c("bis_parc_b2", "bas_parc_b2", "bor_parc_b2", "bez_parc_b2", "blo_parc_b2"), zygo = "DNK_zigotnost", estimator = "MLM")



# values / outcomes (extracted from the object by object$value; can of course be saved into an object in its own right and included in further analyses):

# names = variable names, 

# orig = original variable names

# renamed = renamed variables 

# complete.matrix = matrix used in the algorithm, with renamed variables

# uni.ace (...ae, ade, de, ce, e) = basic information on the fit of the independent pathways models

# com.ace (...ae, ade, de, ce, e) = basic information on the fit of the common pathways model 

# par.uni.ace (...ae, ade, de, ce, e) = squared parameters of the independent pathways models 

# since the outcome is a list, the results are obtained the following way: object$variable$result (e.g. twinanalysis$agreeableness$par.uni.ace)





# example - specifying analysis parameters:



#univariate.fit <- univariate(data = uop29_data, first = c("bis_parc", "bas_parc", "bor_parc", "bez_parc", "blo_parc"), second = c("bis_parc_b2", "bas_parc_b2", "bor_parc_b2", "bez_parc_b2", "blo_parc_b2"), zygo = "DNK_zigotnost", estimator = "MLM")



# values / outcomes (extracted from the object by object$value; can of course be saved into an object in its own right and included in further analyses):

# names = variable names, 

# orig = original variable names

# renamed = renamed variables 

# complete.matrix = matrix used in the algorithm, with renamed variables

# uni.ace (...ae, ade, de, ce, e) = basic information on the fit of the independent pathways models

# com.ace (...ae, ade, de, ce, e) = basic information on the fit of the common pathways model 

# par.uni.ace (...ae, ade, de, ce, e) = squared parameters of the independent pathways models 

# since the outcome is a list, the results are obtained the following way: object$variable$result (e.g. twinanalysis$agreeableness$par.uni.ace)







# values / outcomes (extracted from the object by object$value; can of course be saved into an object in its own right and included in further analyses):

# names = variable names, 

# orig = original variable names

# renamed = renamed variables 

# complete.matrix = matrix used in the algorithm, with renamed variables

# uni.ace (...ae, ade, de, ce, e) = basic information on the fit of the independent pathways models

# com.ace (...ae, ade, de, ce, e) = basic information on the fit of the common pathways model 

# par.uni.ace (...ae, ade, de, ce, e) = squared parameters of the independent pathways models 

# since the outcome is a list, the results are obtained the following way: object$variable$result (e.g. twinanalysis$agreeableness$par.uni.ace)