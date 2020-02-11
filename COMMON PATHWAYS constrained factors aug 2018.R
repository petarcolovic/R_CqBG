# common pathways - "constrained" models with specified factor loadings - preliminary syntax

# arguments: for the most part the same as in other functions, with an important differences: factor loadings are specified as e.g. f1 = c(variables), f2 = c(variables), etc. 

# prior to specification of factor loadings, it is important to enter the number of factors as the "nfaccom" argument


cspecif <- function(data, first, second, nfaccom, zygo,...){
 
  library(lavaan)
  library(semPlot)
  library(stringr)
  library(stringi)
  library(car)
  library(dplyr)
  library(psych)
  
  if (length(first)!=length(second)){
    stop ('Number of variables is different for twin1 and twin2')}
  else
    
  # arguments
    
  data <- as.data.frame(data)
  
  zygosity <- zygo
  
  zygosity <- subset(data, select = c(zygosity))
  
  colnames(zygosity) <- "zygosity"
  
  nfaccom = nfaccom
  
  nvar = length(first)
  
  str.sub.first <- stri_extract_all_words(first)
  
  class(str.sub.first)
  
  str.sub.second <- stri_extract_all_words(second)
  
  class(str.sub.second)
  
  firsttwin <- function(nvar){for(i in 1:nvar)
    (cat("V",i,".twin1", sep = "", fill = TRUE))
  }
  
  secondtwin <- function(nvar){for(i in 1:nvar)
    (cat("V",i,".twin2", sep = "", fill = TRUE))
  }
  
  firstlabels <- capture.output(firsttwin(nvar))
  
  secondlabels <- capture.output(secondtwin(nvar))
  
  for(i in 1:nvar){
    labels.first <- (str_sub(str.sub.first))
    labels.second <- (str_sub(str.sub.second))
    twinvar <- subset(data, select = c(labels.first, labels.second))}
  twinvard <- as.data.frame(twinvar)
  twinvardren <- as.data.frame(twinvard)
  colnames(twinvardren) <- c(firstlabels, secondlabels)
  completed <- cbind(twinvardren, zygosity) 
  
  arguments.list <- list(...) 
  
  length(arguments.list[[1]])
  
  first.list <- list(first)
  
  first.vect <- as.vector(unlist(first.list))
  
  variablenames <- list()
  
  for(i in 1:length(first)){variablenames[[i]]<- paste(str_extract_all(unlist(strsplit(first[[i]], split = "")), unlist(strsplit(second[[i]], split = "")), simplify = TRUE), collapse = "")}
  
  variablenames <- unlist(variablenames)
  
  intersection.list <- vector("list")
  
  for(i in 1:length(arguments.list)){intersection.list[[i]]<- grep(paste(arguments.list[[i]], collapse = "|"), first.vect)}
  
  intersection.array <- array(intersection.list, dim = length(intersection.list))
  
  factorloadings <- function(array){for(i in 1:length(array)){for (j in 1:length(array[[i]])) 
    cat("F",i, ".twin1","=~", "c(l",",l",")*V", as.vector((array[[i]])[j]),".twin1", "\n",
        sep="", fill = FALSE)}} 
  
  loads <- capture.output(factorloadings(intersection.array)) 
  
  loads.pure <- loads[grep("NA", x= loads, invert = TRUE)]
  
  number <- length(loads.pure)
  
  vector.number <- 1:number
  
  vector.letter <- 1:number
  
  vector.letter[] <- "l"
  
  vector.comma <- 1:number
  
  vector.comma[] <-","
  
  half <- paste(vector.letter, vector.number, vector.comma, vector.letter, vector.number, sep="", collapse = NULL)
  
  loads.repl.p <- str_replace_all(string = loads.pure, pattern = "l,l", replacement=half)
  
  loads.repl.p.2 <- str_replace_all(string = loads.repl.p, pattern = "twin1", replacement = "twin2")
  
  Ftwin1 <- loads.repl.p # raw
  
  Ftwin2 <- loads.repl.p.2 #raw
  
  Fboth <- c(Ftwin1, Ftwin2) 
  
  common.variables <- function(nvar){for (j in 1:nvar)(cat("As",j,".1=~ c(a",j,j,",a",j,j,")*V",j,".twin1","\n",
                                                               "As",j,".2=~ c(a",j,j,",a",j,j,")*V",j,".twin2","\n",
                                                               "Cs",j,".1=~ c(c",j,j,",c",j,j,")*V",j,".twin1","\n",
                                                               "Cs",j,".2=~ c(c",j,j,",c",j,j,")*V",j,".twin2","\n",
                                                               "Es",j,".1=~ c(e",j,j,",e",j,j,")*V",j,".twin1","\n",
                                                               "Es",j,".2=~ c(e",j,j,",e",j,j,")*V",j,".twin2","\n",
                                                               
                                                               
                                                               "As",j,".1","~~c(1,.5)*As",j,".2", "\n",
                                                               "Cs",j,".1","~~c(1,1)*Cs",j,".2", "\n",
                                                               "Es",j,".1","~~c(0,0)*Es",j,".2", "\n",
                                                               "V",j,".twin1~~c(0,0)*V",j,".twin1","\n",
                                                               "V",j,".twin2~~c(0,0)*V",j,".twin2", "\n", 
                                                               
                                                               sep="", fill=FALSE))}
  
  common.variables.out <- capture.output(common.variables(nvar))
  
  
  
  common.factors <- function(nfaccom){for (i in 1:nfaccom)(cat("Aop",i,"twin.1","=~","c(a",i,",","a",i,")","*F",i,".twin1","\n",
                                                                   "Cop",i,"twin.1","=~","c(c",i,",","c",i,")","*F",i,".twin1","\n",
                                                                   "Eop",i,"twin.1","=~","c(e",i,",","e",i,")","*F",i,".twin1","\n",
                                                                   
                                                                   "Aop",i,"twin.2","=~","c(a",i,",","a",i,")","*F",i,".twin2","\n",
                                                                   "Cop",i,"twin.2","=~","c(c",i,",","c",i,")","*F",i,".twin2","\n",
                                                                   "Eop",i,"twin.2","=~","c(e",i,",","e",i,")","*F",i,".twin2","\n",
                                                                   
                                                                   "Aop",i,"twin.1","~~c(1,.5)*Aop",i,"twin.2", "\n",
                                                                   "Cop",i,"twin.1","~~c(1,1)*Cop",i,"twin.2", "\n",
                                                                   "Eop",i,"twin.1","~~c(0,0)*Eop",i,"twin.2", "\n",
                                                                   
                                                                   "F",i,".twin1","~~c(0,0)*","F",i,".twin1","\n",
                                                                   "F",i,".twin2","~~c(0,0)*","F",i,".twin2","\n",
                                                                   sep="", fill=FALSE))}
  # model formulation
  
  common.factors.out <- capture.output(common.factors(nfaccom))
  
  cspec.ace.model.mat <- c(common.factors.out, common.variables.out, Fboth)
  
  cspecif.ace.mod <- unique(cspec.ace.model.mat) #
  
  
  cspecif.ade.mod <- cspecif.ace.mod #
  
  cspecif.ade.mod <- gsub("Cop", "Dop", x = cspecif.ade.mod)
  
  cspecif.ade.mod <- gsub("Cs", "Ds", x = cspecif.ade.mod)
  
  cspecif.ade.mod <- gsub("1,1", "1,.25", x = cspecif.ade.mod)
  
  cspecif.ae.mod <- cspecif.ace.mod[grep("Cop|Cs", x= cspecif.ace.mod, invert = TRUE)]
  
  cspecif.de.mod <- cspecif.ade.mod[grep("Aop|As", x= cspecif.ade.mod, invert = TRUE)]
  
  cspecif.ce.mod <- cspecif.ace.mod[grep("Aop|As", x= cspecif.ace.mod, invert = TRUE)]
  
  cspecif.e.mod <- cspecif.ce.mod[grep("Cop|Cs", x= cspecif.ce.mod, invert = TRUE)]
  
  
  
  
  # model fitting
  
  cspecif.ace.fit <- cfa(cspecif.ace.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = FALSE, estimator = "ML")   
  
  cspecif.ade.fit <- cfa(cspecif.ade.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = FALSE, estimator = "ML")   
  
  cspecif.ae.fit <- cfa(cspecif.ae.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = FALSE, estimator = "ML")
  
  cspecif.de.fit <- cfa(cspecif.de.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = FALSE, estimator = "ML")
  
  cspecif.ce.fit <- cfa(cspecif.ce.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = FALSE, estimator = "ML")   
  
  cspecif.e.fit <- cfa(cspecif.e.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = FALSE, estimator = "ML")   
  
  # parameters
  
  # function for the extraction of parameters - common
  
  
  ext.par.com <- function(model){
    
    ord <- lavInspect(model, what = "std", add.labels = TRUE)
    
    parameters <- ord$`1`$lambda
    
    factors <- ord$`1`$beta
    
    factors.transp <- t(factors)
    
    factors.transp.sel <- factors.transp['^F']
    
    factors.transp.data <- as.data.frame(factors)
    
    factors.columns <- select(as.data.frame(factors.transp), starts_with("F")) 
    
    round(factors.columns,3)
    
    factors.columns.data <- as.data.frame(factors.columns)
    
    factors.sel <- select(factors.columns.data, ends_with("twin1"))
    
    factors.rows <- t(factors.sel)
    
    factors.rows.transp.sel <- select(as.data.frame(factors.rows), ends_with("1"))
    
    factors.sel.mat <- t(factors.rows.transp.sel)
    
    factors.sel.mat.data <- as.data.frame(factors.sel.mat)
    
    
    parameters.data <- as.data.frame(parameters)
    
    columns.parameters <- select(parameters.data, ends_with("1"))
    
    columns.parameters.t.data <- as.data.frame(t(columns.parameters))
    
    columns.parameters.data.colnames <- columns.parameters.t.data
    
    rows <- select(columns.parameters.t.data, ends_with("twin1"))
    
    rows.col <- rows
    
    rows.factors <- cbind(rows.col, factors.sel.mat.data)
    
    result <- rows.factors
    
    result.rn <- tibble::rownames_to_column(result) # labels of the variables from the first columns to row names
    
    result.rn.data <- as.data.frame(result.rn) # conversion to data frame
    
    result.sub1 <- result.rn.data[grep("^F", result.rn.data$rowname),] # extraction of a matrix with factors as rows
    
    result.sub2 <- result.rn.data[grep("^Aop|^Cop|^Eop|^Dop", result.rn.data$rowname),] # extraction of a matrix with common genetic factors
    
    result.factors <- as.data.frame(result.sub1) # conversion to data frame
    
    result.gen <- as.data.frame(result.sub2) # conversion to data frame
    
    columns.gen <- colnames(result.gen) # labels of the columns in matrix with common genetic factors as rows
    
    columns.fact <- colnames(result.factors) # labels of the columns in the matrix with factors as rows
    
    columns.gen.sub <- columns.gen[grep((glob2rx("*F*")), x = columns.gen, invert = FALSE)] # characters: labels of genetic and environmental effects
    
    columns.fact.sub <- columns.fact[grep((glob2rx("*V*")), x = columns.fact, invert = FALSE)] # extraction of columns, a.k.a. observed variables from the object where factors are rows
    
    result.gen.factors.sub <- subset(result.gen, select = c("rowname", columns.gen.sub)) # extraction of genetic and environmental common effects
    
    result.factors.sub <- subset(result.factors, select = c("rowname", columns.fact.sub)) # extraction of factors
    
    
    # multiplication
    
    
    genetic <- as.data.frame(result.gen.factors.sub, optional = FALSE) # data frame sa opstim effectsma
    
    row.names(genetic) <- NULL # deleting row names
    
    genetic <- tibble::column_to_rownames(genetic) # first column as row names
    
    factorial <- result.factors.sub # data frame for factor loadings
    
    row.names(factorial) <- NULL # deleting row names
    
    factorial <- tibble::column_to_rownames(factorial) # first column as rownames
    
    contributions = as.matrix(genetic)%*%as.matrix(factorial) # multiplication of matrices with genetic and environmental influences by factor loadings
    
    # final result
    
    common.squared <- contributions^2
    
    specific <- result.rn.data[grep("^As|^Cs|^Es|^Ds", result.rn.data$rowname),] # matrix with genetic influences as rows
    
    specific.colnames <- colnames(specific)
    
    columns.specific.sub <- specific.colnames[grep((glob2rx("*F*")), x = specific.colnames, invert = TRUE)] # characters: labels of genetic and environmental effects
    
    specific.sub <- subset(specific, select = c(columns.specific.sub))
    
    row.names(specific.sub) <- NULL
    
    specific.sub <- tibble::column_to_rownames(specific.sub)
    
    specific.squared <- as.matrix(specific.sub)^2
    
    squared.effects <- rbind(common.squared, specific.squared)
    
    loadings <- factorial
    
    all.together <- rbind(squared.effects, loadings)
    
    colnames(squared.effects) <- variablenames
    
    colnames(loadings) <- variablenames
    
    colnames(all.together) <- variablenames
    
    return(list(effects = squared.effects, factor_loadings = loadings, complete = all.together))}
  
 
# parameter calculation
  
  par.common.ace.fit <- ext.par.com(cspecif.ace.fit)
  
  par.common.ae.fit <- ext.par.com(cspecif.ae.fit)
  
  par.common.ade.fit <- ext.par.com(cspecif.ade.fit)
  
  par.common.de.fit <- ext.par.com(cspecif.de.fit)
  
  par.common.ce.fit <- ext.par.com(cspecif.ce.fit)
  
  par.common.e.fit <- ext.par.com(cspecif.e.fit)
  
  # model comparison
  
  comparison <- anova(cspecif.ace.fit, cspecif.ae.fit, cspecif.ade.fit, cspecif.de.fit, cspecif.ce.fit, cspecif.e.fit)
  
  comp.sub <- subset(comparison, select=c("BIC"))
  
  comp.sub.sort <- sort(comp.sub, decreasing = FALSE)
  
  
  return(list(com.ace = cspecif.ace.fit, com.ae = cspecif.ae.fit, com.ade = cspecif.ade.fit, com.de = cspecif.de.fit,
              com.ce = cspecif.ce.fit, compared = comparison, order.fit = comp.sub.sort, par.c.ace = par.common.ace.fit, par.c.ae = par.common.ae.fit, par.c.ade = par.common.ade.fit, par.c.de = par.common.de.fit, par.c.ce = par.common.ce.fit, par.c.e = par.common.e.fit))}
  
# outputs:

# com.ace (ae, ade, de, ce, e) = basic information on model fit

# par.c.ace (ae, ade, de, ce, e) = parameters (all tables in a single object, can be easily extracted from it)

# compared = model comparison

# order.fit = models ordered by the BIC size (ascending)

###test:

#constr.trial <- cspecif(data = twindata, first = c("bis_1", "bas_1", "fight_1", "flight_1", "freeze_1"), second = c("bis_2", "bas_2", "fight_2", "flight_2", "freeze_2"), zygo = "zygo", nfaccom = 2, f1 = c("bas", "fight", "flight"), f2 = c("bis", "flight", "freeze"))  

#constr.trial$par.c.ae
