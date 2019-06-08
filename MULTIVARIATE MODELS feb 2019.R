# multivariate models: independent and common pathways - preliminary syntax

# the script provides estimates of standard multivariate twin models: independent pathways model, and common pathways model; both contain common and specific genetic and environmental factors, however in common pathways model there is an assumption that all observed variables fall under the same latent construct

# for the common pathways models, the function permits to set multiple common factors, but not to precisely determine their loadings 

# in other words, all common factors will load on all variables

    # this may be obsolete, but it seems that it could be helpful in determination of the number of non-trivial common factors

    # the aforementioned assumption is, of course, to be further explored by simulation studies

# arguments:

# data: dataset (data frame); the function should convert the data into data frame

# first, second: variables for the first / second twin - exact names and order as in the dataset
# the "first" and "second" arguments should be stated as e.g. first = c("V1_tw1", "V2_tw1", "V3_tw1"), second = c("V1_tw2", "V2_tw2", "V3_tw2")

# zygosity: zygosity variable; preferably, the cases with the missing zygosity data (if any) should be excluded prior to the analysis

# nfaccom: number of common factors in the common pathways models

multiv <- function(data, first, second, zygo, nfaccom){

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

    data <- as.data.frame(data)
  
  zygosity <- zygo
  
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
  
  independent.ace <- function(nvar){
    for(i in 1:nvar)(cat("Acom.1","=~ c(a",i,",a",i,")*V",i,".twin1", "\n", 
                         "Acom.2","=~ c(a",i,",a",i,")*V",i,".twin2", "\n",
                         "Ccom.1","=~ c(c",i,",c",i,")*V",i,".twin1", "\n", 
                         "Ccom.2","=~ c(c",i,",c",i,")*V",i,".twin2", "\n",
                         "Ecom.1","=~ c(e",i,",e",i,")*V",i,".twin1", "\n", 
                         "Ecom.2","=~ c(e",i,",e",i,")*V",i,".twin2", "\n",
                         "As",i,".1=~ c(a",i,i,",a",i,i,")*V",i,".twin1","\n",
                         "As",i,".2=~ c(a",i,i,",a",i,i,")*V",i,".twin2","\n",
                         "Cs",i,".1=~ c(c",i,i,",c",i,i,")*V",i,".twin1","\n",
                         "Cs",i,".2=~ c(c",i,i,",c",i,i,")*V",i,".twin2","\n",
                         "Es",i,".1=~ c(e",i,i,",e",i,i,")*V",i,".twin1","\n",
                         "Es",i,".2=~ c(e",i,i,",e",i,i,")*V",i,".twin2","\n",
                         "Acom.1~~c(1,.5)*Acom.2", "\n",
                         "Ccom.1~~c(1,1)*Ccom.2", "\n",
                         "Ecom.1~~c(0,0)*Ecom.2", "\n",
                         "As",i,".1","~~c(1,.5)*As",i,".2", "\n",
                         "Cs",i,".1","~~c(1,1)*Cs",i,".2", "\n",
                         "Es",i,".1","~~c(0,0)*Es",i,".2", "\n",
                         "V",i,".twin1~~c(0,0)*V",i,".twin1","\n",
                         "V",i,".twin2~~c(0,0)*V",i,".twin2", "\n", sep="", fill=FALSE))}
    
  independent.ace.out <- capture.output(independent.ace(nvar))
  
  independent.ace.mod <- unique(independent.ace.out)
  
  independent.ade.mod <- independent.ace.mod
  
  independent.ade.mod <- gsub("Ccom", "Dop", x = independent.ade.mod)
  
  independent.ade.mod <- gsub("Cs", "Ds", x = independent.ade.mod)
  
  independent.ade.mod <- gsub("1,1", "1,.25", x = independent.ade.mod)
  
  independent.ae.mod <- independent.ace.mod[grep("Ccom|Cs", x= independent.ace.mod, invert = TRUE)]
  
  independent.de.mod <- independent.ade.mod[grep("Acom|As", x= independent.ade.mod, invert = TRUE)]
  
  independent.ce.mod <- independent.ace.mod[grep("Acom|As", x= independent.ace.mod, invert = TRUE)]
  
  independent.e.mod <- independent.ce.mod[grep("Ccom|Cs", x= independent.ce.mod, invert = TRUE)]
  
  
  independent.ace.fit <- cfa(independent.ace.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")   
  
  independent.ade.fit <- cfa(independent.ade.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")   
  
  independent.ae.fit <- cfa(independent.ae.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")
  
  independent.de.fit <- cfa(independent.de.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")
  
  independent.ce.fit <- cfa(independent.ce.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")   
  
  independent.e.fit <- cfa(independent.e.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")   
  
  
  common.factors <- function(nfaccom, nvar){
      for(i in 1:nfaccom){for (j in 1:nvar)(cat("Acom",i,"twin.1","=~","c(a",i,",","a",i,")","*F",i,".twin1","\n",
                                                "Ccom",i,"twin.1","=~","c(c",i,",","c",i,")","*F",i,".twin1","\n",
                                                "Ecom",i,"twin.1","=~","c(e",i,",","e",i,")","*F",i,".twin1","\n",
                                                
                                                "F",i,".twin1","=~","c(l",j,",","l",j,")*V",j,".twin1", "\n",
                                                
                                                "Acom",i,"twin.2","=~","c(a",i,",","a",i,")","*F",i,".twin2","\n",
                                                "Ccom",i,"twin.2","=~","c(c",i,",","c",i,")","*F",i,".twin2","\n",
                                                "Ecom",i,"twin.2","=~","c(e",i,",","e",i,")","*F",i,".twin2","\n",
                                               
                                                "F",i,".twin2","=~","c(l",j,",","l",j,")*V",j,".twin2", "\n",
                                                
                                                "As",j,".1=~ c(a",j,j,",a",j,j,")*V",j,".twin1","\n",
                                                "As",j,".2=~ c(a",j,j,",a",j,j,")*V",j,".twin2","\n",
                                                "Cs",j,".1=~ c(c",j,j,",c",j,j,")*V",j,".twin1","\n",
                                                "Cs",j,".2=~ c(c",j,j,",c",j,j,")*V",j,".twin2","\n",
                                                "Es",j,".1=~ c(e",j,j,",e",j,j,")*V",j,".twin1","\n",
                                                "Es",j,".2=~ c(e",j,j,",e",j,j,")*V",j,".twin2","\n",
                                                
                                                "Acom",i,"twin.1","~~c(1,.5)*Acom",i,"twin.2", "\n",
                                                "Ccom",i,"twin.1","~~c(1,1)*Ccom",i,"twin.2", "\n",
                                                "Ecom",i,"twin.1","~~c(0,0)*Ecom",i,"twin.2", "\n",
                                                
                                                "As",j,".1","~~c(1,.5)*As",j,".2", "\n",
                                                "Cs",j,".1","~~c(1,1)*Cs",j,".2", "\n",
                                                "Es",j,".1","~~c(0,0)*Es",j,".2", "\n",
                                                "V",j,".twin1~~c(0,0)*V",j,".twin1","\n",
                                                "V",j,".twin2~~c(0,0)*V",j,".twin2", "\n", 
                                                "F",i,".twin1","~~c(0,0)*","F",i,".twin1","\n",
                                                "F",i,".twin2","~~c(0,0)*","F",i,".twin2","\n",
                                                sep="", fill=FALSE))}}
    
    
    common.out <- capture.output(common.factors(nfaccom, nvar))
    
    common.ace.mod <- unique(common.out)
    
    common.ade.mod <- common.ace.mod
    
    common.ade.mod <- gsub("Ccom", "Dop", x = common.ade.mod)
    
    common.ade.mod <- gsub("Cs", "Ds", x = common.ade.mod)
    
    common.ade.mod <- gsub("1,1", "1,.25", x = common.ade.mod)
    
    common.ae.mod <- common.ace.mod[grep("Ccom|Cs", x= common.ace.mod, invert = TRUE)]
    
    common.de.mod <- common.ade.mod[grep("Acom|As", x= common.ade.mod, invert = TRUE)]
    
    common.ce.mod <- common.ace.mod[grep("Acom|As", x= common.ace.mod, invert = TRUE)]
    
    common.e.mod <- common.ce.mod[grep("Ccom|Cs", x= common.ce.mod, invert = TRUE)]
    
    
    
    common.ace.fit <- cfa(common.ace.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")   
    
    common.ade.fit <- cfa(common.ade.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")   
    
    common.ae.fit <- cfa(common.ae.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")
    
    common.de.fit <- cfa(common.de.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")
    
    common.ce.fit <- cfa(common.ce.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")   
    
    common.e.fit <- cfa(common.e.mod, data = completed, group = "zygosity", orthogonal = TRUE, std.lv = TRUE, missing = "fiml")   
    
 

   # funkcija za ekstrakciju parametara independent
    
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
      
      colnames(rows.squared.data.coln) <- (str.sub.first)
      
      return(rows.squared.data.coln)}
    
   
# parameters: common pathway
    
    
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
      
      result.rn <- tibble::rownames_to_column(result) # labels varijabli iz prve columns u nazive redova
      
      result.rn.data <- as.data.frame(result.rn) # konverzija u data frame
      
      result.sub1 <- result.rn.data[grep("^F", result.rn.data$rowname),] # izdvajanje matrice s factorsma kao rowsma
      
      result.sub2 <- result.rn.data[grep("^Acom|^Ccom|^Ecom|^Dop", result.rn.data$rowname),] # izdvajanje matrice s opstim geneticm uticajima kao rowsma
      
      result.factors <- as.data.frame(result.sub1) # konverzija u data frame
      
      result.gen <- as.data.frame(result.sub2) # konverzija u data frame
      
      columns.gen <- colnames(result.gen) # labels kolona u matrici s opstim efektima kao rowsma
      
      columns.fakt <- colnames(result.factors) # labels kolona u matrici s factorsma kao rowsma
      
      columns.gen.sub <- columns.gen[grep((glob2rx("*F*")), x = columns.gen, invert = FALSE)] # karakteri: labels genetich i sredinskih efekata
      
      columns.fakt.sub <- columns.fakt[grep((glob2rx("*V*")), x = columns.fakt, invert = FALSE)] # izdvajanje kolona, odnosno manifestnih varijabli iz objekta gde su factors rows
      
      result.gen.factors.sub <- subset(result.gen, select = c("rowname", columns.gen.sub)) # izdvajanje genetich i sredinskih opstih efekata
      
      result.factors.sub <- subset(result.factors, select = c("rowname", columns.fakt.sub)) # izdvajanje faktora
      
      
      # multiplication (at last!)
      
      
      genetic <- as.data.frame(result.gen.factors.sub, optional = FALSE) # data frame sa opstim efektima
      
      row.names(genetic) <- NULL # brisanje imena redova
      
      genetic <- tibble::column_to_rownames(genetic) # prva kolona kao naziv reda
      
      factorloadings <- result.factors.sub # data frame za faktorska loadings
      
      row.names(factorloadings) <- NULL # brisanje naziva redova
      
      factorloadings <- tibble::column_to_rownames(factorloadings) # prva kolona kao naziv reda
      
      contributions = as.matrix(genetic)%*%as.matrix(factorloadings) # mnozenje matrica sa geneticm i sredinskim contributionsma
      
      # final result
      
      common.squared <- contributions^2
      
      specific <- result.rn.data[grep("^As|^Cs|^Es|^Ds", result.rn.data$rowname),] # izdvajanje matrice sa specificm geneticm uticajima kao rowsma
      
      specific.colnames <- colnames(specific)
      
      columns.specific.sub <- specific.colnames[grep((glob2rx("*F*")), x = specific.colnames, invert = TRUE)] # karakteri: labels genetich i sredinskih efekata
      
      specific.sub <- subset(specific, select = c(columns.specific.sub))
      
      row.names(specific.sub) <- NULL
      
      specific.sub <- tibble::column_to_rownames(specific.sub)
      
      specific.squared <- as.matrix(specific.sub)^2
      
      squared.effects <- rbind(common.squared, specific.squared) # squared effects, comprising general genetic and environmental effects multiplied by factor loadings, as well as the effects of specific factors
      
      loadings <- factorloadings
      
      all.together <- rbind(squared.effects, loadings) # both squared effects and factor loadings (not appropriate for summing-up, because the impact of factor loadings has already been counted in the common factor effects, so the sum would be larger than 1)
      
      colnames(squared.effects) <- str.sub.first
      
      colnames(loadings) <- str.sub.first
      
      colnames(all.together) <- str.sub.first
      
      return(list(effects = squared.effects, load = loadings, all = all.together))}
    
    
    
    par.independent.ace.fit <- ext.par.ind(independent.ace.fit)
    
    par.independent.ae.fit <- ext.par.ind(independent.ae.fit)
    
    par.independent.ade.fit <- ext.par.ind(independent.ade.fit)
    
    par.independent.de.fit <- ext.par.ind(independent.de.fit)
    
    par.independent.ce.fit <- ext.par.ind(independent.ce.fit)
    
    par.independent.e.fit <- ext.par.ind(independent.e.fit)
    
    
    
    par.common.ace.fit <- ext.par.com(common.ace.fit)
    
    par.common.ae.fit <- ext.par.com(common.ae.fit)
    
    par.common.ade.fit <- ext.par.com(common.ade.fit)
    
    par.common.de.fit <- ext.par.com(common.de.fit)
    
    par.common.ce.fit <- ext.par.com(common.ce.fit)
    
    par.common.e.fit <- ext.par.com(common.e.fit)
    
    
    
    comparison <- as.matrix(anova(independent.ace.fit, independent.ae.fit, independent.ade.fit, independent.de.fit, independent.ce.fit, independent.e.fit, common.ace.fit, common.ae.fit, common.ade.fit, common.de.fit, common.ce.fit, common.e.fit))
    
    comp.sub <- as.matrix(subset(comparison, select=c("BIC")))
    
    comp.sub.sort <- as.matrix(sort(comp.sub, decreasing = FALSE))
    
  return(list(names = colnames(twinvardren), orig = as.data.frame(twinvar), renamed = twinvardren, complete.matrix = completed, ind.ace = independent.ace.fit, ind.ae = independent.ae.fit, ind.ade = independent.ade.fit, ind.de = independent.de.fit,
              ind.ce = independent.ce.fit, ind.e = independent.e.fit, com.ace = common.ace.fit, com.ae = common.ae.fit, com.ade = common.ade.fit, com.de = common.de.fit,
              com.ce = common.ce.fit, comparison = comparison, ordered.fit = comp.sub.sort, par.ind.ace = par.independent.ace.fit, par.ind.ae = par.independent.ae.fit, par.ind.ade = par.independent.ade.fit, par.ind.de = par.independent.de.fit, par.ind.ce = par.independent.ce.fit, par.ind.e = par.independent.e.fit, par.com.ace = par.common.ace.fit$effects, par.com.ae = par.common.ae.fit$effects, par.com.ade = par.common.ade.fit$effects, par.com.de = par.common.de.fit$effects, par.com.ce = par.common.ce.fit$effects, par.com.e = par.common.e.fit$effects, loadings.com.ace = par.common.ace.fit$load, loadings.com.ae = par.common.ae.fit$load, loadings.com.ade = par.common.ade.fit$load, loadings.com.de = par.common.ace.fit$load, loadings.com.ce = par.common.ce.fit$load, loadings.com.e = par.common.e.fit$load, common.all = par.common.ace.fit$all, common.all.ace = par.common.ace.fit$all, common.all.ae = par.common.ae.fit$all, common.all.ade = par.common.ade.fit$all, common.all.de = par.common.de.fit$all, common.all.ce = par.common.ce.fit$all, common.all.e = par.common.e.fit$all))}

# values / outcomes (extracted from the object by object$value; can of course be saved into an object in its own right and included in further analyses):

# names = variable names, 

# orig = original variable names

# renamed = renamed variables 

# complete.matrix = matrix used in the algorithm, with renamed variables

# ind.ace (...ae, ade, de, ce, e) = basic information on the fit of the independent pathways models

# com.ace (...ae, ade, de, ce, e) = basic information on the fit of the common pathways model 

# par.ind.ace (...ae, ade, de, ce, e) = squared parameters of the independent pathways models 

# par.com.ace (...ae, ade, de, ce, e) = squared parameters of the common pathways models (common factor effects: genetic / environmental effects are multiplied by factor loadings and squared; thus factor loadings are non included here)

# loadings.com.ace (...ae, ade, de, ce, e) = factor loadings for the common pathways models

# common.all.ace (...ae, ade, de, ce, e) = table with both squared parameters and factor loadings


# trying out the function

# multiv.trial <- multiv(data = twindata, first = c("bis_1", "bas_1", "fight_1", "flight_1", "freeze_1"), second = c("bis_2", "bas_2", "fight_2", "flight_2", "freeze_2"), zygo = "zygo", nfaccom = 3)


