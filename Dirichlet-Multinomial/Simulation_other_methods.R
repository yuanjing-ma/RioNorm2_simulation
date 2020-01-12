## Apply other DA-tests on Dirichlet-Multinomial based simulation data 
## e.g. DESeq, DESeq2, metagenomeSeq, Omnibus, RAIDA, ANCOM
## Simulation data available upon request

## load in required package
reqpkg = c("cluster", "doParallel", "edgeR", "DESeq", "DESeq2", "stats",
           "foreach", "ggplot2", "grid", "scales", "metagenomeSeq", 
           "phyloseq", "plyr", "reshape2", "ROCR", "RioNorm2", 
           "igraph", "pscl", "MASS", "exactRankTests", "nlme")
for (i in reqpkg) {
  print(i)
  print(packageVersion(i))
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}


## load in simulated dataset
setwd("~/Documents/Microbiome/Corrected_simulation_DirMulti")
load("simulation-data-dirmulti-Mock-corrected.RData")

## simulation parameters
data("GlobalPatterns")
comdelim = "_"
simparams = names(simlist0)
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")

## register for parallel computing
library(foreach)
library(doParallel)
Ncores = detectCores() - 1
cl <- makeCluster(Ncores)
registerDoParallel(cl)


############ RAIDA
library(RAIDA)
RAIDA_test <- function(OTU_table) {
  nsamples = dim(OTU_table)[2]
  n.samples = c(nsamples/2, nsamples/2)
  OTU_table = data.frame(OTU_table)
  res = raida(OTU_table, n.samples, show.ref.features=FALSE, show.all.features=TRUE)
  return(res)
}
RAIDAreslist <- foreach(i = 1: length(simparams), .packages = c("phyloseq", "RAIDA")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix"))
  raida_res = RAIDA_test(OTU_table)
  # rename column "p.adj" to "padj"
  res = cbind(raida_res$p, raida_res$p.adj)
  colnames(res) = c("pvalue", "padj")
  rownames(res) = rownames(raida_res)
  return(res)
}
names(RAIDAreslist) <- simparams
save(RAIDAreslist, file = "RAIDAreslist_dirmulti_Mock_corrected.RData")




########### Omnibus
library(mbzinb)
library(stats)
library(phyloseq)
Omnibus_test <- function(OTU_table) {
  nsamples = dim(OTU_table)[2]
  class = c(rep(0, nsamples/2), rep(1, nsamples/2))
  sample_names = colnames(OTU_table)
  meta = data.frame(class)
  rownames(meta) = sample_names
  # create mbzinb object
  data_class = mbzinb.dataset(data.matrix(OTU_table, rownames.force = NA), meta, taxon = NULL)
  # Omnibus test
  result = mbzinb.test(data_class, group = "class")
  praw = unname(result$results$PValue)
  pval = praw[!is.na(praw)]
  padj = p.adjust(pval, method = "BH", n = length(pval))
  id = names(result$results$PValue)[!is.na(praw)]
  res = data.frame(pval, padj, id)
  colnames(res) = c("pvalue", "padj", "id")
  rownames(res) = id
  return(res)
}
Omnibusreslist <- foreach(i = 1: length(simparams), .packages = c("phyloseq", "mbzinb")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix"))
  res = Omnibus_test(OTU_table)
  return(res)
}
names(Omnibusreslist) <- simparams
save(Omnibusreslist, file = "Omnibusreslist_dirmulti_Mock_corrected.RData")




############ metagenomeSeq
library("metagenomeSeq")
# function for converting phyloseq object to metagenomeSeq object
make_metagenomeSeq <- function(physeq){
  require("metagenomeSeq")
  require("phyloseq")
  if(!taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  OTU = as(otu_table(physeq), "matrix")
  # convert sample data to AnnotatedDataFrame
  ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))
  # define dummy 'feature' data for OTUs
  # use their name helps with extraction and relating to taxonomy later
  TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), row.names = taxa_names(physeq)))
  # create metagenomeSeq object
  MGS = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
  # calculate cumulative sum scaling factor
  MGS = cumNorm(MGS)
  return(MGS)
}
# function for performing relevant test (zero-inflated gaussian)
test_metagenomeSeq <- function(MGS, variable, physeq = NULL){
  # MGS - metagenomeSeq object, produced from function make_metagenomeSeq()
  # variable - variable among the sample_data ('phenoData' in metagenomeSeq)
  # physeq - optional, phyloseq data object that has tax_table
  require("metagenomeSeq")
  require("phyloseq")
  # Create the `mod` variable used in the fitZig test.
  if (inherits(variable, "factor")) {
    # If variable is already a factor, use directly in model.matrix
    mod = model.matrix(~variable)
  } else if (inherits(variable, "matrix")) {
    # If it is a matrix, assume that model.matrix() has been used already
  } else if (inherits(variable, "character")) {
    # If it is a character that specifies a variable in phenoData, use the corresponding variable from MGS
    if (variable %in% colnames(phenoData(MGS)@data)) {
      mod = model.matrix(~phenoData(MGS)@data[, variable])
    } else {
      stop("The phenoData variable name you specified is not present in `phenoData(MGS)`")
    }
  } else {
    stop("Improper specification of the experimental design variable for testing. See `variable` argument")
  }
  # run the Expectation-maximization algorithm 
  # and estimate $f_count$ fits with the zero-inflated Guassian
  fit = fitZig(MGS, mod)
  # specify all OTUs to get the full table from MRfulltable.
  x = MRfulltable(fit, number = nrow(assayData(MGS)$counts))
  # if any OTUs left out, rm those from x. Detected by NA rownames.
  x = x[!is.na(rownames(x)), ]
  # Modify this data.frame by adding the OTUnames. Clip the ':1' added to the OTU names
  rownames(x) <- gsub(":1", "", x = rownames(x), fixed = TRUE)
  x$OTUnames <- as.character(rownames(x))
  if (!is.null(tax_table(physeq, errorIfNULL = FALSE))) {
    # Attach the bacterial taxonomy to the table, if available
    TAX = data.frame(tax_table(physeq))
    TAX$OTUnames <- as.character(rownames(TAX))
    y = merge(x, TAX, by = "OTUnames")
  } else {
    y = x
  }
  # Sort and return
  y = y[order(y$adjPvalue), ]
  return(y)
}
# perform metagenomeSeq differential abundance detection on simulated data
MGSreslist <- foreach(physeq = simlist0, .packages = c("phyloseq", "metagenomeSeq")) %dopar% 
{
  MGSi = make_metagenomeSeq(physeq)
  designfac = factor(gsub("[[:print:]]+\\;", "", sample_names(physeq)))
  y = test_metagenomeSeq(MGSi, designfac)
  rownames(y) <- y[, "OTUnames"]
  # our name convention for B-H adjusted P-value, and OTU-IDs.  
  # Necessary for standardized evaluation of performance
  colnames(y)[colnames(y) == "adjPvalues"] <- "padj"
  colnames(y)[colnames(y) == "OTUnames"] <- "id"
  return(y)
}
names(MGSreslist) <- names(simlist0)
save(MGSreslist, file = "MGSreslist_dirmulti_Mock_corrected.RData")




############ DESeq
library("DESeq")
# function for performing the DESeq test for two classes
# Fisher exact test
deseq_binom_test = function(physeq) {
  designfac = factor(gsub("[[:print:]]+\\;", "", sample_names(physeq)))
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  # Convert to matrix, round up to nearest integer
  x = ceiling(as(otu_table(physeq), "matrix")) + 1
  taxADF = as(data.frame(as(tax_table(physeq), "matrix")), "AnnotatedDataFrame")
  cds = newCountDataSet(x, conditions = designfac, featureData = taxADF)
  cds = estimateSizeFactors(cds)
  sizeFactors(cds)
  cds = estimateDispersions(cds, fitType = c("local"))
  res = nbinomTest(cds, levels(designfac)[1], levels(designfac)[2])
  return(res)
}
# perform DESeq differential abundance detection on simulated data
DESeqreslist <- foreach(physeq = simlist0, .packages = c("phyloseq", "DESeq")) %dopar% {
  return(deseq_binom_test(physeq))
}
names(DESeqreslist) <- names(simlist0)
save(DESeqreslist, file = "DESeqreslist_dirmulti_Mock_corrected.RData")




############ DESeq2
library(DESeq2)
# function for wrapping DESeq2
eval_DESeq2 = function(physeq, IndFilt = NULL) {
  require("DESeq2")
  require("phyloseq")
  # physeq = simlist0[[1]] physeq = rarelist[[7]] Enforce orientation. Samples
  # are columns
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits = 0)
  countData = countData + 1L
  colData = data.frame(sample_data(physeq))
  # Re-order the levels so the NULL set is first.
  colData$postfix <- relevel(factor(colData$postfix), levels(factor(colData$postfix))[2])
  # Create the DESeq data set, dds.
  dds <- DESeqDataSetFromMatrix(countData, colData, design = ~postfix)
  # Run DESeq.
  suppressWarnings(dds <- try(DESeq(dds, quiet = TRUE), silent = TRUE))
  if (inherits(dds, "try-error")) {
    # If the parametric fit failed, try the local.
    suppressWarnings(dds <- try(DESeq(dds, fitType = "local", quiet = TRUE), 
                                silent = TRUE))
    if (inherits(dds, "try-error")) {
      # If local fails, try the mean
      suppressWarnings(dds <- try(DESeq(dds, fitType = "mean", quiet = TRUE), 
                                  silent = TRUE))
    }
    if (inherits(dds, "try-error")) {
      # If still bad, quit with error.
      return(NULL)
    }
  }
  res = results(dds)
  # Independent Filtering
  if (!is.null(IndFilt)) {
    use <- res$baseMean >= IndFilt & !is.na(res$pvalue)
    resFilt <- res[use, ]
    resFilt$padj <- p.adjust(resFilt$pvalue, method = "BH")
    res = as(resFilt, "data.frame")
    res[order(res$padj), ]
  }
  res$id <- rownames(res)
  return(res)
}
# perform DESeq2 differential abundance detection on simulated data
DESeq2reslist <- foreach(physeq = simlist0, .packages = c("phyloseq", "DESeq2")) %dopar% {
  DESeq2_form  = eval_DESeq2(physeq)
  res = data.frame(DESeq2_form)
  return(res)
}
names(DESeq2reslist) <- names(simlist0)
save(DESeq2reslist, file = "DESeq2reslist_dirmulti_Mock_corrected.RData")




############ ANCOM
# ANCOM function 
library(exactRankTests)
library(nlme)
library(ggplot2)
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }
ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}
ANCOMlist <- foreach(i = 1:length(simparams), .packages = c("phyloseq", "exactRankTests", "nlme", "ggplot2")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = as(physeq@otu_table, "matrix") # OTU_Table: Samples x OTUs
  OTU_table = data.frame(Sample.ID = rownames(OTU_table), OTU_table)
  nsamples = dim(OTU_table)[1]
  class = c(rep(0, nsamples/2), rep(1, nsamples/2))
  meta = data.frame(cbind(rownames(OTU_table), class))
  colnames(meta) = c("Sample.ID", "Class")
  comparison_test=ANCOM.main(OTUdat=OTU_table,
                             Vardat = meta,
                             adjusted=FALSE,
                             repeated=F,
                             main.var="Class",
                             adj.formula=NULL,
                             repeat.var=NULL,
                             longitudinal=FALSE,
                             multcorr=2,
                             sig=0.05,
                             prev.cut=0.9)
  res = comparison_test$W.taxa
  return(res)
}
names(ANCOMlist) <- simparams
save(ANCOMlist, file = "ANCOMreslist_dirmulti.RData")


# Stop the cluster
stopCluster(cl)

