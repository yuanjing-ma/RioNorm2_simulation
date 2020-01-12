## Compare ANCOM and RioNorm2 on sebset of Multinomial simulated data

## load in required package
reqpkg = c("cluster", "doParallel", "edgeR", "DESeq", "DESeq2", "stats",
           "foreach", "ggplot2", "grid", "scales", "metagenomeSeq", 
           "phyloseq", "plyr", "reshape2", "ROCR", "RioNorm2", 
           "igraph", "pscl", "MASS")
for (i in reqpkg) {
  print(i)
  print(packageVersion(i))
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}


## load in simulated dataset
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/ANCOM_RioNorm2")
load("../simulation-data-small-effectsize.RData")
# subset of simulated data
subset <- list()
simparams <- c()
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")
for (i in names(simlist0_small_effectsize)){
  params = strsplit(i, '_')[[1]]
  names(params) = simparamslabels
  ns = as.numeric(params["nreads"])
  nsamp = as.numeric(params["nsamples"])
  effectsize = as.numeric(params["EffectSize"])
  if (ns == 5000 & nsamp == 25 & effectsize != 1.5){
    subset <- append(subset, list(simlist0_small_effectsize[[i]]))
    simparams <- c(simparams, i)
  }
}
names(subset) = simparams
remove(simlist0_small_effectsize)

## simulation parameters
data("GlobalPatterns")
comdelim = "_"

# register for parallel computing
library(foreach)
library(doParallel)
Ncores = detectCores()
cl <- makeCluster(Ncores)
registerDoParallel(cl)


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

ANCOMreslist <- foreach(i = 1:length(subset), .packages = c("phyloseq", "exactRankTests", "nlme", "ggplot2")) %dopar% {
  physeq = subset[[i]]
  OTU_table = as(physeq@otu_table, "matrix") # OTU_Table: Samples x OTUs
  OTU_table = data.frame(Sample.ID = rownames(OTU_table), OTU_table)
  nsamples = dim(OTU_table)[1]
  class = c(rep(0, nsamples/2), rep(1, nsamples/2))
  meta = data.frame(cbind(rownames(OTU_table), class))
  colnames(meta) = c("Sample.ID", "Class")
  res = NULL  
  tryCatch({
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
  }, error = function(e){})
  return(res)
}
names(ANCOMreslist) <- simparams
save(ANCOMreslist, file = "ANCOMreslist_small_effectsize_ns_5000_nsamples_25.RData")
stopCluster(cl)



#### Compare ANCOM and RioNorm2
library(plyr)
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/ANCOM_RioNorm2")
load("ANCOMreslist_small_effectsize_ns_5000_nsamples_25.RData")
load("../RioNorm2reslist_small_effectsize.RData")

# extract comparable RioNorm2 results
RioNorm2reslist = RioNorm2reslist[names(ANCOMreslist)]
comdelim = '_'
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")

eval_res_list = function(resi) {
  require("ROCR")
  # Some DESeq2 results (for example) had NA adjusted p-values Replace NA
  # $padj values to highest possible value (1.0)
  resi[is.na(resi[, "padj"]), "padj"] <- 1
  # Evaluate detection performance.
  wh.pred = (as.numeric(resi[, "padj"]) < 0.05)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("[[:print:]]+\\-TP$", rownames(resi))
  NTP = length(wh.TP)
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = nrow(resi) - NTP
  FNs = sum(wh.neg %in% wh.TP)
  Power = TPs/NTP
  # Sensitivity
  Sensitivity = TPs/(TPs + FNs)
  # Specificity
  Specificity = TNs/(TNs + FPs)
  wh.truth = (1:nrow(resi) %in% wh.TP)
  pred <- prediction(as.numeric(wh.pred), factor(wh.truth))
  AUC = performance(pred, "auc")@y.values[[1]]
  # FDR
  FDR = FPs / length(wh.pos)
  return(c(NTP = NTP, FP = FPs, TP = TPs, Power = Power, AUC = AUC, Specificity = Specificity, 
           Sensitivity = Sensitivity, FDR = FDR))
}
make_power_df = function(reslist, comdelim, simparamslabels) {
  require("plyr")
  powerdf = ldply(reslist)
  colnames(powerdf)[1] <- "parameters"
  paramdf = ldply(strsplit(powerdf[, "parameters"], comdelim))
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)}

perflist_RioNorm2 <- lapply(RioNorm2reslist, eval_res_list)
perflist_RioNorm2 <- make_power_df(perflist_RioNorm2, comdelim, simparamslabels)

# extract results from ANCOM
eval_res_list_ANCOM = function(resi) {
  require("ROCR")
  rownames(resi) <- resi$otu.names
  # Evaluate detection performance.
  wh.pred = (as.numeric(resi[, "detected_0.9"]) == TRUE)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("[[:print:]]+\\TP$", rownames(resi))
  NTP = length(wh.TP)
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = nrow(resi) - NTP
  FNs = sum(wh.neg %in% wh.TP)
  Power = TPs/NTP
  # Sensitivity
  Sensitivity = TPs/(TPs + FNs)
  # Specificity
  Specificity = TNs/(TNs + FPs)
  wh.truth = (1:nrow(resi) %in% wh.TP)
  pred <- prediction(as.numeric(wh.pred), factor(wh.truth))
  AUC = performance(pred, "auc")@y.values[[1]]
  # FDR
  FDR = FPs / length(wh.pos)
  return(c(NTP = NTP, FP = FPs, TP = TPs, Power = Power, AUC = AUC, Specificity = Specificity, 
           Sensitivity = Sensitivity, FDR = FDR))
}

perflist_ANCOM <- lapply(ANCOMreslist, eval_res_list_ANCOM)
perflist_ANCOM <- make_power_df(perflist_ANCOM, comdelim, simparamslabels)

perfdflist <- list()
perfdflist <- append(perfdflist, list(perflist_RioNorm2))
perfdflist <- append(perfdflist, list(perflist_ANCOM))
names(perfdflist) <- c("RioNorm2", "ANCOM")

df = ldply(perfdflist)
colnames(df)[1] <- "Approach"
df$Normalization <- "Model/None"
df$Method <- df$Approach
df$nsamples <- as.numeric(df$nsamples)
df$EffectSize <- as.numeric(df$EffectSize)
df$FP <- df$FP 

# calculate means, store as separate dfmean variable
nonrepvars = c("Approach", "Normalization", "Method", simparamslabels[!simparamslabels %in% c("Replicate")])
dfmeansd = ddply(df, nonrepvars[nonrepvars != "SampleType"], function(x, vars) {
  xdf = data.frame(x[1, nonrepvars, drop = FALSE], AUC = mean(x$AUC), sd.AUC = sd(x$AUC))
  # xdf = cbind(xdf, Shannon = mean(x$Shannon), InvSimpson = mean(x$InvSimpson))
  xdf = cbind(xdf, FP = mean(x$FP), sd.FP = sd(x$FP))
  xdf = cbind(xdf, Power = mean(x$Power), sd.Power = sd(x$Power))
  xdf = cbind(xdf, Sensitivity = mean(x$Sensitivity), sd.Sensitivity = sd(x$Sensitivity))
  xdf = cbind(xdf, Specificity = mean(x$Specificity), sd.Specificity = sd(x$Specificity))
  xdf = cbind(xdf, FDR = mean(x$FDR[!is.na(x$FDR)]), sd.FDR = sd(x$FDR[!is.na(x$FDR)]))
  return(xdf)
}, vars = nonrepvars[nonrepvars != "SampleType"])
dfmeansd$replications = 90
dfmeansd <- dfmeansd[, -which(names(dfmeansd) %in% c("SampleType"))]

# Visualization
library(ggplot2)

metric_plot = ggplot(dfmeansd, aes(EffectSize, Power, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) + coord_cartesian(ylim = c(0., 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Power plot - All environments, Multi"))

metric_plot = ggplot(dfmeansd, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0., 0.05, 0.1, 0.15, 0.2)) + coord_cartesian(ylim = c(0., 0.2)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("FDR plot - All environments, Multi"))

