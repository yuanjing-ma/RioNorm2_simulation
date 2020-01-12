## Compare results of different DA-tests
## Organize results and generate figures in the paper and supplementary file 
## Contain codes for all analysis in the paper
## such as (1). robustness analysis of h value 
##         (2). compare common divisors of RAIDA and RioNorm2
##         (3). compare ANCOM and RioNorm2 etc.


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
setwd("~/Documents/Microbiome/Corrected_simulation_DirMulti")
filename = "simulation-data-dirmulti-Mock-corrected.RData"
load(filename)

## simulation parameters
data("GlobalPatterns")
comdelim = "_"
simparams = names(simlist0)
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")
nOTUs_lib = c()
for (i in simparams) {
  OTU_table = t(simlist0[[i]]@otu_table@.Data)
  nOTUs_lib <- c(nOTUs_lib, dim(OTU_table)[1])
}
names(nOTUs_lib) = simparams


#============================================================================================================
#======== Compare DESeq, DESeq2, MetagenomeSeq, RAIDA, RioNorm2, Omnibus ====================================
#============================================================================================================
library(foreach)
library(doParallel)
Ncores = detectCores() - 1
cl <- makeCluster(Ncores)
registerDoParallel(cl)

library(plyr)
load("RioNorm2reslist_dirmulti_Mock_corrected.RData")
load("RAIDAreslist_dirmulti_Mock_corrected.RData")
load("DESeqreslist_dirmulti_Mock_corrected.RData")
load("DESeq2reslist_dirmulti_Mock_corrected.RData")
load("Omnibusreslist_dirmulti_Mock_corrected.RData")
load("MGSreslist_dirmulti_Mock_corrected.RData")

# DESeq
for (i in 1:length(DESeqreslist)) {
  rownames(DESeqreslist[[i]]) = DESeqreslist[[i]][,"id"]
}
# RioNorm2reslist extract only results
RioNorm2reslist <- list()
list_names <- c()
for (i in 1:length(RioNorm2list)) {
  if (is.null(RioNorm2list[[i]]$result)){
    next
  }
  RioNorm2reslist <- append(RioNorm2reslist, list(RioNorm2list[[i]]$result))
  list_names <- c(list_names, names(RioNorm2list)[i])
}
names(RioNorm2reslist) <- list_names
remove(RioNorm2list)

# funtion for creating a data.frame summarizing the difference abundance performance
make_power_df = function(reslist, comdelim, simparamslabels) {
  require("plyr")
  powerdf = ldply(reslist)
  colnames(powerdf)[1] <- "parameters"
  paramdf = ldply(strsplit(powerdf[, "parameters"], comdelim))
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)}

# put every results list into a coherently-names superlist
superlist = list(RioNorm2 = RioNorm2reslist, RAIDA = RAIDAreslist,
                 DESeq = DESeqreslist, DESeq2 = DESeq2reslist,
                 Omnibus = Omnibusreslist, MetagenomeSeq = MGSreslist)

perfdflist = foreach(resultslist = superlist, .packages = c("ROCR")) %dopar% {
  perflist <- list()
  for (i in names(resultslist)) {
    param = strsplit(i, '_')[[1]]
    names(param) = simparamslabels
    resi = resultslist[[i]]
    resi[is.na(resi[, "padj"]), "padj"] <- 1
    # Evaluate detection performance.
    wh.pred = (as.numeric(resi[, "padj"]) < 0.05)
    wh.pos = which(wh.pred)
    wh.neg = which(!wh.pred)
    wh.TP = grep("[[:print:]]+\\-TP$", rownames(resi))
    # total number of TPs: NTP
    nOTUs = as.numeric(nOTUs_lib[i])
    NTP = as.integer(nOTUs * as.numeric(param['nTP']))
    FPs = sum(!wh.pos %in% wh.TP)
    TPs = sum(wh.pos %in% wh.TP)
    TNs = nOTUs - NTP
    FNs = NTP - TPs
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
    perflist <- append(perflist, list(c(NTP = NTP, FP = FPs, TP = TPs, Power = Power, AUC = AUC, Specificity = Specificity, 
                                   Sensitivity = Sensitivity, FDR = FDR)))
  }
  names(perflist) <- names(resultslist)
  perfdf = make_power_df(perflist, comdelim, simparamslabels)
  return(perfdf)
}
names(perfdflist) <- names(superlist)
# Stop the cluster
stopCluster(cl)

df = ldply(perfdflist)
colnames(df)[1] <- "Approach"
df$Normalization <- "Model/None"
df$Method <- df$Approach
df$nsamples <- as.numeric(df$nsamples)
df$EffectSize <- as.numeric(df$EffectSize)
df$TPpercentage <- as.numeric(df$nTP)
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
dfmeansd$replications = 10
dfmeansd <- dfmeansd[, -which(names(dfmeansd) %in% c("SampleType"))]

# Visualization
library(ggplot2)
dfmeansd_0.10 = dfmeansd[dfmeansd[,'nTP']=="0.10",]
dfmeansd_0.15 = dfmeansd[dfmeansd[,'nTP']=="0.15",]

metric_plot = ggplot(dfmeansd_0.15, aes(EffectSize, Power, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1)) + coord_cartesian(ylim = c(0.4, 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Power plot with 15% of true positives - Mock, Multi-Dirich-corrected"))

metric_plot = ggplot(dfmeansd_0.15, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8)) + coord_cartesian(ylim = c(0.0, 0.8)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("FDR plot with 15% of true positives - Mock, Multi-Dirich-corrected"))




#============================================================================================================
#========= Compare RAIDA's and RioNorm2's common divisor using same test ====================================
#============================================================================================================
# Here, we use Mock data as an example
setwd("~/Documents/Microbiome/Simulation_DirMulti/simulation-results")
load("RioNorm2reslist_dirmulti_Mock.RData")
load("RAIDA2stagereslist_dirmulti_Mock.RData")
load("RAIDAreslist_dirmulti.RData")

# extract RAIDA results of mock environment only 
RAIDAres_Mock <- list()
RAIDAres_Mock_names <- c()
for (i in 1:length(RAIDAreslist)) {
  name = names(RAIDAreslist)[i]
  env_name = strsplit(name, "_")[[1]][2]
  if (env_name == "Mock") {
    RAIDAres_Mock <- append(RAIDAres_Mock, list(RAIDAreslist[[i]]))
    RAIDAres_Mock_names <- c(RAIDAres_Mock_names, name)
  }
}
names(RAIDAres_Mock) <- RAIDAres_Mock_names
RAIDAreslist = RAIDAres_Mock

# extract RioNorm2 results
for (i in 1:length(RioNorm2list)) {
  RioNorm2list[[i]] <- RioNorm2list[[i]]$result
}

# extract RAIDA2stage results
for(i in 1:length(RAIDA2stagelist)) {
  RAIDA2stagelist[[i]] <- RAIDA2stagelist[[i]]$result
}


#### Organize results and calculate mean and sd for each eval metric
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

# funtion for creating a data.frame summarizing the difference abundance performance
make_power_df = function(reslist, comdelim, simparamslabels) {
  require("plyr")
  powerdf = ldply(reslist)
  colnames(powerdf)[1] <- "parameters"
  paramdf = ldply(strsplit(powerdf[, "parameters"], comdelim))
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)}


# Parallel to organize all data into a list of data.frames
library(foreach)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)
comdelim = "_"
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")

# put every results list into a coherently-names superlist
superlist = list(RAIDA = RAIDAreslist,
                 RAIDA2StageTest = RAIDA2stagelist,
                 RioNorm2 = RioNorm2list)

perfdflist = foreach(resultslist = superlist) %dopar% {
  perflist <- lapply(resultslist, eval_res_list)
  names(perflist) <- names(resultslist)
  perfdf = make_power_df(perflist, comdelim, simparamslabels)
  return(perfdf)
}
names(perfdflist) <- names(superlist)
save(perfdflist, file = "eval-matrix-compare-RAIDA-RioNorm2-common-divisor.RData")

# Stop the cluster
stopCluster(cl)


df = ldply(perfdflist)
colnames(df)[1] <- "Approach"
df$Normalization <- "Model/None"
df$Method <- df$Approach
df$nsamples <- as.numeric(df$nsamples)
df$EffectSize <- as.numeric(df$EffectSize)
df$TPpercentage <- as.numeric(df$nTP)
df$FP <- df$FP
df[is.na(df[, "FDR"]), "FDR"] <- 0
write.csv(df, "eval_matrix_compare_RAIDA_RioNorm2_common_divisor.csv")

# calculate means, store as separate dfmean variable
nonrepvars = c("Approach", "Normalization", "Method", simparamslabels[!simparamslabels %in% c("Replicate")])
dfmeansd = ddply(df, nonrepvars[nonrepvars != "SampleType"], function(x, vars) {
  xdf = data.frame(x[1, nonrepvars, drop = FALSE], AUC = mean(x$AUC), sd.AUC = sd(x$AUC))
  # xdf = cbind(xdf, Shannon = mean(x$Shannon), InvSimpson = mean(x$InvSimpson))
  xdf = cbind(xdf, FP = mean(x$FP), sd.FP = sd(x$FP))
  xdf = cbind(xdf, Power = mean(x$Power), sd.Power = sd(x$Power))
  xdf = cbind(xdf, Sensitivity = mean(x$Sensitivity), sd.Sensitivity = sd(x$Sensitivity))
  xdf = cbind(xdf, Specificity = mean(x$Specificity), sd.Specificity = sd(x$Specificity))
  xdf = cbind(xdf, FDR = mean(x$FDR), sd.FDR = sd(x$FDR))
  return(xdf)
}, vars = nonrepvars[nonrepvars != "SampleType"])
dfmeansd$replications = 10 
dfmeansd <- dfmeansd[, -which(names(dfmeansd) %in% c("SampleType"))]

# Visualization
library(ggplot2)

dfmeansd_0.05 = dfmeansd[dfmeansd[,'nTP']=="0.05",]
dfmeansd_0.10 = dfmeansd[dfmeansd[,'nTP']=="0.10",]
dfmeansd_0.15 = dfmeansd[dfmeansd[,'nTP']=="0.15",]

# FDR plot
metric_plot = ggplot(dfmeansd_0.05, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + coord_cartesian(ylim = c(0., 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("FDR plot with 5% of true positives"))

# Power plot
metric_plot = ggplot(dfmeansd_0.05, aes(EffectSize, Power, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + coord_cartesian(ylim = c(0., 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Power plot with 5% of true positives"))




#============================================================================================================
#========= Compare ANCOM and RioNorm2  ======================================================================
#============================================================================================================
# Since ANCOM takes long time to run (~20mins on an OTU table with ~1000 features and ~100 samples)
# we only use subset of simulated data
# env = Mock, effect size = 2,3,4,5, median library size = 5000, nTP = 5%, samples per condition = 25

setwd("~/Documents/Microbiome/Simulation_DirMulti/simulation-results")
load("RioNorm2reslist_dirmulti_Mock.RData")
load("ANCOMreslist_dirmulti_mock_sublist.RData")
comdelim = "_"
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")

# extract comparable RioNorm2 results
RioNorm2list = RioNorm2list[names(ANCOMlist)]

eval_res_list = function(resi) {
  resi = resi$result
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

perflist_RioNorm2 <- lapply(RioNorm2list, eval_res_list)
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

perflist_ANCOM <- lapply(ANCOMlist, eval_res_list_ANCOM)
perflist_ANCOM <- make_power_df(perflist_ANCOM, comdelim, simparamslabels)

perfdflist <- list()
perfdflist <- append(perfdflist, list(perflist_RioNorm2))
perfdflist <- append(perfdflist, list(perflist_ANCOM))
names(perfdflist) <- c("RioNorm2", "ANCOM")
save(perfdflist, file = "eval-matrix-compare-ANCOM-RioNorm2.RData")

df = ldply(perfdflist)
colnames(df)[1] <- "Approach"
df$Normalization <- "Model/None"
df$Method <- df$Approach
df$nsamples <- as.numeric(df$nsamples)
df$EffectSize <- as.numeric(df$EffectSize)
df$TPpercentage <- as.numeric(df$nTP)
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
dfmeansd$replications = 10
dfmeansd <- dfmeansd[, -which(names(dfmeansd) %in% c("SampleType"))]

# Visualization
library(ggplot2)

dfmeansd_0.05 = dfmeansd[dfmeansd[,'nTP']=="0.05",]

metric_plot = ggplot(dfmeansd_0.05, aes(EffectSize, Power, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) + coord_cartesian(ylim = c(0.5, 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Compare ANCOM and RioNorm2: Power plot with 5% of TPs"))

metric_plot = ggplot(dfmeansd_0.05, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.02, 0.04)) + coord_cartesian(ylim = c(0., 0.05)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Compare ANCOM and RioNorm2: FDR plot with 5% of TPs"))




#============================================================================================================
#============== Robustness of h value  ======================================================================
#============================================================================================================
setwd("~/Documents/Microbiome/Simulation_DirMulti/simulation-results")
load("RioNorm2reslist_dirmulti_Mock_h_0.02.RData")
h_0.02 = RioNorm2list_h
load("RioNorm2reslist_dirmulti_Mock.RData")
h_0.025 = RioNorm2list
load("RioNorm2reslist_dirmulti_Mock_h_0.03.RData")
h_0.03 = RioNorm2list_h

#### Organize results and calculate mean and sd for each eval metric
eval_res_list = function(resi) {
  resi = resi$result
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

# funtion for creating a data.frame summarizing the difference abundance performance
make_power_df = function(reslist, comdelim, simparamslabels) {
  require("plyr")
  powerdf = ldply(reslist)
  colnames(powerdf)[1] <- "parameters"
  paramdf = ldply(strsplit(powerdf[, "parameters"], comdelim))
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)}


# Parallel to organize all data into a list of data.frames
library(foreach)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)
comdelim = "_"
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")

# put every results list into a coherently-names superlist
superlist = list(RioNorm2_h_0.02 = h_0.02,
                 RioNorm2_h_0.025 = h_0.025,
                 RioNorm2_h_0.03 = h_0.03)

perfdflist = foreach(resultslist = superlist) %dopar% {
  perflist <- lapply(resultslist, eval_res_list)
  names(perflist) <- names(resultslist)
  perfdf = make_power_df(perflist, comdelim, simparamslabels)
  return(perfdf)
}
names(perfdflist) <- names(superlist)
save(perfdflist, file = "eval-matrix-robustness_analysis_of_h_value.RData")

# Stop the cluster
stopCluster(cl)


df = ldply(perfdflist)
colnames(df)[1] <- "Approach"
df$Normalization <- "Model/None"
df$Method <- df$Approach
df$nsamples <- as.numeric(df$nsamples)
df$EffectSize <- as.numeric(df$EffectSize)
df$TPpercentage <- as.numeric(df$nTP)
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

dfmeansd_0.05 = dfmeansd[dfmeansd[,'nTP']=="0.05",]
dfmeansd_0.10 = dfmeansd[dfmeansd[,'nTP']=="0.10",]
dfmeansd_0.15 = dfmeansd[dfmeansd[,'nTP']=="0.15",]

metric_plot = ggplot(dfmeansd_0.05, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0., 0.2)) + coord_cartesian(ylim = c(0., 0.3)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Robustness of h value: FDR plot with 5% of true positives"))




#============================================================================================================
#============== RioNorm2 normalization + ZIP or ZINB ========================================================
#============================================================================================================
setwd("~/Documents/Microbiome/Simulation_DirMulti/simulation-results")
load("RioNormZIPreslist_dirmulti_Mock.RData")
load("RioNorm2reslist_dirmulti_Mock.RData")
load("RioNormZINBreslist_dirmulti_Mock.RData")

for (i in 1:length(RioNorm2list)) {
  RioNorm2list[[i]] <- RioNorm2list[[i]]$result
}

#### Organize results and calculate mean and sd for each eval metric
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

# funtion for creating a data.frame summarizing the difference abundance performance
make_power_df = function(reslist, comdelim, simparamslabels) {
  require("plyr")
  powerdf = ldply(reslist)
  colnames(powerdf)[1] <- "parameters"
  paramdf = ldply(strsplit(powerdf[, "parameters"], comdelim))
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)}


# Parallel to organize all data into a list of data.frames
library(foreach)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)
comdelim = "_"
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")

# put every results list into a coherently-names superlist
superlist = list(RioNorm2 = RioNorm2list,
                 RioNorm_ZIP = RioNormZIPlist,
                 RioNorm_ZINB = RioNormZINBlist)

perfdflist = foreach(resultslist = superlist) %dopar% {
  perflist <- lapply(resultslist, eval_res_list)
  names(perflist) <- names(resultslist)
  perfdf = make_power_df(perflist, comdelim, simparamslabels)
  return(perfdf)
}
names(perfdflist) <- names(superlist)
save(perfdflist, file = "eval-matrix-RioNorm_with_ZIP_ZINB.RData")

# Stop the cluster
stopCluster(cl)

df = ldply(perfdflist)
colnames(df)[1] <- "Approach"
df$Normalization <- "Model/None"
df$Method <- df$Approach
df$nsamples <- as.numeric(df$nsamples)
df$EffectSize <- as.numeric(df$EffectSize)
df$TPpercentage <- as.numeric(df$nTP)
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

dfmeansd_0.05 = dfmeansd[dfmeansd[,'nTP']=="0.05",]
dfmeansd_0.10 = dfmeansd[dfmeansd[,'nTP']=="0.10",]
dfmeansd_0.15 = dfmeansd[dfmeansd[,'nTP']=="0.15",]

metric_plot = ggplot(dfmeansd_0.15, aes(EffectSize, Power, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) + coord_cartesian(ylim = c(0.5, 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("New normalization with other tests: Power plot with 15% of TPs"))

metric_plot = ggplot(dfmeansd_0.15, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4)) + coord_cartesian(ylim = c(0., 0.5)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("New normalization with other tests: FDR plot with 15% of TPs"))



setwd("~/Documents/Microbiome/Corrected_Simulation_DirMulti/")
load("simulation-data-dirmulti-Mock-corrected.RData")
zero_props = c()
for(i in 1:length(simlist0)){
  OTU_table = t(simlist0_small_effectsize[[i]]@otu_table@.Data)
  zeros = sum(apply(OTU_table, 1, function(x) {sum(x == 0)}))
  prop = zeros/(dim(OTU_table)[1]*dim(OTU_table)[2])
  zero_props <- c(zero_props, prop)
}
histogram(zero_props, title = "histogram of zero proportions of each simulated OTU table")

