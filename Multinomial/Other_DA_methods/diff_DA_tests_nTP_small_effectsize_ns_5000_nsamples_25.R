## Performance comparison of different downstream DA-tests 
## using same RioNorm2 common divisors
## Under different percentage of DA-OTUs

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


## load in simulated dataset, will be stored in "simlist0"
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/Other_DA_methods")
filename = "simulation_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData"
load(filename)

## simulation parameters
data("GlobalPatterns")
comdelim = "_"
simparams = names(simlist0)
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")


#============================================================================================================
#=================================== Compare 2stage, ZIP, ZINB and ttest ====================================
#============================================================================================================

library(plyr)
load("RioNorm2reslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")
load("RioNormZIPreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")
load("RioNormZINBreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")
load("RioNormttestreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")

# funtion for creating a data.frame summarizing the difference abundance performance
make_power_df = function(reslist, comdelim, simparamslabels) {
  require("plyr")
  powerdf = ldply(reslist)
  colnames(powerdf)[1] <- "parameters"
  paramdf = ldply(strsplit(powerdf[, "parameters"], comdelim))
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)}

library(foreach)
library(doParallel)
Ncores = detectCores() - 1
cl <- makeCluster(Ncores)
registerDoParallel(cl)

# put every results list into a coherently-names superlist
superlist = list(RioNorm2 = RioNorm2reslist,
                 RioNorm_ZIP = RioNormZIPreslist,
                 RioNorm_ZINB = RioNormZINBreslist,
                 RioNorm_ttest = RioNormttestreslist)

perfdflist = foreach(resultslist = superlist, .packages = c("ROCR")) %dopar% {
  perflist <- list()
  for (i in names(resultslist)) {
    param = strsplit(i, '_')[[1]]
    names(param) = simparamslabels
    resi = resultslist[[i]]$result
    resi[is.na(resi[, "padj"]), "padj"] <- 1
    # Evaluate detection performance.
    wh.pred = (as.numeric(resi[, "padj"]) < 0.05)
    wh.pos = which(wh.pred)
    wh.neg = which(!wh.pred)
    wh.TP = grep("[[:print:]]+\\-TP$", rownames(resi))
    # total number of TPs: NTP
    nOTUs = dim(resi)[1]
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
dfmeansd$replications = 90
dfmeansd <- dfmeansd[, -which(names(dfmeansd) %in% c("SampleType"))]

# Visualization
library(ggplot2)

metric_plot = ggplot(dfmeansd, aes(EffectSize, Power, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nsamples ~ nTP)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1)) + coord_cartesian(ylim = c(0.4, 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Power plot with different DA-tests - All environments, Multi"))

metric_plot = ggplot(dfmeansd, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nsamples ~ nTP)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4)) + coord_cartesian(ylim = c(0.0, 0.4)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("FDR plot with different DA-tests - All environments, Multi"))



