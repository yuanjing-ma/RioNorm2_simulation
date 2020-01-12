## Robustness analysis of h value using
## Multinomial-distribution based simulation data

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


## load in simulated dataset, will be stored in "env_dat"
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/h_robustness_analysis")
load("../simulation-data-small-effectsize.RData")
# subset of data
subset <- list()
simparams <- c()
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")
comdelim = "_"
for (i in names(simlist0_small_effectsize)){
  params = strsplit(i, comdelim)[[1]]
  names(params) = simparamslabels
  ns = as.numeric(params["nreads"])
  nsamples = as.numeric(params["nsamples"])
  effectsize = as.numeric(params["EffectSize"])
  if (ns == 5000 & nsamples == 25 & effectsize != 1.5) {
    subset <- append(subset, list(simlist0_small_effectsize[[i]]))
    simparams <- c(simparams, i)
  }
}
names(subset) <- simparams
remove(simlist0_small_effectsize)


## register for parallel computing
library(foreach)
library(doParallel)
Ncores = detectCores()
cl <- makeCluster(Ncores)
registerDoParallel(cl)


############ RioNorm2
hk_find_simulation <- function(OTU_table, min_avg_counts = 5, q = 0.03){
  nsamples = dim(OTU_table)[2]
  samobs = apply(OTU_table, 1, function(x) sum(x != 0))
  hk_pool = OTU_table[samobs >= nsamples * 0.9,]
  avg_count = apply(hk_pool, 1, mean)
  hk_pool = hk_pool[avg_count >= min_avg_counts,]
  
  nOTUs = dim(hk_pool)[1]
  OTU_ID = rownames(hk_pool)
  # create symmetric distance matrix between selected OTUs
  ratio_var = matrix(0, nrow = nOTUs, ncol = nOTUs)
  for (i in 1:nOTUs){
    for (j in 1:nOTUs){
      mul = (hk_pool[i,]*hk_pool[j,])
      ind = unname(unlist(mul)) != 0
      ratio_var[i,j] = var(unlist(log(hk_pool[i,][ind]/hk_pool[j,][ind])))
    }
  }
  ratio_var[lower.tri(ratio_var, diag = TRUE)] <- 0
  dist = ratio_var[ratio_var > 0]
  
  nodes = data.frame(seq(1, nOTUs, 1), OTU_ID)
  colnames(nodes) = c("order","OTU_ID")
  # Build links dataset
  links = matrix(0, nrow = 1, ncol=3)
  h = quantile(dist, probs = q)
  order = seq(1,nOTUs,1)
  for (i in 1:dim(ratio_var)[1]){
    check = ratio_var[i,]
    ind = order[check<h & check>0]
    if (length(ind)>0){
      dist = check[ind]
      rep = replicate(length(ind),i)
      record = cbind(rep,ind,dist)
      links = rbind(links,record)
    }
  }
  links = links[-1,]
  library(igraph)
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  largest_cliques = largest_cliques(net) # cliques with max number of nodes
  intersect_cliques = largest_cliques[[1]]
  
  # Normalize OTU count data using sum of housekeeping OTUs
  riOTUs = hk_pool[intersect_cliques,]
  size_factor = colSums(riOTUs)
  riOTUs_ID = OTU_ID[intersect_cliques]
  return(list(riOTUs_ID = riOTUs_ID, size_factor = size_factor))
}

#### Parallel computing
hknorm_2stage_test <- function(OTU_table) {
  # OTU table: OTUs x samples
  nsamples = dim(OTU_table)[2]
  class = c(rep(0, nsamples/2), rep(1, nsamples/2))
  
  # Step 1: find RiOTUs
  size_factor = hk_find_simulation(OTU_table, min_avg_counts = 5, q = 0.02)$size_factor
  log_sizefactor = log(size_factor)
  log_sizefactor[is.infinite(log_sizefactor)] = 0 
  
  # Step 2: apply "overdisp_scoretest" function in the package for over-dispersion test
  # output will be two lists of OTUs: overdispersed and non-overdispersed
  scoretest = overdisp_scoretest(OTU_table, class, log_sizefactor)
  ID_nondisp = scoretest$ID_nondisp
  ID_nondisp = intersect(ID_nondisp, rownames(OTU_table))
  ID_disp = scoretest$ID_disp
  ID_disp = intersect(ID_disp, rownames(OTU_table))
  
  # Step 3: apply "ZIP_test" function in the package to test differential abundance for non-overdispersed OTUs
  nondisp_OTU = OTU_table[ID_nondisp,]
  if (length(ID_nondisp) == 1){
    nondisp_OTU = matrix(nondisp_OTU,nrow = 1)
    rownames(nondisp_OTU) = ID_nondisp
  }
  if (dim(nondisp_OTU)[1] != 0){
    nondisp_res_list = ZIP_test(nondisp_OTU, class, log_sizefactor)
    nondisp_res = cbind(nondisp_res_list$padj, nondisp_res_list$pvalue)
    colnames(nondisp_res) = c("padj", "pvalue")
    rownames(nondisp_res) = nondisp_res_list$id
  }
  
  # Step 4: apply "ZINB_test" function in the package to test differential abundance for overdispersed OTUs
  disp_OTU = OTU_table[ID_disp,]
  if (length(ID_disp) == 1){
    disp_OTU = matrix(disp_OTU, nrow = 1)
    rownames(disp_OTU) = ID_disp
  }
  if (dim(disp_OTU)[1] != 0){
    disp_res_list = ZINB_test(disp_OTU, class, log_sizefactor)
    disp_res = cbind(disp_res_list$padj, disp_res_list$pvalue)
    colnames(disp_res) = c("padj", "pvalue")
    rownames(disp_res) = disp_res_list$id
  }
  
  # combine test results from ZIP and ZINB
  if (dim(nondisp_OTU)[1] == 0){
    result = disp_res
  } else if (dim(disp_OTU)[1] == 0){
    result = nondisp_OTU
  } else {
    result = rbind(nondisp_res, disp_res)
  }
  return(list(result = result, size_factor = size_factor))
}

RioNorm2reslist <- foreach(i = 1:length(subset), .packages = c("igraph", "phyloseq", "MASS", "pscl", "RioNorm2")) %dopar% {
  physeq = subset[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix")) # OTU_Table: OTUs x Samples
  combined_res = NULL
  
  tryCatch({
    combined_res = hknorm_2stage_test(OTU_table)
  }, error = function(e){})
  
  return(combined_res)
}

names(RioNorm2reslist) <- simparams
save(RioNorm2reslist, file = "RioNorm2reslist_small_effectsize_h_0.02_ns_5000_nsamples_25.RData")

# Stop the cluster
stopCluster(cl)


#### Compare results using different h values
library(plyr)
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/h_robustness_analysis/")
load("RioNorm2reslist_small_effectsize_h_0.02_ns_5000_nsamples_25.RData")
h_0.02 = RioNorm2reslist
load("RioNorm2reslist_small_effectsize_h_0.025_ns_5000_nsamples_25.RData")
h_0.025 = RioNorm2reslist
load("RioNorm2reslist_small_effectsize_h_0.03_ns_5000_nsamples_25.RData")
h_0.03 = RioNorm2reslist
load("RioNorm2reslist_small_effectsize_h_0.035_ns_5000_nsamples_25.RData")
h_0.035 = RioNorm2reslist
load("RioNorm2reslist_small_effectsize_h_0.04_ns_5000_nsamples_25.RData")
h_0.04 = RioNorm2reslist
remove(RioNorm2reslist)

superlist = list(h_0.02 = h_0.02,
                 h_0.025 = h_0.025,
                 h_0.03 = h_0.03,
                 h_0.035 = h_0.035,
                 h_0.04 = h_0.04)

# pre-process 
for(s in names(superlist)){
  sublist = superlist[[s]]
  for(l in names(sublist)){
    superlist[[s]][[l]] <- superlist[[s]][[l]]$result
  }
}

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
Ncores = 3
cl <- makeCluster(Ncores)
registerDoParallel(cl)
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")
comdelim = "_"

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
    NTP = length(wh.TP)
    FPs = sum(!wh.pos %in% wh.TP)
    TPs = sum(wh.pos %in% wh.TP)
    TNs = nrow(resi) - NTP
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
colnames(df)[1] <- "h_value"
df$Normalization <- "Model/None"
df$Method <- df$Approach
df$nsamples <- as.numeric(df$nsamples)
df$EffectSize <- as.numeric(df$EffectSize)
df$FP <- df$FP 

# calculate means, store as separate dfmean variable
nonrepvars = c("h_value", "Normalization", "Method", simparamslabels[!simparamslabels %in% c("Replicate")])
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
  scale_y_continuous(breaks = c(0.6, 0.8, 1)) + coord_cartesian(ylim = c(0.6, 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Power plot - Multinomial, all environments"))

metric_plot = ggplot(dfmeansd, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.0, 0.2)) + coord_cartesian(ylim = c(0.0, 0.2)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("FDR plot - Multinomial, all environments"))
