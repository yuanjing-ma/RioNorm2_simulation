## RioNorm2 common divisor + t test
## Multinomial-distribution based simulation

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
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/Other_DA_methods/")
filename = "simulation_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData"
load(filename)


## load in common divisors of RioNorm2
load("RioNorm2reslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")
common_divisors <- list()
for (i in 1:length(RioNorm2reslist)){
  cd = RioNorm2reslist[[i]]$size_factor
  common_divisors <- append(common_divisors, list(cd))
}
names(common_divisors) <- names(simlist0)


## simulation parameters
data("GlobalPatterns")
comdelim = "_"
simparams = names(simlist0)
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")


## register for parallel computing
library(foreach)
library(doParallel)
Ncores = detectCores()
cl <- makeCluster(Ncores)
registerDoParallel(cl)


############ RioNorm + ttest
RioNorm_ttest <- function(OTU_table, size_factor) {
  # normalization using size factor 
  for (i in 1:ncol(OTU_table)){
    OTU_table[,i] = OTU_table[,i]/size_factor[i]
  }
  nsamples = dim(OTU_table)[2]

  # two sample t test
  pvalue <- c()
  for (i in 1:nrow(OTU_table)){
    s1 = unname(OTU_table[i,1:(nsamples/2)])
    s2 = unname(OTU_table[i, (nsamples/2 + 1):ncol(OTU_table)])
    result = t.test(s1, s2, var.equal=TRUE, paired=FALSE)
    pvalue <- c(pvalue, result$p.value)
  }
  padj <- p.adjust(pvalue, method = "BH", n = length(pvalue))
  res = cbind(padj, pvalue)
  colnames(res) = c("padj", "pvalue")
  rownames(res) = rownames(OTU_table)
  
  return(list(result = res, size_factor = size_factor))
}

RioNormttestreslist <- foreach(i = names(simlist0), .packages = c("igraph", "phyloseq", "MASS", "pscl", "RioNorm2")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix")) # OTU_Table: OTUs x Samples
  size_factor = common_divisors[[i]]
  combined_res = NULL
  
  tryCatch({
    combined_res = RioNorm_ttest(OTU_table, size_factor)
  }, error = function(e){})
  
  return(combined_res)
}

names(RioNormttestreslist) <- names(simlist0)
save(RioNormttestreslist, file = "RioNormttestreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")

# Stop the cluster
stopCluster(cl)
