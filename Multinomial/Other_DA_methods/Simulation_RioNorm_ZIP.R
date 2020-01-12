## RioNorm2 common divisor + ZIP DA-test
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


############ RioNorm2
RioNorm_ZIP_test <- function(OTU_table, size_factor) {
  # OTU table: OTUs x samples
  nsamples = dim(OTU_table)[2]
  class = c(rep(0, nsamples/2), rep(1, nsamples/2))
  
  log_sizefactor = log(size_factor)
  log_sizefactor[is.infinite(log_sizefactor)] = 0 
  
  # apply "ZIP_test" 
  nondisp_OTU = OTU_table
  nondisp_res_list = ZIP_test(nondisp_OTU, class, log_sizefactor)
  nondisp_res = cbind(nondisp_res_list$padj, nondisp_res_list$pvalue)
  colnames(nondisp_res) = c("padj", "pvalue")
  rownames(nondisp_res) = nondisp_res_list$id
  
  return(list(result = nondisp_res, size_factor = size_factor))
}

RioNormZIPreslist <- foreach(i = names(simlist0), .packages = c("igraph", "phyloseq", "MASS", "pscl", "RioNorm2")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix")) # OTU_Table: OTUs x Samples
  size_factor = common_divisors[[i]]
  combined_res = NULL
  
  tryCatch({
    combined_res = RioNorm_ZIP_test(OTU_table, size_factor)
  }, error = function(e){})
  
  return(combined_res)
}

names(RioNormZIPreslist) <- names(simlist0)
save(RioNormZIPreslist, file = "RioNormZIPreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")

# Stop the cluster
stopCluster(cl)
