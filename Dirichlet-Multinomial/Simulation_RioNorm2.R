## Apply RioNorm2 on Dirichlet-Multinomial distribution based simulation
## Simulated data available upon request

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
setwd("~/Documents/Microbiome/Corrected_Simulation_DirMulti")
filename = "simulation-data-dirmulti-Mock-corrected.RData"
load(filename)
print("RioNorm2")
print(filename)


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
## Find a group of relatively invariant OTUs (riOTUs)
hk_find_simulation <- function(OTU_table, min_avg_counts = 5){
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
  h = quantile(dist, probs = 0.025)
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

## 2-stage differential abundant test
hknorm_2stage_test <- function(OTU_table) {
  # OTU table: OTUs x samples
  nsamples = dim(OTU_table)[2]
  class = c(rep(0, nsamples/2), rep(1, nsamples/2))
  
  # Step 1: find RiOTUs
  size_factor = hk_find_simulation(OTU_table, min_avg_counts = 10)$size_factor
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

RioNorm2list <- foreach(i = 1:length(simlist0), .packages = c("igraph", "phyloseq", "MASS", "pscl", "RioNorm2")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix")) # OTU_Table: OTUs x Samples
  combined_res = NULL
  
  tryCatch({
    combined_res = hknorm_2stage_test(OTU_table)
  }, error = function(e){})
  
  return(combined_res)
}
## save results
names(RioNorm2list) <-names(simlist0)
save(RioNorm2list, file = "RioNorm2reslist_dirmulti_Mock_corrected.RData")
# Stop the cluster
stopCluster(cl)
