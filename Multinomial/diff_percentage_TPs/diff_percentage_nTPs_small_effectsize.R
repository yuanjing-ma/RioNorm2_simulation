## Generate multinomial distribution-based simulation data
## with different percentages of DA-OTUs
## The simulation is based on McMurdie and Holmes(2014).

## load in required package
reqpkg = c("cluster", "doParallel", "foreach", "ggplot2", 
           "grid", "scales", "phyloseq", "plyr", "reshape2", "ROCR")
for (i in reqpkg) {
  print(i)
  print(packageVersion(i))
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}


## set parameters for simulation
set.seed(20180921)
data("GlobalPatterns")
comdelim = "_"
# sample types
sampletypes = levels(get_variable(GlobalPatterns, "SampleType"))
# number of true positives per sample
nTP = c(0.1, 0.15, 0.20, 0.25, 0.3)
# number of OTUs
nOTUs = 1000
# number of samples per class
J = 25
# different values of effect size (should be integers) 
foldeffect = c(2, 3, 4, 5)
#foldeffect = c(1)
# number of reads per sample 
ns = 5000
# minimal num of reads to be considered an OTU as "observed"
minobs = 1
# replication of each parameter setting
reps = 1:10
# combine different simulation parameters
simparams = apply(expand.grid(ns, sampletypes, reps, foldeffect, J, nTP), 1, paste0, collapse = comdelim)
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples", "nTP")


## define template
data("GlobalPatterns")
sampsums = sample_sums(GlobalPatterns)
keepsamples = sample_data(GlobalPatterns)$SampleType %in% sampletypes
# trim OTUs, sort by prevalance and then by abundance value
template = prune_samples(keepsamples, GlobalPatterns)
# list of source templates
templatelist = lapply(sampletypes, function(i, tempall, minobs, nOTUs){
  cat(i, "\n")
  whtemp = (get_variable(tempall, "SampleType") %in% i)
  templatei = prune_samples(whtemp, tempall)
  samobs = apply(otu_table(templatei), 1, function(x,m) sum(x > m), m = minobs)
  otudf = data.frame(prev = samobs, sums = taxa_sums(templatei))
  otudf = otudf[order(-otudf$prev, -otudf$sums), ]
  return(prune_taxa(rownames(otudf)[1:nOTUs], templatei))
}, template, minobs, nOTUs)
names(templatelist) <- sampletypes


## simulate microbiome data from a given sample type template
microbesim = function(postfix = 'sim', template, J, n = 10000){
  # generate J samples with n total reads each all the same
  # OR different if n has length equal to the value of J 
  # postfix is used to distinguish simulated samples in downstream analysis
  require("phyloseq")
  pi = taxa_sums(template)
  if (length(J) != 1){
    stop("Length of J should be 1.")
  }
  if (length(n) != 1 & length(n) != J){
    stop("n should be length 1 or length J.")
  }
  simat = mapply(function(i, x, sample.size){
    phyloseq:::rarefaction_subsample(x, sample.size)
  }, i = 1:J, sample.size = n, MoreArgs = list(x = pi), SIMPLIFY = TRUE)
  simat = t(simat)
  colnames(simat) <- names(pi)
  rownames(simat) <- paste(i, "::", 1:nrow(simat), postfix, sep = "")
  OTU = otu_table(simat, taxa_are_rows = FALSE)
  SDF = data.frame(sample = sample_names(OTU), TableNumber = i, type = "simulated")
  SDF$postfix <- postfix
  rownames(SDF) <- sample_names(OTU)
  SD = sample_data(SDF)
  return(phyloseq(OTU, SD))
}

# function for sampling from the library sizes
sumsim = function(n, sumtemplate, J){
  scaledSums = round(n*(sumtemplate/median(sumtemplate)))
  return(sample(scaledSums, size = J, replace = TRUE))
}

# Register parallel clusters for parallel computing
Ncores = detectCores() 
cl <- makeCluster(Ncores)
registerDoParallel(cl)


# parallelized simulation
simlist <- foreach(i = simparams, .packages = c("phyloseq")) %dopar% {
  # "nreads", "SampleType", "Replicate", "EffectSize", "nsamples"
  params = strsplit(i, comdelim)[[1]]
  names(params) = simparamslabels
  n = sim = sim1 = sim2 = n1 = n2 = NULL
  n = as.numeric(params["nreads"])
  sampletypei = params["SampleType"]
  Ji = as.integer(params["nsamples"])
  templatei = templatelist[[sampletypei]]
  # rarely a simulation has a weird value and fails, catch these with try
  # repeat the simulation call if error
  tryAgain = TRUE
  infiniteloopcounter = 1
  while(tryAgain & infiniteloopcounter < 5){
    n1 = sumsim(n, sampsums, Ji)
    n2 = sumsim(n, sampsums, Ji)
    sim1 = microbesim(paste0(sampletypei, ";grp1"), templatei, Ji, n1)
    sim2 = microbesim(paste0(sampletypei, ";grp2"), templatei, Ji, n2)
    if (is.null(sim1) | is.null(sim2) | is.null(n1) | is.null(n2) | 
        inherits(sim1, "try-error") | inherits(sim2, "try-error")){
      tryAgain = TRUE
      infiniteloopcounter = infiniteloopcounter + 1
    }else{
      tryAgain = FALSE
    }
  }
  if (infiniteloopcounter > 5){
    stop("Consistent error found during simulation. Neet further investigation.")
  }
  sim = merge_phyloseq(sim1, sim2)
  sim = merge_phyloseq(sim, tax_table(GlobalPatterns), phy_tree(GlobalPatterns))
  return(sim)
}
names(simlist) <- simparams

# function for minimal trimming
simpletrim = function(physeq, minobs){
  Ji = nsamples(physeq)
  if (taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  prevalence = apply(as(otu_table(physeq), "matrix"), 2, function(x, minobs){
    return(sum(x>minobs))}, minobs)/Ji
  keepOTUs = prevalence > 0.05 & taxa_sums(physeq) > (0.5 * Ji)
  return(prune_taxa(keepOTUs, physeq))
}
simlist0 <- lapply(simlist, simpletrim, minobs)

# Make artificial differentially abundant OTUs
# Note: We simulate all samples from the same environment template, such that every OTU is a true null
# for different abundance. Then we randomly modify a subset of OTUs with a known effect size
simlistnoeffect = simlist0
simlist0 <- foreach(i = simparams, .packages = c("phyloseq")) %dopar%{
  physeq = simlist0[[i]]
  params = strsplit(i, comdelim)[[1]]
  names(params) <- simparamslabels
  effectsize = as.numeric(params["EffectSize"])
  # number of differentially abundant OTUs 
  perc_nTP = as.numeric(params["nTP"])
  nTP = as.integer(dim(otu_table(physeq))[2] * perc_nTP)
  # randomly select OTUs and add artificial effect
  TPOTUs = sample(taxa_names(physeq), nTP, replace = FALSE)
  # define the samples that will have effect 
  effectsamples = grep(";grp1", sample_names(physeq), fixed = TRUE)
  # apply effect
  otu_table(physeq)[effectsamples, TPOTUs] <- effectsize * otu_table(physeq)[effectsamples, TPOTUs]
  # rename new 'true positive' OTUs with 'TP'
  wh.TP = taxa_names(physeq) %in% TPOTUs
  newname = paste0(taxa_names(physeq)[wh.TP], "-TP")
  colnames(physeq@otu_table)[wh.TP] <- newname
  physeq@phy_tree$tip.label[wh.TP] <- newname
  rownames(physeq@tax_table)[wh.TP] <- newname
  return(physeq)
}
names(simlist0) <- names(simlist)

save(simlist0, file = "simulation_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")
# stop the cluster
stopCluster(cl)
