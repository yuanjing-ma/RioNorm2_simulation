## Apply other DA-tests on Multinomial-based simulation data
## with different percentages of DA-OTUs
## e.g. DESeq, DESeq2, metagenomeSeq, Omnibus, RAIDA

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


## load in simulated dataset, will be stored in "simlist0"
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/diff_percentage_TPs")
filename = "simulation_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData"
load(filename)

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
save(RAIDAreslist, file = "RAIDAreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")



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
save(Omnibusreslist, file = "Omnibusreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")



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
save(MGSreslist, file = "MGSreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")



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
save(DESeqreslist, file = "DESeqreslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")



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
save(DESeq2reslist, file = "DESeq2reslist_diff_percentage_nTPs_small_effectsize_ns_5000_nsamples_25.RData")

# Stop the cluster
stopCluster(cl)

