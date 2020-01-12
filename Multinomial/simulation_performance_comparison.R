## Run different DA-tests on multinomial-based simulation data
## e.g. RioNorm2, MetagenomeSeq, DESeq, DESeq2, RAIDA, Omnibus

## Compare results of different DA-tests
## Organize results and generate figures in the paper and supplementary material


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

## simulation parameters
data("GlobalPatterns")
comdelim = "_"
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")

## register for parallel computing
library(foreach)
library(doParallel)
Ncores = detectCores() - 1
cl <- makeCluster(Ncores)
registerDoParallel(cl)


############ RioNorm2
hknorm_2stage_test <- function(OTU_table) {
  # OTU table: OTUs x samples
  nsamples = dim(OTU_table)[2]
  class = c(rep(0, nsamples/2), rep(1, nsamples/2))
  size_factor = hk_find(OTU_table, min_avg_counts = 5)$size_factor
  log_sizefactor = log(size_factor)
  scoretest = overdisp_scoretest(OTU_table, class, log_sizefactor)
  ID_nondisp = scoretest$ID_nondisp
  ID_disp = scoretest$ID_disp
  nondisp_OTU = OTU_table[ID_nondisp,]
  nondisp_res = ZIP_test(nondisp_OTU, class, log_sizefactor)
  disp_OTU = OTU_table[ID_disp,]
  disp_res = ZINB_test(disp_OTU, class, log_sizefactor)
  combined_res = apply(cbind(disp_res, nondisp_res),1,unlist)
  # if needed, rownames of combined_res can reflect overdispersion or not
  return(combined_res)
}
hknormreslist <- foreach(i = 1:length(simparams), .packages = c("igraph", "phyloseq", "MASS", "pscl", "RioNorm2")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix")) # OTU_Table: OTUs x Samples
  combined_res = hknorm_2stage_test(OTU_table)
  rownames(combined_res) = combined_res[,"id"]
  return(combined_res)
}
names(hknormreslist) <- simlist0
save(hknormreslist, file = "RioNorm2reslist_small_effectsize.RData")




############ RAIDA
library(RAIDA)
RAIDA_test <- function(OTU_table) {
  nsamples = dim(OTU_table)[2]
  n.samples = c(nsamples/2, nsamples/2)
  OTU_table = data.frame(OTU_table)
  res = raida(OTU_table, n.samples)
  return(res)
}
RAIDAreslist <- foreach(i = 1: length(simparams), .packages = c("phyloseq", "RAIDA")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix"))
  res = RAIDA_test(OTU_table)
  return(res)
}
names(RAIDAreslist) <- simparams
save(RAIDAreslist, file = "RAIDAreslist_small_effctsize.RData")
# post-processing
# rename column "p.adj" to "padj"
for (i in 1:length(RAIDAreslist)){
  colnames(RAIDAreslist[[i]])[2] <- c("padj")
  RAIDAreslist[[i]][, "id"] <- rownames(RAIDAreslist[[i]])
}




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
  id = rownames(OTU_table)[!is.na(praw)]
  res = data.frame(pval, padj, id)
  return(res)
}
Omnibusreslist <- foreach(i = 1: length(simparams), .packages = c("phyloseq", "mbzinb")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix"))
  res = Omnibus_test(OTU_table)
  return(res)
}
names(Omnibusreslist) <- simparams
save(Omnibusreslist, file = "Omnibusreslist_small_effectsize.RData")
# post-processing
for (i in 1:length(Omnibusreslist)){
  rownames(Omnibusreslist[[i]]) = Omnibusreslist[[i]][, 'id']
}




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
save(MGSreslist, file = "MGSreslist_small_effectsize.RData")




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
save(DESeqreslist, file = "DESeqreslist_small_effectsize.RData")
# post-processing
for (i in 1:length(DESeqreslist)){
  rownames(DESeqreslist[[i]]) = DESeqreslist[[i]][, 'id']
}




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
  return(eval_DESeq2(physeq))
}
names(DESeq2reslist) <- names(simlist0)
save(DESeq2reslist, file = "DESeq2reslist_small_effectsize.RData")

# post-processing
for (i in 1:length(DESeq2reslist)){
  DESeq2reslist[[i]] = data.frame(DESeq2reslist[[i]])
}


#============================================================================================================
# Following part is separated from above
# Load DA-tests results generated from different DA approaches
library(plyr)
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/")
load("DESeqreslist_small_effectsize.RData")
load("DESeq2reslist_small_effectsize.RData")
load("MGSreslist_small_effectsize.RData")
load("RioNorm2reslist_small_effectsize.RData")
load("Omnibusreslist_small_effectsize.RData")
load("RAIDAreslist_small_effctsize.RData")
load("simulation-data-small-effectsize.RData")
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")
simparams <- names(simlist0_small_effectsize)
comdelim = "_"
nTP = 30

nOTUs_lib = c()
for(i in simparams){
  OTU_table = t(simlist0_small_effectsize[[i]]@otu_table@.Data)
  nOTUs_lib <- c(nOTUs_lib, dim(OTU_table)[1])
}
names(nOTUs_lib) = simparams

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
superlist = list(DESeq = DESeqreslist, RioNorm2 = RioNorm2reslist,
                 DESeq2 = DESeq2reslist, metagenomeSeq = MGSreslist,
                 Omnibus = Omnibusreslist, RAIDA = RAIDAreslist)

library(foreach)
library(doParallel)
Ncores = detectCores() - 1
cl <- makeCluster(Ncores)
registerDoParallel(cl)

perfdflist = foreach(resultslist = superlist, .packages = c("ROCR")) %dopar% {
  perflist <- list()
  for (i in 1: length(resultslist)){
    param = strsplit(names(resultslist)[i], '_')[[1]]
    names(param) = simparamslabels
    resi = resultslist[[i]]
    resi[is.na(resi[, "padj"]), "padj"] <- 1
    wh.pred = (as.numeric(resi[, "padj"]) < 0.05)
    wh.pos = which(wh.pred)
    wh.neg = which(!wh.pred)
    wh.TP = grep("[[:print:]]+\\-TP$", rownames(resi))
    nOTUs = as.numeric(nOTUs_lib[names(resultslist)[i]])
    NTP = 30
    FPs = sum(!wh.pos %in% wh.TP)
    TPs = sum(wh.pos %in% wh.TP)
    TNs = nOTUs - NTP
    FNs = NTP - TPs 
    Power = TPs/NTP
    Sensitivity = TPs/(TPs + FNs)
    Specificity = TNs/(TNs + FPs)
    wh.truth = (1:nrow(resi) %in% wh.TP)
    pred <- prediction(as.numeric(wh.pred), factor(wh.truth))
    AUC = performance(pred, "auc")@y.values[[1]]
    FDR = FPs / length(wh.pos)
    perflist <- append(perflist, list(c(NTP = NTP, FP = FPs, TP = TPs, Power = Power, AUC = AUC, Specificity = Specificity, 
                                        Sensitivity = Sensitivity, FDR = FDR)))
  }
  if (is.null(names(resultslist))) {
    names(perflist) <- simparams
  } else {
    names(perflist) <- names(resultslist)
  }
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
df$FP <- df$FP
df[is.na(df[, "FDR"]), "FDR"] <- 0
# extract effect size >= 2 
df = df[df[, 'EffectSize'] >= 2, ]

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
dfmeansd$replications = 90 
dfmeansd <- dfmeansd[, -which(names(dfmeansd) %in% c("SampleType"))]
#write.csv(dfmeansd, "summarized_table_params_without_ZIP_ZINB_small_effectsize.csv")

# Visualization
library(ggplot2)
nreadskeep = c(5000, 10000, 50000)

metric_plot = ggplot(dfmeansd, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6)) + coord_cartesian(ylim = c(0.0, 0.6)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot)

metric_plot = ggplot(dfmeansd, aes(EffectSize, Power, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.6, 0.8, 1.0)) + coord_cartesian(ylim = c(0.5, 1.0)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot)
