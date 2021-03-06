## compare common divisors of RioNorm2 and RAIDA
## using multinomial-distribution based simulation data

# Compute ratios
ratios <- function(count.data, divisor){
  c.data <- as.matrix(count.data);
  rownames(c.data) <- rownames(count.data);
  r.data <- t(apply(c.data, 1, "/", c.data[divisor,]));
  colnames(r.data) <- colnames(count.data);
  return(r.data);
}

# Randomly generate missing (0) ratios using non-zero ratios
get.missing.value <- function(dat, n.missing){
  nonzero.dat <- dat[dat>0];
  est.par <- fitdistr(log(nonzero.dat), "normal")$estimate;
  repeat{
    missing.value <- exp(rnorm(n.missing, est.par[1], est.par[2]));
    if (sum(missing.value==0) == 0){
      break;
    }
  }
  return(missing.value);
}

# Get ratios of features in samples for each condition
get.ratios <- function(data1, data2){
  #Find features with nonzero counts in all the samples for each condition.
  da1 <- data1[!apply(data1, 1, function(y) any(y==0)),];
  da2 <- data2[!apply(data2, 1, function(y) any(y==0)),];
  #Find features with nonzero counts in all the samples for both conditions.
  common.nonzero.rows <- intersect(rownames(da1), rownames(da2));
  n.common.nonzero.rows <- length(common.nonzero.rows);
  #If no features with nonzero counts in all the samples for both conditions,
  #use a feature with nonzero counts in all the samples in each condition as
  #a common divisor to compute ratios. Otherwise, use the feature with nonzero
  #counts in all the samples for both conditions.
  if (n.common.nonzero.rows < 1){
    nonzero.rows1 <- rownames(da1);
    nonzero.rows2 <- rownames(da2);
    n.nonzero.rows1 <- length(nonzero.rows1);
    n.nonzero.rows2 <- length(nonzero.rows2);
    divisor1 <- nonzero.rows1[ceiling(n.nonzero.rows1/2)];
    divisor2 <- nonzero.rows2[ceiling(n.nonzero.rows2/2)];
    rdata1 <- ratios(data1, divisor1);
    rdata2 <- ratios(data2, divisor2);
  } else {
    divisor1 <- common.nonzero.rows[ceiling(n.common.nonzero.rows/2)];
    divisor2 <- divisor1;
    rdata1 <- ratios(data1, divisor1);
    rdata2 <- ratios(data2, divisor2);
  }
  return(list(r1 = rdata1, r2 = rdata2, divisor = c(divisor1, divisor2)));
}

# Cluster using minimax linkage with h=0.05 as a default
clustering <- function(dist, cut=0.05){
  pc <- protoclust(dist);
  group.index <- cutree(pc, h=cut);
  return(group.index);
}

# Combine counts of features in a cluster
combine.mcfs <- function(dat1, dat2, mcfs){
  #Add counts of features in a cluster and store them in the first feature in the cluster.
  mcfs.minus.first.element <- mcfs[-1];
  tmp.dat1 <- dat1[-which(rownames(dat1) %in% mcfs.minus.first.element), ];
  tmp.dat2 <- dat2[-which(rownames(dat2) %in% mcfs.minus.first.element), ];
  tmp.dat1[mcfs[1], ] <- colSums(dat1[mcfs, ]);
  tmp.dat2[mcfs[1], ] <- colSums(dat2[mcfs, ]);
  return(list(x1 = tmp.dat1, x2 = tmp.dat2, divisor = mcfs[1]));
}

# Compare means estimated by EM between two conditions
comparison.bw <- function(ldata, p.adj, ep){
  n.lib1 <- dim(ldata$x1)[2];
  n.lib2 <- dim(ldata$x2)[2];
  dat1 <- as.matrix(ldata$x1);
  dat2 <- as.matrix(ldata$x2);
  n.zero.in.divisor1 <- length(dat1[ldata$divisor, ][dat1[ldata$divisor, ] == 0]);
  #If a divisor contains zeros in its elements, estimate them.
  if(n.zero.in.divisor1 > 0){
    tmp.divisor <- which(apply(dat1, 1, function(y) sum(y==0))==0)[1];
    missing.rvalues <- get.missing.value(dat1[ldata$divisor,]/dat1[tmp.divisor,], n.zero.in.divisor1);
    missing.locations <- which(dat1[ldata$divisor,] == 0);
    est.missing.c1 <- round(missing.rvalues*dat1[tmp.divisor, missing.locations]);
    dat1[ldata$divisor, missing.locations] <- ifelse(est.missing.c1==0, 1, est.missing.c1);
  }
  n.zero.in.divisor2 <- length(dat2[ldata$divisor, ][dat2[ldata$divisor, ] == 0]);
  if(n.zero.in.divisor2 > 0){
    tmp.divisor <- which(apply(dat2, 1, function(y) sum(y==0))==0)[1];
    missing.rvalues <- get.missing.value(dat2[ldata$divisor,]/dat2[tmp.divisor,], n.zero.in.divisor2);
    missing.locations <- which(dat2[ldata$divisor,] == 0);
    est.missing.c2 <- round(missing.rvalues*dat2[tmp.divisor, missing.locations]);
    dat2[ldata$divisor, missing.locations] <- ifelse(est.missing.c2==0, 1, est.missing.c2);
  }
  #Compute ratios using a common divisor across different conditions.
  dat1 = as.data.frame(dat1)
  dat2 = as.data.frame(dat2)
  
  rdat1 <- ratios(dat1, ldata$divisor);
  rdat2 <- ratios(dat2, ldata$divisor);
  #Find epsilon.
  x1.eps <- quantile(rdat1[rdat1>0], ep);
  x2.eps <- quantile(rdat2[rdat2>0], ep);
  #Compute ratios using a common divisor across different conditions after adding epsilon to count data.
  rdat1 <- ratios(dat1+x1.eps, ldata$divisor)[-which(rownames(rdat1) == ldata$divisor), ];
  rdat2 <- ratios(dat2+x2.eps, ldata$divisor)[-which(rownames(rdat2) == ldata$divisor), ];
  
  n.comparisons <- dim(rdat1)[1];
  min.n.lib <- min(n.lib1, n.lib2);
  #Estimate parameters using EM algorithm.
  param1 <- t(apply(rdat1, 1, function(y) em(y, x1.eps)));
  param2 <- t(apply(rdat2, 1, function(y) em(y, x2.eps)));
  #Number of ratios from a lognormal
  n1 <- round((1-param1[,1])*n.lib1);
  n2 <- round((1-param2[,1])*n.lib2);
  dfg <- n1+n2-2; # degree of freedom
  dfg <- ifelse(dfg<1, 1, dfg);
  var.pool <- pooled.var(param1[,3], n1, param2[,3], n2);
  var.pool <- ifelse(dfg<1, 1e-6, var.pool);
  #If the number of libraries is greater than 50, use a two sample t-test.
  #Otherwise, use a moderated t-test.
  if(min.n.lib <= 50){
    eg <- log(var.pool)-digamma(dfg/2)+log(dfg/2);
    mean.eg <- mean(eg);
    f.prime <- (eg-mean.eg)^2*n.comparisons/(n.comparisons-1)-trigamma(dfg/2);
    mean.f.prime <- mean(f.prime);
    if (mean.f.prime <= 0){
      df.prior <- 1e10;
      var.prior <- exp(mean.eg);
    } else{
      df.prior <- trigammaInverse(mean(f.prime))*2;
      var.prior <- exp(mean.eg+digamma(df.prior/2)-log(df.prior/2));
    }
    mod.var <- (df.prior*var.prior + dfg*var.pool)/(df.prior+dfg);
    t.test.rslt <- 2*pt(abs((param1[,2]-param2[,2])/(sqrt(mod.var*(1/n1+1/n2)))), df.prior+dfg, lower.tail=FALSE);
  } else{
    t.test.rslt <- 2*pt(abs((param1[,2]-param2[,2])/(sqrt(var.pool*(1/n1+1/n2)))), dfg, lower.tail=FALSE);
  }
  #Options for multiple test correction
  if(p.adj == "q"){
    t.test.rslt.adj <- qvalue(t.test.rslt)$qvalues;
  } else{
    t.test.rslt.adj <- p.adjust(t.test.rslt, method=p.adj);
  }
  rslt <- data.frame(p = t.test.rslt, p.adj = t.test.rslt.adj, par1 = param1, par2 = param2, mod.var = mod.var);
  ref.counts = list(dat1[ldata$divisor,], dat2[ldata$divisor,]);
  return(list(results=rslt, c.ref=ref.counts, eps=c(x1.eps, x2.eps), priors=c(var.prior, df.prior)));
}

# Remove some samples if no features has non-zero counts for all samples
# and separate the features containing <= 2 non-zero counts from the others
preprocess <- function(dat){
  n.zeros <- apply(dat, 1, function(y) sum(y==0));
  min.n.zeros <- min(n.zeros);
  min.zero.rows <- which(n.zeros == min.n.zeros);
  tmp.dat <- NULL;
  if(min.n.zeros != 0){
    if(length(min.zero.rows) > 1){
      row.sums <- rowSums(dat[c(min.zero.rows),]);
      max.row.sum <- max(row.sums);
      min.zero.row <- min.zero.rows[which(row.sums == max.row.sum)];
      min.n.zeros.loc <- colnames(dat[which(dat[min.zero.row,] == 0)]);
    } else{
      min.zero.row <- min.zero.rows;
      min.n.zeros.loc <- colnames(dat[which(dat[min.zero.row,] == 0)]);
    }
    tmp.dat <- dat[,-which(colnames(dat) %in% min.n.zeros.loc)];
    n.zeros <- apply(tmp.dat, 1, function(y) sum(y==0));
  } else{
    tmp.dat <- dat;
  }
  n.lib <- dim(tmp.dat)[2];
  max.n.zeros <- max(n.zeros);
  max.zero.lib <- n.lib - 3;
  max.n.zeros.loc <- NULL;
  #If a feature has less than 3 nonzero ratios in all the samples in each condition,
  #temporarily remove it from the data.
  if(max.n.zeros > max.zero.lib){
    max.n.zeros.loc <- which(n.zeros > max.zero.lib);
  }
  return(list(x=tmp.dat, rm.rows=names(max.n.zeros.loc)));
}

# Equation of the Bhattacharrya distance for two normal distributions
Db <- function(mu, va){
  mu.va <- cbind(mu, va);
  b.dist.f <- function(y) { 0.25*(log(0.25*(y[3]/y[4] + y[4]/y[3] + 2)) + (y[1]-y[2])^2/(y[3]+y[4])); }
  b.dist <- apply(mu.va, 1, function(y) b.dist.f(y));
  return(b.dist);
}

# Compute the distance between features using the Bhattacharrya distance
get.feature.distance <- function(count.dat, dat, divisor){
  r.eps <- min(dat[dat>0]);
  tmp.dat <- ratios(count.dat+r.eps, divisor)[-which(rownames(count.dat) == divisor), ];
  n.comparisons <- dim(tmp.dat)[1];
  param <- matrix(c(rep(0, n.comparisons*3)), n.comparisons);
  for (i in 1:n.comparisons){
    param[i,] <- em(tmp.dat[i,], r.eps);
  }
  mu.pairs <- t(combn(param[,2], 2));
  var.pairs <- t(combn(param[,3], 2));
  b.dist <- Db(mu.pairs, var.pairs);
  dist.mat <- matrix(numeric(n.comparisons^2), n.comparisons);
  dist.mat[lower.tri(dist.mat)] <- b.dist;
  rownames(dist.mat) <- rownames(tmp.dat);
  colnames(dist.mat) <- rownames(tmp.dat);
  return(as.dist(dist.mat));
}

# Get the common features in clusters across two conditions
get.common.features <- function(cl1, cl2){
  ct1 <- table(cl1);
  ct2 <- table(cl2);
  n.ct1 <- length(ct1);
  n.ct2 <- length(ct2);
  k <- 1;
  cf <- list();
  for(i in 1:n.ct1){
    for(j in 1:n.ct2){
      tmp.cf <- intersect(names(cl1[which(cl1 == i)]), names(cl2[which(cl2 == j)]));
      if(length(tmp.cf) > 0){
        cf[[k]] <- tmp.cf;
        k <- k+1;
      }
    }
  }
  return(cf[sapply(cf, length)>1]);
}  

# Compute pooled variance for the two sample t-test
pooled.var <- function(s1, n1, s2, n2){
  return(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2));
}

# Maximize a log-likelihood function with bound constraints using optim()
opt1 = function(r, z, eta, mu, sigma, eps){
  h = dnorm(log(r), mu, sigma);
  fn = function(x){
    return(sum(ifelse(r<eps, z*log(eta), 0) + ifelse(z==1, 0, (1-z)*log(1-eta)+(1-z)*log(h))));
  }
  
  grad.fn = function(x){
    dl.deta = sum(ifelse(r<eps, z/eta, 0) - ifelse(z==1, 0, (1-z)/(1-eta)));
    dl.dmu = sum(ifelse(z==1, 0, (1-z)*(log(r)-mu)/sigma^2));
    dl.dsigma = sum(ifelse(z==1, 0, (1-z)*((log(r)-mu)^2/sigma^2 - 1/sigma)));
    return(c(dl.deta, dl.dmu, dl.dsigma));
  }
  
  return( optim(c(eta, mu, sigma), fn, grad.fn, method = "L-BFGS-B", lower = c(0, -20, 0), upper = c(1, 20, 100))$par );
}

# Error function
erf <- function(x){ 2*pnorm(x*sqrt(2))-1; }

# Q function
Q <- function(r, eta, mu, sigma, eps){
  n.r <- length(r);
  z <- NULL;
  for(i in 1:n.r){
    if(r[i]<eps){
      z[i] <- eta/(eta + (1-eta)/2*(1+erf((log(r[i])-mu)/(sqrt(2)*sigma))));
    } else{
      z[i] <- 0;
    }
  }
  return(z);
}

# EM Algorithm
em <- function(r, eps){
  mu.H <- s.H <- eta.H <- NULL;
  if(length(r[r>eps]) == 0){
    mu.H <- log(eps);
    s.H <- 1e-6;
    eta.H <- 1;
  } else{
    mu.hat <- mean(log(r[r>eps]));
    s.hat <- sd(log((r[r>eps])));
    s.hat <- ifelse((s.hat==0 | is.na(s.hat)), 1e-6, s.hat);
    eta.hat <- length(r[r<=eps])/length(r);
    n.iter <- 0; max.iter <- 1000; ep <- 1e-6;
    mu.H <- mu.hat; s.H <- s.hat; eta.H <- eta.hat;
    deta <- 10; dmu <- 10; ds <- 10;
    while ((abs(deta) > ep | abs(dmu) > ep | abs(ds) > ep) & n.iter < max.iter){
      z.hat <- Q(r, eta.hat, mu.hat, s.hat, eps);
      param <- opt1(r, z.hat, eta.hat, mu.hat, s.hat, eps);
      eta.hat <- param[1];
      mu.hat <- param[2];
      s.hat <- param[3];
      deta <- eta.H - eta.hat;
      dmu <- mu.H - mu.hat;
      ds <- s.H - s.hat;
      eta.H <- eta.hat;
      mu.H <- mu.hat;
      s.H <- s.hat;
      n.iter <- n.iter + 1;
    }
  }
  return(c(eta.H, mu.H, s.H));
}

# Compare features that contains less than 3 non-zero samples in each condition
comparison.bw.rm <- function(dat1, dat2, ref, eps, prior.param){
  rm.rdat1 <- t(apply(as.matrix(dat1)+eps[1], 1, function(y) y/as.matrix(ref[[1]]+eps[1])));
  colnames(rm.rdat1) <- colnames(dat1);
  rm.rdat2 <- t(apply(as.matrix(dat2)+eps[2], 1, function(y) y/as.matrix(ref[[2]]+eps[2])));
  colnames(rm.rdat2) <- colnames(dat2);
  rm.param1 <- t(apply(rm.rdat1, 1, function(y) em(y, eps[1])));
  rm.param2 <- t(apply(rm.rdat2, 1, function(y) em(y, eps[2])));
  rm.n1 <- round((1-rm.param1[,1])*dim(rm.rdat1)[2]);
  rm.n2 <- round((1-rm.param2[,1])*dim(rm.rdat2)[2]);
  rm.dfg <- rm.n1+rm.n2-2;
  rm.var.pool <- pooled.var(rm.param1[,3], rm.n1, rm.param2[,3], rm.n2);
  rm.var.pool <- ifelse(rm.dfg<1, 1e-6, rm.var.pool);
  rm.dfg <- ifelse(rm.dfg<1, 1, rm.dfg);
  rm.mod.var <- (prior.param[1]*prior.param[2] + rm.dfg*rm.var.pool)/(prior.param[2]+rm.dfg);
  n.rm.dat <- dim(dat1)[1];
  rm.t.test.rslt <- as.numeric(rep(0, n.rm.dat));
  for(i in 1:n.rm.dat){
    #If a feature has more than 1 nonzero counts in all the samples in each condition,
    #use a moderated two sample t-test. Otherwise, use a moderated one sample t-test.
    if(rm.n1[i]==0|rm.n2[i]==0){
      if(rm.n1[i]==0){
        if(rm.n2[i]<2){
          rm.t.test.rslt[i] <- pt((rm.param2[i,2]-log(eps[1]))/sqrt(rm.param2[i,3]), 1e-6, lower.tail=FALSE);
        } else{
          rm.t.test.rslt[i] <- pt((rm.param2[i,2]-log(eps[1]))/(rm.param2[i,3]/sqrt(rm.n2[i])), rm.n2[i]-1, lower.tail=FALSE);
        }
      } else{
        if(rm.n1[i]<2){
          rm.t.test.rslt[i] <- pt((rm.param1[i,2]-log(eps[2]))/rm.param1[i,3], 1e-6, lower.tail=FALSE);
        } else{
          rm.t.test.rslt[i] <- pt((rm.param1[i,2]-log(eps[2]))/(rm.param1[i,3]/sqrt(rm.n1[i])), rm.n1[i]-1, lower.tail=FALSE);
        }
      }
    } else{
      rm.t.test.rslt[i] <- 2*pt(abs((rm.param1[i,2]-rm.param2[i,2])/(sqrt(rm.mod.var[i]*(1/rm.n1[i]+1/rm.n2[i])))), prior.param[2]+rm.dfg[i], lower.tail=FALSE);
    }
  }
  names(rm.t.test.rslt) <- rownames(rm.rdat1);
  rm.rslt <- data.frame(p = rm.t.test.rslt, par1 = rm.param1, par2 = rm.param2, mod.var = rm.mod.var);
  return(rm.rslt);
}

# RAIDA main function
raida_common_divisor <- function(c.data, n.lib, show.ref.features=FALSE, show.all.features=FALSE, mtcm="BH", zsc=0.00){
  if(!is.data.frame(c.data)) stop("Input data must be a data frame!", call.=FALSE);
  if(dim(c.data)[2] != sum(n.lib)) stop("The length of columns for the input data must be equal to the number of all samples!", call.=FALSE);
  #Total number of zeros across conditions
  tn.zeros <- apply(c.data, 1, function(y) sum(y==0));
  zero.features <- names(tn.zeros[which(tn.zeros == sum(n.lib))]);
  #If a feature doesn't have nonzero counts in all the samples for both conditions, remove it.
  if(length(zero.features) > 0){
    c.data <- c.data[-which(rownames(c.data) %in% zero.features),];
  }
  n.features <- dim(c.data)[1];
  # Split data into two sets: features with more than 2 nonzero counts in all samples for each condition
  # and features with less than or equal to 2 nonzero counts in all the samples for each condition.
  n.nonzeros.feature = apply(c.data, 1, function(y) sum(y!=0));
  informative.features = names(n.nonzeros.feature[which(n.nonzeros.feature > 2)]);
  trimmed.c.data = c.data[informative.features,];
  pop1 = trimmed.c.data[,1:n.lib[1]];
  pop2 = trimmed.c.data[,(n.lib[1]+1):(n.lib[1]+n.lib[2])]
  
  pre.pop1 <- NULL; pre.pop2 <- NULL;
  pre.pop1 <- preprocess(pop1);
  pre.pop2 <- preprocess(pop2);
  removed.rows <- union(pre.pop1$rm.rows, pre.pop2$rm.rows);
  removed.pop1 <- NULL; removed.pop2 <- NULL;
  if(is.null(removed.rows) == FALSE){
    removed.pop1 <- pre.pop1$x[removed.rows,];
    removed.pop2 <- pre.pop2$x[removed.rows,];
    pop1 <- NULL; pop2 <- NULL;
    pop1 <- pre.pop1$x[-which(rownames(pre.pop1$x) %in% removed.rows), ];
    pop2 <- pre.pop2$x[-which(rownames(pre.pop2$x) %in% removed.rows), ];
  } else{
    pop1 <- pre.pop1$x;
    pop2 <- pre.pop2$x;
  }
  pre.pop1 <- NULL; pre.pop2 <- NULL;
  #Get ratios
  rdata <- get.ratios(pop1, pop2);
  #Get the Bhattacharyya distance
  dist1 <- get.feature.distance(pop1, rdata$r1, rdata$divisor[1]);
  #Cluster with minimax linkage using the Bhattacharyya distance
  cluster.within1 <- clustering(dist1, cut=0.05);
  dist2 <- get.feature.distance(pop2, rdata$r2, rdata$divisor[2]);
  cluster.within2 <- clustering(dist2, cut=0.05);
  #Get common features across two conditions
  common.features <- get.common.features(cluster.within1, cluster.within2);
  n.common.features <- sapply(common.features, length);
  names(n.common.features)<-seq(1:length(n.common.features));
  n.common.features <- sort(n.common.features, decreasing=TRUE);
  max.n.common.feature <- max(n.common.features);
  #Get possible common divisors
  candidates <- as.numeric(names(which(n.common.features > round(max.n.common.feature/2))));
  if(length(candidates) < 4){
    candidates <- as.numeric(names(n.common.features)[1:min(4, length(n.common.features))]);
  }
  #Compute the number of DAFs and identify the common divisor that gives the smallest number of DAFs.
  n.old.daf <- n.features;
  for(i in 1:length(candidates)){
    comp.bw.rslt <- list();
    aggregated.data <- combine.mcfs(pop1, pop2, common.features[[candidates[i]]]);
    comp.bw.rslt <- comparison.bw(aggregated.data, mtcm, zsc);
    tmp.rslt <- comp.bw.rslt$results;
    n.new.daf <- sum(tmp.rslt$p.adj < 0.05);
    if(n.new.daf < n.old.daf){
      rslt <- tmp.rslt;
      n.old.daf <- n.new.daf;
      ref.fea <- common.features[[candidates[i]]];
      n.ref.fea <- length(ref.fea);
      ref.fea.count <- comp.bw.rslt$c.ref;
      eps <- comp.bw.rslt$eps;
      priors <- comp.bw.rslt$priors;
    }
  }
  if(n.ref.fea < 3){
    warning(paste("The number of features in the final common divisor is ", n.ref.fea, sep=""));
  }
  ref.fea.values = c.data[ref.fea,]
  common_divisors = apply(ref.fea.values, 2, sum)
  return(list(common_divisors = common_divisors, ref.fea = ref.fea))
}


library(RAIDA)
library(RioNorm2)
library(foreach)
library(doParallel)
Ncores = detectCores() - 1 
cl <- makeCluster(Ncores)
registerDoParallel(cl)

setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/RAIDA_RioNorm2_common_divisor/")
load("../simulation-data-small-effectsize.RData")
simparams <- names(simlist0_small_effectsize)
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")

# effect size 2, 3, 4, 5
simlist0 <- list()
names_simlist0 <- c()
for (i in simparams) {
  params = strsplit(i, '_')[[1]]
  names(params) = simparamslabels
  effectsize = as.numeric(params['EffectSize'])
  if (effectsize != 1.5) {
    simlist0 <- append(simlist0, list(simlist0_small_effectsize[[i]]))
    names_simlist0 <- c(names_simlist0, i)
  }
}
names(simlist0) <- names_simlist0
remove(simlist0_small_effectsize)

RAIDA_common_divisor_compute <- function(OTU_table) {
  nsamples = dim(OTU_table)[2]
  n.samples = c(nsamples/2, nsamples/2)
  OTU_table = data.frame(OTU_table)
  res = raida_common_divisor(OTU_table, n.samples)
  return(res)
}

RAIDA_common_divisor_list <- foreach(i = 1: length(simlist0), .packages = c("phyloseq", "RAIDA")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix"))
  raida_res = RAIDA_common_divisor_compute(OTU_table)
  return(raida_res)
}
names(RAIDA_common_divisor_list) <- names(simlist0)
save(RAIDA_common_divisor_list, file = "RAIDA_common_divisors_small_effectsize.RData")


## RAIDA common divisor + 2-stage DA-test 
RAIDA_2stage_test <- function(OTU_table, size_factor) {
  # OTU table: OTUs x samples
  nsamples = dim(OTU_table)[2]
  class = c(rep(0, nsamples/2), rep(1, nsamples/2))
  
  log_sizefactor = log(size_factor)
  log_sizefactor[is.infinite(log_sizefactor)] = 0 
  
  scoretest = overdisp_scoretest(OTU_table, class, log_sizefactor)
  ID_nondisp = scoretest$ID_nondisp
  ID_nondisp = intersect(ID_nondisp, rownames(OTU_table))
  ID_disp = scoretest$ID_disp
  ID_disp = intersect(ID_disp, rownames(OTU_table))
  
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

RAIDA2stagereslist <- foreach(i = names(simlist0), .packages = c("igraph", "phyloseq", "MASS", "pscl", "RAIDA")) %dopar% {
  physeq = simlist0[[i]]
  OTU_table = t(as(physeq@otu_table, "matrix")) # OTU_Table: OTUs x Samples
  size_factor = RAIDA_common_divisor_list[[i]]$common_divisors
  combined_res = NULL
  
  tryCatch({
    combined_res = RAIDA_2stage_test(OTU_table, size_factor)
  }, error = function(e){})
  
  return(combined_res)
}

names(RAIDA2stagereslist) <-names(simlist0)
save(RAIDA2stagereslist, file = "RAIDA2stagereslist_small_effectsize.RData")


#### Compare RAIDA+2stage and RioNorm2
library(plyr)
setwd("~/Documents/Microbiome/Norm_test_microbiome/simulation_data_result/small_effect_size/RAIDA_RioNorm2_common_divisor")
load("RAIDA2stagereslist_small_effectsize.RData")
load("../RioNorm2reslist_small_effectsize.RData")
load("../RAIDAreslist_small_effctsize.RData")

# pre-process RioNorm2reslist 
RioNorm2_sub <- list()
names <- c()
comdelim = '_'
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")
for (i in names(RioNorm2reslist)){
  params = strsplit(i, comdelim)[[1]]
  names(params) = simparamslabels
  effectsize = params['EffectSize']
  if (effectsize != 1.5){
    RioNorm2_sub <- append(RioNorm2_sub, list(RioNorm2reslist[[i]]))
    names <- c(names, i)
  }
}
names(RioNorm2_sub) <- names
RioNorm2reslist = RioNorm2_sub
remove(RioNorm2_sub)

# pre-process RAIDA
RAIDA_sub <- list()
names <- c()
comdelim = '_'
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")
for (i in names(RAIDAreslist)){
  params = strsplit(i, comdelim)[[1]]
  names(params) = simparamslabels
  effectsize = params['EffectSize']
  if (effectsize != 1.5){
    res = RAIDAreslist[[i]]
    res[, 'padj'] = res[, 'p.adj']
    RAIDA_sub <- append(RAIDA_sub, list(res))
    names <- c(names, i)
  }
}
names(RAIDA_sub) <- names
RAIDAreslist = RAIDA_sub
remove(RAIDA_sub)

# pre-process RAIDA + 2stage test
for (i in names(RAIDA2stagereslist)){
  RAIDA2stagereslist[[i]] = RAIDA2stagereslist[[i]]$result
}

make_power_df = function(reslist, comdelim, simparamslabels) {
  require("plyr")
  powerdf = ldply(reslist)
  colnames(powerdf)[1] <- "parameters"
  paramdf = ldply(strsplit(powerdf[, "parameters"], comdelim))
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)}

# put every results list into a coherently-names superlist
# superlist = list(RioNorm2 = RioNorm2reslist, 
#                  RAIDA2stage = RAIDA2stagereslist,
#                  RAIDA = RAIDAreslist)
superlist = list(RioNorm2 = RioNorm2reslist,
                 RAIDA2stage = RAIDA2stagereslist)

# number of total OTUs in each simulated table
nOTUs_lib = c()
for (i in names(RioNorm2reslist)){
  nOTUs_lib = c(nOTUs_lib, dim(RioNorm2reslist[[i]])[1])
}
names(nOTUs_lib) = names(RioNorm2reslist)


library(foreach)
library(doParallel)
Ncores = 3
cl <- makeCluster(Ncores)
registerDoParallel(cl)

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
    NTP = 30
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
  scale_y_continuous(breaks = c(0.6, 0.8, 1)) + coord_cartesian(ylim = c(0.6, 1)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("Power plot - Multinomial, all environments"))

metric_plot = ggplot(dfmeansd, aes(EffectSize, FDR, linetype = Method, color = Method, size=Method)) + geom_path(size = 1) +
  facet_grid(nreads ~ nsamples)  +
  scale_x_continuous(breaks = c(2, 3, 4, 5)) + coord_cartesian(xlim = c(0, 5)) + 
  scale_y_continuous(breaks = c(0.0, 0.2)) + coord_cartesian(ylim = c(0.0, 0.2)) + 
  theme(legend.position = "bottom", legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))
print(metric_plot + ggtitle("FDR plot - Multinomial, all environments"))
