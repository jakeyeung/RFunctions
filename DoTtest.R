# 2016-02-02
# Perform T-test on samples given 
# sample means, variances, and sample sizes

TwoSampVar <- function(var1, var2, n1, n2){
  # get the denominator for Ttest unequal variance 
  return((var1 / n1) + (var2 / n2))
}

DegFree <- function(var1, var2, n1, n2){
  degfree <- ((var1 / n1 + var2 / n2) ^ 2) / ((var1 / n1) ^ 2 / (n1 - 1) + (var2 / n2) ^ 2 / (n2 - 1))
  return(degfree)
}

PvalFromTDist <- function(t.stat, degfree, two.sided = TRUE){
  if (t.stat > 0) t.stat <- t.stat * -1  # take the negative part so our pvals are always less than 1
  if (two.sided){
    pval <- pt(t.stat, degfree) * 2
  } else {
    pval <- pt(t.stat, degfree)
  }
  return(pval)
}

DoTtest <- function(mean1, mean2, var1, var2, n1, n2){
  # https://en.wikipedia.org/wiki/Student%27s_t-test#Equal_or_unequal_sample_sizes.2C_unequal_variances
  t.stat <- (mean1 - mean2) / sqrt(TwoSampVar(var1, var2, n1, n2))
  degfree <- DegFree(var1, var2, n1, n2)
  pval <- PvalFromTDist(t.stat, degfree, two.sided=TRUE)
  return(pval)
}

TtestDat <- function(dat, mean.colname = "exprs", var.colname = "exprs.var.log2", n = 3, except.time = 0, nsd.i = 1){
  # do ttest on dataframe
  # 3 samples per time point except at time 0 where the NSD has 4 samples
  # nsd.i: index of the NSD sample, either 1 or 2. SD will be 2 or 1 respectively 
  
  # check dat is right format
  if (nrow(dat) != 2) warning("Warning, expect dat to have 2 rows")
  if (dat$trtmnt[[nsd.i]] != "NSD"){
    warning(paste("Expected NSD to be at position:", nsd.i))
  }
  
  time <- dat$time[[1]]
  gene <- dat$gene[[1]]
  
  n2 <- n  # 
  if (time %% 24 != except.time){
    n1 <- n + 1  # this time has exceptionally 1 + n (4 normally)
  } else {
    n1 <- n
  }
  
  mean1 <- dat[[mean.colname]][1]; mean2 <- dat[[mean.colname]][2]
  var1 <- dat[[var.colname]][1]; var2 <- dat[[var.colname]][2]
  
  pval <- DoTtest(mean1, mean2, var1, var2, n1, n2)
  
  return(data.frame(nsd.minus.sd = mean1 - mean2, pval = pval))
}