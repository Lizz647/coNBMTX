# compute the zero proportion
zero_prop <- function(x){
  length(which(x==0))/length(x)
}

# compute normalize factor for mtx and mgx
Normalize_factor <- function(mtx,mgx,method="TSS",meta){
  if(method == "TSS"){
    sample_sum <- apply(mtx,2,sum)
    factors_mtx = sample_sum / mean(sample_sum)
    sample_sum <- apply(mgx,2,sum)
    factors_mgx = sample_sum / mean(sample_sum)
    factors = list(factors_mtx, factors_mgx)
  } else if (method == "RLE"){
    # library(DESeq2)  # the following code produce the exact same result with DEseq2
    geom_means <- apply(mtx, 1, geom_mean)
    factors_mtx <- apply(mtx, 2, get_median, geom_means = geom_means)
    factors_mtx = factors_mtx/mean(factors_mtx)

    geom_means <- apply(mgx, 1, geom_mean)
    factors_mgx <- apply(mgx, 2, get_median, geom_means = geom_means)
    factors_mgx = factors_mgx/mean(factors_mgx)
    factors = list(factors_mtx, factors_mgx)

  } else if (method == "CSS"){ # not accessible since no taxa-level info
    library(metagenomeSeq)

    #mtx
    feature_anot <- data.frame(feature=rownames(mtx),Taxonomy=rownames(mtx))

    rownames(feature_anot) = feature_anot[,1]

    # metarow = rownames(meta)
    # metarow = paste("S",metarow,sep="")
    # print(metarow)
    # print(colnames(mtx))

    rownames(meta) = colnames(mtx)

    phenotypeData = AnnotatedDataFrame(meta)
    featuredata = AnnotatedDataFrame(feature_anot)

    rownames(mgx) = rownames(mtx)
    obj_mtx = newMRexperiment(mtx,phenoData=phenotypeData,featureData=featuredata)
    p = cumNormStatFast(obj_mtx)
    obj_mtx = cumNorm(obj_mtx, p = p)

    #mgx
    obj_mgx = newMRexperiment(mgx,phenoData=phenotypeData,featureData=featuredata)
    p = cumNormStatFast(obj_mgx)
    obj_mgx = cumNorm(obj_mgx, p = p)

    factors_mtx = normFactors(obj_mtx)/mean(normFactors(obj_mtx))
    factors_mgx = normFactors(obj_mgx)/mean(normFactors(obj_mgx))
    factors = list(factors_mtx, factors_mgx)

  } else if (method == "TMM"){
    library(edgeR) # need more considerations
    sample_sum <- apply(mtx,2,sum)
    factors_mtx = sample_sum / mean(sample_sum)
    # adjust lib size
    normalize_factors = edgeR::calcNormFactors(mtx)
    factors_mtx = factors_mtx*normalize_factors

    sample_sum <- apply(mtx,2,sum)
    factors_mgx = sample_sum / mean(sample_sum)
    # adjust lib size
    normalize_factors = edgeR::calcNormFactors(mgx)
    factors_mgx = factors_mgx*normalize_factors

    factors = list(factors_mtx, factors_mgx)

  }else if(method == "GMPR"){
    library(GMPR)
    library(DESeq2)
    require(vegan)
    require(GUniFrac)

    factors_mtx = GMPR(mtx)
    factors_mgx = GMPR(mgx)

    factors = list(factors_mtx,factors_mgx)

  } else {
    print("No such method!")
  }
}

# used in function: Normalize_factor
geom_mean <- function(x){
  exp(1/length(x)*sum(log(x)))
}
get_median <- function(x,geom_means){
  index <- which(geom_means==0)
  x = x[-index]; geom_means = geom_means[-index]
  median(x/geom_means)
}

# get dispersion trend
GetTrendReg <- function(trenddataframe_mtx,genewise_phi_mtx){
  print(trenddataframe_mtx)
  reg <- glm(dispersion~mu,data=trenddataframe_mtx, family = Gamma(link = "identity"), start = c(0.5,0.5))
  ## fit iteratively
  flag = 1000000; maxit = 10; index=c()
  genewise_phi_mtx_copy <- genewise_phi_mtx
  while(flag>1e-3 & maxit > 0){
    coef_ss <- sum(as.numeric(reg$coefficients)^2)
    if(is.null(index)){
      MyPredict = predict(reg)
      ratio = genewise_phi_mtx/MyPredict
    } else {
      genewise_phi_mtx_copy[index] = 0
      res = predict(reg)
      index_left = setdiff(1:length(genewise_phi_mtx_copy),index)
      MyPredict[index_left] = res
      ratio = genewise_phi_mtx_copy/MyPredict
    }
    index <- which(ratio>10 | ratio<1e-2)
    if(length(index)<=0){
      flag = 0
      break
    }

    trenddataframe_mtx_new <- trenddataframe_mtx[-index,]
    reg <- glm(dispersion~mu,data=trenddataframe_mtx_new, family = Gamma(link = "identity"), start = c(0.5,0.5))

    new_coef_ss <- sum(as.numeric(reg$coefficients)^2)
    flag <- abs(new_coef_ss - coef_ss)
    maxit = maxit - 1
  }
  return(list(reg=reg,index=index))
}

est_trend_sigma <- function(reg,genewise_phi_mtx,trend_outlier,m,p){

  log_disp = log(genewise_phi_mtx)
  MyPredict = predict(reg)
  log_mean = log(MyPredict)

  if(length(trend_outlier)<=0){
    slr = mad(log_disp-log_mean)
  } else {
    slr = mad(log_disp[-trend_outlier]-log_mean)
  }

  # adjust for sampling dist
  adj = trigamma((m-p)/2)

  est_sigma_trend <- slr - adj
  est_sigma_trend
}
