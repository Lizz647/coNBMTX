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
    # require(Biobase)
    # if (!require("metagenomeSeq", quietly = TRUE)){
    #   BiocManager::install("metagenomeSeq")
    #   require("metagenomeSeq")
    # }
    # require(metagenomeSeq)

    #mtx
    feature_anot <- data.frame(feature=rownames(mtx),Taxonomy=rownames(mtx))

    rownames(feature_anot) = feature_anot[,1]

    # metarow = rownames(meta)
    # metarow = paste("S",metarow,sep="")
    # print(metarow)
    # print(colnames(mtx))

    rownames(meta) = colnames(mtx)

    phenotypeData = Biobase::AnnotatedDataFrame(meta)
    featuredata = Biobase::AnnotatedDataFrame(feature_anot)

    rownames(mgx) = rownames(mtx)
    obj_mtx = metagenomeSeq::newMRexperiment(mtx,phenoData=phenotypeData,featureData=featuredata)
    p = metagenomeSeq::cumNormStatFast(obj_mtx)
    obj_mtx = metagenomeSeq::cumNorm(obj_mtx, p = p)

    #mgx
    obj_mgx = metagenomeSeq::newMRexperiment(mgx,phenoData=phenotypeData,featureData=featuredata)
    p = metagenomeSeq::cumNormStatFast(obj_mgx)
    obj_mgx = metagenomeSeq::cumNorm(obj_mgx, p = p)

    factors_mtx = metagenomeSeq::normFactors(obj_mtx)/mean(metagenomeSeq::normFactors(obj_mtx))
    factors_mgx = metagenomeSeq::normFactors(obj_mgx)/mean(metagenomeSeq::normFactors(obj_mgx))
    factors = list(factors_mtx, factors_mgx)

  } else if (method == "TMM"){

    # require(edgeR)

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
  stats::median(x/geom_means)
}

# get dispersion trend
GetTrendReg <- function(trenddataframe_mtx,genewise_phi_mtx){
  print(trenddataframe_mtx)
  reg <- stats::glm(dispersion~mu,data=trenddataframe_mtx, family = stats::Gamma(link = "identity"), start = c(0.5,0.5))
  ## fit iteratively
  flag = 1000000; maxit = 10; index=c()
  genewise_phi_mtx_copy <- genewise_phi_mtx
  while(flag>1e-3 & maxit > 0){
    coef_ss <- sum(as.numeric(reg$coefficients)^2)
    if(is.null(index)){
      MyPredict = stats::predict(reg)
      ratio = genewise_phi_mtx/MyPredict
    } else {
      genewise_phi_mtx_copy[index] = 0
      res = stats::predict(reg)
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
    reg <- stats::glm(dispersion~mu,data=trenddataframe_mtx_new, family = stats::Gamma(link = "identity"), start = c(0.5,0.5))

    new_coef_ss <- sum(as.numeric(reg$coefficients)^2)
    flag <- abs(new_coef_ss - coef_ss)
    maxit = maxit - 1
  }
  return(list(reg=reg,index=index))
}

est_trend_sigma <- function(reg,genewise_phi_mtx,trend_outlier,m,p){

  log_disp = log(genewise_phi_mtx)
  MyPredict = stats::predict(reg)
  log_mean = log(MyPredict)

  if(length(trend_outlier)<=0){
    slr = stats::mad(log_disp-log_mean)
  } else {
    slr = stats::mad(log_disp[-trend_outlier]-log_mean)
  }

  # adjust for sampling dist
  adj = trigamma((m-p)/2)

  est_sigma_trend <- slr - adj
  est_sigma_trend
}

GetLackGroupGenes <- function(simu_dna,simu_rna,design_mat){
  lackgroup_gene = c()
  for(i in 1:dim(simu_rna)[1]){
    Y = design_mat
    rnazero <- which(simu_rna[i,]==0)
    dnazero <- which(simu_dna[i,]==0)
    index <- which(simu_dna[i,]!=0)
    del_index <- dnazero
    if(length(del_index)>0){
      Y = Y[-del_index,]
    }
    if(is.null(dim(Y))){
      lackgroup_gene = c(lackgroup_gene,i)
    } else if(length(unique(Y[,2]))==1){
      lackgroup_gene = c(lackgroup_gene,i)
    }
  }
  return(lackgroup_gene)
}
