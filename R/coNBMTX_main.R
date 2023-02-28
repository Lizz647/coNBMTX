#' Conditional negative binomial regression for sample-paired mgx and mtx data
#'
#' @param datamat A list contains mtx, mgx and metadata.
#' @param group_formula Specify the group information.
#' @param normal_method Specify the way to estimate the scaling factor.
#'
#' @return A list contains the analysis result of each gene.
#' @export
#'
#' @examples
#' mtxtable <- matrix(rnbinom(5000,size=10,mu=50),nrow=100)
#' mgxtable <- matrix(rnbinom(5000,size=5,mu=30),nrow=100)
#' metatable <- data.frame(ID=paste("S",seq(1:50),sep=""),group=c(rep("a",25),rep("b",25)))
#' datamat = list(mtx=mtxtable,mgx=mgxtable,metadata=metatable)
#' group_formula = "~group"
#' output = coNBMTX(datamat=datamat,group_formula=group_formula)
#'

coNBMTX <- function(datamat,group_formula,normal_method="CSS"){

  simu_rna <- datamat$mtx
  simu_dna <- datamat$mgx
  simu_meta <- datamat$metadata

  # check col and row names
  check_rname = rownames(simu_rna)
  if(is.null(check_rname)){
    rownames(simu_dna) = paste("G",seq(1:dim(simu_dna)[1]),sep="")
    rownames(simu_rna) = paste("G",seq(1:dim(simu_rna)[1]),sep="")
  }

  check_cname = colnames(simu_rna)
  if(is.null(check_cname)){
    colnames(simu_dna) = paste("S",seq(1:dim(simu_dna)[2]),sep="")
    colnames(simu_rna) = paste("S",seq(1:dim(simu_rna)[2]),sep="")
    rownames(simu_meta) = paste("S",seq(1:dim(simu_rna)[2]),sep="")
  }

  # DE_genes = simu_datamat$DE_gene

  message("Preprocessing...")

  # Design Metrix: ~ diagnosis + site_name
  group_formula = substr(group_formula,2,nchar(group_formula))
  groups = strsplit(group_formula,'\\+')
  assign(groups[[1]][1],groups[[1]][1])

  groupind = which(colnames(simu_meta)==groups[[1]][1])

  # we add this part later
  # simu_meta$diagnosis <- relevel(as.factor(simu_meta$diagnosis),"nonIBD")
  design_mat <- stats::model.matrix(~get(groups[[1]][1]),data=simu_meta)

  # prepare corresponding data to return
  genewise_phi_mtx.full = rep(-1,dim(simu_rna)[1])
  genewise_phi_mtx_nc.full = rep(-1,dim(simu_rna)[1])
  post_phi_mtx.full = rep(-1,dim(simu_rna)[1])
  pvals.full = rep(1,dim(simu_rna)[1])
  padj.full = rep(1,dim(simu_rna)[1])
  padj.bf.full = rep(1,dim(simu_rna)[1])
  chivars.full = rep(NA,dim(simu_rna)[1])
  FCs.full = matrix(-1,nrow = dim(simu_rna)[1], ncol = dim(design_mat)[2])
  FCs_shrink.full = matrix(-1,nrow = dim(simu_rna)[1], ncol = dim(design_mat)[2])

  pvals_list = list()
  padj_list = list()
  padj.bf_list = list()
  FCs_H0_list = list()
  chivars_list = list()

  # just deal with the situation of two groups
  lackgroup_gene = GetLackGroupGenes(simu_dna = simu_dna,simu_rna = simu_rna, design_mat = design_mat)

  # filter low expression genes
  zero_mgx_gene <- apply(simu_dna,1,zero_prop)
  zero_mtx_gene <- apply(simu_rna,1,zero_prop)
  low_mgx_gene <- which(zero_mgx_gene>=0.95)
  low_mtx_gene <- which(zero_mtx_gene>=0.95)
  low_gene <- union(low_mgx_gene,low_mtx_gene)

  low_gene = union(low_gene,lackgroup_gene)

  left_gene <- setdiff(1:dim(simu_rna)[1],low_gene)

  if(length(low_gene)>0){
    simu_dna <- simu_dna[-low_gene,]
    simu_rna <- simu_rna[-low_gene,]
  }

  mtx_count <- simu_rna
  mgx_count <- simu_dna
  meta <- simu_meta

  ###########################
  # get normalization factors
  ###########################

  message("Caculating normalize factors...")

  Normalize_factors <- Normalize_factor(mtx_count,mgx_count,method=normal_method,
                                        meta=meta)
  size_factor_mtx <- Normalize_factors[[1]]
  size_factor_mgx <- Normalize_factors[[2]]


  ##############################
  # estimate genewise-dispersion
  ##############################

  # # example:
  # y <- as.numeric(mtx_count[5,])
  # mu <- as.numeric(ini_mu_mtx[5,])
  # calc.nb.loglik(phi,y,mu,design_mat)
  # res <- estimate_phi(y,mu,design_mat)

  message("Estimating genewise-dispersion...")

  genewise_phi_mtx <- c()
  genewise_phi_mtx_nc <- c()

  ## non-cond estimation of mtx

  ###############################################
  # "Estimating non-conditional genewise-dispersion..."
  ###############################################

  message("Estimating non-conditional genewise-dispersion...")

  ini_mu_mtx_nc <- t(apply(mtx_count,1,GetNbPredict, meta = meta, groups = groups))

  for(i in 1:dim(mtx_count)[1]){
    y <- as.numeric(mtx_count[i,])
    mu <- as.numeric(ini_mu_mtx_nc[i,])
    mu = mu * size_factor_mtx
    phi <- estimate_phi(y,mu,design_mat)
    genewise_phi_mtx_nc <- c(genewise_phi_mtx_nc,phi)
  }

  genewise_phi_mtx_nc.full[left_gene] = genewise_phi_mtx_nc

  ###############################################
  # "Estimating conditional genewise-dispersion..."
  ###############################################

  message("Estimating conditional genewise-dispersion...")

  # approach: normalize the mtx matrix using the mgx, then estimate rij

  normalized_mtx <- matrix(0, nrow = dim(simu_rna)[1], ncol = dim(simu_rna)[2])

  for(i in 1:dim(simu_rna)[1]){

    if(i%%1000==0){
      mes = paste("normalized mtx for ",i,"th gene",sep="")
      message(mes)
    }

    rnazero <- which(simu_rna[i,]==0)
    dnazero <- which(simu_dna[i,]==0)
    index <- which(simu_dna[i,]!=0)
    del_index <- dnazero
    if(length(del_index)>0){
      normalized_mtx[i,del_index] = -1
    }
    normalized_mtx[i,index] = as.numeric(simu_rna[i,index] * size_factor_mgx[index] /
                                           (size_factor_mtx[index]*simu_dna[i,index]))
    val = as.numeric(stats::quantile(normalized_mtx[i,index],c(0.975)))
    del_index2 <- which(normalized_mtx[i,index]/val>20)
    if(length(del_index2)>0){
      normalized_mtx[i,index][del_index2] = -1
    }
  }

  # ini.rate <- t(apply(normalized_mtx,1,Compute_mean,meta=simu_meta, groups=groups))

  ini.rate <- t(apply(normalized_mtx,1,GetLMPredict_normalized,meta=simu_meta, groups=groups))

  ini_mu_mtx <- matrix(0,nrow = dim(simu_rna)[1], ncol = dim(simu_rna)[2])

  start <- Sys.time()

  for(i in 1:dim(simu_rna)[1]){

    if(i%%1000==0){
      mes = paste("mtx disp for ",i,"th gene",sep="")
      message(mes)
    }

    y <- as.numeric(simu_rna[i,])
    rate <- as.numeric(ini.rate[i,])

    del_index <- which(normalized_mtx[i,]==-1)
    index <- setdiff(1:dim(simu_rna)[2],del_index)

    mu = as.numeric(rate*size_factor_mtx*simu_dna[i,]/(size_factor_mgx))

    ini_mu_mtx[i,] <- mu
    ini_mu_mtx[i,del_index] = -1

    phi <- try(estimate_phi(y[index],mu[index],design_mat[index,]))

    if('try-error' %in% class(phi)){
      genewise_phi_mtx <- c(genewise_phi_mtx,NA)
    } else {
      genewise_phi_mtx <- c(genewise_phi_mtx,phi)
    }
  }

  genewise_phi_mtx.full[left_gene] = genewise_phi_mtx

  #####################
  # check dispersions
  #####################

  # Use a global approach
  # filter extreme low disp
  lowdisp_ind <- which(abs(genewise_phi_mtx.full)<0.001)
  if(length(lowdisp_ind)>0){
    genewise_phi_mtx.full[lowdisp_ind] = NA
  }
  # check high disp/nc_disp
  ratios <- genewise_phi_mtx.full/genewise_phi_mtx_nc.full
  if(length(which(ratios>0.9))){
    badmodel_ind <- setdiff(which(ratios>0.9),low_gene)
    zeros_badmodel <- zero_mtx_gene[badmodel_ind]
  }
  cor_list <- c()
  for(i in badmodel_ind){
    cor_list = c(cor_list,stats::cor(as.numeric(datamat$mtx[i,]),as.numeric(datamat$mgx[i,])))
  }

  ## filter many zeros
  many_zeros <- which(zeros_badmodel>0.60 & cor_list<0.4)
  few_zeros <- which(zeros_badmodel<=0.60 & cor_list<0.050)
  # genewise_phi_mtx.full[badmodel_ind[c(many_zeros,few_zeros)]] = NA

  genewise_phi_mtx = genewise_phi_mtx.full[left_gene]

  newnaind = c(lowdisp_ind,badmodel_ind[c(many_zeros,few_zeros)])

  # get NA index
  naind = which(is.na(genewise_phi_mtx))

  #############################
  # dispersion trend estimation
  #############################

  mtxmean <- apply(mtx_count,1,mean)
  mgxmean <- apply(mgx_count,1,mean)

  if(length(naind)>0){
    trenddataframe_mtx = data.frame(dispersion=genewise_phi_mtx[-naind],mu=1/mtxmean[-naind])

    trendres <- GetTrendReg(trenddataframe_mtx,genewise_phi_mtx[-naind])
    reg <- trendres$reg
    trend_outlier <- trendres$index

    ## estimating sigma in trend
    est_sigma_trend <- est_trend_sigma(reg,genewise_phi_mtx[-naind],trend_outlier,m=dim(simu_rna)[2],p=2)

    # find outliers
    reciprocal_mu <- 1/mtxmean[-naind]
    trend_predict <- stats::predict(reg,newdata = data.frame(mu=reciprocal_mu))
    trend_disp_outlier <- which(abs(log(genewise_phi_mtx[-naind])-log(trend_predict))>2*est_sigma_trend)
    # length(trend_disp_outlier)
  } else {
    trenddataframe_mtx = data.frame(dispersion=genewise_phi_mtx,mu=1/mtxmean)

    trendres <- GetTrendReg(trenddataframe_mtx,genewise_phi_mtx)
    reg <- trendres$reg
    trend_outlier <- trendres$index

    ## estimating sigma in trend
    est_sigma_trend <- est_trend_sigma(reg,genewise_phi_mtx,trend_outlier,m=dim(simu_rna)[2],p=2)

    # find outliers
    reciprocal_mu <- 1/mtxmean
    trend_predict <- stats::predict(reg,newdata = data.frame(mu=reciprocal_mu))
    trend_disp_outlier <- which(abs(log(genewise_phi_mtx)-log(trend_predict))>2*est_sigma_trend)
  }

  ###############################
  # Estimate posterior dispersion
  ###############################

  message("Estimating posterior dispersion...")

  cond_disp_est = T
  post_phi_mtx <- c()

  # recall that ----------------------
  # mtxmean <- apply(mtx_count,1,mean)
  # mgxmean <- apply(mtx_count,1,mean)

  if(length(naind)>0){
    trend <- stats::predict(reg, newdata = data.frame(mu=1/(mtxmean[-naind])))
  } else {
    trend <- stats::predict(reg, newdata = data.frame(mu=1/(mtxmean)))
  }

  for(i in 1:dim(mtx_count)[1]){

    if(i%%1000==0){
      mes = paste("mtx MAP disp for ",i,"th gene",sep="")
      message(mes)
    }

    flag = 0
    y <- as.numeric(mtx_count[i,])

    if(cond_disp_est){
      mu = as.numeric(ini_mu_mtx[i,])
      index <- setdiff(1:dim(mtx_count)[2],which(ini_mu_mtx[i,]<0))
    } else {
      mu <- as.numeric(ini_mu_mtx[i,])*size_factor_mtx
      index <- setdiff(1:dim(mtx_count)[2],which(ini_mu_mtx[i,]<0))
    }

    if(i %in% trend_disp_outlier){

      post_phi_mtx <- c(post_phi_mtx, genewise_phi_mtx[i])

    } else {

      est_sigma_DNA = 3 # nonsense
      flag <- 2
      phi <- try(estimate_post_phi(y[index],mu[index],design_mat[index,],
                                   trend[i],est_sigma_trend,flag))

      if('try-error' %in% class(phi)){
        post_phi_mtx <- c(post_phi_mtx, NA)
      } else {
        post_phi_mtx <- c(post_phi_mtx, phi)
      }

      # post_phi_mtx <- c(post_phi_mtx, phi)

    }
  }

  post_phi_mtx[naind] = NA
  post_phi_mtx.full[left_gene] = post_phi_mtx

  # if genewise!=NA, post=NA -> post = genewise
  diffind = setdiff(which(is.na(post_phi_mtx)),which(is.na(genewise_phi_mtx)))
  post_phi_mtx[diffind] = genewise_phi_mtx[diffind]

  ######################
  # Estimate fold change
  ######################

  message("Estimating fold changes...")
  FCshrink = T

  FCs <- rep(0,dim(design_mat)[2])
  for(i in 1:dim(mtx_count)[1]){

    if(i%%1000==0){
      mes <- paste("est for FCs for ",i,"th gene",sep = "")
      message(mes)
    }

    del_index = which(normalized_mtx[i,]==-1)
    y <- as.numeric(mtx_count[i,])
    yg <- as.numeric(mgx_count[i,])
    phi <- post_phi_mtx[i]

    if(is.na(phi)){
      FC = rep(NA,dim(design_mat)[2])
      FCs <- rbind(FCs, FC)
    } else {
      FC <- try(estimate_FC_zero(y,yg,design_mat,phi,size_factor_mtx, size_factor_mgx,del_index))
      if('try-error' %in% class(FC)){
        FC = rep(NA,dim(design_mat)[2])
        FCs <- rbind(FCs, FC)
      } else {
        FCs <- rbind(FCs, FC)
      }
    }

    if(i==dim(mtx_count)[1]){
      FCs <- FCs[-1,]
    }
  }

  FCs.full[left_gene,] = FCs

  sigma_list  = c()

  if(FCshrink){

    FCs_noNA = FCs[which(!is.na(FCs[,dim(design_mat)[2]])),]

    for(j in 2:dim(design_mat)[2]){
      sigma <- GetFCsigma(FCs_noNA,j)
      sigma_list = c(sigma_list,sigma)
    }

    # shrinkage estimate
    FCs_shrink <- rep(0,dim(design_mat)[2])
    for(i in 1:dim(mtx_count)[1]){

      if(i%%1000==0){
        mes <- paste("shrinkage est for FCs for ",i,"th gene",sep = "")
        message(mes)
      }

      del_index = which(normalized_mtx[i,]==-1)
      y <- as.numeric(mtx_count[i,])
      yg <- as.numeric(mgx_count[i,])
      phi <- post_phi_mtx[i]

      if(is.na(phi)){
        FC = rep(NA,dim(design_mat)[2])
        FCs_shrink <- rbind(FCs_shrink, FC)
      } else {
        FC <- try(estimate_FC_shrink_zero(y,yg,design_mat,phi,FC_sigma=sigma_list,size_factor_mtx, size_factor_mgx,del_index))
        if('try-error' %in% class(FC)){
          FC = rep(NA,dim(design_mat)[2])
          FCs_shrink <- rbind(FCs_shrink, FC)
        } else {
          FCs_shrink <- rbind(FCs_shrink, FC)
        }
      }

      if(i==dim(mtx_count)[1]){
        FCs_shrink <- FCs_shrink[-1,]
      }
    }
    FCs = FCs_shrink
  }


  FCs.full[left_gene,] = FCs

  ###########
  # Inference
  ###########

  message("Inference...")

  X = design_mat
  validind = c()

  # do inference for all coefficients.
  for(j in 2:dim(X)[2]){

    FCs_H0 <- rep(0,dim(design_mat)[2]-1)
    for(i in 1:dim(mtx_count)[1]){

      if(i%%1000==0){
        mes = paste("mtx FCs_H0 for ",i,"th gene",sep="")
        message(mes)
      }

      del_index = which(normalized_mtx[i,]==-1)
      y <- as.numeric(mtx_count[i,])
      yg <- as.numeric(mgx_count[i,])
      phi <- post_phi_mtx[i]

      if(is.na(phi)){
        FC_H0 = rep(NA,dim(design_mat)[2]-1)
        FCs_H0 <- rbind(FCs_H0, FC_H0)
      } else {
        FC_H0 <- try(estimate_H0_FC_zero(y,yg,design_mat,phi,size_factor_mtx, size_factor_mgx,j,del_index))
        if('try-error' %in% class(FC_H0)){
          FC_H0 = rep(NA,dim(design_mat)[2]-1)
          FCs_H0 <- rbind(FCs_H0, FC_H0)
        } else {
          FCs_H0 <- rbind(FCs_H0, FC_H0)
        }
      }

      # phi = phi*1.2

      if(i==dim(mtx_count)[1]){
        FCs_H0 <- FCs_H0[-1,]
      }
    }

    FCs_H0_list[[j]] = FCs_H0

    chivars = c()
    for(i in 1:dim(mtx_count)[1]){

      if(i%%1000==0){
        mes = paste("LRT for ",i,"th gene",sep="")
        message(mes)
      }

      del_index = which(normalized_mtx[i,]==-1)
      y <- as.numeric(mtx_count[i,])
      yg <- as.numeric(mgx_count[i,])
      phi <- post_phi_mtx[i]

      if(is.na(phi)){
        chivar = NA
      } else {
        if(is.na(FCs[i,dim(X)[2]])){
          chivar = NA
        } else {
          if(is.null(dim(FCs_H0))){
            if(is.na(FCs_H0[i])){
              chivar = NA
            } else {
              chivar <- LRT_zero(y,yg,design_mat,phi,size_factor_mtx, size_factor_mgx,FCs[i,],FCs_H0[i],j,del_index)
            }
          } else {
            if(is.na(FCs_H0[i,2])){
              chivar = NA
            } else {
              chivar <- LRT_zero(y,yg,design_mat,phi,size_factor_mtx, size_factor_mgx,FCs[i,],FCs_H0[i,],j,del_index)
            }
          }
        }
      }
      chivars <- c(chivars,chivar)
    }

    # get pvals
    pvals = c()
    for(i in 1:length(chivars)){
      pval <- 1 - stats::pchisq(chivars[i],df=1)
      pvals = c(pvals,pval)
    }
    pvals.full[left_gene] = pvals

    validind = which(!is.na(pvals))
    padj = pvals
    padj.bf = pvals

    chivars.full[left_gene] = chivars

    padj[validind] <- stats::p.adjust(pvals[validind],method = "BH")
    padj.full[left_gene] = padj

    padj.bf[validind] <- stats::p.adjust(pvals[validind],method = "bonferroni")
    padj.bf.full[left_gene] = padj.bf

    pvals_list[[j-1]] = pvals.full
    padj_list[[j-1]] = padj.full
    padj.bf_list[[j-1]] = padj.bf.full
    chivars_list[[j-1]] = chivars.full
  }

  ########################################
  # Construct the component of the output
  ########################################
  # LeftGene
  LeftGene = left_gene

  # GeneName
  GeneName=rownames(simu_rna)

  # NAind for each covariates
  NA_list = list()
  for(i in 1:dim(design_mat)[2]-1){
    if(i!=0){
      NA_list[[i]] = which(is.na(chivars_list[[i]][left_gene]))
    }
  }

  # AvgDNA
  AvgDNA = mgxmean

  # AvgRNA
  AvgRNA = mtxmean

  # GroupName
  GroupName = colnames(design_mat)[-1]

  # LogFCs
  LogFCs=FCs.full[left_gene,]

  # LogFCs
  Chi_list = list()
  for(i in 1:dim(design_mat)[2]-1){
    if(i!=0){
      Chi_list[[i]] = chivars_list[[i]][left_gene]
    }
  }

  # P.vals
  P.Vals = list()
  for(i in 1:dim(design_mat)[2]-1){
    if(i!=0){
      P.Vals[[i]] = pvals_list[[i]][left_gene]
    }
  }

  adj.P.Vals = list()
  for(i in 1:dim(design_mat)[2]-1){
    if(i!=0){
      adj.P.Vals[[i]] = padj.bf_list[[i]][left_gene]
    }
  }

  message("Done!")

  output = list(LeftGene=LeftGene,GeneName=GeneName,NA_list=NA_list,AvgDNA=AvgDNA,AvgRNA=AvgRNA,
                Dispersions=post_phi_mtx.full[left_gene],nonconDispersions=genewise_phi_mtx_nc.full[left_gene],
                GroupName=GroupName,LogFCs=LogFCs,
                chivars_list=Chi_list,P.Vals=P.Vals,adj.P.Vals=adj.P.Vals)

  return(output)
}

#' Get a table for the analysis result of coNBMTX
#'
#' @param output The output of coNBMTX.
#'
#' @return A dataframe contains the analysis result of all genes.
#' @export
#'
#' @examples
#' mtxtable <- matrix(rnbinom(5000,size=10,mu=50),nrow=100)
#' mgxtable <- matrix(rnbinom(5000,size=5,mu=30),nrow=100)
#' metatable <- data.frame(ID=paste("S",seq(1:50),sep=""),group=c(rep("a",25),rep("b",25)))
#' datamat = list(mtx=mtxtable,mgx=mgxtable,metadata=metatable)
#' group_formula = "~group"
#' output = coNBMTX(datamat=datamat,group_formula=group_formula)
#' Tab = GetTab(output)
#'
GetTab <- function(output){
  Tab = data.frame(GeneName=output$GeneName,AvgDNA=output$AvgDNA,AvgRNA=output$AvgRNA,
                   Disps=output$Dispersions,nc_Disps=output$nonconDispersions)
  for(i in 1:length(output$GroupName)+1){
    if(i==1){ID = "Intercept"}
    else{
      ID = paste("LogFC",output$GroupName[i-1],sep="")
    }
    Tab[ID] = output$LogFCs[,i]
  }
  for(i in 1:length(output$GroupName)){
    ID = paste("Chi",output$GroupName[i],sep="")
    Tab[ID] = output$chivars_list[[i]]
  }
  for(i in 1:length(output$GroupName)){
    ID = paste("P.vals",output$GroupName[i],sep="")
    Tab[ID] = output$P.Vals[[i]]
  }
  for(i in 1:length(output$GroupName)){
    ID = paste("adj.P.vals",output$GroupName[i],sep="")
    Tab[ID] = output$adj.P.Vals[[i]]
  }
  return(Tab)
}

#' Get the top table for the analysis result of coNBMTX
#'
#' @param output The output of coNBMTX.
#' @param groupname The name of the group you want to check.
#'
#' @return A dataframe contains the analysis result of top 10 genes.
#' @export
#'
#' @examples
#' mtxtable <- matrix(rnbinom(5000,size=10,mu=50),nrow=100)
#' mgxtable <- matrix(rnbinom(5000,size=5,mu=30),nrow=100)
#' metatable <- data.frame(ID=paste("S",seq(1:50),sep=""),group=c(rep("a",25),rep("b",25)))
#' datamat = list(mtx=mtxtable,mgx=mgxtable,metadata=metatable)
#' group_formula = "~group"
#' output = coNBMTX(datamat=datamat,group_formula=group_formula)
#' topTab = TopTab(output,groupname="a)
#'
TopTab <- function(output,groupname){
  Tab = data.frame(GeneName=output$GeneName,AvgDNA=output$AvgDNA,AvgRNA=output$AvgRNA,
                   Disps=output$Dispersions,nc_Disps=output$nonconDispersions)
  for(i in 1:length(output$GroupName)+1){
    if(i==1){ID = "Intercept"}
    else{
      ID = paste("LogFC",output$GroupName[i-1],sep="")
    }
    Tab[ID] = output$LogFCs[,i]
  }
  for(i in 1:length(output$GroupName)){
    ID = paste("Chi",output$GroupName[i],sep="")
    Tab[ID] = output$chivars_list[[i]]
  }
  for(i in 1:length(output$GroupName)){
    ID = paste("P.vals",output$GroupName[i],sep="")
    Tab[ID] = output$P.Vals[[i]]
  }
  for(i in 1:length(output$GroupName)){
    ID = paste("adj.P.vals",output$GroupName[i],sep="")
    Tab[ID] = output$adj.P.Vals[[i]]
  }

  colind = which(stringr::str_detect(colnames(Tab),groupname) & stringr::str_detect(colnames(Tab),"P.vals") & !stringr::str_detect(colnames(Tab),"adj.P.vals"))
  orderedTab = Tab[order(Tab[,colind]),]
  return(orderedTab[1:10,])
}
