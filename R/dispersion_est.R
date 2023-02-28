# get gene-wise dispersion estimate
## MLE + Cox-Reid
## Notice: We can refine this crude estimate by refine initial mu !!!!!

estimate_phi <- function(y,mu,X){
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mean(y), n)
  }
  min_mu = min(mu)
  ini.phi = (stats::var(y)/mean(y)-1)/mean(y)

  flex.vec = stats::optim(par = ini.phi,fn=calc.nb.loglik, method = "L-BFGS-B", lower=1e-10,
                   upper = 3*(stats::var(y)/min_mu-1)/min_mu, y = y, mu = mu, X=X)

  return(flex.vec$par)
}

calc.nb.loglik <- function(phi, y, mu, X){
  n = length(y)
  func1 <- sum(lgamma(1/phi+y))
  func2 <- -n*lgamma(1/phi)
  func3 <- n*1/phi*log(1/phi)
  func4 <- -sum((1/phi+y)*log(1/phi+mu))

  # Cox-Reid
  W = matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    W[i,i] = 1/(1/mu[i]+phi)
  }

  return(-sum(func1, func2, func3, func4))
}

GetNbPredict <- function(x, meta, groups){
  meta$gene <- as.numeric(x)
  glmnb_gene <- MASS::glm.nb(gene~get(groups[[1]][1]),data = meta)
  ind = which(colnames(meta)==groups[[1]][1])
  # newdata = data.frame(diagnosis=meta[,ind])
  # print(dim(newdata)
  ini_mu <- stats::predict(glmnb_gene, type="response")
  ini_mu
}

GetLMPredict_normalized <- function(x, meta, groups){

  n = length(x)
  ini_mu <- rep(0,length(x))
  del_index <- which(x==-1)
  if(length(del_index)>0){
    x = x[-del_index]
    meta = meta[-del_index,]
    ini_mu[del_index] = -1
  }
  meta$gene <- as.numeric(x)

  # print(groups)
  # print(meta$diagnosis)

  lm_gene <- stats::lm(gene~get(groups[[1]][1]),data = meta)
  # newdata = data.frame(diagnosis=meta$diagnosis)
  ini_mu[setdiff(1:n,del_index)] <- stats::predict(lm_gene, type="response")
  ini_mu
}

###############################
# Estimate posterior dispersion
###############################

estimate_post_phi <- function(y,mu,X,trend,est_sigma_trend, flag){
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mean(y), n)
  }
  min_mu = min(mu)
  ini.phi = (stats::var(y)/mean(y)-1)/mean(y)

  loglik = calc.nb.post.loglik.noDNA

  # if(flag==1){
  #   loglik = calc.nb.post.loglik.noTrend
  # } else if (flag==2) {
  #   loglik = calc.nb.post.loglik.noDNA
  # } else if (flag==3) {
  #   loglik = calc.nb.post.loglik
  # } else {
  #   print("Wrong Flag!")
  # }

  flex.vec = stats::optim(par = ini.phi,fn=loglik, method = "L-BFGS-B", lower=1e-10,
                   upper = 3*(stats::var(y)/min_mu-1)/min_mu, y = y, mu = mu, X=X,
                   trend = trend,est_sigma_trend = est_sigma_trend)

  return(flex.vec$par)
}

calc.nb.post.loglik.noDNA <- function(phi, y, mu, X, trend, est_sigma_trend){
  n = length(y)
  func1 <- sum(lgamma(1/phi+y))
  func2 <- -n*lgamma(1/phi)
  func3 <- n*1/phi*log(1/phi)
  func4 <- -sum((1/phi+y)*log(1/phi+mu))

  # Cox-Reid
  W = matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    W[i,i] = 1/(1/mu[i]+phi)
  }

  # similarly, we do not use it n
  # adj = -1/2*log(det(t(X) %*% W %*% X))

  ## posterior adjustment
  # trend
  tr_adj = - (log(phi)-log(trend))^2/(2*est_sigma_trend^2)

  return(-sum(func1, func2, func3, func4, tr_adj))
}


