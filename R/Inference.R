estimate_H0_FC_zero <- function(y,yg,X,post_phi_mtx, size_factor_mtx, size_factor_mgx, j,del_index){

  # get useful data
  ini.FCs <- rep(0,dim(X)[2]-1)
  ini.gamma = rep(0,length(y))

  index = setdiff(1:length(y),del_index)

  # delete data with mgx=0 and mtx!=0
  if(length(del_index)>0){
    y <- y[-del_index]
    yg <- yg[-del_index]
    X <- X[-del_index,]
    size_factor_mtx <- size_factor_mtx[-del_index]
    size_factor_mgx <- size_factor_mgx[-del_index]
  }

  # print(ini.FCs)
  opt.vec <- optim(par = ini.FCs,fn=calc.nb.FC.H0.loglik, method = "L-BFGS-B", lower=-40,
                   upper = 40, y = y, yg = yg, X=X,  phi = post_phi_mtx,
                   size_factor_mtx = size_factor_mtx, size_factor_mgx = size_factor_mgx, j=j)
  FCs <- opt.vec$par

  return(FCs)
}

calc.nb.FC.H0.loglik <- function(FCs, y, yg, X, phi, size_factor_mtx, size_factor_mgx,j){
  n = length(y)

  X = X[,-j]
  gamma <- exp(as.matrix(X) %*% as.matrix(FCs))
  mu = (size_factor_mtx/size_factor_mgx)*yg*as.numeric(gamma)

  func1 <- -sum((1/phi+y)*log(1/phi+mu+0.0001))
  func2 <- sum(y*log(mu+0.0001))

  return(-sum(func1, func2))
}

LRT_zero <- function(y,yg,X,phi,size_factor_mtx,size_factor_mgx,FCs,FCs_H0,j,del_index){

  # delete data with mgx=0 and mtx!=0
  if(length(del_index)>0){
    y <- y[-del_index]
    yg <- yg[-del_index]
    X <- X[-del_index,]
    size_factor_mtx <- size_factor_mtx[-del_index]
    size_factor_mgx <- size_factor_mgx[-del_index]
  }

  gamma1 <- exp(as.matrix(X) %*% as.matrix(FCs))
  mu1 = (size_factor_mtx/size_factor_mgx)*yg*as.numeric(gamma1)

  X1 = X[,-j]
  gamma2 <- exp(as.matrix(X1) %*% as.matrix(FCs_H0))
  mu2 = (size_factor_mtx/size_factor_mgx)*yg*as.numeric(gamma2)

  chivar = sum(-2*(1/phi*log(1/phi+mu1)-1/phi*log(1/phi+mu2)+y*log(mu2*(1/phi+mu1)+0.0001)-y*log(mu1*(1/phi+mu2)+0.0001)))
  chivar
}

# shrinkage estimation (DESeq2)

GetFCsigma <- function(FCs,j){
  extreme_FCs = which(abs(FCs[,j])>=5)
  # print(extreme_FCs)
  if(length(extreme_FCs)>0){
    FCs_moderate <- FCs[-extreme_FCs,j]
  } else {
    FCs_moderate <- FCs[,j]
  }
  high = quantile(FCs_moderate,0.95)
  low = quantile(FCs_moderate,0.05)

  med <- (abs(high)+abs(low))/2
  sigma = med/(qnorm(0.95,mean=0,sd=1))
  sigma
}
