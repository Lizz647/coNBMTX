estimate_FC_zero <- function(y,yg,X,post_phi_mtx, size_factor_mtx, size_factor_mgx,del_index){

  # get useful data
  ini.FCs <- rep(0,dim(X)[2])

  if(length(del_index)>0){
    y <- y[-del_index]
    yg <- yg[-del_index]
    X <- X[-del_index,]
    size_factor_mtx <- size_factor_mtx[-del_index]
    size_factor_mgx <- size_factor_mgx[-del_index]
  }

  opt.vec <- stats::optim(par = ini.FCs,fn=calc.nb.FC.loglik, method = "L-BFGS-B", lower=-40,
                   upper = 40, y = y, yg = yg, X=X,  phi = post_phi_mtx,
                   size_factor_mtx = size_factor_mtx, size_factor_mgx = size_factor_mgx)
  FCs <- opt.vec$par

  return(FCs)
}

calc.nb.FC.loglik <- function(FCs, y, yg, X, phi, size_factor_mtx, size_factor_mgx){
  n = length(y)

  gamma <- exp(X %*% FCs)
  mu = (size_factor_mtx/size_factor_mgx)*yg*gamma

  func1 <- -sum((1/phi+y)*log(1/phi+mu+0.0001))
  func2 <- sum(y*log(mu+0.0001))

  return(-sum(func1, func2))
}

estimate_FC_shrink_zero <- function(y,yg,X,post_phi_mtx,FC_sigma,size_factor_mtx, size_factor_mgx,del_index){

  ini.FCs <- rep(0,dim(X)[2])
  ini.gamma = rep(0,length(y))

  index = setdiff(1:length(y),del_index)

  ini.gamma[index] = y[index]/((size_factor_mtx[index]/size_factor_mgx[index])*yg[index])

  if(length(del_index)>=1){
    ini.gamma = ini.gamma[-del_index]
    if(length(ini.gamma)<=3){
      return(rep(NA,dim(X)[2]))
    }
  }

  if(length(del_index)>0){
    y <- y[-del_index]
    yg <- yg[-del_index]
    X <- X[-del_index,]
    size_factor_mtx <- size_factor_mtx[-del_index]
    size_factor_mgx <- size_factor_mgx[-del_index]
  }

  opt.vec <- stats::optim(par = ini.FCs,fn=calc.nb.FC.loglik.shrink, method = "L-BFGS-B", lower=-40,
                          upper = 40, y = y, yg = yg, X = X, phi = post_phi_mtx, FC_sigma = FC_sigma,
                          size_factor_mtx = size_factor_mtx, size_factor_mgx = size_factor_mgx)
  FCs <- opt.vec$par

  return(FCs)
}

calc.nb.FC.loglik.shrink <- function(FCs, y, yg, X, phi, FC_sigma, size_factor_mtx, size_factor_mgx){
  n = length(y)


  gamma <- exp(X %*% FCs)
  mu = (size_factor_mtx/size_factor_mgx)*yg*gamma

  func1 <- -sum((1/phi+y)*log(1/phi+mu+0.0001))
  func2 <- sum(y*log(mu+0.0001))

  # shirnkage
  adj = 0
  for(j in 2:dim(X)[2]){
    adj = adj - (FCs[j])^2/(2*FC_sigma[j-1]^2)
  }

  return(-sum(func1, func2, adj))
}
