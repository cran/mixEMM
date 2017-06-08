#' A mixed-effects model for analyzing cluster-level non-ignorable missing data
#' 
#' This function fits a mixed-effects model for clustered data with cluster-level missing values in the outcome.
#' 
#' The model consists of two parts, the outcome model and the missing-data model. The outcome model 
#' is a mixed-effects model,
#' \deqn{\mathbf{y}_{i} =  \mathbf{X}_{i}\boldsymbol{\alpha}+\mathbf{Z}_{i}\boldsymbol{b}_{i}+\mathbf{e}_{i},}
#' where \eqn{\mathbf{y}_{i}} is the outcome for the i-th cluster, \eqn{\mathbf{X}_{i}} is the covariate matrix,
#' \eqn{\boldsymbol{\alpha}} is the fixed-effects, \eqn{\mathbf{Z}_{i}} is the design matrix for
#' the random-effects \eqn{\mathbf{b}_i}, and \eqn{\mathbf{e}_{i}} is the error term. 
#' 
#' The non-ignorable batch-level (or cluster-level) abundance-dependent missing-data model (BADMM) can be written as
#' \deqn{\textrm{Pr}\left(M_{i}=1|\mathbf{y}_{i}\right)= \mathrm{exp}\left(-\gamma_{0} - \gamma \bar\mathbf{y}_{i}
#' \right),}
#' where \eqn{M_{i}} is the missing indicator for the i-th cluster, and \eqn{\bar\mathbf{y}_{i}} is the average of \eqn{\mathbf{y}_{i}}. 
#' If \eqn{M_{i}=1}, the outcome of the i-th cluster  
#' \eqn{\mathbf{y}_{i}} would be missing altogether. 
#' The estimation of the mixEMM model is implemented via an ECM algorithm. If \eqn{\gamma \neq 0}, i.e., 
#' the missingness depends on the outcome, the missing-data mechanism is missing not at random (MNAR), 
#' otherwise it is missing completely at random (MCAR) for the current model. The parameter \eqn{\gamma} can 
#' be estimated by borrowing information across outcomes and finding the common missing-data patterns
#' in the high-dimensional data. For example, by estimating the relationship
#' the observed average value of \eqn{\bar\mathbf{y}_{i}} and the missing rate, or the parameter can be
#' selected by the log-likelihood profile (see the Reference).
#' 
#' @param Ym is an N by p outcome data from N clusters/batches/experiments; p is the number of samples within each cluster.
#' The first sample within each cluster is assumed to be a reference sample with different error variance. 
#' Missing values are coded as NAs.
#' @param Xm is a covariate array of dimension N by k by p, where k is the number of covariates.
#' @param Zm is a design array for random-effects, with a dimension of N by h by p, where h is the number of variables with random effects.
#' @param gamma is the parameter for the missing-data mechanism. The missingness of the outcome in cluster i 
#' depends on the mean of the outcome. The missing probability is modelled as exp(-gamma0 - gamma*mean(y)). The parameter gamma can 
#' be estimated by borrowing information across outcomes and finding the common missing-data patterns
#' in the high-dimensional data. For example, by estimating the relationship
#' the observed average value of \eqn{\bar\mathbf{y}_{i}} and the missing rate, or the parameter can be
#' selected by the log-likelihood profile (see the Reference).
#' If gamma = 0, the missingness is ignorable. The parameter gamma0 does not affect the estimation of the EM algorithm,
#' and is mostly determined by the missing rate. So it is set as 0 in the estimation here.
#' @param maxIter the maximum number of iterations in the estimation of the EM algorithm.
#' @param tol the tolerance level for the absolute change in the observed-data log-likelihood function.
#' 
#' @return A list containing
#' \item{alpha.hat}{the estimated fixed-effects.}
#' \item{alpha.se}{the standard errors for the estimated fixed-effects.}
#' \item{sigma0.hat, sigma2.hat}{the estimated sample error variances. It returns 
#'  the variances for the first (reference) sample and the other samples within each cluster/batch.}
#' \item{D}{the estimated covariance matrix for the random-effects.}
#' \item{RE}{the estimated random-effects.}
#' \item{loglikelihood}{the observed-data log-likelihood values.}
#' 
#' @references Chen, L. S., Wang, J., Wang, X., & Wang, P. (2017). A mixed-effects model for incomplete data from 
#' labeling-based quantitative proteomics experiments. The Annals of Applied Statistics, 11(1), 114-138. \doi{10.1214/16-AOAS994}
#' 
#' @export
#' @importFrom stats cov
#' @examples
#' data(sim_dat)
#'
#' Z = sim_dat$X[, 1, , drop = FALSE]
#' fit0 = mixEMM(Ym = sim_dat$Ym, Xm = sim_dat$X, Zm = Z, gamma = 0.14)


mixEMM <- function(Ym, Xm, Zm, gamma, maxIter = 100, tol = 0.001) {
  
  ep = 0.1
  N = nrow(Ym)
  p = ncol(Ym)
  K = dim(Xm)[2]
  h = dim(Zm)[2]
  
  ### for each cluster, assume missing happens to either all or none samples  
  idxObs = which(rowSums(is.na(Ym)) != p)
  idxMis = which(rowSums(is.na(Ym)) == p)
  
  I.m = matrix(1, nrow = p, ncol = p)
  I.v = matrix(1, nrow = p, ncol = 1)
  
  ############ initial values
  #### variance from available cases
  SIG.obs <- cov(Ym[idxObs, ], use = "pairwise")  
  delta2.hat = max(mean(SIG.obs[upper.tri(SIG.obs)]), ep)
  D = diag(delta2.hat, h, h)
  
  sigma2.hat = max(mean(diag(SIG.obs)[-1]) - delta2.hat, ep)
  sigma0.hat = max(diag(SIG.obs)[1] - delta2.hat, ep)
  R = diag(c(sigma0.hat, rep(sigma2.hat, p - 1)))
  
  Ym.hat = Ym
  
  ZmO = array(Zm[idxObs, , ], dim = c(length(idxObs), h, p))
  b0 = matrix(0, h, length(idxObs))
  Zi = matrix(1, p, h)
  SS = R + Zi %*% D %*% t(Zi)
  
  #### using available cases
  alpha.result = get.alpha(Ym[idxObs, ], Xm[idxObs, , ], ZmO, SS, b = b0) 
  alpha.hat = matrix(alpha.result$est, ncol = 1)
  alpha.var = matrix(alpha.result$var, nrow = nrow(alpha.hat))
  
  ############# constant
  gamma0 = 0
  diff = 999
  iter = 0
  nllik = NULL
  
  dellis = delta2.hat
  siglis = c(sigma0.hat, sigma2.hat)
  alphalis = alpha.hat
  lliklis = NULL
  
  
  ############# START the ECM algorithm
  
  while (iter < maxIter & diff > tol) {
    iter = iter + 1
    
    ###################### E-step:
    
    
    Ebi = matrix(NA, h, N)
    Varbi = array(NA, dim = c(h, h, N))
    
    Eei = matrix(NA, nrow = N, ncol = p)
    Varei = array(NA, dim = c(p, p, N))
    
    
    for (i in 1:N) {
      
      pi = sum(!is.na(Ym[i, ]))
      if (pi != 0) {
        #### E(bi|yobs, theta.hat) E(ei|yobs, theta.hat)
        
        idxo = which(!is.na(Ym[i, ]))
        pi = length(idxo)
        Yi = t(matrix(Ym[i, idxo], ncol = pi))
        Xi = t(matrix(Xm[i, , idxo], ncol = pi))
        Zi = t(matrix(Zm[i, , idxo], ncol = pi))
        Sigma.i = Zi %*% D %*% t(Zi) + matrix(R[idxo, idxo], ncol = pi)
        Wi = solve(Sigma.i)
        
        Ebi[, i] = D %*% t(Zi) %*% Wi %*% (Yi - Xi %*% alpha.hat)
        Eei[i, idxo] = as.vector(Yi - Xi %*% alpha.hat - Zi %*% 
                                   Ebi[, i])
        
        Varbi[, , i] = D - D %*% t(Zi) %*% Wi %*% Zi %*% D
        Varei[idxo, idxo, i] = Zi %*% Varbi[, , i] %*% t(Zi)
        
      } else {
        
        #### E(Yi|Mi=1, theta.hat) E(bi|Mi=1, theta.hat) E(ei|Mi=1, theta.hat)
        
        Xi = t(matrix(Xm[i, , ], ncol = p))
        Zi = t(matrix(Zm[i, , ], ncol = p))
        Sigma.i = Zi %*% D %*% t(Zi) + R
        Ebi[, i] = -gamma * rowMeans(D %*% t(Zi))
        Eei[i, ] = -gamma * colMeans(R)
        Varei[, , i] = R
        Varbi[, , i] = D
        Ym.hat[i, ] = t(Xi %*% alpha.hat - gamma * rowMeans(Sigma.i))
      }
    }
    
    ################################### M step
    
    
    D = (Ebi %*% t(Ebi) + apply(Varbi, c(1, 2), sum))/N
    
    ### R is from last iteration
    alpha.result = get.alpha(Ym.hat, Xm, Zm, R, Ebi)  
    alpha.hat = matrix(alpha.result$est, ncol = 1)
    for (i in 1:N) {
      idxo = which(!is.na(Ym.hat[i, ]))
      pi = length(idxo)
      Yi = t(matrix(Ym.hat[i, idxo], ncol = pi))
      Xi = t(matrix(Xm[i, , idxo], ncol = pi))
      Zi = t(matrix(Zm[i, , idxo], ncol = pi))
      Eei[i, idxo] = as.vector(Yi - Xi %*% alpha.hat - Zi %*% Ebi[, i])
      ### update the residue using current alpha
    }  
    
    df0 = sum(!is.na(Eei[, 1]))
    df1 = sum(!is.na(Eei[, -1]))
    sigma0.hat = (sum(Varei[1, 1, ], na.rm = T) + sum(Eei[, 1]^2, na.rm = T))/df0
    sigma2.hat = (sum(diag(apply(Varei[-1, -1, ], c(1, 2), sum, na.rm = T))) + 
                    sum(Eei[, -1]^2, na.rm = T))/df1
    R = diag(c(sigma0.hat, rep(sigma2.hat, p - 1)))
    
    dellis = c(dellis, D)
    siglis = cbind(siglis, c(sigma0.hat, sigma2.hat))
    alphalis = cbind(alphalis, alpha.hat)
    
    nllikt = target.f(Ym = Ym, Xm = Xm, Zm = Zm, alpha = alpha.hat, 
                      R = R, D = D, gamma0 = gamma0, gamma = gamma)
    nllik = c(nllik, nllikt$result)
    
    if (length(nllik) >= 2) 
      diff = abs((nllik[iter - 1] - nllik[iter]))
    if (is.na(diff)) 
      diff = 999
  }
  
  V = matrix(0, K, K)
  for (i in idxObs) {
    pi = sum(!is.na(Ym.hat[i, ]))
    idxo = which(!is.na(Ym.hat[i, ]))
    ## Zi is pi by h
    Zi = t(matrix(Zm[i, , idxo], ncol = pi))  
    Sigma.i = Zi %*% D %*% t(Zi) + as.matrix(R[idxo, idxo])
    Wi = solve(Sigma.i)
    Xii = Xm[i, , idxo]
    V = V + Xii %*% Wi %*% t(Xii)
  }
  alpha.var = solve(V)
  
  return(list(RE = as.vector(Ebi), loglikelihood = nllik[length(nllik)], alpha.hat = as.vector(alpha.hat), 
              alpha.se = sqrt(diag(alpha.var)), 
              sigma2.hat = sigma2.hat, sigma0.hat = sigma0.hat, D = D))
}


################################# the target (observed-data) log-likelihood function

target.f = function(Ym, Xm, Zm, alpha, R, D, gamma0, gamma) {
  
  
  result = 0
  p = ncol(Ym)
  N = nrow(Ym)
  K = dim(Xm)[2]
  
  D.inv = solve(D)
  
  
  I.v = matrix(1, nrow = p, ncol = 1)
  
  I1T = I2T = 0
  for (i in 1:N) {
    pi = sum(!is.na(Ym[i, ]))
    if (pi > 0) {
      idxo = which(!is.na(Ym[i, ]))
      Yi = t(matrix(Ym[i, idxo], ncol = pi))
      Xi = t(matrix(Xm[i, , idxo], ncol = pi))  ## Xi is pi by k
      Zi = t(matrix(Zm[i, , idxo], ncol = pi))  ## Zi is pi by h
      Sigma.i = Zi %*% D %*% t(Zi) + as.matrix(R[idxo, idxo])
      
      Wi = solve(Sigma.i)
      
      ## I2 = log(1-exp(gamma0-gamma*mean(Yi))) ## I3 do not change
      I1 = 0.5 * log(det(Wi)) - 0.5 * t(Yi - Xi %*% alpha) %*% Wi %*% 
        (Yi - Xi %*% alpha)
      
      I1T = I1T + I1
      
    } else {
      ## for the missing experiments, calculate I4, I5, I6
      Xi = t(matrix(Xm[i, , ], ncol = p))
      Zi = t(matrix(Zm[i, , ], ncol = p))
      
      Sigma.i = Zi %*% D %*% t(Zi) + R
      # Yi= Xi%*%alpha- gamma/p*rowSums(Sigma.i) I2 = gamma0-gamma*mean(Yi)
      I2 = gamma0 - gamma * mean(Xi %*% alpha) + gamma^2 * mean(Sigma.i)/2
      I2T = I2T + I2
      
    }
  }
  result = result + I1T + I2T
  
  return(list(result = result))
}



get.alpha <- function(Y, X, Z, R, b) {
  #### Ym: nxp Xm: nxKxp SigmaHat: pxp
  
  K = dim(X)[2]
  N = nrow(Y)
  xwx = matrix(0, K, K)
  xwy = matrix(0, K, 1)
  xwg = matrix(0, K, 1)
  
  for (i in 1:N) {
    pi = sum(!is.na(Y[i, ]))
    
    idxo = which(!is.na(Y[i, ]))
    Yi = t(matrix(Y[i, idxo], ncol = pi))
    Xi = t(matrix(X[i, , idxo], ncol = pi))  ## Xi is pi by k
    Zi = t(matrix(Z[i, , idxo], ncol = pi))  ## Zi is pi by h
    bi = matrix(b[, i], ncol = 1)
    Ri.inv = solve(R[idxo, idxo])
    
    xwx = xwx + t(Xi) %*% Ri.inv %*% Xi
    xwy = xwy + t(Xi) %*% Ri.inv %*% (Yi - Zi %*% bi)
  }
  var.alpha = solve(xwx)
  ooe = var.alpha %*% xwy
  
  return(list(est = ooe, var = var.alpha))
}

