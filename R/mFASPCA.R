#' Estimation of univariate functional principal components under CBS
#'
#' This function gives the estimation of the FPCA under CBS which is a slightly adapted version of
#' .PACE-function in the package \code{\link[funData]{funData}}. It is an internal function, and it is called by PACE_S.
#'
#' @param X A vector of xValues
#' @param Y A matrix of observed functions (by row).
#' @param W A weight matrix
#' @param Y.pred A matrix of functions (by row) to be approximated using the
#'               functional principal components.
#' @param nbasis An integer, giving the number of B-spline basis to use.
#'               Defaults to \code{10}.
#' @param pve    A value between 0 and 1, giving the percentage of variance
#'               explained in the data by the functional principal components. This value is
#'               used to choose the number of principal components. Defaults to \code{0.99}
#' @param npc    The number of principal components to be estimated. Defaults to
#'               \code{NULL}. If given, this overrides \code{pve}.
#' @param makePD Logical, should positive definiteness be enforced for the
#'   covariance estimate? Defaults to \code{FALSE}.
#' @param cov.weight.type  The type of weighting used for the smooth covariance
#'   estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively,
#'   \code{"counts"} (corresponds to \code{fpca.sc} in \strong{refund}) weights the pointwise estimates of the covariance function
#'   by the number of observation points.
#'
#' @return  \item{fit}{The approximation of \code{Y.pred} (if \code{NULL}, the
#'   approximation of \code{Y}) based on the functional principal components.}
#'   \item{scores}{A matrix containing the estimated scores (observations by
#'   row).}
#'   \item{mu}{The estimated mean function.}
#'   \item{efunctions}{A matrix
#'   containing the estimated eigenfunctions (by row).}
#'   \item{evalues}{The estimated eigenvalues.}
#'   \item{npc}{The number of principal comopnents that
#'   were calculated.}
#'   \item{sigma2}{The estimated variance of the measurement
#'   error.}
#'   \item{estVar}{The estimated smooth variance function of the data.}
#'
#' @references  Happ, C. (2021). Multivariate functional principal component
#'   analysis for data observed on different dimensional domains, R package version
#'   Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
#'   Multilevel functional principal component analysis. Annals of Applied
#'   Statistics, 3, 458--488. Yao, F., Mueller, H.-G., and Wang, J.-L. (2005).
#'   Functional data analysis for sparse longitudinal data. Journal of the
#'   American Statistical Association, 100, 577--590.
.PACE_S <- function(X, Y, W, Y.pred = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{
  if (is.null(Y.pred))
    Y.pred = Y
  D = NCOL(Y)
  if(D != length(X)) # check if number of observation points in X & Y are identical
    stop("different number of (potential) observation points differs in X and Y!")
  I = NROW(Y)
  I.pred = NROW(Y.pred)
  d.vec = rep(X, each = I) # use given X-values for estimation of mu
  # weight matrix
  W <- (1/2)*(W+t(W))
  W <- W/sum(W)# weight matrix
  Y.weight = W%*%Y
  gam0 = mgcv::gam(as.vector(Y.weight) ~ s(d.vec, k = nbasis))
  mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec = X))
  Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
  cov.sum = cov.count = cov.mean = matrix(0, D, D)
  for (i in seq_len(I)) {
    obs.points = which(!is.na(Y.weight[i, ]))
    cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
    cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
  }
  G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)

  diag.G0 = diag(G.0)
  diag(G.0) = NA
  row.vec = rep(X, each = D) # use given X-values
  col.vec = rep(X, D) # use given X-values
  cov.weights <- switch(cov.weight.type,
                        none = rep(1, D^2),
                        counts = as.vector(cov.count),
                        stop("cov.weight.type ", cov.weight.type, " unknown in smooth covariance estimation"))

  npc.0 = matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0)~te(row.vec, col.vec, k = nbasis), weights = cov.weights),
                                   newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
  npc.0 = (npc.0 + t(npc.0))/2
  # no extra-option (useSymm) as in fpca.sc-method
  if (makePD) { # see fpca.sc
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE,
                            do2eigen = TRUE, trace = options()$verbose)
      as.matrix(tmp$mat)
    }
  }

  ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Chapter 8)
  w <- funData::.intWeights(X, method = "trapezoidal")
  Wsqrt <- diag(sqrt(w))
  Winvsqrt <- diag(1/(sqrt(w)))
  V <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values
  ###
  evalues = replace(evalues, which(evalues <= 0), 0)

  npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
  efunctions = matrix(Winvsqrt%*%eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values[seq_len(npc)]  # use correct matrix for eigenvalue problem
  cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  ### numerical integration for estimation of sigma2
  T.len <- X[D] - X[1] # total interval length
  T1.min <- min(which(X >= X[1] + 0.25*T.len)) # left bound of narrower interval T1
  T1.max <- max(which(X <= X[D] - 0.25*T.len)) # right bound of narrower interval T1
  DIAG = (diag.G0 - diag(cov.hat))[T1.min :T1.max] # function values
  # weights
  w <- funData::.intWeights(X[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(1/(X[T1.max]-X[T1.min]) * sum(DIAG*w, na.rm = TRUE), 0) #max(1/T.len * sum(DIAG*w), 0)
  ####
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  fit = matrix(0, nrow = I.pred, ncol = D)
  scores = matrix(NA, nrow = I.pred, ncol = npc)
  # no calculation of confidence bands, no variance matrix
  for (i.subj in seq_len(I.pred)) {
    obs.points = which(!is.na(Y.pred[i.subj, ]))
    if (sigma2 == 0 & length(obs.points) < npc) {
      stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
    }
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points),
                  ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
    scores[i.subj, ] = ZtZ_sD.inv %*% crossprod(Zcur, Y.tilde[i.subj, obs.points])
    fit[i.subj, ] = t(as.matrix(mu)) + tcrossprod(scores[i.subj, ], efunctions)
  }
  ret.objects = c("fit", "scores", "mu", "efunctions", "evalues",
                  "npc", "sigma2") # add sigma2 to output
  ret = lapply(seq_len(length(ret.objects)), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  ret$estVar <- diag(cov.hat)
  return(ret)
}

#' Estimation of univariate functional principal components under CBS. This is an slightly adapted version of the
#' PACE function in the package \code{\link[funData]{funData}}. This is an internal function
#' which is called by the mfasPCA function.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}} or
#'   \code{\link[funData]{irregFunData}} containing the functional data
#'   observed, for which the functional principal component analysis is
#'   calculated. If the data is sampled irregularly (i.e. of class
#'   \code{\link[funData]{irregFunData}}), \code{funDataObject} is transformed
#'   to a \code{\link[funData]{funData}} object first.
#' @param W A weight matrix
#' @param predData An object of class \code{\link[funData]{funData}}, for which
#'   estimated trajectories based on a truncated Karhunen-Loeve representation
#'   should be estimated. Defaults to \code{NULL}, which implies prediction for
#'   the given data.
#' @param nbasis An integer, representing the number of  B-spline basis
#'   functions used for estimation of the mean function and bivariate smoothing
#'   of the covariance surface. Defaults to \code{10} (cf.
#'   \code{fpca.sc} in \strong{refund}).
#' @param pve A numeric value between 0 and 1, the proportion of variance
#'   explained: used to choose the number of principal components. Defaults to
#'   \code{0.99} (cf. \code{fpca.sc} in \strong{refund}).
#' @param npc An integer, giving a prespecified value for the number of
#'   principal components. Defaults to \code{NULL}. If given, this overrides
#'   \code{pve} (cf. \code{fpca.sc} in \strong{refund}).
#' @param makePD Logical: should positive definiteness be enforced for the
#'   covariance surface estimate? Defaults to \code{FALSE} (cf.
#'   \code{fpca.sc} in \strong{refund}).
#' @param cov.weight.type The type of weighting used for the smooth covariance
#'   estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively,
#'   \code{"counts"} (corresponds to \code{fpca.sc} in \strong{refund}) weights the
#'   pointwise estimates of the covariance function by the number of observation
#'   points.
#'
#' @return \item{mu}{A \code{\link[funData]{funData}} object with one
#'   observation, corresponding to the mean function.}
#'   \item{values}{A vector
#'   containing the estimated eigenvalues.}
#'   \item{functions}{A \code{\link[funData]{funData}} object containing the estimated functional
#'   principal components.}
#'   \item{scores}{An matrix of estimated scores for the
#'   observations in \code{funDataObject}. Each row corresponds to the scores of
#'   one observation.}
#'   \item{fit}{A \code{\link[funData]{funData}} object
#'   containing the estimated trajectories based on the truncated Karhunen-Loeve
#'   representation and the estimated scores and functional principal components
#'   for \code{predData} (if this is not \code{NULL}) or \code{funDataObject}
#'   (if \code{predData} is \code{NULL}).}
#'   \item{npc}{The number of functional
#'   principal components: either the supplied \code{npc}, or the minimum number
#'   of basis functions needed to explain proportion \code{pve} of the variance
#'   in the observed curves.}
#'   \item{sigma2}{The estimated measurement error variance.}
#'   \item{estVar}{The estimated smooth
#'   variance function of the data.}
#'@references  Happ, C. (2021). Multivariate functional principal component
#'   analysis for data observed on different dimensional domains, R package version
#'   Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
#'   Multilevel functional principal component analysis. Annals of Applied
#'   Statistics, 3, 458--488. Yao, F., Mueller, H.-G., and Wang, J.-L. (2005).
#'   Functional data analysis for sparse longitudinal data. Journal of the
#'   American Statistical Association, 100, 577--590.
PACE_S <- function(funDataObject, W, predData = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{

  # check inputs
  if(! class(funDataObject) %in% c("funData", "irregFunData"))
    stop("Parameter 'funDataObject' must be a funData or irregFunData object.")
  if(funData::dimSupp(funDataObject) != 1)
    stop("PACE_S: Implemented only for funData objects with one-dimensional support.")
  if(methods::is(funDataObject, "irregFunData")) # for irregular functional data, use funData representation
    funDataObject <- funData::as.funData(funDataObject)

  if(is.null(predData))
    Y.pred = NULL # use only funDataObject
  else
  {
    if(!isTRUE(all.equal(funDataObject@argvals, predData@argvals)))
      stop("PACE_S: funDataObject and predData must be defined on the same domains!")

    Y.pred = predData@X
  }

  if(!all(is.numeric(nbasis), length(nbasis) == 1, nbasis > 0))
    stop("Parameter 'nbasis' must be passed as a number > 0.")

  if(!all(is.numeric(pve), length(pve) == 1, 0 <= pve, pve <= 1))
    stop("Parameter 'pve' must be passed as a number between 0 and 1.")

  if(!is.null(npc) & !all(is.numeric(npc), length(npc) == 1, npc > 0))
    stop("Parameter 'npc' must be either NULL or passed as a number > 0.")

  if(!is.logical(makePD))
    stop("Parameter 'makePD' must be passed as a logical.")

  if(!is.character(cov.weight.type))
    stop("Parameter 'cov.weight.type' must be passed as a character.")


  res <- .PACE_S(X = funDataObject@argvals[[1]], funDataObject@X, W ,Y.pred = Y.pred,
                 nbasis = nbasis, pve = pve, npc = npc, makePD = makePD,
                 cov.weight.type = cov.weight.type)

  return(list(mu = funData::funData(funDataObject@argvals, matrix(res$mu, nrow = 1)),
              values = res$evalues,
              functions = funData::funData(funDataObject@argvals, t(res$efunctions)),
              scores = res$scores,
              fit = funData::funData(funDataObject@argvals, res$fit),
              npc = res$npc,
              sigma2 = res$sigma2,
              estVar = funData::funData(funDataObject@argvals, matrix(res$estVar, nrow = 1))
  ))
}

#' This function calculates a multivariate functional areal spatial functional principal component analysis
#' (mfasPCA) based on i.i.d. multivariate areal spatial functional data.
#'
#' @param mData An object of of \code{\link[multifunData]{multifunData}} containing the multivariate
#'              areal spatial functional data observed.
#' @param W     A spatial weight matrix.
#' @param J     The number of multivariate functional principal components to be calculated.
#' @param fit   Logical: If fit is TRUE, then the fitting of the multivariate areal spatial functional data
#'              is performed based on the calculated fasPCs.
#'
#' @return  \item{values}{Eigenvalues}
#' \item{functions}{An object of \code{\link[multifunData]{multifunData}}. }
#' \item{scores}{ A vector of eigenvalues.}
#' \item{vectors}{A matrix of eigenvectors.}
#' \item{argvals}{A vector of numerics, defining a grid on the interval for which the basis functions are computed}
#' \item{meanFunction}{The estimated mean function.}
#' \item{fit}{fit = TRUE returns an object of \code{\link[multifunData]{multifunData}} containing
#'            the approximation of {X_i} using the mfasPCs with J components.}
#' \item{MaxP}{The total number of univariate basis function.}
#' @export
#'
#' @examples

mfasPCA <- function(mData,W,J,fit = TRUE)
{
  # number of components
  p <- length(mData)
  # number of observations
  N <- funData::nObs(mData)
  # creating empty lists
  mylist <- list()
  Function_pr <- list()
  moy <- list()
  # Weight matrix
  W <- (1/2)*(W+t(W))
  W <- W/sum(W)
  # Univariate FPCA
  for (j in seq_len(p)) {
    FPCA_univ        <- PACE_S(funDataObject=mData[[j]], W, predData = NULL, nbasis = 15, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
    mylist[[j]]      <- FPCA_univ[["scores"]]    # put all scores in the list
    Function_pr[[j]] <- FPCA_univ[["functions"]]
    moy[[j]]         <- FPCA_univ[["mu"]]
  }
  Xi <- do.call("cbind",mylist) # combine all scores
  Pmax <- vapply(mylist, function(x){dim(x)[2]}, FUN.VALUE = 0)
  # Calculation of mean function
  m <- funData::multiFunData(do.call("list", moy))
  # Multivariate FPCA
  if(J > sum(Pmax))
  {
    J <- sum(Pmax)
    warning("Function mfasPCA: total number of univariate basis functions is smaller than given J. J was set to ", sum(Pmax), ".")
  }
  Z <- t(Xi)%*%W%*%Xi
  # Calculation of PCA on scores
  e <- eigen(Z)
  values    <- Re(e$values[seq_len(J)])
  vectors   <- Re(e$vectors[,seq_len(J)])
  # Calculate eigenfunctions
  PmaxCum   <- cumsum(c(0,Pmax))
  Fp        <- list()
  for(j in seq_len(p)) {
    Tt      <- mData[[j]]@argvals[[1]]
    vec     <- t(vectors)
    V       <- as.matrix(vec[ ,PmaxCum[j]+seq_len(Pmax[j])])
    FUNC    <- Function_pr[[j]]
    Fp[[j]] <- MFPCA::univExpansion(type = "default", scores = V, functions = FUNC)
  }
  multiFuns <- funData::multiFunData(do.call("list",Fp))
  # Calculate scores
  scores    <- Xi %*% vectors

  # Calculate the fitted multivariate data
  if(fit){
    univExp <- list()
    for(j in seq_len(p)) {
      univExp[[j]] = MFPCA::univExpansion(type = "default", scores = scores, functions = multiFuns[[j]])
    }

    Xfit <- m + funData::multiFunData(do.call("list", univExp))
    names(Xfit) <- names(mData)
  }
  # Create the list of the output
  res <- list(values = values,
              functions = multiFuns,
              scores = scores,
              vectors = vectors[seq_len(J),seq_len(J)],
              argval = lapply(1:p, function(i) mData[[i]]@argvals[[1]]),
              meanFunction=m,
              fit = Xfit,
              MaxP = Pmax, mylistX=mylist,Function_proj=Function_pr)

  return(res)
}
