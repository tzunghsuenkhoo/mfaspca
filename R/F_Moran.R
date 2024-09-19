
#' Title
#'
#' @param mData An object of of \code{\link[multifunData]{funData}} containing the multivariate
#'              areal spatial functional data observed.
#' @param W     A spatial weight matrix
#' @param mylist
#' @param Function_pr
#' @param k
#' @param l
#'
#' @return
#' @export
#'
#' @examples
FMORAN_I <- function(mData,W,mylist,Function_pr, k,l)
{
  # number of components
  p <- length(mData)
  # number of observations
  N <- funData::nObs(mData)
  # calculating the functional Moran I
  Fpm   <- list()
  Tt    <- mData[[1]]@argvals[[1]]
  Fpms1 <- matrix(0, nrow=length(Tt),ncol=p)
  Fpms2 <- matrix(0, nrow=length(Tt),ncol=p)
  ## functional Moran
  for(j in seq_len(p)) {
    Tt       <- mData[[j]]@argvals[[1]]
    vec      <- mylist[[j]]
    FUNC     <- Function_pr[[j]]
    Fpm[[j]] <- MFPCA::univExpansion(type = "default", scores =vec, functions = FUNC)
    for(i in c(Tt)) {
      Fpms1[i,j] = (t(Fpm[[j]]@X)[i,]%*%W%*%Fpm[[j]]@X[,i])
      Fpms2[i,j] = (t(Fpm[[j]]@X)[i,]%*%Fpm[[j]]@X[,i])
    }
  }

  Xim1    <- matrix(apply(Fpms1, 1, sum),nrow=1)
  Xim2    <- matrix(apply(Fpms2, 1, sum),nrow=1)
  Xim     <- matrix(c(Xim1/Xim2),nrow=1)
  # multivariate functional Moran
  fmmoram <- funData::funData(argvals = Fpm[[1]]@argvals,X = Xim)

  # bivariate functional Moran
  Fpmb   <- list()
  Fpmsb1 <- matrix(0, nrow=length(Tt),ncol=2)
  Fpmsb2 <- matrix(0, nrow=length(Tt),ncol=2)
  # functional Moran

  Tt    <- mData[[k]]@argvals[[1]]
  veck  <- mylist[[k]]
  FUNCk <- Function_pr[[k]]
  vecl  <- mylist[[l]]
  FUNCl <- Function_pr[[l]]
  Fpmb[[1]] = MFPCA::univExpansion(type = "default", scores =veck, functions = FUNCk)
  Fpmb[[2]] = MFPCA::univExpansion(type = "default", scores =vecl, functions = FUNCl)
  for(i in c(Tt)) {
    Fpmsb1[i,1] =(t(Fpmb[[1]]@X)[i,]%*%W%*%Fpmb[[1]]@X[,i])
    Fpmsb2[i,1] =(t(Fpmb[[1]]@X)[i,]%*%Fpmb[[1]]@X[,i])
    Fpmsb1[i,2] =(t(Fpmb[[2]]@X)[i,]%*%W%*%Fpmb[[2]]@X[,i])
    Fpmsb2[i,2] =(t(Fpmb[[2]]@X)[i,]%*%Fpmb[[2]]@X[,i])
  }

  Ximb1 <- matrix(apply(Fpmsb1, 1, sum),nrow=1)
  Ximb2 <- matrix(apply(Fpmsb2, 1, sum),nrow=1)
  Ximb <- matrix(c(Ximb1/Ximb2),nrow=1)
  fmmoramb <- funData::funData(argvals = Fpmb[[1]]@argvals,X = Ximb) #multivariate functional Moran

  #create the list of the output
  res <- list(FMMoran=fmmoram, BVMoran=fmmoramb)
  return(res)
}
