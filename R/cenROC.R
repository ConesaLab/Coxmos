#' Estimation of the time-dependent ROC curve for right censored survival data
#'
#' @description  This function computes the time-dependent ROC curve for right censored survival data using the cumulative sensitivity and dynamic specificity definitions.
#'  The ROC curves can be either empirical (non-smoothed) or smoothed with/wtihout boundary correction. It also calculates the time-dependent area under the ROC curve (AUC).
#'  Edited by Pedro Salguero to remove the PLOT argument.
#' @usage cenROC(Y, M, censor, t, U = NULL, h = NULL, bw = "NR", method = "tra",
#'     ktype = "normal", ktype1 = "normal", B = 0, alpha = 0.05, plot = "TRUE")
#' @param Y The numeric vector of event-times or observed times.
#' @param M The numeric vector of marker values for which the time-dependent ROC curves is computed.
#' @param censor The censoring indicator, \code{1} if event, \code{0} otherwise.
#' @param t A scaler time point at which the time-dependent ROC curve is computed.
#' @param U The vector of grid points where the ROC curve is estimated. The default is a sequence of \code{151} numbers between \code{0} and \code{1}.
#' @param bw A character string specifying the bandwidth estimation method for the ROC itself. The possible options are "\code{NR}" for the normal reference, the plug-in "\code{PI}" and the cross-validation "\code{CV}". The default is the "\code{NR}" normal reference method. The user can also introduce a numerical value.
#' @param h A scaler for the bandwidth of Beran's weight calculaions. The defualt is the value obtained by using the method of Sheather and Jones (1991).
#' @param method The method of ROC curve estimation. The possible options are "\code{emp}" emperical metod; "\code{untra}" smooth without boundary correction and "\code{tra}" is smooth ROC curve estimation with boundary correction. The default is the "\code{tra}" smooth ROC curve estimate with boundary correction.
#' @param ktype A character string giving the type kernel distribution to be used for smoothing the ROC curve: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @param ktype1 A character string specifying the desired kernel needed for Beran weight calculation. The possible options are "\code{normal}", "\code{epanechnikov}", "\code{tricube}", "\code{boxcar}", "\code{triangular}", or "\code{quartic}". The defaults is "\code{normal}" kernel density.
#' @param B The number of bootstrap samples to be used for variance estimation. The default is \code{0}, no variance estimation.
#' @param alpha The significance level. The default is \code{0.05}.
#' @param plot The logical parameter to see the ROC curve plot. The default is \code{TRUE}.
#' @details The empirical (non-smoothed) ROC estimate and the smoothed ROC estimate with/without boundary correction can be obtained using this function.
#' The smoothed ROC curve estimators require selecting two bandwidth parametrs: one for Beran’s weight calculation and one for smoothing the ROC curve.
#' For the latter, three data-driven methods: the normal reference "\code{NR}", the plug-in "\code{PI}" and the cross-validation "\code{CV}" were implemented.
#' To select the bandwidth parameter needed for Beran’s weight calculation, by default, the plug-in method of Sheather and Jones (1991) is used but it is also possible introduce a numeric value.
#' See Beyene and El Ghouch (2020) for details.
#' @return Returns the following items:
#' @return    \code{ROC     } The vector of estimated ROC values. These will be numeric numbers between zero
#' @return    \code{        } and one.
#' @return    \code{U       } The vector of grid points used.
#' @return    \code{AUC      } A data frame of dimension \eqn{1 \times 4}. The columns are: AUC, standard error of AUC, the lower
#' @return    \code{         }               and upper limits of bootstrap CI.
#' @return    \code{bw       } The computed value of bandwidth. For the empirical method this is always \code{NA}.
#' @return    \code{Dt      } The vector of estimated event status.
#' @return    \code{M       } The vector of Marker values.
#' @importFrom stats pnorm qnorm quantile approx bw.SJ integrate sd
#' @examples library(cenROC)
#'
#' data(mayo)
#' cenROC(Y=mayo$time, M=mayo$mayoscore5, censor=mayo$censor, t=365*6)$AUC
#'
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @references Sheather, S. J. and Jones, M. C. (1991). A Reliable data-based bandwidth selection method for kernel density estimation. \emph{Journal of the Royal Statistical Society}. Series B (Methodological) 53(3): 683–690.

cenROC <- function(Y, M, censor, t, U = NULL, h = NULL, bw = "NR",  method = "tra", ktype = "normal", ktype1 = "normal", B = 0, alpha = 0.05, plot = F)
{
  if (is.null(U)) {U <- seq(0, 1, length.out = 151)}
  if (!is.vector(Y, mode = "numeric") |
      !is.vector(M, mode = "numeric") |
      !is.vector(censor, mode = "numeric"))
    print("Error! all numeric vectors Y, M and censor should be specified")
  else{
    Dt <- Csurv(Y = Y, M = M, censor = censor, t = t, h = h, kernel = ktype1)$positive
    estim <- RocFun(U = U, D = Dt, M = M, method = method, bw = bw, ktype = ktype)
    ROC <- estim$roc
    AUCc <- 1 - estim$auc
    AUC <- data.frame(AUC = round(AUCc, 4) , sd = NA, LCL = NA, UCL = NA)
  }

  if (B > 0){
    data <- data.frame(Y=Y, M=M, cen=censor)
    aucb <- NULL
    rocb <- matrix(NA, nrow = length(U), ncol = B)
    for (i in 1:B){
      bootsample <- sample(1:nrow(data), nrow(data), replace=TRUE)
      dat <- data[bootsample, ]
      Dt <- Csurv(Y = dat$Y, M = dat$M, censor = dat$cen, t = t, h = h, kernel = ktype1)$positive
      estim <- RocFun(U = U, D = Dt, M = dat$M, method = method, bw = bw, ktype = ktype)
      aucb[i] <- 1 - estim$auc
      rocb[, i] <- estim$roc
    }
    SP <- unname(quantile(aucb, p=c(alpha/2, 1-alpha/2), na.rm = TRUE))
    AUC <- data.frame(AUC = round(AUCc, 4) , sd = round(sd(aucb, na.rm = TRUE), 4), LCL = round(SP[1], 4), UCL = round(SP[2], 4))
    qroc <- unname(apply(rocb, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE))
    ROC <- data.frame(ROC = ROC, LCL = qroc[1, ], UCL = qroc[2, ])
  }
  return(list(AUC = AUC, ROC = ROC, U = U, Dt = Dt, M = M, bw = estim$bw))
}

#' Survival probability conditional to the observed data estimation for right censored data.
#'
#' @param Y The numeric vector of event-times or observed times.
#' @param M The numeric vector of marker values for which we want to compute the time-dependent ROC curves.
#' @param censor The censoring indicator, \code{1} if event, \code{0} otherwise.
#' @param t A scaler time point at which we want to compute the time-dependent ROC curve.
#' @param h A scaler for the bandwidth of Beran's weight calculaions. The defualt is using the method of Sheather and Jones (1991).
#' @param kernel A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", , "\code{tricube}", "\code{boxcar}", "\code{triangular}", or "\code{quartic}". The defaults is "\code{normal}" kernel density.
#' @return Return a vectors:
#' @return \code{positive    }    \code{P(T<t|Y,censor,M)}.
#' @return \code{negative    }     \code{P(T>t|Y,censor,M)}.
#' @references Beyene, K. M. and El Ghouch A. (2019). Smoothed time-dependent ROC curves for right-censored survival data. <\url{https://dial.uclouvain.be/pr/boreal/object/boreal:219643}>.
#' @references Li, Liang, Bo Hu and Tom Greene (2018).  A simple method to estimate the time-dependent receiver operating characteristic curve and the area under the curve with right censored data, Statistical Methods in Medical Research, 27(8): 2264-2278.
#' @references Pablo Martínez-Camblor and Gustavo F. Bayón and Sonia Pérez-Fernández (2016). Cumulative/dynamic roc curve estimation, Journal of Statistical Computation and Simulation, 86(17): 3582-3594.
#' @keywords internal

Csurv <- function(Y, M, censor, t, h = NULL, kernel="normal") {
  if (is.null(h)) {
    h <- bw.SJ(M, method = "dpi")
  }
  if(kernel=="normal"){
    kernel <- "gaussian"
  }
  n <- length(M)
  positive <- rep(NA, n)
  for (i in 1:n) {
    if (Y[i] > t) {
      positive[i] <- 0

    } else {
      if (censor[i] == 1) {
        positive[i] <- 1

      } else {
        St <- Beran(time = Y, status = censor, covariate = M, x = M[i], y = t, kernel = kernel, bw = h)
        Sy <- Beran(time = Y, status = censor, covariate = M, x = M[i], y = Y[i], kernel = kernel, bw = h)
        if (Sy == 0) {
          positive[i] <- 1

        } else {
          positive[i] <- 1 - St / Sy

        }
      }
    }
  }
  negative <- 1 - positive

  return(list(positive = positive, negative = negative))

}

#'  ROC estimation function
#'
#' @param U The vector of grid points where the ROC curve is estimated.
#' @param D The event indicator.
#' @param M The numeric vector of marker values for which the time-dependent ROC curves is computed.
#' @param bw The bandwidth parameter for smoothing the ROC function. The possible options are \code{NR} normal reference method; \code{PI} plug-in method and \code{CV} cross-validation method. The default is the \code{NR} normal reference method.
#' @param method is the method of ROC curve estimation. The possible options are \code{emp} emperical metod; \code{untra} smooth without boundary correction and \code{tra} is smooth ROC curve estimation with boundary correction.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}".
#' @author Beyene K. Mehari and El Ghouch Anouar
#' @references Beyene, K. M. and El Ghouch A. (2019). Smoothed time-dependent ROC curves for right-censored survival data. <\url{https://dial.uclouvain.be/pr/boreal/object/boreal:219643}>.
#' @keywords internal

RocFun <- function(U, D, M, bw = "NR", method, ktype) {
  oM <- order(M)
  D <- (D[oM])
  nD <- length(D)
  sumD <- sum(D)
  Z <- 1 - cumsum(1 - D) / (nD - sumD)
  AUC <- sum(D * Z) / sumD
  if (method == "emp") {
    difmat <- (outer(U, Z, "-"))
    resul <- (difmat >= 0)
    roc1 <- sweep(resul, 2, D, "*")
    roc <- apply(roc1, 1, sum) / sumD
    bw1 <- NA
  }
  else if (method == "untra") {
    Zt <- Z
    Ut <- U
    Ztt <- Zt[D != 0]
    wt <- D[D != 0]
    bw1 <- wbw(X = Ztt, wt = wt, bw = bw, ktype = ktype)$bw
    difmat <- (outer(Ut, Ztt, "-")) / bw1
    resul <- kfunc(ktype = ktype, difmat = difmat)
    w <- wt / sum(wt)
    roc1 <- sweep(resul, 2, w, "*")
    roc <- apply(roc1, 1, sum)
  }
  else if (method == "tra") {
    mul <- nD / (nD + 1)
    Zt <- qnorm(mul * Z + (1 / nD ^ 2))
    Ut <- qnorm(mul * U + (1 / nD ^ 2))
    Ztt <- Zt[D != 0]
    wt <- D[D != 0]
    bw1 <- wbw(X = Ztt, wt = wt, bw = bw, ktype = ktype)$bw
    difmat <- (outer(Ut, Ztt, "-")) / bw1
    resul <- kfunc(ktype = ktype, difmat = difmat)
    w <- wt / sum(wt)
    roc1 <- sweep(resul, 2, w, "*")
    roc <- apply(roc1, 1, sum)
  }
  else{
    stop("The specified method is not correct.")
  }
  return(list(roc = roc, auc = AUC, bw = bw1))
}

#' Function to select the bandwidth parameter needed for smoothing the time-dependent ROC curve.
#'
#' @description This function computes the data-driven bandwidth value for smoothing the ROC curve.
#'              It contains three methods: the normal refrence, the plug-in and the cross-validation methods.
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param bw A character string specifying the bandwidth selection method. The possible options are "\code{NR}" for the normal reference, the plug-in "\code{PI}" and cross-validation "\code{CV}".
#' @param ktype A character string indicating the type of kernel function: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". Default is "\code{normal}" kernel.
#' @return Returns the estimated value for the bandwith parameter.
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @keywords internal

wbw <- function(X, wt, bw = "NR", ktype = "normal")
{
  if (is.numeric(bw)) {
    bwv <- bw
  } else if (bw == "CV") {
    bwv <- CV(X = X, wt = wt, ktype = ktype)$bw
  } else if (bw == "NR") { #default
    bwv <- NR(X = X, wt = wt, ktype = ktype)$bw
  } else if (bw == "PI") {
    bwv <- PI(X = X, wt = wt, ktype = ktype)$bw
  } else{
    print("Please check your bandwidth options")
  }
  return(list(bw = bwv))
}

#' Function to evaluate the matrix of data vector minus the grid points divided by the bandwidth value.
#'
#' @param difmat A numeric matrix of sample data (X) minus evaluation points (x0) divided by bandwidth value (bw).
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @return Returns the matrix resulting from evaluating \code{difmat}.
#' @keywords internal

kfunc <- function(ktype = "normal", difmat)
{
  if (ktype == "normal")
  {
    estim <- kfunction(ktype = "normal", X = difmat)
  }
  else if (ktype == "epanechnikov")
  {
    estim <- difmat
    low <- (difmat <= -1)
    up <- (difmat >= 1)
    btwn <- (difmat > -1 & difmat < 1)
    estim[low] <- 0
    estim[up] <- 1
    value <- estim[btwn]
    estim[btwn] <- kfunction(ktype = "epanechnikov", X = value)
  }
  else if (ktype == "biweight")
  {
    estim <- difmat
    low <- (difmat <= -1)
    up <- (difmat >= 1)
    btwn <- (difmat > -1 & difmat < 1)
    estim[low] <- 0
    estim[up] <- 1
    value <- estim[btwn]
    estim[btwn] <- kfunction(ktype = "biweight", X = value)
  }
  else if (ktype == "triweight")
  {
    estim <- difmat
    low <- (difmat <= -1)
    up <- (difmat >= 1)
    btwn <- (difmat > -1 & difmat < 1)
    estim[low] <- 0
    estim[up] <- 1
    value <- estim[btwn]
    estim[btwn] <- kfunction(ktype = "triweight", X = value)
  }
  return(estim)
}

#' The normal reference bandwidth selection for weighted data
#'
#' @description This function computes the data-driven bandwidth for smoothing the ROC (or distribution) function using the NR method of Beyene and El Ghouch (2020). This is an extension of the classical (unweighted) normal reference bandwith selection method to the case of weighted data.
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @return Returns the computed value for the bandwith parameter.
#' @details See Beyene and El Ghouch (2020) for details.
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @examples library(cenROC)
#'
#' X <- rnorm(100) # random data vector
#' wt <- runif(100) # weight vector
#'
#' ## Normal reference bandwidth selection
#' NR(X = X, wt = wt)$bw
#'


NR <- function (X, wt, ktype="normal") {
  nx <- length(X)
  mul <- (nx * sum(wt * wt)) / ((sum(wt)) ^ 2)
  stdv <- sqrt(wvar(X = X, wt = wt))
  IQR <- wIQR(X = X, wt = wt)
  sigma <- min(stdv, IQR / 1.349)
  c <- (4 * sqrt(pi) * (muro(ktype)$ro) / ((muro(ktype = ktype)$mu2) ^ 2)) ^ (1 / 3)
  wbw <- (c * sigma) * (mul ^ (1 / 3)) * nx ^ (-1 / 3)

  return(list(bw = wbw))
}


#' The plug-in bandwidth selection for weighted data
#'
#' @description This function computes the data-driven bandwidth for smoothing the ROC (or distribution) function using the PI method of Beyene and El Ghouch (2020). This is an extension of the classical (unweighted) direct plug-in bandwith selection method to the case of weighted data.
#'
#' @param X The numeric vector of random variable.
#' @param wt The non-negative weight vector.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @return Returns the computed value for the bandwith parameter.
#' @details See Beyene and El Ghouch (2020) for details.
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @examples library(cenROC)
#'
#' X <- rnorm(100) # random data vector
#' wt <- runif(100) # weight vector
#'
#' ## Plug-in bandwidth selection
#' PI(X = X, wt = wt)$bw
#'


PI <- function(X, wt, ktype="normal")
{
  ### bandwidth estimation ###################
  band <- function (X, wt, psi, ktype) {
    n <- length(X)
    mul <- (n * sum(wt * wt)) / ((sum(wt)) ^ 2)
    #### Estimates of E(wt^2)/((E(wt))^2)
    co <- ((muro(ktype = ktype)$ro) / ((muro(ktype = ktype)$mu2) ^ 2)) ^ (1 / 3)
    wbww <- (co) * (mul ^ (1 / 3)) * (n ^ (-1 / 3)) * ((-psi) ^ (-1 / 3))

    return(wbww)
  }

  ####### Estimation of psi  ################
  psi <- function(r, g) {
    n <- length(X)
    w <- (n * wt) / (sum(wt)) #### Estimate for wt/E(wt)
    ww <- outer(w, w, "*")
    aux <- outer(X, X, "-") / g
    aux <- (ww) * (dnorkernel(r, aux))
    result <- (sum(aux)) * (((g) ^ (-r - 1)) * (n ^ (-2)))
    return(result)
  }
  g0 <- NR(X, wt, ktype = ktype)$bw ### Intial bandwidt estimatioh using normal reference
  psi2 <- psi(2, g0) ### Plug-in estimate
  PIbw <- band(X, wt, psi2, ktype = ktype)

  return(list(bw = PIbw))
}


#############################################################################
## This function is modified from kerdiest R package CVbw function ##########
##  of Quintela-del-Rio and Estevez-Perez (2015)                   ##########
#############################################################################

#' The cross-validation bandwidth selection for weighted data
#'
#' @description This function computes the data-driven bandwidth for smoothing the ROC (or distribution) function using the CV method of Beyene and El Ghouch (2020). This is an extension of the classical (unweighted) cross-validation bandwith selection method to the case of weighted data.
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @return Returns the computed value for the bandwith parameter.
#' @details Bowman et al (1998) proposed the cross-validation bandwidth selection method for unweighted kernal smoothed distribution function. This method is implemented in the \code{R} package \code{kerdiest}.
#' We adapted this for the case of weighted data by incorporating the weight variable into the cross-validation function of Bowman's method. See Beyene and El Ghouch (2020) for details.
#'
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @references Bowman A., Hall P. and Trvan T.(1998). Bandwidth selection for the smoothing of distribution functions. \emph{Biometrika} 85:799-808.
#' @references Quintela-del-Rio, A. and Estevez-Perez, G. (2015). \code{kerdiest:} Nonparametric kernel estimation of the distribution function, bandwidth selection and estimation of related functions. \code{R} package version 1.2.
#' @examples
#' \dontrun{library(cenROC)
#'
#' X <- rnorm(100) # random data vector
#' wt <- runif(100) # weight vector
#'
#' ## Cross-validation bandwidth selection
#' CV(X = X, wt = wt)$bw
#'
#' }


CV <- function(X, wt, ktype = "normal")
{
  mul <- length(wt) / sum(wt)
  prob_quantile <- 0
  ss <- quantile(X, c(prob_quantile, 1 - prob_quantile))
  y <- seq(ss[1], ss[2], length.out = 100)
  seq_bws = seq((max(X) - min(X)) / 200, (max(X) - min(X)) / 2, length.out = 50)
  n_bws <- length(seq_bws)
  CVfunction <- numeric(length = n_bws)
  for (i in 1:n_bws)
  {
    integrand <- apply((t((mul * wt) * (t(outer(y, X, "-") >= 0))) - t(ker_dis_i(ktype = ktype, y = y, X = X, wt = wt, bw = seq_bws[i]))) ^ 2, 1, mean)
    CVfunction[i] <- integ(x = y, fx = integrand, method = "simps")
  }
  i0 <- which.min(CVfunction)
  CVbw_val <- seq_bws[i0]
  wbw <- CVbw_val
  return(list(bw = wbw))
}

#' Kernel distribution function
#'
#' @param X A numeric vector of sample data.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}".
#' @return Returns a vector resulting from evaluating X.
#' @keywords internal

kfunction <- function(ktype, X) {
  if (ktype == "normal") {
    result <- pnorm(X)

  }
  else if (ktype == "epanechnikov") {
    result <- (0.75 * X * (1 - (X ^ 2) / 3) + 0.5)

  }
  else if (ktype == "biweight") {
    result <- ((15 / 16) * X - (5 / 8) * X ^ 3 + (3 / 16) * X ^ 5 + 0.5)

  }
  else if (ktype == "triweight") {
    result <- ((35 / 32) * X - (35 / 32) * X ^ 3 + (21 / 32) * X ^ 5 - (5 / 32) * X ^ 7 + 0.5)
  }
  return(result)
}

#' Weighted inter-quartile range estimation
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @keywords internal

wIQR <- function(X, wt) {
  (wquantile(X = X, wt = wt, p = 0.75) - wquantile(X = X, wt = wt, p = 0.25))
}

#' The value of squared integral x^2 k(x) dx and integral x k(x) K(x) dx
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}".
#' @keywords internal

muro <- function(ktype)
{
  if (ktype == "normal") {
    ro <- 2 * 0.28209
    mu2 <- 1
  } else if (ktype == "epanechnikov") {
    ro <- 2 * 0.12857
    mu2 <- 1 / 5
  } else if (ktype == "biweight") {
    ro <- 2 * 0.10823
    mu2 <- 1 / 7
  } else if (ktype == "triweight") {
    ro <- 2 * 0.095183
    mu2 <- 1 / 9
  }

  return(list(ro = ro, mu2 = mu2))
}

#' Derivative of normal distribution
#'
#' @param X The numeric data vector.
#' @param ord The order of derivative.
#' @keywords internal

dnorkernel <- function(ord, X)
{
  if (ord == 2)
    # second derivative
    result <- (1 / (sqrt(2 * pi))) * exp(-(X ^ 2) / 2) * ((X ^ 2) - 1)
  else if (ord == 4)
    # fourth derivative
    result <- (1 / (sqrt(2 * pi))) * exp(-(X ^ 2) / 2) * (3 - (6 * (X ^ 2)) + X ^ 4)
  else if (ord == 6)
    # sixth derivative
    result <- (1 / (sqrt(2 * pi))) * exp(-(X ^ 2) / 2) * (X ^ 6 - (15 * (X ^ 4)) + (45 * (X ^ 2)) - 15)
  else if (ord == 8)
    # eighth derivative
    result <- (1 / (sqrt(2 * pi))) * exp(-(X ^ 2) / 2) * (X ^ 8 - (28 * (X ^ 6)) + (210 * (X ^ 4)) - (420 * (X ^ 2)) + 105)
  return(result)
}

#' Distribution function without the ith observation
#'
#' @param X The numeric data vector.
#' @param y The vector where the kernel estimation is computed.
#' @param wt The non-negative weight vector.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}".
#' @param bw A numeric bandwidth value.
#' @return Returns the estimated value for the bandwith parameter.
#' @author Kassu Mehari Beyene  and Anouar El Ghouch
#' @keywords internal

ker_dis_i <- function(X, y, wt, ktype, bw)
{
  n <- length(X);
  AUX <- matrix(0, n, n);
  zero <- rep(0, n);
  ww <- outer(wt, zero, "-");
  diag(ww) <- 0;
  den <- apply(ww, 2, sum);
  resu <- matrix(0, n, length(y));
  for (j in 1:length(y))
  {
    AUX <- matrix(rep.int(outer(y[j], X, "-"), n), nrow = n, byrow = TRUE) / bw;
    aux <- kfunc(ktype = ktype, difmat = AUX );
    aux1 <- t(wt * t(aux));
    diag(aux1) <- 0;
    resu[, j] <- (apply(aux1, 1, sum)) / den;
  }
  return(resu)
}

#' Numerical Integral function using Simpson's rule
#'
#' @param x The numeric data vector.
#' @param fx The function.
#' @param n.pts Number of points.
#' @param method The character string specifying method of numerical integration. The possible options are \code{trap} for trapezoidal rule and \code{simps} for simpson'r rule.
#' @keywords internal

integ <- function(x, fx, method, n.pts = 256) {
  n = length(x)
  if (method == "simps") {
    if (class(fx) == "function")
      fx = fx(x)
    if (n != length(fx))
      stop("Unequal input vector lengths")
    if (n.pts < 64)
      n.pts = 64
    ap = approx(x, fx, n = 2 * n.pts + 1)
    h = diff(ap$x)[1]
    integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] + ap$y[2 * (1:n.pts) + 1]) / 3
    value = sum(integral)
  }
  if (method == "trap") {
    if (!is.numeric(x) | !is.numeric(fx))
    {
      stop('The variable of integration "x" or "fx" is not numeric.')
    }
    if (length(x) != length(fx))
    {
      stop("The lengths of the variable of integration and the integrand do not match.")
    }
    # integrate using the trapezoidal rule
    integral <- 0.5 * sum((x[2:(n)] - x[1:(n - 1)]) * (fx[1:(n - 1)] + fx[2:n]))
    value <- integral
  }
  return(value)
}

#' Weighted quartile estimation
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param p The percentile value. The defult is 0.5.
#' @keywords internal

wquantile <- function(X, wt, p = 0.5)
{
  if (!is.numeric(wt) || length(X) != length(wt))
    stop("X and wt must be numeric and equal-length vectors")
  if (!is.numeric(p) || any(p < 0 | p > 1))
    stop("Quartiles must be 0<=p<=1")
  if (min(wt) < 0)
    stop("Weights must be non-negative numbers")
  ord <- order(X)
  X <- X[ord]
  cusumw <- cumsum(wt[ord])
  sumW <- sum(wt)
  plist <- cusumw / sumW
  qua <- withCallingHandlers(approx(plist, X, p)$y, warning=function(w){invokeRestart("muffleWarning")})

  return(qua)
}

#' Weighted variance estimation
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param na.rm The character indicator wether to consider missing value(s) or not. The defult is FALSE.
#' @keywords internal

wvar <- function(X, wt, na.rm = FALSE) {
  if (na.rm) {
    wt <- wt[i <- !is.na(X)]
    X <- X[i]
  }
  wsum <- sum(wt)
  wmean = sum(wt * X) / wsum
  varr = sum(wt * (X - wmean) ^ 2) / (wsum)
  return(varr)
}
