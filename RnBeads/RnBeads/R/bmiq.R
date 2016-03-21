########################################################################################################################
## bmiq.R
## created: 2013-04-01
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Beta-mixture quantile normalization method, as implemented by Andrew Teschendorff and improved by Steve Horvath.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

### Author: Andrew Teschendorff
### Date: 4th Oct 2012

#' BMIQ
#'
#' Performs Beta-mixture quantile normalization, adjusting for type II bias in Infinium 450K data.
#'
#' @param beta.v   \code{double} vector consisting of beta values. Missing values (\code{NA}s) cannot be handled, so
#'                 these must be removed or imputed prior to running BMIQ. Before normalization, beta values that are
#'                 exactly 0 and exactly 1 are replaced by the minimum positive and maximum value below 1, respectively.
#' @param design.v \code{integer} vector of length \code{length(beta.v)}, containing the values \code{1} and \code{2}
#'                 only. These values specify probe design type.
#' @param doH      Flag indicating if normalization for hemimethylated type II probes is to be performed.
#' @param nfit     Number of probes of a given design to use for the fitting. Smaller values will make BMIQ faster at
#'                 the expense of accuracy. Values between 10000 and 50000 seem to work well.
#' @param th1.v    Thresholds "type 1" to use for the initialization of the EM algorithm. These values should represent
#'                 best guesses for calling type I probes hemi-methylated and methylated, and are refined in further
#'                 steps by the algorithm.
#' @param th2.v    Thresholds "type 2" to used for the initialization of the EM algorithm. These values should represent
#'                 best guesses for calling type II probes hemi-methylated and methylated, and are refined in further
#'                 steps by the EM algorithm. If this is \code{NULL} (default), the thresholds are estimated based on
#'                 \code{th1.v} and a modified PBC correction method.
#' @param niter    Maximum number of EM iterations to be performed.
#' @param tol      Tolerance threshold for EM algorithm.
#' @return List with the following elements:
#'         \describe{
#'           \item{\code{"all"}}{The normalised beta-profile for the sample.}
#'           \item{\code{"class1"}}{Methylation state assigned to the type I probes.}
#'           \item{\code{"class2"}}{Methylation state assigned to the type II probes.}
#'           \item{\code{"av1"}}{Mean beta values for the \code{nL} classes for type I probes.}
#'           \item{\code{"av2"}}{Mean beta values for the \code{nL} classes for type II probes.}
#'           \item{\code{"hf"}}{Hubble dilation factor.}
#'           \item{\code{"th1"}}{Estimated thresholds used for type I probes.}
#'           \item{\code{"th2"}}{Estimated thresholds used for type II probes.}
#'         }
#' 
#' @author Andrew Teschendorff and Steve Horvath; with minor modifications by Yassen Assenov
#' @export
BMIQ <- function(beta.v,design.v,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001){

	rnb.require("RPMM")

	type1.idx <- which(design.v==1L)
	type2.idx <- which(design.v==2L)
	beta1.v <- beta.v[type1.idx]
	beta2.v <- beta.v[type2.idx]

	## Replace zeros by minimum non-zero
	zero.idx <- which(beta1.v==0)
	suppressWarnings(beta1.v[zero.idx] <- min(beta1.v[-zero.idx]))
	zero.idx <- which(beta2.v==0)
	suppressWarnings(beta2.v[zero.idx] <- min(beta2.v[-zero.idx]))
	rm(zero.idx)

	## Estimates initial weight matrix and trains the model
	fit.blc <- function(betam.v, th.v, maxiter = niter, toler = tol) {
		w.m <- matrix(0,nrow=length(betam.v),ncol=3L)
		w.m[which(betam.v <= th.v[1]),1] <- 1
		w.m[which(th.v[1] < betam.v & betam.v <= th.v[2]),2] <- 1
		w.m[which(betam.v > th.v[2]),3] <- 1
		blc2(matrix(betam.v,ncol=1),w=w.m,maxiter=maxiter,tol=toler,verbose=FALSE)
	}

	### Fit type I
	betam.v <- sample(beta1.v,nfit)
	em1.o <- fit.blc(betam.v, th1.v)
	subsetclass.v <- apply(em1.o$w,1,which.max)
	i2 <- which(subsetclass.v==2L)
	subsetth1.v <- c(
		mean(c(max(betam.v[subsetclass.v==1L]),min(betam.v[i2]))),
		mean(c(max(betam.v[i2]),min(betam.v[subsetclass.v==3L]))))
	class1.v <- rep(2L,length(beta1.v))
	class1.v[which(beta1.v < subsetth1.v[1])] <- 1L
	class1.v[which(beta1.v > subsetth1.v[2])] <- 3L

	## Estimate modes
	d1U.o <- density(beta1.v[class1.v==1L])
	d1M.o <- density(beta1.v[class1.v==3L])
	mod1U <- d1U.o$x[which.max(d1U.o$y)]
	mod1M <- d1M.o$x[which.max(d1M.o$y)]
	
	d2U.o <- density(beta2.v[which(beta2.v<0.4)])
	d2M.o <- density(beta2.v[which(beta2.v>0.6)])
	mod2U <- d2U.o$x[which.max(d2U.o$y)]
	mod2M <- d2M.o$x[which.max(d2M.o$y)]

	### Deal with type II fit
	if (is.null(th2.v)) {
		th2.v <- c(subsetth1.v[1] + (mod2U-mod1U), subsetth1.v[2] + (mod2M-mod1M))
	}
	betam.v <- sample(beta2.v, nfit)
	em2.o <- fit.blc(betam.v,th2.v)
	
	### For type II probes: assign to state (unmethylated, hemi- or fully methylated)
	subsetclass.v <- apply(em2.o$w,1,which.max)
	i2 <- which(subsetclass.v==2L)
	subsetth2.v <- c(
		mean(c(max(betam.v[subsetclass.v==1L]),min(betam.v[i2]))),
		mean(c(max(betam.v[i2]),min(betam.v[subsetclass.v==3L]))))
	class2.v <- rep(2L,length(beta2.v))
	class2.v[which(beta2.v < subsetth2.v[1])] <- 1L
	class2.v[which(beta2.v > subsetth2.v[2])] <- 3L

	classAV1.v <- as.vector(em1.o$mu)
	classAV2.v <- as.vector(em2.o$mu)

	### Start normalising type II probes
	nbeta2.v <- beta2.v
	## Select unmethylated probes
	lt <- 1L
	selU.idx <- which(class2.v == lt)
	selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])]
	selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])]
	## Probability according to type II distribution -> corresponding quantile in type I distribution
	p.v <- pbeta(beta2.v[selUR.idx], em2.o$a[lt,1], em2.o$b[lt,1], lower.tail = FALSE)
	nbeta2.v[selUR.idx] <- qbeta(p.v, em1.o$a[lt,1], em1.o$b[lt,1], lower.tail = FALSE)
	p.v <- pbeta(beta2.v[selUL.idx], em2.o$a[lt,1], em2.o$b[lt,1], lower.tail = TRUE)
	nbeta2.v[selUL.idx] <- qbeta(p.v, em1.o$a[lt,1], em1.o$b[lt,1], lower.tail = TRUE)

	## Select methylated probes
	lt <- 3L
	selM.idx <- which(class2.v == lt)
	selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])]
	selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])]
	## Probability according to type II distribution -> corresponding quantile in type I distribution
	p.v <- pbeta(beta2.v[selMR.idx], em2.o$a[lt,1], em2.o$b[lt,1], lower.tail = FALSE)
	nbeta2.v[selMR.idx] <- qbeta(p.v, em1.o$a[lt,1], em1.o$b[lt,1], lower.tail = FALSE)

	if (doH) {
		## Select hemimethylated probes and include ML probes
		lt <- 2L
		selH.idx <- c(which(class2.v == lt), selML.idx)
		minH <- min(beta2.v[selH.idx])
		maxH <- max(beta2.v[selH.idx])

		## Set new maximum and minimum of hemimethylated probes
		nmaxH <- min(nbeta2.v[selMR.idx]) - min(beta2.v[selMR.idx]) + maxH
		nminH <- max(nbeta2.v[selU.idx]) - max(beta2.v[selU.idx]) + minH
		
		## Postulate a conformal transformation (shift + dilation)
		## new_beta_H(i) = a + hf * (beta_H(i) - minH)
		hf <- (nmaxH - nminH) / (maxH - minH)
		nbeta2.v[selH.idx] <- nminH + hf * (beta2.v[selH.idx] - minH)
	}

	pnbeta.v <- beta.v
	pnbeta.v[type1.idx] <- beta1.v
	pnbeta.v[type2.idx] <- nbeta2.v

	return(list(all = pnbeta.v, class1 = class1.v, class2 = class2.v, av1 = classAV1.v, av2 = classAV2.v,
			hf = hf, th1 = subsetth1.v, th2 = th2.v))
}

########################################################################################################################

## @author Steve Horvath
betaEst2 <- function(y, w, weights) {
	yobs = which(!is.na(y))
	if (length(yobs) <= 1) {
		return(c(1, 1))
	}
	y <- y[yobs]
	w <- w[yobs]
	weights <- weights[yobs]
	N <- sum(weights * w)
	p <- sum(weights * w * y ) / N
	v <- sum(weights * w * y * y) / N - p * p
	logab = log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p) / v - 1))
	if (sum(yobs) == 2) {
		return(exp(logab))
	}
	opt <- try(optim(logab, RPMM::betaObjf, ydata = y, wdata = w, weights = weights, method = "Nelder-Mead",
			control = list(maxit = 50)), silent = TRUE)
	if (inherits(opt, "try-error")) {
		return(c(1, 1))
	}
	exp(opt$par)
}

########################################################################################################################

## @author Steve Horvath
blc2 <- function(Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = TRUE) {
	Ymn = min(Y[Y > 0], na.rm = TRUE)
	Ymx = max(Y[Y < 1], na.rm = TRUE)
	Y = pmax(Y, Ymn/2)
	Y = pmin(Y, 1 - (1 - Ymx)/2)
	Yobs = !is.na(Y)
	J = dim(Y)[2]
	K = dim(w)[2]
	n = dim(w)[1]
	if (n != dim(Y)[1]) 
		stop("Dimensions of w and Y do not agree")
	if (is.null(weights)) 
		weights = rep(1, n)
	mu = a = b = matrix(Inf, K, J)
	crit = Inf
	for (i in 1:maxiter) {
		warn0 = options()$warn
		options(warn = -1)
		eta = apply(weights * w, 2, sum)/sum(weights)
		mu0 = mu
		for (k in 1:K) {
			for (j in 1:J) {
				ab = betaEst2(Y[, j], w[, k], weights)
				a[k, j] = ab[1]
				b[k, j] = ab[2]
				mu[k, j] = ab[1]/sum(ab)
			}
		}
		ww = array(0, dim = c(n, J, K))
		for (k in 1:K) {
			for (j in 1:J) {
				ww[Yobs[, j], j, k] = dbeta(Y[Yobs[, j], j], 
					a[k, j], b[k, j], log = TRUE)
			}
		}
		options(warn = warn0)
		w = apply(ww, c(1, 3), sum, na.rm = TRUE)
		wmax = apply(w, 1, max)
		for (k in 1:K) w[, k] = w[, k] - wmax
		w = t(eta * t(exp(w)))
		like = apply(w, 1, sum)
		w = (1/like) * w
		llike = weights * (log(like) + wmax)
		crit = max(abs(mu - mu0))
		if (verbose) 
			print(crit)
		if (crit < tol) 
			break
	}
	return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}
