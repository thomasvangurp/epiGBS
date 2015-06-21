########################################################################################################################
## differentialMethylation.R
## created: 2012-12-06
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for determining differential methylation between groups.
########################################################################################################################

P.VAL.CUT <- 0.05
#parameters for density scatterplots on the site level
DENS.SCATTER.SPARSE.POINTS.PERC <- 0.01 #percentage of points to plot in the sparsely populated regions
DENS.SCATTER.SPARSE.POINTS.MAX <- 1e4 #maximum number of points to plot in the sparsely populated regions
DENS.SCATTER.SUBSAMPLE.THRES <- 2e6 #threshold to induce subsampling

### tTestP
###
### wrapper for t.test that catches errors
### @author Fabian Mueller
### @param ... parameters to be passed on to t.test
### @return vector of p-values resulting from t.test. Returns NA if error occurred
### @seealso \code{\link{t.test}}
tTestP <- function(...) {
	tryCatch(suppressWarnings(t.test(...)$p.value), error = function(e) { NA })
}

#' rowWelchP
#'
#' performs a two-sided Welch's t-test (unequal variances, unequal sample sizes) on each row of a matrix X with the indices inds.1 vs indices inds.g2 as group assignments.
#' @author Fabian Mueller
#' @param X Matrix on which the test is performed for every row
#' @param inds.g1 column indices of group 1 members
#' @param inds.g2 column indices of group 2 members
#' @param na.rm Should NAs be removed (logical)
#' @param alternative Testing alternative. Must be one of "two.sided" (default),"less","greater" or "all".
#' 		  in case of "all" a data frome with corresping alternative variables is returned. 
#' 		  Otherwise the result is a vector.
#' @return vector (or data.frame if alternative=="all") of p-values resulting from the Welch's t-test
#' @export
#' @note Requires \code{matrixStats} package
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' meth.mat <- meth(rnb.set.example)
#' sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
#' p.vals <- rowWelchP(meth.mat,sample.groups[[1]],sample.groups[[2]])
#' }
rowWelchP <- function(X,inds.g1,inds.g2=-inds.g1,na.rm=FALSE,alternative="two.sided"){
	if (!(alternative %in% c("two.sided","less","greater","all"))) {
		stop("invalid value for testing alternative")
	}
	X.1 <- X[,inds.g1]
	X.2 <- X[,inds.g2]
	if (na.rm){
		n.1 <- rowSums(!is.na(X.1), FALSE)
		n.2 <- rowSums(!is.na(X.2), FALSE)
	} else {
		n.1 <- length(inds.g1)
		n.2 <- length(inds.g2)
	}
	rm.1 <- rowMeans(X.1, na.rm = na.rm)
	rm.2 <- rowMeans(X.2, na.rm = na.rm)
	rv.1 <- rowVars(X.1, na.rm = na.rm)
	rv.2 <- rowVars(X.2, na.rm = na.rm)
	rq.1 <- rv.1/n.1
	rq.2 <- rv.2/n.2
	t.stat <- (rm.1 - rm.2)/sqrt(rq.1 + rq.2)
	rdf <- (rq.1 + rq.2)^2/(rq.1^2/(n.1-1) + rq.2^2/(n.2-1)) #degrees of freedom
	rp <- rep(NA,nrow(X))
	if (alternative == "two.sided" || alternative == "all") {
		rp.2s <- 2*pt(-abs(t.stat),rdf)
	}
	if (alternative == "less" || alternative == "all") {
		rp.l <- pt(t.stat,rdf)
	}
	if (alternative == "greater" || alternative == "all") {
		rp.g <- pt(t.stat,rdf,lower.tail=FALSE)
	}
	if (alternative == "two.sided") rp <- rp.2s
	if (alternative == "greater")   rp <- rp.g
	if (alternative == "less")      rp <- rp.l
	if (alternative == "all")		rp <- data.frame(less=rp.l,greater=rp.g,two.sided=rp.2s)
	return(rp)
}

#' rowPairedTP
#'
#' performs a two-sided t-test for paired samples on each row of a matrix X with the indices inds.1 vs indices inds.g2 as group assignments.
#' @author Fabian Mueller
#' @param X Matrix on which the test is performed for every row
#' @param inds.g1 column indices of group 1 members. \code{length(inds.g1)==length(inds.g2)} has to hold true.
#' @param inds.g2 column indices of group 2 members. \code{length(inds.g1)==length(inds.g2)} has to hold true.
#' @param alternative Testing alternative. Must be one of "two.sided" (default),"less","greater" or "all".
#' 		  in case of "all" a data frome with corresping alternative variables is returned. 
#' 		  Otherwise the result is a vector.
#' @return vector (or data.frame if alternative=="all") of p-values from a paired t-test
#' @export
#' @note Requires \code{matrixStats} package
rowPairedTP <- function(X,inds.g1,inds.g2=-inds.g1,alternative="two.sided"){
	if (!(alternative %in% c("two.sided","less","greater","all"))) {
		stop("invalid value for testing alternative")
	}
	X.1 <- X[,inds.g1]
	X.2 <- X[,inds.g2]
	if (ncol(X.1)!=ncol(X.2)){
		stop("unequal number of indices for the two groups")
	}
	X.d <-  X.1 - X.2
	n <- rowSums(!is.na(X.d), FALSE)
	d.bar <- rowMeans(X.d, na.rm = TRUE)
	s.d <- sqrt(rowVars(X.d, na.rm = TRUE))

	t.stat <- sqrt(n)*(d.bar)/s.d
	
	rp <- rep(NA,nrow(X))
	if (alternative == "two.sided" || alternative == "all") {
		rp.2s <- 2*pt(-abs(t.stat),n-1)
	}
	if (alternative == "less" || alternative == "all") {
		rp.l <- pt(t.stat,n-1)
	}
	if (alternative == "greater" || alternative == "all") {
		rp.g <- pt(t.stat,n-1,lower.tail=FALSE)
	}
	if (alternative == "two.sided") rp <- rp.2s
	if (alternative == "greater")   rp <- rp.g
	if (alternative == "less")      rp <- rp.l
	if (alternative == "all")		rp <- data.frame(less=rp.l,greater=rp.g,two.sided=rp.2s)
	return(rp)
}

#' rowOneSampleTP
#'
#' performs a two-sided t-test for paired samples on each row of a matrix X with the indices inds.1 vs indices inds.g2 as group assignments.
#' @author Fabian Mueller
#' @param X Matrix on which the test is performed for every row
#' @param mu The mean that is tested against
#' @param alternative Testing alternative. Must be one of "two.sided" (default),"less","greater" or "all".
#' 		  in case of "all" a data frome with corresping alternative variables is returned. 
#' 		  Otherwise the result is a vector.
#' @return vector (or data.frame if alternative=="all") of p-values from a paired t-test
#' @export
#' @note Requires \code{matrixStats} package
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' meth.mat <- meth(rnb.set.example)
#' p.vals <- rowOneSampleTP(meth.mat,mu=0,alternative="greater")
#' }
rowOneSampleTP <- function(X,mu=0,alternative="two.sided"){
	if (!(alternative %in% c("two.sided","less","greater","all"))) {
		stop("invalid value for testing alternative")
	}

	n <- rowSums(!is.na(X), FALSE)
	x.bar <- rowMeans(X, na.rm = TRUE)
	s.d <- sqrt(rowVars(X, na.rm = TRUE))
	
	t.stat <- sqrt(n)*(x.bar - mu)/s.d
	
	rp <- rep(NA,nrow(X))
	if (alternative == "two.sided" || alternative == "all") {
		rp.2s <- 2*pt(-abs(t.stat),n-1)
	}
	if (alternative == "less" || alternative == "all") {
		rp.l <- pt(t.stat,n-1)
	}
	if (alternative == "greater" || alternative == "all") {
		rp.g <- pt(t.stat,n-1,lower.tail=FALSE)
	}
	if (alternative == "two.sided") rp <- rp.2s
	if (alternative == "greater")   rp <- rp.g
	if (alternative == "less")      rp <- rp.l
	if (alternative == "all")		rp <- data.frame(less=rp.l,greater=rp.g,two.sided=rp.2s)
	return(rp)
}

#' limmaP
#'
#' applies hierarchical modeling anlalogous to differential expression employed in the \code{limma} package and returns
#' p-values for differential methylation
#' @author Fabian Mueller
#' @param X Matrix on which the test is performed for every row
#' @param inds.g1 column indices of group 1 members
#' @param inds.g2 column indices of group 2 members
#' @param adjustment.table a \code{data.frame} containing variables to adjust for in the testing
#' @param fun.conversion conversion function to transform the beta values into M values. By default, it is the logit function with adjustment
#' 						 for infinity values. See \code{\link{rnb.beta2mval}} for details.
#' @param paired should a paired analysis model be used. If so, the first index in \code{inds.g1} must correspond to the first
#' 				 index in  \code{inds.g2} and so on.
#' @return vector of p-values resulting from limma's differential analysis
#' @note Requires \code{limma} package
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' meth.mat <- meth(rnb.set.example)
#' sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
#' p.vals <- limmaP(meth.mat,sample.groups[[1]],sample.groups[[2]])
#' }
limmaP <- function(X,inds.g1,inds.g2=-inds.g1,adjustment.table=NULL,fun.conversion=rnb.beta2mval,paired=FALSE){
	# suppressPackageStartupMessages(require(limma))
	if (is.logical(inds.g1)) inds.g1 <- which(inds.g1)
	if (is.logical(inds.g2)) inds.g2 <- which(inds.g2)
	n.g1 <- length(inds.g1)
	n.g2 <- length(inds.g2)
	n <- n.g1 + n.g2
	if (!is.null(adjustment.table)){
		if (!(is.data.frame(adjustment.table) && nrow(adjustment.table)==n && (!any(is.na(adjustment.table))))) {
			stop("invalid value for adjustment.table")
		}
		m <- ncol(adjustment.table)
		if (m == 0) {
			adjustment.table <- NULL
		} else {
			colnames(adjustment.table) <- paste0("x",1:m,"x")
		}
	}

	ind.vec <- c(inds.g1,inds.g2)
	if (length(ind.vec) < 2) stop("need at least two samples indices to compare")
	X.m <- fun.conversion(X[,ind.vec,drop=FALSE])

	## Set a covariate defining the two groups
	df <- data.frame(xg = factor(rep(c("group1","group2"), c(n.g1,n.g2)), levels=c("group1","group2")))
	if (!is.null(adjustment.table)){
		## Add covariates to adjust for
		df <- cbind(df,adjustment.table)
	}
	if (paired){
		## Add a covariate for pairing
		if (n.g1 != n.g2) {
			stop("Could not conduct paired limma analysis: unequal groupsizes")
		}
		df$xp <- as.factor(rep(1:n.g1,2))
	}

	formula.text <- paste0(c("~0",colnames(df)),collapse="+")
	design.m <- model.matrix(as.formula(formula.text),data=df)
	colnames(design.m) <- make.names(colnames(design.m),unique=TRUE)
	colnames(design.m)[1:2] <- c("group1","group2")
	fit <- limma::lmFit(X.m,design.m)
	contrasts.m <- makeContrasts(group1vs2=group1-group2,levels=design.m)
	fit <- limma::contrasts.fit(fit,contrasts.m)
	fit <- limma::eBayes(fit)
	return(fit$p.value[,"group1vs2"])
}

#' computeDiffTab.site
#'
#' computes a difference table containing multiple difference measures,
#' In the simple version the difference in means,
#' quotients in means and a p-value for the comparison of two groups in a table are computed.
#' This is computed for each row of the input table. The extended version contains additional columns
#' @rdname computeDiffTab.site
#' @author Fabian Mueller
#' @aliases computeDiffTab.site
#' @aliases computeDiffTab.default.site
#' @aliases computeDiffTab.extended.site
#' @param X Matrix on which the difference measures are calculated for every row
#' @param inds.g1 column indices of group 1 members
#' @param inds.g2 column indices of group 2 members
#' @param diff.method Method to determine p-values for differential methylation. Currently supported are 
#' 				"ttest" for a two-sided Welch t-test, "refFreeEWAS" for adjusting for cell mixtures,
#' 				and "limma" for p-values resulting from linear modeling of the transformed beta values (M-values)
#'				and using techniques from expression microarray analysis employed in the \code{limma} package.
#' @param paired should a paired a analysis be performed. If \code{TRUE} then inds.g1 and inds.g2 should have exactly the same length and should be
#' 			     order, such that the first element of inds.g1 corresponds to the first element of inds.g2 and so on.
#' @param adjustment.table a table of variables to be adjusted for in the differential methylation test. Currently this is only supported for
#'        \code{diff.method=="limma"}
#' @param eps Epsilon for computing quotients (avoid division by 0 by adding this value to denominator and enumerator before calculating the quotient)
#' @param covg coverage information (should be NULL for disabled or of equal dimensions as X)
#' @param covg.thres a coverage threshold
#' @return a dataframe containing the following variables:
#' \item{mean.g1}{Mean of group 1}
#' \item{mean.g2}{Mean of group 2}
#' \item{mean.diff}{Difference in means}
#' \item{mean.quot.log2}{log2 of the quotient of means}
#' \item{diffmeth.p.val}{P-value (as determined by \code{diff.method})}
#' \item{max.g1/max.g2}{[extended version only] Group maxima}
#' \item{min.g1/min.g2}{[extended version only] Group minima}
#' \item{sd.g1/sd.g2}{[extended version only] Group standard deviations}
#' \item{min.diff}{[extended version only] Minimum of 0 and single linkage difference between the groups}
#' \item{diffmeth.p.adj.fdr}{[extended version only] FDR adjusted p-values}
#' \item{num.na.g1/num.na.g2}{[extended version only] number of NA methylation values for groups 1 and 2 respectively}
#' \item{mean.covg.g1/mean.covg.g2}{[extended version with coverage information only] mean coverage of groups 1 and 2 respectively}
#' \item{min.covg.g1/min.covg.g2}{[extended version with coverage information only] minimum coverage of groups 1 and 2 respectively}
#' \item{max.covg.g1/max.covg.g2}{[extended version with coverage information only] maximum coverage of groups 1 and 2 respectively}
#' \item{covg.thresh.nsamples.g1/2}{[extended version with coverage information only] number of samples in group 1 and 2 respectively exceeding the
#' 									coverage threshold for this site.}
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' meth.mat <- meth(rnb.set.example)
#' sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
#' dm <- computeDiffTab.extended.site(meth.mat,sample.groups[[1]],sample.groups[[2]])
#' summary(dm)
#' }
computeDiffTab.default.site <- function(X,inds.g1,inds.g2,
		diff.method=rnb.getOption("differential.site.test.method"),
		paired=FALSE,adjustment.table=NULL,eps=0.01){
	if (!(diff.method %in% c("ttest","limma","refFreeEWAS"))) {
		stop("Invalid method for differential site methylation test method")
	}
	# require(matrixStats)
	tab.g1 <- X[,inds.g1]
	tab.g2 <- X[,inds.g2]
	if(length(inds.g1)<2) {
		logger.info("Group 1 has less than 2 members")
		tab.g1 <- as.matrix(tab.g1)
	}
	if(length(inds.g2)<2) {
		logger.info("Group 2 has less than 2 members")
		tab.g2 <- as.matrix(tab.g2)
	}

	if (!is.null(adjustment.table)){
		if (!is.element(diff.method,c("limma","refFreeEWAS"))){
			logger.warning("Adjust for covariates is currently not supported for the selected p-value method --> skipping covariate adjustment")
		} else {
			logger.info(paste0("Adjusting for covariates: ",paste(colnames(adjustment.table),collapse=",")))
		}
	}

	mean.g1 <- rowMeans(tab.g1, na.rm = TRUE)
	mean.g2 <- rowMeans(tab.g2, na.rm = TRUE)

	if (paired){
		mean.diff <- rowMeans(tab.g1 - tab.g2, na.rm = TRUE)
		mean.quot <- rowMeans((tab.g1+eps)/(tab.g2+eps), na.rm = TRUE)
	} else {
		mean.diff <- mean.g1 - mean.g2
		mean.quot <- (mean.g1+eps)/(mean.g2+eps)
	}
	mean.quot.log2 <- log2(mean.quot)

	p.vals <- rep(as.double(NA),nrow(X))
	do.p.vals <- ncol(tab.g1) > 1 || ncol(tab.g2) > 1
	if (do.p.vals) {
		if (diff.method == "limma"){
			logger.info("Conducting differential analysis using limma")
			tryCatch(
				p.vals <- limmaP(X,inds.g1,inds.g2,adjustment.table=adjustment.table,paired=paired),
				error = function(ee) {
					logger.warning(c("Could not compute p-values using limma:",ee$message))
				}
			)
		} else if (diff.method == "refFreeEWAS"){
			p.vals <- refFreeEWASP(X, inds.g1, inds.g2, adjustment.table=adjustment.table, paired=paired)
		} else if (paired){
			logger.info("Conducting differential analysis using paired Student t-test")
			p.vals <- rowPairedTP(X,inds.g1,inds.g2)
		} else if (length(inds.g1)>1 && length(inds.g2)>1) {
			logger.info("Conducting differential analysis using two-sided Welch t-test")
			p.vals <- rowWelchP(X,inds.g1,inds.g2,na.rm=TRUE)
		} else if(length(inds.g1)>1) {
			logger.info("Conducting differential analysis using two-sided Welch t-test")
			p.vals <- rowOneSampleTP(tab.g1,mu=tab.g2[,1])
		} else if(length(inds.g2)>1) {
			logger.info("Conducting differential analysis using two-sided Welch t-test")
			p.vals <- rowOneSampleTP(tab.g2,mu=tab.g1[,1])
		}
	} else {
		logger.warning("Skipping p-value computation due to insufficient sample numbers")
	}
	tt <- data.frame(mean.g1=mean.g1,mean.g2=mean.g2,mean.diff=mean.diff,mean.quot.log2=mean.quot.log2,diffmeth.p.val=p.vals)
	return(tt)
}
#' @rdname computeDiffTab.site
#' @export
computeDiffTab.extended.site <- function(X,inds.g1,inds.g2,
		diff.method=rnb.getOption("differential.site.test.method"),
		paired=FALSE,adjustment.table=NULL,
		eps=0.01,covg=NULL,covg.thres=rnb.getOption("filtering.coverage.threshold")){
	# require(matrixStats)
	tt.basic <- computeDiffTab.default.site(
		X,inds.g1=inds.g1,inds.g2=inds.g2,
		diff.method=diff.method,paired=paired,
		adjustment.table=adjustment.table,eps=eps
	)
	tab.g1 <- X[,inds.g1]
	tab.g2 <- X[,inds.g2]
	if (length(inds.g1)>1){
		max.g1 <- rowMaxs(tab.g1,na.rm=TRUE)
		min.g1 <- rowMins(tab.g1,na.rm=TRUE)
		sd.g1  <- rowSds(tab.g1,na.rm=TRUE)
		num.na.g1 <- rowSums(is.na(tab.g1))
	} else {
		max.g1  <- tab.g1
		min.g1  <- tab.g1
		sd.g1   <- NA
		num.na.g1 <- rep(0,length(tab.g1))
		num.na.g1[is.na(tab.g1)] <- 1
	}
	if (length(inds.g2)>1){
		max.g2 <- rowMaxs(tab.g2,na.rm=TRUE)
		min.g2 <- rowMins(tab.g2,na.rm=TRUE)
		sd.g2  <- rowSds(tab.g2,na.rm=TRUE)
		num.na.g2 <- rowSums(is.na(tab.g2))
	} else {
		max.g2  <- tab.g2
		min.g2  <- tab.g2
		sd.g2   <- NA
		num.na.g2 <- rep(0,length(tab.g2))
		num.na.g2[is.na(tab.g2)] <- 1
	}

	min.diff <- ifelse(max.g1 < min.g2, max.g1 - min.g2,ifelse(max.g2 < min.g1, min.g1 - max.g2, 0))
	p.vals.t.na.adj <- tt.basic$diffmeth.p.val
	p.vals.is.na <- is.na(tt.basic$diffmeth.p.val)
	if (!all(p.vals.is.na)){
		if (any(p.vals.is.na)){
			logger.info(c(sum(p.vals.is.na),"p-values are NA. They are treated as 1 in FDR adjustment"))
			p.vals.t.na.adj[is.na(p.vals.t.na.adj)] <- 1
		}
		p.vals.adj.t <- p.adjust(p.vals.t.na.adj, method = "fdr")
	} else {
		p.vals.adj.t <- rep(NA,nrow(tt.basic))
	}
	

	tt.ext <- data.frame(max.g1=max.g1,min.g1=min.g1,sd.g1=sd.g1,max.g2=max.g2,min.g2=min.g2,sd.g2=sd.g2,min.diff=min.diff,diffmeth.p.adj.fdr=p.vals.adj.t,
						 num.na.g1=num.na.g1,num.na.g2=num.na.g2)
	#coverage information
	if (!is.null(covg) & all(dim(covg)==dim(X))){
		covg[is.na(covg)] <- 0 #set NA to 0 coverage
		
		tab.covg.g1 <- covg[,inds.g1]
		tab.covg.g2 <- covg[,inds.g2]
		if(length(inds.g1)<2) {
			logger.info("Group 1 has less than 2 members")
			tab.covg.g1 <- as.matrix(tab.covg.g1)
		}
		if(length(inds.g2)<2) {
			logger.info("Group 2 has less than 2 members")
			tab.covg.g2 <- as.matrix(tab.covg.g2)
		}
		mean.covg.g1 <- rowMeans(tab.covg.g1)
		mean.covg.g2 <- rowMeans(tab.covg.g2)
		min.covg.g1 <- rowMins(tab.covg.g1)
		min.covg.g2 <- rowMins(tab.covg.g2)
		max.covg.g1 <- rowMaxs(tab.covg.g1)
		max.covg.g2 <- rowMaxs(tab.covg.g2)
		covg.thresh.nsamples.g1 <- rowSums(tab.covg.g1>=covg.thres)
		covg.thresh.nsamples.g2 <- rowSums(tab.covg.g2>=covg.thres)

		tt.ext.covg	<- data.frame(mean.covg.g1=mean.covg.g1,mean.covg.g2=mean.covg.g2,
								  min.covg.g1=min.covg.g1,min.covg.g2=min.covg.g2,
								  max.covg.g1=max.covg.g1,max.covg.g2=max.covg.g2,
								  covg.thresh.nsamples.g1=covg.thresh.nsamples.g1,covg.thresh.nsamples.g2=covg.thresh.nsamples.g2)
		tt.ext <- cbind(tt.ext,tt.ext.covg)	
	}

	tt <- cbind(tt.basic,tt.ext)
	return(tt)
}

#' combineTestPvalsMeth
#'
#' combine p-values of multiple tests using (a generalization of) Fisher's method. The parameter setting here is taylored to DNA methylation, but can be adapted. 
#' Reference: Makambi, K. (2003). Weighted inverse chi-square method for correlated significance tests. Journal of Applied Statistics, 30(2), 225-234.
#' @author Fabian Mueller, Christoph Bock
#' @aliases combineTestPvalsMeth
#' @param pvalues p-values to combine
#' @param testWeights weights for the individual tests
#' @param correlated are the individual tests correlated
#' @param methExpectedTestCorrelation expected correlation. Empirically approximated to the default value of 0.8 for DNA-methylation
#' @return the combined p-value
#' @export
#' @examples
#' \dontrun{
#' p.vals <- 10^-c(0,1,5)
#' combineTestPvalsMeth(p.vals)
#' }
combineTestPvalsMeth <- function(pvalues,testWeights=NULL,correlated=FALSE,methExpectedTestCorrelation = 0.8) {
	if (is.null(pvalues)){
		return(NA)	
	}
	if (!is.numeric(pvalues)){
		logger.warning(c("Non numeric value for pvalues in combination:",pvalues))
		return(NA)
	} 
	
	if (!is.null(testWeights)) {
		# check if weights are valid
		if (length(pvalues) != length(testWeights)) stop("Number of items in <pvalues> and in <testWeights> must be identical if weights are to be used")
		if (sum(is.na(testWeights))>0) stop("NA values are not permitted for the test weights")
		if (sum(testWeights<0)>0) stop("Weights must be positive")
		# standardize weights
		testWeights = testWeights/sum(testWeights)
	} else {
		# use equal weighting
		testWeights = rep(1/length(pvalues),length(pvalues))
	}
	pvalues[is.na(pvalues)] = 1
	
	if (length(pvalues) < 1){
		return(NA)	
	} else if (length(pvalues) == 1) {
		return(pvalues[1])
	} else if (length(pvalues) > sqrt(.Machine$integer.max)){
		logger.info(c("Too many p-values to combine --> using subsampling"))
		nn <- trunc(sqrt(.Machine$integer.max))
		ss <- sample(length(pvalues),nn)
		pvalues <- pvalues[ss]
		testWeights <- testWeights[ss]
		testWeights <- testWeights/sum(testWeights)
	}
	
	##estimating the correlation
	#s <- -2*log(pvalues)
	#q.t <- var(s)
	#to.root <- 10.028 - 4*q.t/3
	#rho <- -2.167 + sqrt(to.root)
	#rho <- max(rho,0)
	
	if (correlated==FALSE & is.null(testWeights)) {
		# use Fisher's classical method for combining p-values in the case of independence
		tcombined = sum(-2*log(pvalues))
		return(pchisq(tcombined,2*length(pvalues),lower.tail=FALSE))
	} else {
		# use the method proposed in Makambi (2003) Journal of Applied Statistics
		r = ifelse(correlated,methExpectedTestCorrelation,0)
		m = length(pvalues)
		M.Fm = sum(-2*log(pvalues)*testWeights)
		var.M.Fm = 4*sum(testWeights^2)

		ij.pairs <- expand.grid(1:m,1:m) #create all combinations of 2 indices
		ij.pairs <- ij.pairs[ij.pairs[,1]!=ij.pairs[,2],] #remove those pairs where i==j
		tw.i <- testWeights[ij.pairs[,1]]
		tw.j <- testWeights[ij.pairs[,2]]
		vv <- tw.i * tw.j * (3.25*r+0.75*r^2)
		var.M.Fm <- var.M.Fm + sum(vv)
		
		nu = 8/var.M.Fm
		tcombined = nu*M.Fm/2 # chi square test generalization of Makambi	   	   
		return(pchisq(tcombined,nu,lower.tail=FALSE))
	}
}

#' computeDiffTab.region
#'
#' computes a difference table containing multiple difference measures,
#' In the simple version the mean of the difference in means,
#' the mean quotient in means and a combination of p-values on the site level are computed.
#' This is computed for each row of the input table. The extended version contains additional columns
#' @rdname computeDiffTab.region
#' @author Fabian Mueller
#' @aliases computeDiffTab.region
#' @aliases computeDiffTab.default.region
#' @param dmtp differential methylation table on the site level (as obtained from \code{\link{computeDiffTab.default.site}})
#' @param regions2sites a list containing for each region the indices of the corresponding sites in the site differential methylation table
#' @param includeCovg flag indicating whether to include coverage information
#' @return a dataframe containing the following variables for a given genomic region:
#' \item{mean.mean.g1,mean.mean.g2}{mean of mean methylation levels for group 1 and 2 across all sites in a region}
#' \item{mean.mean.diff}{Mean difference in means across all sites in a region}
#' \item{mean.mean.quot.log2}{Mean quotient in means across all sites in a region}
#' \item{comb.p.val}{Combined p-value using a generalization of Fisher's method. See \code{\link{combineTestPvalsMeth}} for details.}
#' \item{comb.p.adj.fdr}{FDR adjusted combined p-value}
#' \item{num.sites}{number of sites that were considered for a region}
#' \item{mean.num.na.g1/2}{mean number (accross all considered sites) of samples that contained an NA for group 1 and 2 respectively}
#' \item{mean.mean.covg.g1/2}{Mean value of mean coverage values (across all samples in a group) across all sites in a region}
#' \item{mean.nsamples.covg.thresh.g1/2}{mean number (accross all considered sites) of samples that have a coverage larger than the specified threshold
#' 		(see \code{\link{computeDiffTab.default.site}} for details) for group 1 and 2 respectively}
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' meth.mat <- meth(rnb.set.example)
#' sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
#' dm.sites <- computeDiffTab.extended.site(meth.mat,sample.groups[[1]],sample.groups[[2]])
#' map.regions.to.sites <- regionMapping(rnb.set.example,"promoters")
#' dm.promoters <- computeDiffTab.default.region(dm.sites,map.regions.to.sites)
#' }
computeDiffTab.default.region <- function(dmtp,regions2sites,includeCovg=FALSE){
	if (length(regions2sites)<1){
		stop("regions2sites argument must have length > 0")
	}
	col.id.g1 <- "mean.g1"
	col.id.g2 <- "mean.g2"
	col.id.diff <- "mean.diff"
	col.id.quot <- "mean.quot.log2"
	col.id.p  <- "diffmeth.p.val"
	col.id.num.na.g1 <- "num.na.g1"
	col.id.mean.covg.g1 <- "mean.covg.g1"
	col.id.covg.thresh.ns.g1 <- "covg.thresh.nsamples.g1"
	col.id.num.na.g2 <- "num.na.g2"
	col.id.mean.covg.g2 <- "mean.covg.g2"
	col.id.covg.thresh.ns.g2 <- "covg.thresh.nsamples.g2"
	n.regs.with.sites <- length(regions2sites)
	mean.g1 <- rep(NA,n.regs.with.sites)
	mean.g2 <- rep(NA,n.regs.with.sites)
	diff <- rep(NA,n.regs.with.sites)
	quot <- rep(NA,n.regs.with.sites)
	p.vals <- rep(NA,n.regs.with.sites)
	num.sites <- rep(NA,n.regs.with.sites)
	mean.num.na.g1 <- rep(NA,n.regs.with.sites)
	mean.num.na.g2 <- rep(NA,n.regs.with.sites)
	mean.mean.covg.g1 <- rep(NA,n.regs.with.sites)
	mean.mean.covg.g2 <- rep(NA,n.regs.with.sites)
	mean.nsamples.covg.thresh.g1 <- rep(NA,n.regs.with.sites)
	mean.nsamples.covg.thresh.g2 <- rep(NA,n.regs.with.sites)
	col.vec <- c(col.id.g1, col.id.g2, col.id.diff, col.id.quot, col.id.p, col.id.num.na.g1, col.id.num.na.g2)
	if (includeCovg) {
		col.vec <- c(col.vec,c(col.id.mean.covg.g1, col.id.mean.covg.g2, col.id.covg.thresh.ns.g1, col.id.covg.thresh.ns.g2))
	}
	dmt4fastProc <- dmtp[,col.vec] #not looking up the columns by name, but by index reduces runtime
	
	if(parallel.isEnabled()){
		dm <- foreach(i=1:n.regs.with.sites, .combine='rbind',.multicombine=TRUE,.maxcombine=200) %dopar% {
			pids <- regions2sites[[i]]
			subtab <- dmt4fastProc[pids,]#these lookups take up most of the time
			
			mean.g1   <- mean(subtab[,1],na.rm=TRUE)
			mean.g2   <- mean(subtab[,2],na.rm=TRUE)
			diff      <- mean(subtab[,3],na.rm=TRUE)
			quot      <- mean(subtab[,4],na.rm=TRUE)
			num.sites <- length(pids)
			mean.num.na.g1 <- mean(subtab[,6])
			mean.num.na.g2 <- mean(subtab[,7])
			if (includeCovg){
				mean.mean.covg.g1 <- mean(subtab[,8])
				mean.mean.covg.g2 <- mean(subtab[,9])
				mean.nsamples.covg.thresh.g1 <- mean(subtab[,10])
				mean.nsamples.covg.thresh.g2 <- mean(subtab[,11])
			} else {
				mean.mean.covg.g1 <- NA
				mean.mean.covg.g2 <- NA
				mean.nsamples.covg.thresh.g1 <- NA
				mean.nsamples.covg.thresh.g2 <- NA
			}
			
			res <- combineTestPvalsMeth(na.omit(subtab[,5]),correlated=TRUE)
			p.vals <- NA
			if (length(res)>0) p.vals <- res
			c(mean.g1, mean.g2, diff, quot, num.sites, mean.num.na.g1, mean.num.na.g2,
			  mean.mean.covg.g1, mean.mean.covg.g2, mean.nsamples.covg.thresh.g1, mean.nsamples.covg.thresh.g2,
			  p.vals)
		}
		mean.g1                      <- dm[, 1]
		mean.g2                      <- dm[, 2]
		diff                         <- dm[, 3]
		quot                         <- dm[, 4]
		num.sites                    <- dm[, 5]
		mean.num.na.g1               <- dm[, 6]
		mean.num.na.g2               <- dm[, 7]
		mean.mean.covg.g1            <- dm[, 8]
		mean.mean.covg.g2            <- dm[, 9]
		mean.nsamples.covg.thresh.g1 <- dm[,10]
		mean.nsamples.covg.thresh.g2 <- dm[,11]
		p.vals                       <- dm[,12]
		
	} else {
		dummy <- sapply(1:n.regs.with.sites,FUN=function(i){
			pids <- regions2sites[[i]]
			subtab <- dmt4fastProc[pids,]#these lookups take up most of the time
			
			mean.g1[i] <<- mean(subtab[,1],na.rm=TRUE)
			mean.g2[i] <<- mean(subtab[,2],na.rm=TRUE)
			diff[i]    <<- mean(subtab[,3],na.rm=TRUE)
			quot[i]    <<- mean(subtab[,4],na.rm=TRUE)
			
			num.sites[i] <<- length(pids)
			mean.num.na.g1[i] <<- mean(subtab[,6])
			mean.num.na.g2[i] <<- mean(subtab[,7])
			if (includeCovg){
				mean.mean.covg.g1[i] <<- mean(subtab[,8])
				mean.mean.covg.g2[i] <<- mean(subtab[,9])
				mean.nsamples.covg.thresh.g1[i] <<- mean(subtab[,10])
				mean.nsamples.covg.thresh.g2[i] <<- mean(subtab[,11])
			}
					
			res <- combineTestPvalsMeth(na.omit(subtab[,5]),correlated=TRUE)
			if (length(res)>0) p.vals[i]  <<- res
			return(TRUE)
		})
	}
	
	p.vals.na.adj <- p.vals
	p.vals.is.na <- is.na(p.vals)
	if (any(p.vals.is.na)){
		logger.info(c(sum(p.vals.is.na),"p-values are NA. They are treated as 1 in FDR adjustment"))
		p.vals.na.adj[is.na(p.vals.na.adj)] <- 1
	}
	p.vals.adj <- p.adjust(p.vals.na.adj, method = "fdr")

	tt <- data.frame(mean.mean.g1=mean.g1,mean.mean.g2=mean.g2,mean.mean.diff=diff,mean.mean.quot.log2=quot,comb.p.val=p.vals,comb.p.adj.fdr=p.vals.adj,
					 num.sites=num.sites,
					 mean.num.na.g1=mean.num.na.g1,mean.num.na.g2=mean.num.na.g2)
	if (includeCovg){
		tt <- cbind(tt,data.frame(mean.mean.covg.g1=mean.mean.covg.g1,mean.mean.covg.g2=mean.mean.covg.g2,
								  mean.nsamples.covg.thresh.g1=mean.nsamples.covg.thresh.g1,mean.nsamples.covg.thresh.g2=mean.nsamples.covg.thresh.g2))
	}
	rownames(tt) <- names(regions2sites)
	return(tt)
}

### combinedRanking.tab
###
### computes the combined ranking for each row as the maximum rank among all columns
### @author Fabian Mueller
### @aliases combinedRanking.tab
### @param tt differential methylation table
### @param rerank if \code{TRUE} then an additional ranking will be performed on the combined rank in order to obtain values in [1,nrow(tt)]
### @return a vector containing the combined ranking
combinedRanking.tab <- function(tt,rerank=FALSE){
	rank.mat <- c()
	for (i in 1:ncol(tt)){
		rrs <- rank(tt[,i],na.last="keep",ties.method="min")
		if (!all(is.na(rrs))) {
			rank.mat <- cbind(rank.mat,rrs)
		}
	}
	if (is.null(rank.mat)){
		logger.warning("Could not compute combined ranking: To few non-NA columns specified")
		return(rep(NA,nrow(tt)))
	}
	res <- rowMaxs(rank.mat,na.rm=FALSE)
	res[res==-Inf] <- NA
	if (rerank) res <- rank(res,na.last="keep",ties.method="min")
	return(res)
}

### extractRankingCols
###
### extracts and transforms the relevant columns from a differential methylation table for ranking on the site and region level
### respectively
### @author Fabian Mueller
### @rdname extractRankingCols
### @aliases extractRankingCols
### @aliases extractRankingCols.site
### @aliases extractRankingCols.region
### @param tt differential methylation table
### @return a matrix containing the transformed and extracted values:
### \item{}{difference in mean methylation (negative absolute value) (mean of differences in means on the region level)}
### \item{}{(mean) quotient in mean methylation (negative absolute value of the logarithm) (mean of quotients in means on the region level)}
### \item{}{p-value from t-test. (combination of p-values using an extension of Fisher's method on the region level)}
extractRankingCols.site <- function(tt){
	return(cbind(-abs(tt$mean.diff),  -abs(tt$mean.quot.log2),  tt$diffmeth.p.val))
}
### @rdname extractRankingCols
extractRankingCols.region <- function(tt){
	return(cbind(-abs(tt$mean.mean.diff),  -abs(tt$mean.mean.quot.log2),  tt$comb.p.val))
}

### doPerm
###
### perform ONE permutation test
### @author Fabian Mueller
### @aliases doPerm
### @param b beta value matrix
### @param all.inds column indices in \code{b} to be used for permutation test
### @param n.inds.g1 number of columns in group 1
### @return vector of ranks resulting from the permutation test
doPerm <- function(b,all.inds,n.inds.g1,...){
	perm.inds.inds.g1 <- sample(1:length(all.inds),n.inds.g1)
	perm.inds.g1 <- all.inds[perm.inds.inds.g1]
	perm.inds.inds.g2 <- setdiff(1:length(all.inds),perm.inds.inds.g1)
	perm.inds.g2 <- all.inds[perm.inds.inds.g2]
	dm <- computeDiffTab.default.site(b,inds.g1=perm.inds.g1,inds.g2=perm.inds.g2,...)
	dm4ranking <- extractRankingCols.site(dm)
	perm.ranks.cur <- combinedRanking.tab(dm4ranking,rerank=FALSE)
	return(perm.ranks.cur)
}


### groupPermutationP.site
###
### computes a p-value for permuting the two sample groups by calculating combined ranks for each permutaion and scoring how many of them yield a better rank for each site.
### @author Fabian Mueller
### @aliases groupPermutationP.site
### @param b beta value matrix
### @param inds.1 column indices in \code{b} of group 1 members
### @param inds.g2 column indices in \code{b} of group 2 members
### @return a vector of p-values
groupPermutationP.site <- function(b,inds.g1,inds.g2,n.perm=500,...){
	report.interval <- 10
	all.inds <- c(inds.g1,inds.g2)
	n.inds.g1 <- length(inds.g1)
	dm <- computeDiffTab.default.site(b,inds.g1=inds.g1,inds.g2=inds.g2,...)
	dm4ranking <- extractRankingCols.site(dm)
	ranking.org <- combinedRanking.tab(dm4ranking,rerank=FALSE)
	perm.ranks <- matrix(NA,ncol=n.perm,nrow=nrow(b))
	for (i in 1:n.perm){
		if (i %% report.interval == 0) logger.status(c("Reached permutation",i))
		perm.ranks[,i] <- doPerm(b,all.inds,n.inds.g1)
	}
	p.perm <- rowSums(perm.ranks <  ranking.org)/n.perm #the m < v operator (m is a matrix, v is a vector) in R works in an equivalent way to the following procedure: m is regareded as the concatetened column vector. v is repeated. elementwise comparison. m is reassemled into a matrix.
	return(p.perm)
}

### @rdname groupPermutationP.site
groupPermutationP.site.parallel <- function(b,inds.g1,inds.g2,n.perm=500,...){
	all.inds <- c(inds.g1,inds.g2)
	n.inds.g1 <- length(inds.g1)
	dm <- computeDiffTab.default.site(b,inds.g1=inds.g1,inds.g2=inds.g2,...)
	dm4ranking <- extractRankingCols.site(dm)
	ranking.org <- combinedRanking.tab(dm4ranking,rerank=FALSE)
	perm.ranks <- foreach(i=1:n.perm,.combine='cbind') %dopar% doPerm(b,all.inds,n.inds.g1,...)
	p.perm <- rowSums(perm.ranks <  ranking.org)/n.perm #the m < v operator (m is a matrix, v is a vector) in R works in an equivalent way to the following procedure: m is regareded as the concatetened column vector. v is repeated. elementwise comparison. m is reassemled into a matrix.
	return(p.perm)
}

### computeDiffMeth.bin.site
###
### computes a differential methylation in the binary case (2 groups) on the site level.
### @author Fabian Mueller
### @aliases computeDiffMeth.bin.site
### @param b beta value matrix
### @param inds.1 column indices in \code{b} of group 1 members
### @param inds.g2 column indices in \code{b} of group 2 members
### @param n.perm number of permutations to be performed for the ranking permutaion tests. Set to values < 1 to disable permutation tests
### @return A data.frame containing differential methylation information with the variables from \code{\link{computeDiffTab.extended.site}} and additionally 
### \item{combinedRank}{the combined rank obtained from the the differential methylation information. As the the worst rank among all columns selected for the ranking.
### 					   \code{\link{extractRankingCols.site}} determines which these are.}
### \item{rankPermP}{[optional] p-value obtained from permuation tests of sample group assignments}
computeDiffMeth.bin.site <- function(b,inds.g1,inds.g2,n.perm=0,...){
	#sanity checks
	if (length(union(inds.g1,inds.g2)) != (length(inds.g1)+length(inds.g2))){
		logger.error("Overlapping sample sets in differential methylation analysis")
	}
	logger.start("Computing Differential Methylation Table (Site Level)")
	diffmeth.tab <- computeDiffTab.extended.site(b,inds.g1=inds.g1,inds.g2=inds.g2,...)
	diffmethTab4ranks <- extractRankingCols.site(diffmeth.tab)
	combRank <- combinedRanking.tab(diffmethTab4ranks,rerank=FALSE)
	diffmeth.tab$combinedRank <- combRank
	logger.completed()
	if (n.perm > 0){
		logger.start("Conducting Permutation Tests")
		if (n.perm > 2000) {
			logger.warning("The number of permutation tests conducted exceeds 2000. 
							Depending on the system's resources this could lead to errors")
		}
		if(parallel.isEnabled()) {
			logger.info("Using multicore")
			diffmeth.tab$rankPermP <- groupPermutationP.site.parallel(b,inds.g1,inds.g2,n.perm=n.perm)
		} else {
			logger.info("Using single core")
			diffmeth.tab$rankPermP <- groupPermutationP.site(b,inds.g1,inds.g2,n.perm=n.perm)
		}
		logger.completed()
	}
	return(diffmeth.tab)
}

### computeDiffMeth.bin.region
###
### computes a differential methylation in the binary case (2 groups) on the region level.
### @author Fabian Mueller
### @aliases computeDiffMeth.bin.region
### @param dmtp differential methylation table on the site level (as obtained from \code{\link{computeDiffMeth.bin.site}})
### @param inds.1 column indices in \code{b} of group 1 members
### @param inds.g2 column indices in \code{b} of group 2 members
### @param n.perm number of permutations to be performed for the ranking permutaion tests. Set to values < 1 to disable permutation tests
### @return blubb
computeDiffMeth.bin.region <- function(rnbSet,dmtp,inds.g1,inds.g2,region.types=rnb.region.types(assembly(rnbSet))){
	#sanity checks
	if (length(union(inds.g1,inds.g2)) != (length(inds.g1)+length(inds.g2))){
		logger.error("Overlapping sample sets in differential methylation analysis")
	}
	logger.start('Computing Differential Methylation Tables (Region Level)')
	diffmeth.tabs <- list()
	for (rt in region.types){
		regions2sites <- regionMapping(rnbSet,rt)
		regions2sites.is.all.na <- sapply(regions2sites,FUN=function(x){all(is.na(x))})
		if (any(regions2sites.is.all.na)) {
			stop(paste("Region mapping of RnBSet from sites to regions is inconsistent (",rt,")"))
		}
#		regions2sites <- regions2sites[!regions2sites.is.all.na]
#		regions2sites <- lapply(regions2sites,FUN=function(x){na.omit(x)})
#		attr(regions2sites, "omitted.regions") <- which(regions2sites.is.all.na)
		dmtr <- computeDiffTab.default.region(dmtp,regions2sites,includeCovg=!is.null(covg(rnbSet)))
		dmtr4ranks <- extractRankingCols.region(dmtr)
		combRank <- combinedRanking.tab(dmtr4ranks,rerank=FALSE)
		dmtr$combinedRank <- combRank
		diffmeth.tabs <- c(diffmeth.tabs,list(dmtr))
		logger.status(c("Computed table for", rt))
	}
	names(diffmeth.tabs) <- region.types
	logger.completed()
	return(diffmeth.tabs)
}

#' exportDMRs2regionFile
#'
#' export differentially methylated regions to region file (standard bed). The output is in BED6 format where the score corresponds to 
#' to the combined rank (rank==1 would receive a score of 1000 and a combined rank equal to the number of regions a score of 0)
#' @author Fabian Mueller
#' @aliases exportDMRs2regionFile
#' @param rnbSet the RnBSet object for which the DMRs were computed.
#' @param diffmeth DiffMeth object. See \code{\link{rnb.execute.computeDiffMeth}} for details.
#' @param dest destination file name
#' @param comp.name name of the comparison
#' @param region.type region type.
#' @param rank.cut rank cutoff. If \code{NULL} (default), all regions are processed.
#' @param rerank flag indicating whether the ranks should be reranked or whether \code{rank.cut} refers to the absolute rank
#' @return \code{NULL}
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' exportDMRs2regionFile(rnb.set.example,dm,tempfile(),get.comparisons(dm)[1],"promoters")
#' }
exportDMRs2regionFile <- function(rnbSet,diffmeth,dest,comp.name,region.type,rank.cut=NULL,rerank=FALSE){
	annot <- annotation(rnbSet, type=region.type, add.names=FALSE)
#	dmt <- diffmeth$region[[comp]][[region.type]]
	dmt <- get.table(diffmeth,comp.name,region.type)
	is.dmr <- rep(TRUE,dim(dmt)[1])
	if (!is.null(rank.cut)){
		rrs <- dmt[,"combinedRank"]
		if (rerank){
			rrs <- rank(rrs,na.last="keep",ties.method="min")
		}
		is.dmr <- rrs<=rank.cut
		dmt <- dmt[is.dmr,]
		annot <- annot[is.dmr,]
	}
	n.regs <- nrow(dmt)
	score <- round((n.regs-dmt[,"combinedRank"]-1)/n.regs * 1000,0)
	outtab <- data.frame(chr=annot$Chromosome,chromStart=annot$Start,chromEnd=annot$End,name=annot$ID,score=score,strand=annot$Strand)
	write.table(outtab,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=dest)
	return(NULL)
}

#' auto.select.rank.cut
#'
#' automatically select a rank cutoff for given ranks and p-values
#' current implementation: sort the p-values according to rank. select as rank cutoff
#' the rank for which the worst (i.e. max) p-value in the top list is still smaller than
#' the best (i.e. min) p-value of the group of worst-ranking p-values of equal size as the top-list
#' @author Fabian Mueller
#' @param p vector of p-values
#' @param r vector of ranks
#' @param alpha the percentile to select the top and bottom part of the list
#' @return the maximum rank fulfilling the criterion
auto.select.rank.cut <- function(p,r,alpha=0.1){
	res <- 0
	lp <- length(p)
	j <- 1L:lp
	oa <- order(r)
	od <- rev(oa)
	p.oa <- p[oa]
	cmin.d <- cummin(p[od])
	cquant <- cummax(p.oa)
	inds.better.than.tail <- which(cquant<cmin.d)
	if (length(inds.better.than.tail) > 0){
		L <- max(inds.better.than.tail)
		res <- r[oa][L]
	}
	return(res)
}

### create.diffMeth.bin.dens.dmp.scatter
###
### Helper function for addReportPlots.diffMeth.bin.site.scatter().
### creates a plot that based on the categorization of Differentially Methylated Probes (DMPs) 
### @author Fabian Mueller
### @param df2p differential methylation table. Must contain the columns "mean.g1","mean.g2","isDMP"
### @param grp1.name name as it appears on the x axis
### @param grp2.name name as it appears on the y axis
### @return ggplot2 plot
create.diffMeth.bin.dens.dmp.scatter <- function(df2p,grp1.name,grp2.name){
	n.points <- nrow(df2p)
	#plot order: plot DMPs last
	df2p$plotOrder <- NA
	is.dmp <- df2p$isDMP
	is.dmp[is.na(is.dmp)] <- FALSE
	num.not.dmp <- sum(!is.dmp)
	df2p$plotOrder[!is.dmp] <- seq_len(num.not.dmp)
	df2p$plotOrder[is.dmp] <- seq((num.not.dmp+1),n.points)
	
	df2p <- na.omit(df2p[,c("mean.g1","mean.g2","isDMP","plotOrder")])
	
	df2p$color <- NA
	if (sum(!df2p$isDMP)>1){
		colors.nodmp <- DENS.COLORS.LOW[1]
		tryCatch(
			colors.nodmp <- densCols(x=df2p[!df2p$isDMP,"mean.g1"],y=df2p[!df2p$isDMP,"mean.g2"],colramp = colorRampPalette(c(DENS.COLORS.LOW[1],DENS.COLORS.HIGH[1]))),
			error=function(ee){
				logger.warning(c("Could not assess density colors using densCols:",ee$message))
			}
		)
		df2p[!df2p$isDMP,"color"] <- colors.nodmp
	} else if (sum(!df2p$isDMP)==1){
		df2p[!df2p$isDMP,"color"] <- DENS.COLORS.LOW[1]
	}
	if (sum(df2p$isDMP)>1){
		colors.dmp <- DENS.COLORS.LOW[2]
		tryCatch(
			colors.dmp   <- densCols(x=df2p[ df2p$isDMP,"mean.g1"],y=df2p[ df2p$isDMP,"mean.g2"],colramp = colorRampPalette(c(DENS.COLORS.LOW[2],DENS.COLORS.HIGH[2]))),
			error=function(ee){
				logger.warning(c("Could not assess density colors using densCols:",ee$message))
			}
		)
		df2p[df2p$isDMP,"color"] <- colors.dmp
	} else if (sum(df2p$isDMP)==1){
		df2p[df2p$isDMP,"color"] <- DENS.COLORS.LOW[2]
	}
	
	pp <- ggplot(df2p) + aes(mean.g1,mean.g2) +
			labs(x=paste("mean.beta",grp1.name,sep="."),y=paste("mean.beta",grp2.name,sep=".")) + coord_fixed() +
			geom_point(aes(color=color,order=plotOrder)) + scale_color_identity()
	
	return(pp)
}

### addReportPlots.diffMeth.bin.site.scatter
###
### adds report scatterplots for differential methylation for the site level binary case to a report.
### @author Fabian Mueller
### @aliases addReportPlots.diffMeth.bin.site.scatter
### @param report the report to be modified
### @param dmt differential methylation table as created by \code{computeDiffMeth.bin.site}
### @param cmpName Comparison name as it will appear in the filename and figure selection box
### @param diffSiteRankCut vector of combined ranking cutoffs for classifying a site as differentially methylated
### @param grp1.name name of group 1 in the compoarison (for labelling in the plots)
### @param grp2.name name of group 2 in the compoarison (for labelling in the plots)
### @return list of report plot objects added
addReportPlots.diffMeth.bin.site.scatter <- function(report,dmt,cmpName,diffSiteRankCut,autoRankCut=NULL,grp1.name="Group1",grp2.name="Group2",
		rerank=TRUE,thres.p.val=0.05){
	df2p <- dmt #data frame to plot
	figPlots <- list()
	
	sparse.points <- DENS.SCATTER.SPARSE.POINTS.PERC
	if (DENS.SCATTER.SPARSE.POINTS.MAX < sparse.points*nrow(df2p)){
		sparse.points <- DENS.SCATTER.SPARSE.POINTS.MAX
	}
	dens.subsample <- FALSE
	if (nrow(df2p) > dens.subsample){
		dens.subsample <- DENS.SCATTER.SUBSAMPLE.THRES
	}

	#scatterplot based on adjusted p-value significance
	if (is.element("diffmeth.p.adj.fdr",colnames(dmt))){
		df2p$isDMP <- df2p[,"diffmeth.p.adj.fdr"] < P.VAL.CUT
		pp <- create.densityScatter(df2p[,c("mean.g1","mean.g2")],is.special=df2p$isDMP,
					dens.subsample=dens.subsample,sparse.points=sparse.points,add.text.cor=TRUE) +
				labs(x=paste("mean.beta",grp1.name,sep="."),y=paste("mean.beta",grp2.name,sep=".")) + coord_fixed()
		cur.cut.name <- "fdrAdjPval"
		figName <- paste("diffMeth_site",cmpName,cur.cut.name,sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots <- c(figPlots,list(report.plot))
	}
	
	#scatterplot based on rank cutoff significance
	rrs <- dmt[,"combinedRank"]
	if (rerank)	rrs <- rank(rrs,na.last="keep",ties.method="min")
	for (i in 1:length(diffSiteRankCut)){
		rc <- diffSiteRankCut[i]
		cur.cut.name <- paste("rc",i,sep="")
		df2p$isDMP <- rrs < rc
		
		pp <- create.densityScatter(df2p[,c("mean.g1","mean.g2")],is.special=df2p$isDMP,
					dens.subsample=dens.subsample,sparse.points=sparse.points,add.text.cor=TRUE) +
				labs(x=paste("mean.beta",grp1.name,sep="."),y=paste("mean.beta",grp2.name,sep=".")) + coord_fixed()
		figName <- paste("diffMeth_site",cmpName,cur.cut.name,sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots <- c(figPlots,list(report.plot))
	}
	if (is.integer(autoRankCut)){
		df2p$isDMP <- dmt[,"combinedRank"] <= autoRankCut
		pp <- create.densityScatter(df2p[,c("mean.g1","mean.g2")],is.special=df2p$isDMP,
					dens.subsample=dens.subsample,sparse.points=sparse.points,add.text.cor=TRUE) +
				labs(x=paste("mean.beta",grp1.name,sep="."),y=paste("mean.beta",grp2.name,sep=".")) + coord_fixed()
		figName <- paste("diffMeth_site",cmpName,"rcAuto",sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots <- c(figPlots,list(report.plot))
	}
	pp <- create.hex.summary.plot(df2p,q="combinedRank") + coord_fixed() +
			labs(x=paste("mean.beta",grp1.name,sep="."),y=paste("mean.beta",grp2.name,sep="."),fill = "median combined rank")
	figName <- paste("diffMeth_site",cmpName,"rankGradient",sep="_")
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))

	if (is.element("rankPermP",colnames(dmt))){
		df2p$isDMP <- dmt[,"rankPermP"] < thres.p.val
		pp <- create.densityScatter(df2p[,c("mean.g1","mean.g2")],is.special=df2p$isDMP,
					dens.subsample=dens.subsample,sparse.points=sparse.points,add.text.cor=TRUE) +
				labs(x=paste("mean.beta",grp1.name,sep="."),y=paste("mean.beta",grp2.name,sep=".")) + coord_fixed()

		figName <- paste("diffMeth_site",cmpName,"permutationP",sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots <- c(figPlots,list(report.plot))
	}

	return(figPlots)
}

### addReportPlots.diffMeth.bin.site.volcano
###
### adds report volcano plots for differential methylation for the site level binary case to a report.
### @author Fabian Mueller
### @aliases addReportPlots.diffMeth.bin.site.volcano
### @param report the report to be modified
### @param dmt differential methylation table as created by \code{computeDiffMeth.bin.site}
### @param cmpName Comparison name as it will appear in the filename and figure selection box
### @param grp1.name name of group 1 in the compoarison (for labelling in the plots)
### @param grp2.name name of group 2 in the compoarison (for labelling in the plots)
### @return list of report plot objects added
addReportPlots.diffMeth.bin.site.volcano <- function(report,dmt,cmpName,grp1.name="Group1",grp2.name="Group2"){
	df2p <- dmt #data frame to plot
	figPlots <- list()
	dont.plot.p.val <- all(is.na(df2p[,"diffmeth.p.val"]))
	
	figName <- paste("diffMeth_site_volcano",cmpName,"diff","pVal",sep="_")
	if (!dont.plot.p.val){
		pp <- ggplot(df2p) + aes(mean.diff,-log10(diffmeth.p.val),color=log10(combinedRank)) + scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
		  				 geom_point()#(alpha=0.3)
	 } else {
		pp <- rnb.message.plot("No p-value available")
	 }
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_site_volcano",cmpName,"quot","pVal",sep="_")
	if (!dont.plot.p.val){
		pp <- ggplot(df2p) + aes(mean.quot.log2,-log10(diffmeth.p.val),color=log10(combinedRank)) + scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
				geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))
	} else {
		pp <- rnb.message.plot("No p-value available")
	}
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_site_volcano",cmpName,"diff","pValAdj",sep="_")
	pp <- ggplot(df2p) + aes(mean.diff,-log10(diffmeth.p.adj.fdr),color=log10(combinedRank)) + 
		  scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
		  geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_site_volcano",cmpName,"quot","pValAdj",sep="_")
	pp <- ggplot(df2p) + aes(mean.quot.log2,-log10(diffmeth.p.adj.fdr),color=log10(combinedRank)) + 
		  scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
		  geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	#mean vs quotient plot
	figName <- paste("diffMeth_site_volcano",cmpName,"diff","quotSig",sep="_")
	pp <- ggplot(df2p) + aes(mean.diff,mean.quot.log2,color=log10(combinedRank)) + scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
			geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_site_volcano",cmpName,"quot","quotSig",sep="_")
	pp <- rnb.message.plot("Quotient--Quotient scatterplot not available")
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	return(figPlots)
}

### create.diffMeth.bin.dens.dmr.scatter
###
### Helper function for addReportPlots.diffMeth.bin.region.scatter().
### creates a plot that based on the categorization of Differentially Methylated Regions (DMRs) 
### @author Fabian Mueller
### @param df2p differential methylation table. Must contain the columns "mean.mean.g1","mean.mean.g2","isDMR"
### @param grp1.name name as it appears on the x axis
### @param grp2.name name as it appears on the y axis
### @return ggplot2 plot
create.diffMeth.bin.dens.dmr.scatter <- function(df2p,grp1.name,grp2.name){
	n.points <- nrow(df2p)
	#plot order: plot DMRs last
	df2p$plotOrder <- NA
	is.dmr <- df2p$isDMR
	is.dmr[is.na(is.dmr)] <- FALSE
	num.not.dmr <- sum(!is.dmr)
	df2p$plotOrder[!is.dmr] <- seq_len(num.not.dmr)
	df2p$plotOrder[is.dmr] <- seq((num.not.dmr+1),n.points)
	
	df2p <- na.omit(df2p[,c("mean.mean.g1","mean.mean.g2","isDMR","plotOrder")])
	
	df2p$color <- NA
	if (sum(!df2p$isDMR)>1){
		colors.nodmr <- DENS.COLORS.LOW[1]
		tryCatch(
			colors.nodmr <- densCols(x=df2p[!df2p$isDMR,"mean.mean.g1"],y=df2p[!df2p$isDMR,"mean.mean.g2"],colramp = colorRampPalette(c(DENS.COLORS.LOW[1],DENS.COLORS.HIGH[1]))),
			error=function(ee){
				logger.warning(c("Could not assess density colors using densCols:",ee$message))
			}
		)
		df2p[!df2p$isDMR,"color"] <- colors.nodmr
	} else if (sum(!df2p$isDMR)==1){
		df2p[!df2p$isDMR,"color"] <- DENS.COLORS.LOW[1]
	}
	if (sum(df2p$isDMR)>1){
		colors.dmr <- DENS.COLORS.LOW[2]
		tryCatch(
			colors.dmr   <- densCols(x=df2p[ df2p$isDMR,"mean.mean.g1"],y=df2p[ df2p$isDMR,"mean.mean.g2"],colramp = colorRampPalette(c(DENS.COLORS.LOW[2],DENS.COLORS.HIGH[2]))),
			error=function(ee){
				logger.warning(c("Could not assess density colors using densCols:",ee$message))
			}
		)
		df2p[df2p$isDMR,"color"] <- colors.dmr
			
	} else if (sum(df2p$isDMR)==1){
		df2p[df2p$isDMR,"color"] <- DENS.COLORS.LOW[2]
	}
	pp <- ggplot(df2p) + aes(mean.mean.g1,mean.mean.g2) +
			labs(x=paste("mean.mean.beta",grp1.name,sep="."),y=paste("mean.mean.beta",grp2.name,sep=".")) + coord_fixed() +
			geom_point(aes(color=color,order=plotOrder)) + scale_color_identity()
	
	return(pp)
}

### addReportPlots.diffMeth.bin.region.scatter
###
### adds report plots for differential methylation for the region level binary case to a report.
### @author Fabian Mueller
### @aliases addReportPlots.diffMeth.bin.region.scatter
### @param report the report to be modified
### @param dmt differential methylation table as created by \code{computeDiffMeth.bin.region}
### @param cmpName Comparison name as it will appear in the filename and figure selection box
### @param regName Region type name as it will appear in the filename and figure selection box
### @param diffSiteRankCut vector of combined ranking cutoffs for classifying a site as differentially methylated
### @param grp1.name name of group 1 in the compoarison (for labelling in the plots)
### @param grp2.name name of group 2 in the compoarison (for labelling in the plots)
### @return list of report plot objects added
addReportPlots.diffMeth.bin.region.scatter <- function(report,dmt,cmpName,regName,diffRegionRankCut,autoRankCut=NULL,grp1.name="Group1",grp2.name="Group2",rerank=TRUE){
	df2p <- dmt #data frame to plot
	figPlots <- list()
	
	#scatterplot based on adjusted p-value significance
	if (is.element("comb.p.adj.fdr",colnames(dmt))){
		df2p$isDMR <- df2p[,"comb.p.adj.fdr"] < P.VAL.CUT

		pp <- create.densityScatter(df2p[,c("mean.mean.g1","mean.mean.g2")],is.special=df2p$isDMR,add.text.cor=TRUE) +
				labs(x=paste("mean.mean.beta",grp1.name,sep="."),y=paste("mean.mean.beta",grp2.name,sep=".")) + coord_fixed()
		cur.cut.name <- "fdrAdjPval"
		figName <- paste("diffMeth_region",cmpName,regName,cur.cut.name,sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots <- c(figPlots,list(report.plot))
	}
	
	rrs <- dmt[,"combinedRank"]
	if (rerank)	rrs <- rank(rrs,na.last="keep",ties.method="min")
	for (i in 1:length(diffRegionRankCut)){
		rc <- diffRegionRankCut[i]
		cur.cut.name <- paste("rc",i,sep="")
		df2p$isDMR <- rrs < rc
		
		pp <- create.densityScatter(df2p[,c("mean.mean.g1","mean.mean.g2")],is.special=df2p$isDMR,add.text.cor=TRUE) +
				labs(x=paste("mean.mean.beta",grp1.name,sep="."),y=paste("mean.mean.beta",grp2.name,sep=".")) + coord_fixed()
		
		figName <- paste("diffMeth_region",cmpName,regName,cur.cut.name,sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots <- c(figPlots,list(report.plot))
	}
	
	if (is.integer(autoRankCut)){
		df2p$isDMR <- dmt[,"combinedRank"] <= autoRankCut
		pp <- create.densityScatter(df2p[,c("mean.mean.g1","mean.mean.g2")],is.special=df2p$isDMR,add.text.cor=TRUE) +
				labs(x=paste("mean.mean.beta",grp1.name,sep="."),y=paste("mean.mean.beta",grp2.name,sep=".")) + coord_fixed()
		figName <- paste("diffMeth_region",cmpName,regName,"rcAuto",sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots <- c(figPlots,list(report.plot))
	}

	figName <- paste("diffMeth_region",cmpName,regName,"rankGradient",sep="_")
	pp <- create.hex.summary.plot(df2p,q="combinedRank") + coord_fixed() +
		labs(x=paste("mean.mean.beta",grp1.name,sep="."),y=paste("mean.mean.beta",grp2.name,sep="."),fill = "median combined rank")
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	return(figPlots)
}

### addReportPlots.diffMeth.bin.region.volcano
###
### adds report volcano plots for differential methylation for the region level binary case to a report.
### @author Fabian Mueller
### @aliases addReportPlots.diffMeth.bin.region.volcano
### @param report the report to be modified
### @param dmt differential methylation table as created by \code{computeDiffMeth.bin.region}
### @param cmpName Comparison name as it will appear in the filename and figure selection box
### @param regName Region type name as it will appear in the filename and figure selection box
### @param grp1.name name of group 1 in the compoarison (for labelling in the plots)
### @param grp2.name name of group 2 in the compoarison (for labelling in the plots)
### @return list of report plot objects added
addReportPlots.diffMeth.bin.region.volcano <- function(report,dmt,cmpName,regName,grp1.name="Group1",grp2.name="Group2"){
	df2p <- dmt #data frame to plot
	figPlots <- list()
	dont.plot.p.val <- all(is.na(df2p[,"comb.p.val"]))
	
	figName <- paste("diffMeth_region_volcano",cmpName,regName,"diff","pVal",sep="_")
	if (!dont.plot.p.val){
		pp <- ggplot(df2p) + aes_string("mean.mean.diff","-log10(comb.p.val)",color="log10(combinedRank)") +
			scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
			geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))#(alpha=0.3)
	} else {
		pp <- rnb.message.plot("No p-value available")
	}
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_region_volcano",cmpName,regName,"diff","pValAdj",sep="_")
	pp <- ggplot(df2p) + aes_string("mean.mean.diff","-log10(comb.p.adj.fdr)",color="log10(combinedRank)") +
		scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
		geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))#(alpha=0.3)
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_region_volcano",cmpName,regName,"quot","pVal",sep="_")
	if (!dont.plot.p.val){
		pp <- ggplot(df2p) + aes_string("mean.mean.quot.log2","-log10(comb.p.val)",color="log10(combinedRank)") +
			scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
			geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))#(alpha=0.3)
	} else {
		pp <- rnb.message.plot("No p-value available")
	}
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_region_volcano",cmpName,regName,"quot","pValAdj",sep="_")
	pp <- ggplot(df2p) + aes_string("mean.mean.quot.log2","-log10(comb.p.adj.fdr)",color="log10(combinedRank)") +
		scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
		geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))#(alpha=0.3)
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_region_volcano",cmpName,regName,"diff","quotSig",sep="_")
	pp <- ggplot(df2p) + aes_string("mean.mean.diff","mean.mean.quot.log2",color="log10(combinedRank)") +
		scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
		geom_point(aes(order=plyr::desc(rank(combinedRank,ties.method="first",na.last=TRUE))))#(alpha=0.3)
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	figName <- paste("diffMeth_region_volcano",cmpName,regName,"quot","quotSig",sep="_")
	pp <- rnb.message.plot("Quotient--Quotient scatterplot not available")
	report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
	report.plot <- off(report.plot,handle.errors=TRUE)
	figPlots <- c(figPlots,list(report.plot))
	
	return(figPlots)
}

### a more robust version of summary.GOHyperGResult (from GOstats)
robustHyperGResultSummary <- function(hgr,maxPval=0.01,htmlLinks=FALSE){
	require(Category)
	AMIGO_URL <- "http://amigo.geneontology.org/amigo/term/%s"
	GOenv <- function(what) {
		annotate::getAnnMap(what, "GO", load=TRUE, type=c("db", "env"))
	}
	goIds <- sigCategories(hgr, maxPval)
	pvals <- pvalues(hgr)[goIds]
	odds <- oddsRatios(hgr)[goIds]
	ec <- expectedCounts(hgr)[goIds]
	cc <- geneCounts(hgr)[goIds]
	ss <- universeCounts(hgr)[goIds]
	goTerms <- sapply(BiocGenerics::mget(goIds, GOenv("TERM"), ifnotfound=NA), Term)
	if (htmlLinks) {
		goTerms <- paste0('<a href="', sprintf(AMIGO_URL, goIds), '">', goTerms, '</a>')
	}
	tt <- data.frame(GOMFID=goIds,Pvalue=pvals,OddsRatio=odds,ExpCount=ec,Count=cc,Size=ss, Term=goTerms)
	return(tt)
}

### addReportPlots.diffMeth.enrich.GO.wordcloud
###
### adds a wordcloud for differential methylation GO annotation to the report
### @author Fabian Mueller
### @param report the report to be modified
### @param hgr GOHyperGREsult object
### @param figName name for the figure as it will appear in the filename
### @return report plot object
addReportPlots.diffMeth.enrich.GO.wordcloud <- function(report,hgr,figName){
	report.plot <- createReportPlot(figName, report,create.pdf=FALSE,high.png=200)
	if (!is.null(hgr)) {	
		tt <- robustHyperGResultSummary(hgr,maxPval=1+1e-10,htmlLinks=FALSE)
		color.scale <- gplots::colorpanel(n=10,low=rnb.getOption("colors.gradient")[2],high=rnb.getOption("colors.gradient")[1])
		suppressWarnings(wordcloud::wordcloud(tt$Term,-log(tt$Pvalue)+1e-10, scale=c(3,.04),min.freq=0,max.words=Inf,
				random.order=T, rot.per=0, colors=color.scale, vfont=c("sans serif","plain")))
	} else {
		print(rnb.message.plot("Enrichment not available"))
	}
	report.plot <- off(report.plot)
	return(report.plot)
}

### rnb.section.diffMeth.introduction
###
### add information to the report for site level analysis
### @author Fabian Mueller
### @aliases rnb.section.diffMeth.introduction
### @param diffmeth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
### @param report report object to which the content is added
### @return the updated report object
rnb.section.diffMeth.introduction <- function(diffmeth,report){
	sectionText <- c("Differential methylation analysis was conducted on site and region level according to the ",
		"sample groups specified in the analysis.")
	report <- rnb.add.section(report, "Introduction: Differential Methylation of Sample Groups", sectionText)
	report <- rnb.add.section(report, "Comparisons", "The following comparisons were made:", level = 2)
	rnb.add.list(report, as.list(get.comparisons(diffmeth)))

	if (.hasSlot(diffmeth,"comparison.info")) { #.hasSlot ensure backwards compatibility
		#construct comparison summary table
		cmp.info <- diffmeth@comparison.info
		is.paired <- vapply(cmp.info,FUN=function(x){x$paired},logical(1))
		is.adj.sva <- vapply(cmp.info,FUN=function(x){x$adj.sva},logical(1))
		is.adj.celltype <- vapply(cmp.info,FUN=function(x){x$adj.celltype},logical(1))
		adj.vars <- lapply(cmp.info,FUN=function(x){
			if (is.null(x$adjustment.table)){
				return(c())
			} else {
				return(colnames(x$adjustment.table))
			}
		})
		has.adj.var <- vapply(adj.vars,FUN=function(x){length(x)>0},logical(1))
		summary.tab <- data.frame(comparison=names(cmp.info))
		if (any(is.paired)){
			summary.tab <- data.frame(summary.tab,paired=is.paired)
		}
		if (any(is.adj.sva)){
			summary.tab <- data.frame(summary.tab,SVAdjust=is.adj.sva)
		}
		if (any(has.adj.var)){
			adj.vars.text <- sapply(adj.vars,FUN=function(x){
				paste(x,collapse=",")
			})
			#write the covariate table to the report directories
			covar.tab.handler <- function(i){
				ccn <- names(cmp.info)[i]
				x <- cmp.info[[i]]
				adj.tab <- x$adjustment.table
				if (is.null(adj.tab)){
					return("")
				}
				if (!is.null(rownames(adj.tab))){
					adj.tab <- data.frame(sample=rownames(adj.tab),adj.tab)
				}
				fname <- paste("adjTable_",i,".csv",sep="")
				fname.rel <- rnb.write.table(
					adj.tab, fname, fpath=rnb.get.directory(report, "data", absolute = TRUE),
					format="csv", gz=FALSE, row.names = FALSE, quote=FALSE
				)
				txt <- paste(c("<a href=\"", rnb.get.directory(report, "data"), "/", fname.rel,"\">","csv","</a>"),collapse="")
				return(txt)
			}
			covar.tab.links <- sapply(1:length(cmp.info),covar.tab.handler)
			summary.tab <- data.frame(summary.tab, adjustment=adj.vars.text, covariateTable=covar.tab.links)
		}
		if (ncol(summary.tab) > 1){
			rownames(summary.tab) <- 1:nrow(summary.tab)
			txt <- c(
				"The table below summarizes information on the comparisons.")
			rnb.add.paragraph(report, txt)
			rnb.add.table(report, summary.tab, first.col.header=TRUE)
		}
	}

	#include information on the p-value method
	site.test.method <- get.site.test.method(diffmeth)
	txt <- c(
		"In the following anlyses, p-values on the site level were computed using the <code>",site.test.method,"</code> method. "
	)
	if (site.test.method == "limma"){
		txt <- c(txt,
			"I.e. hierarchical linear models from the <a href=http://bioinf.wehi.edu.au/limma/>limma</a> package were employed ",
			"and fitted using an empirical Bayes approach on derived M-values."
		)
	} else if (site.test.method == "ttest"){
		txt <- c(txt,
			"I.e. a two-sided Welch t-test was employed."
		)
	} else if (site.test.method == "refFreeEWAS"){
		refText <- c(
			"Houseman, E. A., Molitor, J., and Marsit, C. J. (2014). Reference-Free Cell Mixture Adjustments in Analysis of DNA Methylation Data. ",
			"<i>Bioinformatics</i>, <a href=http://dx.doi.org/doi:10.1093/bioinformatics/btu029>doi:10.1093/bioinformatics/btu029</a>"
		)
		report <- rnb.add.reference(report, refText)

		txt <- c(txt,
			"I.e. hierarchical linear models from the <a href=http://bioinf.wehi.edu.au/limma/>limma</a> package were employed ",
			"and fitted using an empirical Bayes approach on derived M-values.",
			"In this analysis tissue heterogeneity was accounted for by Houseman's reference free EWAS method ",
			rnb.get.reference(report, refText), ".")
	}
	report <- rnb.add.section(report, "P-values", txt, level = 2)

	logger.status("Added introductory section")
	return(report)
}


### rnb.section.replicate.concordance
###
### add information to the report for sample replicat analysis, such as scatterplots
### @author Fabian Mueller
### @param replicateList a list containing replicates as returned by \link{rnb.sample.replicates}
### @param types the vector of site and region types on which the analysis should be conducted
### @param report report object to which the content is added
### @return the updated report object
rnb.section.replicate.concordance <- function(rnbSet,replicateList,types,report){
	if (length(replicateList)<1){
		stop("no valid replicates")
	}
	logger.start("Adding Replicate Level Information")
	sectionText <- c("Sample replicates were compared. This section shows pairwise scatterplots for each sample ",
		"replicate group on both site and region level.")
	report <- rnb.add.section(report, "Analysis of Sample Replicates", sectionText)
	#scatterplots
	logger.start("Adding Scatterplots")
	rep.n <- c()
	rep.v <- c()
	rep.cor <- c()
	addedPlots <- list()
	for (k in 1:length(types)){
		rep.cor.t <- c()
		tt <- types[k]
		cur.type.fname <- paste("type",k,sep="")
		for (i in 1:length(replicateList)){
			rr <- names(replicateList)[i]
			rep.inds <- replicateList[[i]]
			if (length(rep.inds)>1){
				comps <- combn(rep.inds,2) #generate pairwise comparisons --> matrix with 2 rows and (length(rep.inds) choose 2) columns
				for (j in 1:ncol(comps)){
					dd <- data.frame(meth(rnbSet,type=tt)[,comps[,j]])
					cur.cor <- cor(dd[,1],dd[,2],use="pairwise.complete.obs")
					rep.cor.t <- c(rep.cor.t,cur.cor)
					s1 <- colnames(dd)[1]
					s2 <- colnames(dd)[2]
					cur.rep.cmp.fname <- paste("rep",paste(i,j,sep="c"),sep="")
					rep.n <- c(rep.n,cur.rep.cmp.fname)
					rep.v <- c(rep.v,paste(s1," vs. ",s2," (",rr,")",sep=""))
					figName <- paste("replicateScatter",cur.rep.cmp.fname,cur.type.fname,sep="_")
					
					pp <- create.densityScatter(dd[,c(s1,s2)],is.special=NULL,add.text.cor=TRUE) + coord_fixed()
					
					report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
					report.plot <- off(report.plot,handle.errors=TRUE)
					addedPlots <- c(addedPlots,list(report.plot))
				}
			}
		}
		rep.cor <- cbind(rep.cor,rep.cor.t)
		logger.status(c("Processed", tt))
	}
	names(types) <- paste("type",1:length(types),sep="")
	rep.n <- unique(rep.n)
	rep.v <- unique(rep.v)
	names(rep.v) <- rep.n
	setting.names <- list(
			'replicate' = rep.v,
			'site/region' = types)
	description <- 'Scatterplot for replicate methylation comparison.
					The transparency corresponds to point density.
					The 1% of the points in the sparsest populated plot regions are drawn explicitly.'
	report <- rnb.add.figure(report, description, addedPlots, setting.names)
	
	tt <- data.frame(rep.cor)
	colnames(tt) <- types
	rownames(tt) <- rep.v
	tt <- round(tt,4)
	
	txt <- c("The following table contains pearson correlation coefficients:")
	rnb.add.paragraph(report, txt)
	rnb.add.table(report,tt)
	
	logger.completed()
	logger.completed()
	return(report)
}

### rnb.section.diffMeth.site
###
### add information to the report for site level analysis
### @author Fabian Mueller
### @aliases rnb.section.diffMeth.site
### @param diffmeth RnBDiffMeth  object. See \code{\link{RnBDiffMeth-class}} for details.
### @param report report object to which the content is added
### @return the updated report object
rnb.section.diffMeth.site <- function(rnbSet,diffmeth,report,gzTable=FALSE){
	if (length(get.comparisons(diffmeth))<1){
		stop("no valid comparisons")
	}
	diffSiteRankCut <- c(100,1000,10000,100000) #the cutoffs for determining a site as differentially methylated according to combined rank
	logger.start("Adding Site Level Information")
	sectionText <- paste("Differential methylation on the site level was computed based on a variety of metrics. 
						  Of particular interest for the following plots and analyses are the following quantities for each site:
						  a) the difference in mean methylation levels of the two groups being compared, b) the quotient in mean methylation and
						  c) a statistical test (t-test or limma depending on the settings) assessing whether the methylation values in the two groups originate from distinct distributions.
						  Additionally each site was assigned a rank based on each of these three criteria. A combined rank is computed as the maximum (i.e. worst)
						  rank among the three ranks. The smaller the combined rank for a site, the more evidence for differential methylation it exhibits.
						  This section includes scatterplots of the site group means as well as volcano plots
						  of each pairwise comparison colored according to the combined ranks or p-values of a given site.")
	report <- rnb.add.section(report, 'Site Level', sectionText)

	logger.start("Selection of rank cutoffs")
	rank.cuts.auto <- lapply(1:length(get.comparisons(diffmeth)),FUN=function(i){
		cc <- names(get.comparisons(diffmeth))[i]
		ccc <- get.comparisons(diffmeth)[cc]
		dmt <- get.table(diffmeth,ccc,"sites",return.data.frame=TRUE)
		res <- auto.select.rank.cut(dmt$diffmeth.p.adj.fdr,dmt$combinedRank,alpha=0.1)
		return(as.integer(res))
	})
	txt <- paste("The following rank cutfoffs have been automatically selected for the analysis of differentially",
						 "methylated sites:")
	rnb.add.paragraph(report, txt)
	tt <- data.frame(unlist(rank.cuts.auto))
	colnames(tt) <- c("Rank Cutoff")
	rownames(tt) <- get.comparisons(diffmeth)
	rnb.add.table(report,tt)
	logger.completed()
	
	#scatterplots
	logger.start("Adding scatterplots")
	rnb.cleanMem()
	addedPlots <- list()
	if(parallel.isEnabled()){
		addedPlots <- foreach(i=1:length(get.comparisons(diffmeth)),.combine="c") %dopar% {
			cc <- names(get.comparisons(diffmeth))[i]
			ccc <- get.comparisons(diffmeth)[cc]
			dmt <- get.table(diffmeth,ccc,"sites",return.data.frame=TRUE)
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
			grp.names <- get.comparison.grouplabels(diffmeth)[ccc,]
			auto.rank.cut <- rank.cuts.auto[[i]]
			#DEBUG MESSAGE
			logger.info(paste("processing...",cc,"[",i,"] {",ccc,"}{",ccn,"}"))
			res <- addReportPlots.diffMeth.bin.site.scatter(report,dmt,ccn,diffSiteRankCut=diffSiteRankCut,
					autoRankCut=auto.rank.cut,grp1.name=grp.names[1],grp2.name=grp.names[2])
			rnb.cleanMem()
			#DEBUG MESSAGE
			logger.info(paste("done...",cc,"[",i,"] {",ccc,"}{",ccn,"}"))
			res
		}
	} else {
		for (i in 1:length(get.comparisons(diffmeth))){
			cc <- names(get.comparisons(diffmeth))[i]
			ccc <- get.comparisons(diffmeth)[cc]
			dmt <- get.table(diffmeth,ccc,"sites",return.data.frame=TRUE)
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
			grp.names <- get.comparison.grouplabels(diffmeth)[ccc,]
			auto.rank.cut <- rank.cuts.auto[[i]]
			addedPlots <- c(addedPlots,addReportPlots.diffMeth.bin.site.scatter(report,dmt,ccn,diffSiteRankCut=diffSiteRankCut,
							autoRankCut=auto.rank.cut,grp1.name=grp.names[1],grp2.name=grp.names[2]))
			rnb.cleanMem()
		}
	}


	comps <- get.comparisons(diffmeth)
	diffMethType = c(paste("FDR adjusted p-value &lt;",P.VAL.CUT),
					 paste("combined rank among the ",diffSiteRankCut," best ranking sites",sep=""),
					 "automatically selected rank cutoff",
					 "combined rank based gradient","rank permutation test p-value")
	names(diffMethType) = c("fdrAdjPval",paste("rc",1:length(diffSiteRankCut),sep=""),"rcAuto","rankGradient","permutationP")
	setting.names <- list(
		'comparison' = comps,
		'differential methylation measure' = diffMethType)
	description <- c('Scatterplot for differential methylation (sites). If the selected criterion is not <code>rankGradient</code>:
		The transparency corresponds to point density. If the number of points exceeds ',DENS.SCATTER.SUBSAMPLE.THRES,
		' then the number of points for density estimation is reduced to that number by random sampling.',
		'The',round(DENS.SCATTER.SPARSE.POINTS.PERC*100),
		'% of the points in the sparsest populated plot regions are drawn explicitly (up to a maximum of ',DENS.SCATTER.SPARSE.POINTS.MAX,
		" points).",
		'Additionally, the colored points represent differentially methylated sites (according to the selected criterion). 
		If the selected criterion is <code>rankGradient</code>: median combined ranks accross hexagonal bins are shown
		as a gradient according to the color legend.')
	report <- rnb.add.figure(report, description, addedPlots, setting.names)
	logger.completed()
	
	#volcano plots
	logger.start("Adding volcano plots")
	rnb.cleanMem()
	addedPlots <- list()
	if(parallel.isEnabled()){
		addedPlots <- foreach(i=1:length(get.comparisons(diffmeth)),.combine="c") %dopar% {
			cc <- names(get.comparisons(diffmeth))[i]
			ccc <- get.comparisons(diffmeth)[cc]
			dmt <- get.table(diffmeth,ccc,"sites",return.data.frame=TRUE)
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
			grp.names <- get.comparison.grouplabels(diffmeth)[ccc,]
			res <- addReportPlots.diffMeth.bin.site.volcano(report,dmt,ccn,
					grp1.name=grp.names[1],grp2.name=grp.names[2])
			rnb.cleanMem()
			res
		}
	} else {
		for (i in 1:length(get.comparisons(diffmeth))){
			cc <- names(get.comparisons(diffmeth))[i]
			ccc <- get.comparisons(diffmeth)[cc]
			dmt <- get.table(diffmeth,ccc,"sites",return.data.frame=TRUE)
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
			grp.names <- get.comparison.grouplabels(diffmeth)[ccc,]
			addedPlots <- c(addedPlots,addReportPlots.diffMeth.bin.site.volcano(report,dmt,ccn,
							grp1.name=grp.names[1],grp2.name=grp.names[2]))
			rnb.cleanMem()
		}
	}

	diff.measure <- c("diff"="Difference","quot"="Quotient")
	signif.measure <- c("pVal"="p-value","pValAdj"="adjusted p-value","quotSig"="Quotient (only meaningful if 'Difference' is selected above)")
	setting.names <- list(
		'comparison' = comps,
		'difference metric' = diff.measure,
		'significance metric' = signif.measure)
	description <- 'Volcano plot for differential methylation quantified by various metrics. Color scale according to
					combined ranking.'
	report <- rnb.add.figure(report, description, addedPlots, setting.names)
	logger.completed()
	logger.start("Adding tables")
	includeCovg <- !is.null(covg(rnbSet))
	sectionText <- "A tabular overview of measures for differential methylation on the site level for the individual comparisons are provided in this section.
					  Below, a brief explanation of the different columns can be found:
					<ul>
						<li>id: site id</li>
						<li>Chromosome: chromosome of the site</li>
						<li>Start: start coordinate of the site</li>
						<li>Strand: strand of the site</li>
						<li>mean.g1,mean.g2: (where g1 and g2 is replaced by the respective group names in the table) mean methylation in each of the two groups</li>
						<li>mean.diff: difference in methylation means between the two groups: mean.g1-mean.g2. In case of paired analysis, it is the mean of the pairwise differences.</li>
						<li>mean.quot.log2: log2 of the quotient in methylation: log2((mean.g1+epsilon)/(mean.g2+epsilon)), where epsilon:=0.01. In case of paired analysis, it is the mean of the pairwise quotients.</li>
						<li>diffmeth.p.val: p-value obtained from a two-sided Welch t-test or alternatively from linear models employed in the limma package
						(which type of p-value is computed is specified in the differential.site.test.method option). In case of paired analysis, the paired Student's t-test is applied.</li>
						<li>max.g1,max.g2: maximum methylation level in group 1 and 2 respectively</li>
						<li>min.g1,min.g2: minimum methylation level in group 1 and 2 respectively</li>
						<li>sd.g1,sd.g2: standard deviation of methylation levels</li>
						<li>min.diff: Minimum of 0 and the smallest pairwise difference between samples of the two groups</li>
						<li>diffmeth.p.adj.fdr: FDR adjusted p-value of all sites</li>
						<li>combinedRank: mean.diff, mean.quot.log2 and diffmeth.p.val are ranked for all sites. This aggregates them using the maximum, i.e. worst rank of a site among the three measures</li>
						<li>num.na.g1,num.na.g2: number of NA methylation values for groups 1 and 2 respectively</li>"
	
	if (includeCovg){
		ss <- paste("<li>mean.covg.g1,mean.covg.g2: mean coverage of groups 1 and 2 respectively. In case of Infinium array methylation data, coverage is defined as combined beadcount.</li>
			   <li>min.covg.g1,min.covg.g2: minimum coverage of groups 1 and 2 respectively</li>
			   <li>max.covg.g1,max.covg.g2: maximum coverage of groups 1 and 2 respectively</li>
			   <li>covg.thresh.nsamples.g1,covg.thresh.nsamples.g2: number of samples in group 1 and 2 respectively exceeding the coverage threshold (",
	   			get.covg.thres(diffmeth),") for this site.</li>",sep="")
		sectionText <- paste(sectionText,ss,sep="")
	}
	sectionText <- paste(sectionText,"</ul>The tables for the individual comparisons can be found here:\n<ul>\n",sep="")
	annot.cols <- c("Chromosome","Start","Strand")
	sites.info <- annotation(rnbSet,type="sites",add.names=FALSE)[, annot.cols]
	#add cg identifier for infinium datasets
	if (!is.element("RnBiseqSet",class(rnbSet)) && !is.null(rownames(sites.info))){
		sites.info <- data.frame(cgid=rownames(sites.info),sites.info,stringsAsFactors=FALSE)
	}
	grp.names <- get.comparison.grouplabels(diffmeth)
	for (i in 1:length(comps)){
		cc <- comps[i]
		
		annot.vec <- c("mean.g1","mean.g2","mean.diff","mean.quot.log2",
				"diffmeth.p.val","max.g1","min.g1","sd.g1",
				"max.g2","min.g2","sd.g2",
				"min.diff","diffmeth.p.adj.fdr","combinedRank",
				"num.na.g1","num.na.g2")
		if (includeCovg){
			annot.vec <- c(annot.vec,c("mean.covg.g1","mean.covg.g2",
							"min.covg.g1","min.covg.g2","max.covg.g1","max.covg.g2",
							"covg.thresh.nsamples.g1","covg.thresh.nsamples.g2"))
		}
		dmt <- get.table(diffmeth,cc,"sites",return.data.frame=TRUE)[,annot.vec]
		
		g1n <- grp.names[i,1]
		g2n <- grp.names[i,2]
		colname.vec <- c(paste("mean",g1n,sep="."),paste("mean",g2n,sep="."),"mean.diff","mean.quot.log2",
				"diffmeth.p.val",paste("max",g1n,sep="."),paste("min",g1n,sep="."),paste("sd",g1n,sep="."),
				paste("max",g2n,sep="."),paste("min",g2n,sep="."),paste("sd",g2n,sep="."),
				"min.diff","diffmeth.p.adj.fdr","combinedRank",
				paste("num.na",g1n,sep="."),paste("num.na",g2n,sep="."))
		if (includeCovg){
			colname.vec <- c(colname.vec,c(paste("mean.covg",g1n,sep="."),paste("mean.covg",g2n,sep="."),
										   paste("min.covg",g1n,sep="."),paste("min.covg",g2n,sep="."),
										   paste("max.covg",g1n,sep="."),paste("max.covg",g2n,sep="."),
										   paste("nsamples.covg",paste("thres",get.covg.thres(diffmeth),sep=""),g1n,sep="."),
										   paste("nsamples.covg",paste("thres",get.covg.thres(diffmeth),sep=""),g2n,sep=".")))
		}
		colnames(dmt) <- colname.vec
		dmt <- cbind(rownames(dmt),sites.info,dmt)
		colnames(dmt)[1] <- "id"
		rownames(dmt) <- NULL
		ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
		fname <- paste("diffMethTable_site_",ccn,".csv",sep="")
		fname <- rnb.write.table(dmt,fname,fpath=rnb.get.directory(report, "data", absolute = TRUE),format="csv",gz=gzTable,row.names = FALSE,quote=FALSE)
		txt <- paste(c("<a href=\"", rnb.get.directory(report, "data"), "/", fname,"\">",cc,"</a>"),collapse="")
		sectionText <- paste(sectionText,"<li>",txt,"</li>\n",sep="")
	}
	sectionText <- paste(sectionText,"</ul>",sep="")
	report <- rnb.add.section(report, "Differential Methylation Tables", sectionText, level = 2)
	logger.completed()
	
	logger.completed()
	return(report)
}

### rnb.section.diffMeth.region
###
### add information to the report for region level analysis
### @author Fabian Mueller
### @aliases rnb.section.diffMeth.region
### @param diffmeth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
### @param report report object to which the content is added
### @param dm.enrich If enrichment analysis reports are desired this argument should not be \code{NULL} (which is the default value
###                  it should be an object of type \code{DiffMeth.enrich} (see \code{performEnrichment.diffMeth()} for details)
### @return the updated report object
rnb.section.diffMeth.region <- function(rnbSet,diffmeth,report,dm.enrich=NULL,gzTable=FALSE){
	if (length(get.comparisons(diffmeth))<1){
		stop("no valid comparisons")
	}
	if (length(get.region.types(diffmeth))<1){
		stop("no valid region types")
	}

	diffRegionRankCut <- c(100,500,1000) #the cutoffs for determining a site as differentially methylated according to combined rank
	logger.start("Adding Region Level Information")
	refText <- c("Makambi, K. (2003) Weighted inverse chi-square method for correlated significance tests. ",
		"<i>Journal of Applied Statistics</i>, <b>30</b>(2), 225234")
	report <- rnb.add.reference(report, refText)
	
	sectionText <- c("Differential methylation on the region level was computed based on a variety of metrics. ", 
		"Of particular interest for the following plots and analyses are the following quantities for each region: ",
		"the mean difference in means across all sites in a region of the two groups being compared and the mean of quotients ",
		"in mean methylation as well as a combined p-value calculated from all site p-values in the region ",
		rnb.get.reference(report, refText), ". ",
		"Additionally each region was assigned a rank based on each of these three criteria. ",
		"A combined rank is computed as the maximum (i.e. worst) value among the three ranks. The smaller the combined rank for a region, the more evidence for differential methylation it exhibits. ",
		"Regions were defined based on the region types specified in the analysis. ",
		"This section includes scatterplots of the region group means as well as volcano plots of each pairwise comparison ",
		"colored according to the combined rank of a given region.")
	report <- rnb.add.section(report, "Region Level", sectionText)

	comps <- get.comparisons(diffmeth)
	reg.types <- get.region.types(diffmeth)

	logger.start("Selection of rank cutoffs")
	rank.cuts.auto <- lapply(1:length(comps),FUN=function(i){
		lapply(1:length(reg.types),FUN=function(j){
			dmt <- get.table(diffmeth,comps[i],reg.types[j],return.data.frame=TRUE)
			res <- auto.select.rank.cut(dmt$comb.p.adj.fdr,dmt$combinedRank,alpha=0.1)
			return(as.integer(res))
		})
	})
	txt <- paste("The following rank cutfoffs have been automatically selected for the analysis of differentially",
						 "methylated regions:")
	rnb.add.paragraph(report, txt)

	tt <- data.frame(matrix(unlist(rank.cuts.auto),ncol=length(reg.types),nrow=length(comps),byrow=TRUE))
	colnames(tt) <- reg.types
	rownames(tt) <- comps
	rnb.add.table(report,tt)
	logger.completed()

	#scatterplots
	logger.start("Adding scatterplots")
	rnb.cleanMem()
	addedPlots <- list()
	grp.labels <- get.comparison.grouplabels(diffmeth)
	if(parallel.isEnabled()){
		#generate pairs of comparison, region combinations with indices
		iis <- 1:length(comps)
		jjs <- 1:length(reg.types)
		pps <- expand.grid(iis,jjs)
		
		addedPlots <- foreach(k=1:nrow(pps),.combine="c") %dopar% {
			i <- pps[k,1]
			j <- pps[k,2]
			cc <- names(comps)[i]
			ccc <- comps[cc]
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
			rr <- reg.types[j]
			rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
			auto.rank.cut <- rank.cuts.auto[[i]][[j]]
			dmt <- get.table(diffmeth,ccc,rr,return.data.frame=TRUE)
			res <- addReportPlots.diffMeth.bin.region.scatter(report,dmt,ccn,rrn,diffRegionRankCut=diffRegionRankCut,
					autoRankCut=auto.rank.cut,grp1.name=grp.labels[ccc,1],grp2.name=grp.labels[ccc,2])
			rnb.cleanMem()
			res
		}
		
	} else {
		for (i in 1:length(comps)){
			cc <- names(comps)[i]
			ccc <- comps[cc]
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
			for (j in 1:length(reg.types)){
				rr <- reg.types[j]
				rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
				auto.rank.cut <- rank.cuts.auto[[i]][[j]]
				dmt <- get.table(diffmeth,ccc,rr,return.data.frame=TRUE)
				addedPlots <- c(addedPlots,addReportPlots.diffMeth.bin.region.scatter(report,dmt,ccn,rrn,diffRegionRankCut=diffRegionRankCut,
								autoRankCut=auto.rank.cut,grp1.name=grp.labels[ccc,1],grp2.name=grp.labels[ccc,2]))
				rnb.cleanMem()
			}
		}
	}
	
	names(reg.types) <- ifelse(is.valid.fname(reg.types),reg.types,paste("reg",1:length(reg.types),sep=""))
	diffMethType = c(paste("FDR adjusted p-value <",P.VAL.CUT),
					 paste("combined rank among the ",diffRegionRankCut," best ranking regions",sep=""),
					 "automatically selected rank cutoff","combined rank based gradient")
	names(diffMethType) = c("fdrAdjPval",paste("rc",1:length(diffRegionRankCut),sep=""),"rcAuto","rankGradient")
	setting.names <- list(
		'comparison' = comps,
		'regions' = reg.types ,
		'differential methylation measure' = diffMethType)
	description <- 'Scatterplot for differential methylation (regions). If the selected criterion is not <code>rankGradient</code>:
		The transparency corresponds to point density. The 1% of the points in the sparsest populated plot regions are drawn explicitly.
		Additionally, the colored points represent differentially methylated regions (according to the selected criterion). 
		If the selected criterion is <code>rankGradient</code>: median combined ranks accross hexagonal bins are shown
		as a gradient according to the color legend.'
	report <- rnb.add.figure(report, description, addedPlots, setting.names)
	logger.completed()
	
	#volcano plots
	logger.start("Adding volcano plots")
	rnb.cleanMem()
	addedPlots <- list()
	if(parallel.isEnabled()){
		#generate pairs of comparison, region combinations with indices
		iis <- 1:length(comps)
		jjs <- 1:length(reg.types)
		pps <- expand.grid(iis,jjs)
		
		addedPlots <- foreach(k=1:nrow(pps),.combine="c") %dopar% {
			i <- pps[k,1]
			j <- pps[k,2]
			cc <- names(comps)[i]
			ccc <- comps[cc]
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
			rr <- reg.types[j]
			rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
			dmt <- get.table(diffmeth,ccc,rr,return.data.frame=TRUE)
			res <- addReportPlots.diffMeth.bin.region.volcano(report,dmt,ccn,rrn,
					grp1.name=grp.labels[ccc,1],grp2.name=grp.labels[ccc,2])
			rnb.cleanMem()
			res
		}
	} else {
		for (i in 1:length(comps)){
			cc <- names(comps)[i]
			ccc <- comps[cc]
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
			for (j in 1:length(reg.types)){
				rr <- reg.types[j]
				rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
				dmt <- get.table(diffmeth,ccc,rr,return.data.frame=TRUE)
				addedPlots <- c(addedPlots,addReportPlots.diffMeth.bin.region.volcano(report,dmt,ccn,rrn,
								grp1.name=grp.labels[ccc,1],grp2.name=grp.labels[ccc,2]))
				rnb.cleanMem()
			}
		}
	}
	diff.measure <- c("diff"="Difference","quot"="Quotient")
	signif.measure <- c("pVal"="combined p-value","pValAdj"="adjusted combined p-value","quotSig"="Quotient (only meaningful if 'Difference' is selected above)")
	setting.names <- list(
		'comparison' = comps,
		'regions' = reg.types ,
		'difference metric' = diff.measure,
		'significance metric' = signif.measure)
	description <- 'Volcano plot for differential methylation quantified by various metrics. Color scale according to
		combined ranking.'
	report <- rnb.add.figure(report, description, addedPlots, setting.names)
	logger.completed()

	logger.start("Adding tables")
	includeCovg <- !is.null(covg(rnbSet))
	sectionText <- c("A tabular overview of measures for differential methylation on the region level for the ",
		"individual comparisons are provided in this section.")
	report <- rnb.add.section(report, "Differential Methylation Tables", sectionText, level = 2)
	sectionText <- list(
		"id: region id",
		"Chromosome: chromosome of the region",
		"Start: Start coordinate of the region",
		"End: End coordinate of the region",
		"[symbol]: associated gene symbol to the given region [only valid for gene associated regions]",
		"[entrezID]: Entrez ID of the gene associated with the region [only valid for gene associated regions]",
		"mean.mean.g1,mean.mean.g2: (where g1 and g2 is replaced by the respective group names in the table) mean of mean methylation levels for group 1 and 2 across all sites in a region",
		"mean.mean.diff: Mean difference in means across all sites in a region",
		"mean.mean.quot.log2: log2 of the mean quotient in means across all sites in a region",
		c("comb.p.val: Combined p-value aggregating p-values of all sites in the region using a generalization of Fisher's method ", rnb.get.reference(report, refText)),
		"comb.p.adj.fdr: FDR adjusted combined p-value",
		c("combinedRank: mean.mean.diff, mean.mean.quot.log2 and comb.p.val are ranked for all regions. ",
				"This column aggregates them using the maximum, i.e. worst rank of a site among the three measures"),
		"num.sites: number of sites associated with the region",
		"mean.num.na.g1,mean.num.na.g2: Mean number of NA methylation values accross all sites in group 1 and group 2 respectively"
		)
	if (includeCovg){
		sectionText <- c(sectionText,list(
		"mean.mean.covg.g1,mean.mean.covg.g2: Mean value of mean coverage values (across all samples in a group) across all sites in a region",
		c("mean.nsamples.covg.thresh.g1,mean.nsamples.covg.thresh.g2: mean number of samples (accross all considered sites) that have a coverage larger than ",
		  get.covg.thres(diffmeth)," for the site in group 1 and group 2 respectively")
		))
	}
	rnb.add.list(report, sectionText)
	sectionText <- "The tables for the individual comparisons can be found here:"
	rnb.add.paragraph(report, sectionText)
	#create data tables
	region.info.cols <- c("Chromosome","Start","End","symbol","entrezID")#additional information to be included in the output table
	reg.type.infos <- lapply(reg.types,FUN=function(rn){
		aa <- annotation(rnbSet,type=rn,add.names=FALSE)
		return(aa)
	})
	names(reg.type.infos) <- reg.types
	file.tab <- do.call("rbind",lapply(1:length(comps),FUN=function(ic){
		cc <- comps[ic]
		sapply(1:length(reg.types),FUN=function(ir){
			rr <- reg.types[ir]
			reg.info <- reg.type.infos[[rr]]
			region.info.cols.cur <- intersect(region.info.cols,colnames(reg.info))
			reg.info <- reg.info[,region.info.cols.cur]
			
			annot.vec <- c("mean.mean.g1","mean.mean.g2",
					"mean.mean.diff","mean.mean.quot.log2",
					"comb.p.val","comb.p.adj.fdr","combinedRank",
					"num.sites","mean.num.na.g1","mean.num.na.g2")
			if (includeCovg){
				annot.vec <- c(annot.vec,c("mean.mean.covg.g1","mean.mean.covg.g2",
										   "mean.nsamples.covg.thresh.g1","mean.nsamples.covg.thresh.g2"))
			}
			
			dmt <- get.table(diffmeth,cc,rr,return.data.frame=TRUE)[,annot.vec]
			g1n <- grp.labels[ic,1]
			g2n <- grp.labels[ic,2]
			colname.vec <- c(paste("mean.mean",g1n,sep="."),paste("mean.mean",g2n,sep="."),
					"mean.mean.diff","mean.mean.quot.log2",
					"comb.p.val","comb.p.adj.fdr","combinedRank",
					"num.sites",paste("mean.num.na",g1n,sep="."),paste("mean.num.na",g2n,sep="."))
			if (includeCovg){
				colname.vec <- c(colname.vec,c(paste("mean.mean.covg",g1n,sep="."),paste("mean.mean.covg",g2n,sep="."),
											   paste("mean.nsamples.covg",paste("thres",get.covg.thres(diffmeth),sep=""),g1n,sep="."),
											   paste("mean.nsamples.covg",paste("thres",get.covg.thres(diffmeth),sep=""),g2n,sep=".")))
			}
			colnames(dmt) <- colname.vec
			dmt <- cbind("id"=rownames(reg.info),reg.info,dmt)
			
			ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",ic,sep=""))
			rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",ir,sep=""))
			fname <- paste("diffMethTable_region_",ccn,"_",rrn,".csv",sep="")
			fname <- rnb.write.table(dmt,fname,fpath=rnb.get.directory(report, "data", absolute = TRUE),format="csv",gz=gzTable,row.names=FALSE,quote=FALSE)
			txt <- paste(c("<a href=\"", rnb.get.directory(report, "data"), "/", fname,"\">","csv","</a>"),collapse="")
			return(txt)
		})
	}))
	rownames(file.tab) <- comps
	colnames(file.tab) <- reg.types
	rnb.add.table(report,file.tab)
	logger.completed()
	
	sectionText <- "No Enrichment Analysis was conducted"
	if (class(dm.enrich)=="DiffMeth.enrich" & length(dm.enrich$region)>0){
		sectionText <- "Enrichment Analysis was conducted. The wordclouds and tables below contains significant GO terms as determined by a hypergeometric test."
	}
	report <- rnb.add.section(report, 'Enrichment Analysis', sectionText, level = 2)
	if (class(dm.enrich)=="DiffMeth.enrich" && length(dm.enrich$region)>0){
		logger.start("Adding enrichment analysis results")
		require(annotate)
		comps <- names(dm.enrich$region)
		names(comps) <- paste("comp",1:length(comps),sep="")
		ontol <- names(dm.enrich$region[[1]])
		names(ontol) <- paste("ontol",1:length(ontol),sep="")
		names(reg.types) <- paste("reg",1:length(reg.types),sep="")
		hyper.hypo <- c("hypermethylated","hypomethylated")
		names(hyper.hypo) <- c("hyper","hypo")
		#reorder to make genes the default selection if present
		genes.ind <- which(reg.types=="genes")
		if (length(genes.ind)>0){
			reg.types <- reg.types[c(genes.ind,setdiff(1:length(reg.types),genes.ind))]
		}
		rank.cuts <- paste("combined rank among the ",diffRegionRankCut," best ranking regions",sep="")
		rank.cuts <- c(rank.cuts,paste("automatically selected rank cutoff"))
		rank.cuts.names.dm.enrich <-  paste("rankCut_",c(diffRegionRankCut,"autoSelect"),sep="")
		names(rank.cuts) <- c(paste("rc",1:length(diffRegionRankCut),sep=""),"rcAuto")
		names(rank.cuts.names.dm.enrich) <- c(paste("rc",1:length(diffRegionRankCut),sep=""),"rcAuto")
		setting.names <- list(
				'comparison' = comps,
				'Hypermethylation/hypomethylation' = hyper.hypo,
				'ontology' = ontol,
				'regions' = reg.types,
				'differential methylation measure' = rank.cuts)

		colnames2round <- c("Pvalue","OddsRatio","ExpCount")
		do.enrichment.table <- function(ccn,hhn,oon,rrn,rcn){
			require(Category)
			cc <- comps[ccn]
			hh <- hyper.hypo[hhn]
			oo <- ontol[oon]
			rr <- reg.types[rrn]
			rc <- rank.cuts.names.dm.enrich[rcn]

			ee <- dm.enrich$region[[cc]][[oo]][[rr]][[rc]][[hhn]]
			kk <- paste(c(ccn,hhn,oon,rrn,rcn),collapse="_")
			if (!is.null(ee)){
				if (length(sigCategories(ee))>0){
					tt <- robustHyperGResultSummary(ee,htmlLinks=TRUE)
					tt[,colnames2round] <- round(tt[,colnames2round],4)
				}
				else {
					tt <- data.frame("NA"=NA)
				}
			} else {
				tt <- data.frame("NA"=NA)
			}
			return(tt)
		}
		do.enrichment.wordcloud <- function(ccn,hhn,oon,rrn,rcn){
			cc <- comps[ccn]
			hh <- hyper.hypo[hhn]
			oo <- ontol[oon]
			rr <- reg.types[rrn]
			rc <- rank.cuts.names.dm.enrich[rcn]

			ee <- dm.enrich$region[[cc]][[oo]][[rr]][[rc]][[hhn]]
							
			kk <- paste(c(ccn,hhn,oon,rrn,rcn),collapse="_")
			figName <- paste("enrichGOwordcloud_",kk,sep="")
			report.plot <- addReportPlots.diffMeth.enrich.GO.wordcloud(report,ee,figName)
			return(report.plot)
		}

		#generate tuples of parameter combinations
		pps <- expand.grid(names(comps),names(hyper.hypo),names(ontol),names(reg.types),names(rank.cuts),stringsAsFactors=FALSE)
		kks <- paste(pps[,1],pps[,2],pps[,3],pps[,4],pps[,5],sep="_")

		tabs2write <- lapply(1:nrow(pps),FUN=function(k){
			do.enrichment.table(pps[k,1],pps[k,2],pps[k,3],pps[k,4],pps[k,5])

		})
		names(tabs2write) <- kks
		addedPlots <- lapply(1:nrow(pps),FUN=function(k){
			do.enrichment.wordcloud(pps[k,1],pps[k,2],pps[k,3],pps[k,4],pps[k,5])

		})
		names(addedPlots) <- kks

		description <- "Wordclouds for GO enrichment terms."
		report <- rnb.add.figure(report, description, addedPlots, setting.names)
		report <- rnb.add.tables(report, tabs2write, setting.names, row.names = FALSE)
		logger.completed()
	}
	
	logger.completed()
	return(report)
}
#' get.adjustment.variables
#'
#' Given indices for two groups of samples for comparison, this function
#' retrieves \code{data.frame} containing the variables to be adjusted for
#' @author Fabian Mueller
#' @param rnbSet RnBSet object
#' @param inds.g1 sample indices in \code{rnbSet} of group 1 members
#' @param inds.g2 sample indices in \code{rnbSet} of group 2 members
#' @param colnames.adj column names in \code{pheno(rnbSet)} to retrieve
#' @param colname.target column names in \code{pheno(rnbSet)} of the target variable. Only important if \code{adjust.sva==TRUE}
#' @param adjust.sva flag indicating whether the resulting table should also contain surrogate variables (SVs) for the given target variable.
#' @param adjust.celltype flag indicating whether the resulting table should also contain estimated celltype contributions.
#' 				See \code{\link{rnb.execute.ct.estimation}} for details.
#' @return a \code{data.frame} containing one column for each selected variable from the phenotypic data
#'         each row corresponds to a sample in the union of samples of the wto groups with the first
#'         \code{length(inds.g1)} rows corresponding to group 1 and the remaining rows corresponding to group 2
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
#' get.adjustment.variables(rnb.set.example,sample.groups[[1]],sample.groups[[2]],"Cell_Line")
#' }
get.adjustment.variables <- function(rnbSet, inds.g1, inds.g2=-inds.g1, colnames.adj=c(), colname.target="",
		adjust.sva=FALSE, adjust.celltype=FALSE){
	if (!all(colnames.adj %in% colnames(pheno(rnbSet)))){
		stop("invalid adjustment columns specified. please make sure they are a subset of the sample annotaiton table")
	}
	if (length(inds.g1)>=nrow(pheno(rnbSet))){
		stop("invalid group indexing. You must specify 2 groups")
	}
	if (is.element(colname.target,colnames.adj)){
		logger.info(c(colname.target, "is declared as part of the adjustment covariates. --> removed it from the list"))
		colnames.adj <- setdiff(colnames.adj,colname.target)
	}
	sample.names <- c(samples(rnbSet)[inds.g1],samples(rnbSet)[inds.g2])
	res <- data.frame(matrix(nrow=length(sample.names),ncol=0))
	rownames(res) <- sample.names
	if (length(colnames.adj)>0){
		res <- rbind(pheno(rnbSet)[inds.g1,colnames.adj,drop=FALSE],pheno(rnbSet)[inds.g2,colnames.adj,drop=FALSE])
	}
	inds.g1.after <- 1:length(inds.g1)
	inds.g2.after <- (length(inds.g1)+1):nrow(res)
	for (cc in colnames.adj){
		#remove variables which contain NAs
		contains.na <- any(is.na(res[,cc]))
		#remove constant variables
		is.const <- length(unique(res[,cc])) == 1
		#remove unique variables
		is.uniq <- !is.numeric(res[,cc]) && length(unique(res[,cc])) == nrow(res)
		#remove variables in which the values coincide completely with the grouping
		coinc.group <- length(unique(res[inds.g1.after,cc]))==1 && length(unique(res[inds.g2.after,cc]))==1
		#convert character to factor
		if (contains.na || is.const || is.uniq || coinc.group) {
			res <- res[,setdiff(colnames(res),cc),drop=FALSE]
		} else {
			if (is.character(res[,cc])){
				res[,cc] <- factor(res[,cc])
			}
		}
	}
	#SVs
	if (adjust.sva){
		if (has.covariates.sva(rnbSet, colname.target)){
			sv.tab <- get.covariates.sva(rnbSet, colname.target)
			if (is.null(colnames(sv.tab))) colnames(sv.tab) <- paste0("sv",1:ncol(sv.tab))
			res.sv.tab <- rbind(sv.tab[inds.g1,,drop=FALSE],sv.tab[inds.g2,,drop=FALSE])
			res <- cbind(res,res.sv.tab)
		} else {
			logger.warning(c("Could not retrieve surrogate variables for target '",colname.target,"'"))
		}
	}
	# Celltypes
	if (adjust.celltype){
		if (has.covariates.ct(rnbSet)){
			if (!isTRUE(colname.target == attr(get.covariates.ct(rnbSet), "column"))){
				ct.tab <- get.covariates.ct(rnbSet)
				if (is.null(colnames(ct.tab))) colnames(ct.tab) <- paste0("celltype",1:ncol(sv.tab))
				res.ct.tab <- rbind(ct.tab[inds.g1,,drop=FALSE],ct.tab[inds.g2,,drop=FALSE])
				res <- cbind(res,res.ct.tab)
			}
		} else {
			logger.warning(c("Could not retrieve celltype contributions for target '",colname.target,"'"))
		}
	}
	return(res)
}

#' get.comparison.info
#'
#' retrieve the comparison information for an RnBSet object
#' @author Fabian Mueller
#' @param x \code{RnBSet} object
#' @param pheno.cols column names of the pheno slot in \code{x} on which the dataset should be partitioned. Those columns are required to be factors or logical.
#' 				     In case of factors, each group in turn will be compared to all other groups
#' @param region.types which region types should be processed for differential methylation
#' @param pheno.cols.all.pairwise integer or character vector specifying the colomns of \code{pheno(x)} on which all pairwise comparisons should be conducted.
#' 		 A value of \code{NULL} indicates no columns.
#' @param columns.pairs argument passed on to \code{rnb.sample.groups}. See its documentation for details.
#' @param columns.adj Column names or indices in the table of phenotypic information to be used for confounder adjustment in the
#'        differential methylation analysis.
#' @param adjust.sva flag indicating whether the adjustment table should also contain surrogate variables (SVs) for the given target variable.
#' @param adjust.celltype flag indicating whether the resulting table should also contain estimated celltype contributions.
#' 				See \code{\link{rnb.execute.ct.estimation}} for details.
#' @param pheno.cols.adjust.sva Target variables for SVA adjustment. Only important if \code{adjust.sva==TRUE}. Only the intersection of
#'			\code{pheno.cols} and \code{pheno.cols.adjust.sva} is considered for SVA adjustment.
#' @param adjust.na.rm Flag indicating whether NAs in the adjustment table should be removed.
#' @return a list containing one element for each comparison to be conducted. Each element is again a list containing:
#' \describe{
#'   \item{\code{comparison}}{the name of the comparison}
#'   \item{\code{pheno.colname}}{the column name of the sample annotation table the comparison is derived from}
#'   \item{\code{group.names}}{the names of the two groups being compared}
#'   \item{\code{group.inds}}{the sample indices of the samples belonging to the two groups}
#'   \item{\code{paired}}{flag indicating whether paired analysis is conducted}
#'   \item{\code{adj.sva}}{flag indicating whether adjustment for SVA is conducted}
#'   \item{\code{adj.celltype}}{flag indicating whether adjustment for cell type is conducted}
#'   \item{\code{adjustment.table}}{the covariate adjustment table. \code{NULL} if the comparison is not adjusted}
#'   \item{\code{region.types}}{the region types applicable to the analysis}
#' }
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' cmp.info <- get.comparison.info(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' cmp.info[[1]]
#' }
get.comparison.info <- function(x, pheno.cols=rnb.getOption("differential.comparison.columns"),
		region.types=rnb.region.types.for.analysis(x),
		pheno.cols.all.pairwise=rnb.getOption("differential.comparison.columns.all.pairwise"), 
		columns.pairs=rnb.getOption("columns.pairing"), columns.adj=rnb.getOption("covariate.adjustment.columns"),
		adjust.sva=rnb.getOption("differential.adjustment.sva"), pheno.cols.adjust.sva=rnb.getOption("inference.targets.sva"),
		adjust.celltype=rnb.getOption("differential.adjustment.celltype"),
		adjust.na.rm=TRUE){
	#example
	#rr <- get.comparison.info(rnb.set, pheno.cols=rnb.getOption("differential.comparison.columns"))
	if (is.null(pheno.cols)){
		pheno.cols <- colnames(pheno(x))
	}
	group.info <- rnb.sample.groups(x, columns = pheno.cols, columns.pairs = columns.pairs)
	if (length(group.info) < 1){
		logger.warning("No valid grouping information found. NULL returned")
		return(NULL)
	}
	group.info.all.pairwise <- rep(FALSE,length(group.info))
	names(group.info.all.pairwise) <- names(group.info)
	if (!is.null(pheno.cols.all.pairwise)){
		if (is.integer(pheno.cols.all.pairwise)){
			ccns <- colnames(pheno(x))
			if (any(pheno.cols.all.pairwise<1) | any(pheno.cols.all.pairwise>length(ccns))){
				logger.warning("Invalid value for pheno.cols.all.pairwise. NULL returned")
				return(NULL)
			}
			pheno.cols.all.pairwise <- ccns[pheno.cols.all.pairwise]
		}

		pheno.cols.all.pairwise <- intersect(pheno.cols.all.pairwise,names(group.info.all.pairwise))
		if (length(pheno.cols.all.pairwise)<1) {
			logger.warning("No pariwise comparison specifier could be matched to the comparisons conducted")
		} else {
			logger.info(c("Conducting all pairwise comparisons for columns:",paste(pheno.cols.all.pairwise,collapse=",")))
			logger.info("All pairwise comparisons are performed on a subset of sample groupings. Caution: this could increase the runtime significantly
						 due to combinatorial explosion.")
			group.info.all.pairwise[pheno.cols.all.pairwise] <- TRUE
		}
	}
	sva.targets <- intersect(pheno.cols,pheno.cols.adjust.sva)

	#helper function for excluding samples for which the adjustment table contains NAs
	rm.na.from.adj.tab <- function(rnb.set,adj.var.tab,inds1,inds2,cmp.name){
		res <- list(
			adj.var.tab=adj.var.tab,
			inds1 = inds1,
			inds2 = inds2,
			samples.excluded = c()
		)
		if (any(is.na(adj.var.tab))){
			sample.name2ind <- 1:length(samples(rnb.set))
			names(sample.name2ind) <- samples(rnb.set)
			samples.na <- apply(adj.var.tab,1,FUN=function(x){any(is.na(x))})
			names(samples.na) <- rownames(adj.var.tab)
			sample.names.na <- rownames(adj.var.tab)[samples.na]
			sample.inds.na <- sample.name2ind[sample.names.na]

			res$inds1 <- setdiff(inds1,sample.inds.na)
			res$inds2 <- setdiff(inds2,sample.inds.na)
			res$adj.var.tab <- na.omit(adj.var.tab)
			res$samples.excluded <- sample.names.na
			warn.msg <- paste0("The following samples were excluded from comparison '",cmp.name,
				"' because they contain missing values in the adjustment table:", paste(sample.names.na,collapse=",")
			)
			logger.warning(warn.msg)
		}
		return(res)
	}

	res <- list()
	for (i in 1:length(group.info)){
		groups <- group.info[[i]]
		is.paired <- attr(group.info,"paired")[i]
		cc <- names(group.info)[i]
		
		adj.cols <- setdiff(columns.adj,cc)
		adj.sva <- cc %in% sva.targets && adjust.sva && has.covariates.sva(x, cc)
		## FIXME: Cell Type adjustment is disabled when running diff. methylation on the same column
		## It must be disabled in all traits for which adjustment cannot be made + we need to announce this to the user
		adj.celltype <- adjust.celltype && has.covariates.ct(x) && !isTRUE(cc == attr(get.covariates.ct(x), "column"))
		get.adj.tab <- length(adj.cols)>0 || adj.sva || adj.celltype

		if (length(groups)>2){
			if (group.info.all.pairwise[i]) {
				#perform all pairwise comparisons
				pps <- combn(1:length(groups),2)
				for (j in 1:ncol(pps)){
					pp1 <- pps[1,j]
					pp2 <- pps[2,j]
					ll1 <- names(groups)[pp1]
					ll2 <- names(groups)[pp2]
					grps <- c(ll1,ll2)
					comparison <- paste(paste(grps,collapse=" vs. ")," (based on ",cc,")",sep="")

					inds1 <- groups[[pp1]]
					inds2 <- groups[[pp2]]
					adj.var.tab <- NULL
					if (get.adj.tab) {
						adj.var.tab <- get.adjustment.variables(x,inds1,inds2,adj.cols, colname.target=cc, adjust.sva=adj.sva, adjust.celltype=adj.celltype)
					}
					if (adjust.na.rm & !is.null(adj.var.tab)) {
						adj.clean.helper.object <- rm.na.from.adj.tab(x,adj.var.tab,inds1,inds2,comparison)
						inds1 <- adj.clean.helper.object$inds1
						inds2 <- adj.clean.helper.object$inds2
						adj.var.tab <- adj.clean.helper.object$adj.var.tab
					}
					res.cur <- list(
						comparison=comparison,
						pheno.colname=cc,
						group.names=grps,
						group.inds=list(group1=inds1,group2=inds2),
						paired=is.paired,
						adj.sva=adj.sva,
						adj.celltype=adj.celltype,
						adjustment.table=adj.var.tab,
						region.types=region.types
					)
					res.append <- list(res.cur)
					names(res.append) <- res.cur$comparison
					res <- c(res,res.append)
				}
			} else {
				for (j in 1:length(groups)){
					ll <- names(groups)[j]
					grps <- c(ll,paste("non.",ll,sep=""))
					if (is.paired){
						logger.warning(c("Paired analysis is not supported annotations with more than 2 categories and comparing one group vs. all others. ",
									"--> Using unpaired analysis.",
									"Consider enabling the differential.comparison.columns.all.pairwise option or reducing the number of groups ",
									"in this column to 2."))
						is.paired <- FALSE
					}
					comparison <- paste(paste(grps,collapse=" vs. ")," (based on ",cc,")",sep="")
					
					inds1 <- groups[[j]]
					inds2 <- setdiff(unlist(groups),inds1)
					adj.var.tab <- NULL
					if (get.adj.tab) {
						adj.var.tab <- get.adjustment.variables(x,inds1,inds2,adj.cols, colname.target=cc, adjust.sva=adj.sva, adjust.celltype=adj.celltype)
					}
					if (adjust.na.rm & !is.null(adj.var.tab)) {
						adj.clean.helper.object <- rm.na.from.adj.tab(x,adj.var.tab,inds1,inds2,comparison)
						inds1 <- adj.clean.helper.object$inds1
						inds2 <- adj.clean.helper.object$inds2
						adj.var.tab <- adj.clean.helper.object$adj.var.tab
					}
					res.cur <- list(
						comparison=comparison,
						pheno.colname=cc,
						group.names=grps,
						group.inds=list(group1=inds1,group2=inds2),
						paired=is.paired,
						adj.sva=adj.sva,
						adj.celltype=adj.celltype,
						adjustment.table=adj.var.tab,
						region.types=region.types
					)
					res.append <- list(res.cur)
					names(res.append) <- res.cur$comparison
					res <- c(res,res.append)
				}
			}
		} else { # length(groups) == 2
			ll1 <- names(groups)[1]
			ll2 <- names(groups)[2]
			grps <- c(ll1,ll2)
			comparison <- paste0(paste(grps,collapse=" vs. ")," (based on ",cc,")")
			inds1 <- groups[[1]]
			inds2 <- groups[[2]]
			adj.var.tab <- NULL
			if (get.adj.tab) {
				adj.var.tab <- get.adjustment.variables(x,inds1,inds2,adj.cols, colname.target=cc, adjust.sva=adj.sva, adjust.celltype=adj.celltype)
			}
			if (adjust.na.rm & !is.null(adj.var.tab)) {
				adj.clean.helper.object <- rm.na.from.adj.tab(x,adj.var.tab,inds1,inds2,comparison)
				inds1 <- adj.clean.helper.object$inds1
				inds2 <- adj.clean.helper.object$inds2
				adj.var.tab <- adj.clean.helper.object$adj.var.tab
			}
			res.cur <- list(
				comparison=comparison,
				pheno.colname=cc,
				group.names=grps,
				group.inds=list(group1=inds1,group2=inds2),
				paired=is.paired,
				adj.sva=adj.sva,
				adj.celltype=adj.celltype,
				adjustment.table=adj.var.tab,
				region.types=region.types
			)
			res.append <- list(res.cur)
			names(res.append) <- res.cur$comparison
			res <- c(res,res.append)
		}
	}
	return(res)
}

#' rnb.execute.computeDiffMeth
#'
#' computes differential methylation
#' @author Fabian Mueller
#' @aliases rnb.execute.computeDiffMeth
#' @param x RnBSet object
#' @param pheno.cols column names of the pheno slot in \code{x} on which the dataset should be partitioned. Those columns are required to be factors or logical.
#' 				     In case of factors, each group in turn will be compared to all other groups
#' @param region.types which region types should be processed for differential methylation
#' @param covg.thres coverage threshold for computing the summary statistics. See \code{\link{computeDiffTab.extended.site}} for details.
#' @param pheno.cols.all.pairwise integer or character vector specifying the colomns of \code{pheno(x)} on which all pairwise comparisons should be conducted.
#' 		 A value of \code{NULL} (default) indicates no columns.
#' @param columns.pairs argument passed on to \code{rnb.sample.groups}. See its documentation for details.
#' @param columns.adj Column names or indices in the table of phenotypic information to be used for confounder adjustment in the
#'        differential methylation analysis.
#' @param adjust.sva flag indicating whether the adjustment table should also contain surrogate variables (SVs) for the given target variable.
#' @param adjust.celltype flag indicating whether the resulting table should also contain estimated celltype contributions.
#' 				See \code{\link{rnb.execute.ct.estimation}} for details.
#' @param pheno.cols.adjust.sva Column names or indices in the table of phenotypic information to be used for SVA adjustment in the
#'        differential methylation analysis.
#' @param disk.dump Flag indicating whether the resulting differential methylation object should be file backed, ie.e the matrices dumped to disk
#' @param disk.dump.dir disk location for file backing of the resulting differential methylation object. Only meaningful if \code{disk.dump=TRUE}.
#' 						must be a character specifying an NON-EXISTING valid directory.
#' @param ... arguments passed on to binary differential methylation calling. See \code{\link{computeDiffTab.extended.site}} for details.
#' @return an \code{\linkS4class{RnBDiffMeth}} object. See class description for details.
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' get.comparisons(dm)
#' }
rnb.execute.computeDiffMeth <- function(x,pheno.cols,region.types=rnb.region.types.for.analysis(x), covg.thres=rnb.getOption("filtering.coverage.threshold"),
		pheno.cols.all.pairwise=rnb.getOption("differential.comparison.columns.all.pairwise"), columns.pairs=rnb.getOption("columns.pairing"),
		columns.adj=rnb.getOption("covariate.adjustment.columns"),
		adjust.sva=rnb.getOption("differential.adjustment.sva"), pheno.cols.adjust.sva=rnb.getOption("inference.targets.sva"),
		adjust.celltype=rnb.getOption("differential.adjustment.celltype"),
		disk.dump=rnb.getOption("disk.dump.big.matrices"),disk.dump.dir=tempfile(pattern="diffMethTables_"),
		...){

	logger.start("Retrieving comparison info")
	cmp.info <- get.comparison.info(x, pheno.cols=pheno.cols, region.types=region.types, 
		pheno.cols.all.pairwise=pheno.cols.all.pairwise, columns.pairs=columns.pairs, columns.adj=columns.adj,
		adjust.sva=adjust.sva, pheno.cols.adjust.sva=pheno.cols.adjust.sva, adjust.celltype=adjust.celltype)
	logger.completed()
	if (is.null(cmp.info)) {
		return(NULL)
	}

	diff.method <- rnb.getOption("differential.site.test.method")
	dot.args <- list(...)
	if (is.element("diff.method",names(dot.args))){
		diff.method <- dot.args[["diff.method"]]
	}
	logger.start("Computing differential methylation tables")

	diffmeth <- new("RnBDiffMeth",site.test.method=diff.method,disk.dump=disk.dump,disk.path=disk.dump.dir)
	
	for (i in 1:length(cmp.info)){
		cmp.info.cur <- cmp.info[[i]]
		logger.start(c("Comparing:",cmp.info.cur$comparison))
		if (cmp.info.cur$paired){
			logger.status("Conducting PAIRED analysis")
		}

		dm <- computeDiffMeth.bin.site(
				meth(x),inds.g1=cmp.info.cur$group.inds$group1,inds.g2=cmp.info.cur$group.inds$group2,
				covg=covg(x),covg.thres=covg.thres,
				paired=cmp.info.cur$paired, adjustment.table=cmp.info.cur$adjustment.table,
				...
		)
		diffmeth <- addDiffMethTable(diffmeth,dm,cmp.info.cur$comparison,"sites",cmp.info.cur$group.names)
		if (length(cmp.info.cur$region.types)>0){
			dmr <- computeDiffMeth.bin.region(x,dm,
				cmp.info.cur$group.inds$group1,cmp.info.cur$group.inds$group2,
				region.types=cmp.info.cur$region.types
			)			
			for (rt in cmp.info.cur$region.types){
				diffmeth <- addDiffMethTable(diffmeth,dmr[[rt]],cmp.info.cur$comparison, 
					rt, cmp.info.cur$group.names
				)
			}
		}
		logger.completed()
	}

	diffmeth <- addComparisonInfo(diffmeth,cmp.info)
	logger.completed()
	return(diffmeth)
}

