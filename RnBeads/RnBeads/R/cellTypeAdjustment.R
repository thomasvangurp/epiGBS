########################################################################################################################
## cellTypeAdjustment.R
## created: 2014-04-01
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## A set of routines for cell type heterogeneity adjustment.
########################################################################################################################

#' refFreeEWASP
#'
#' Applies the reference-free cell-type heterogeneity adjustment model from [1] and 
#' returns corrected p-values
#' 
#' @param X Matrix on which the test is performed for every row
#' @param inds.g1 column indices of group 1 members
#' @param inds.g2 column indices of group 2 members
#' @param adjustment.table a \code{data.frame} containing variables to adjust for in the testing
#' @param paired should a paired analysis model be used. If so, the first index in \code{inds.g1} must correspond to the first
#' 				 index in  \code{inds.g2} and so on.
#' @param nboot The number of bootstrapping resamples 
#' @param ignore.na in this case all \code{NA} containing rows are removed 
#' @param rescale.residual rescale the residual matrix as z-scores
#' 
#' @return vector of p-values for the "adjusted" regression coefficients from the Reference-free EWAS model
#' 
#' @note Requires the package \pkg{RefFreeEWAS}.
#' 
#' @references 1. Houseman, E. Andres, John Molitor, and Carmen J. Marsit. "Reference-Free Cell Mixture Adjustments in Analysis of DNA Methylation Data." 
#' Bioinformatics (2014): btu029.
#' 
#' @author Pavlo Lutsik
#' @export
refFreeEWASP <- function(
		X,
		inds.g1,
		inds.g2=-inds.g1, 
		adjustment.table=NULL, 
		paired=FALSE, 
		nboot=100,
		ignore.na=TRUE,
		rescale.residual=TRUE) {

	if (!suppressPackageStartupMessages(require(RefFreeEWAS))) {
		stop("missing required package RefFreeEWAS")
	}
	if (paired && packageVersion("RefFreeEWAS")<'1.3'){
		rnb.warning("RefFreeEWAS version >=1.3 is required for paired analysis, pairing will be disregarded")
		paired<-FALSE
	}
	## TODO: Validate parameter values

	rnb.logger.start("Fitting the reference-free EWAS model")

	if (is.logical(inds.g1)) inds.g1 <- which(inds.g1)
	if (is.logical(inds.g2)) inds.g2 <- which(inds.g2)
	n.g1 <- length(inds.g1)
	n.g2 <- length(inds.g2)

	if(ignore.na){
		ncgs<-nrow(X)
		notna.rows<-which(rowSums(is.na(X))==0)
		X<-X[notna.rows,]
		nnas<-length(notna.rows)
	}else{
		nnas<-0
	}
	
	design<-matrix(1L, ncol(X), 2)
	design[inds.g2,2]<-0L
	colnames(design)<-c("(Icept)", "group.f")

	if(!is.null(adjustment.table)){
		design<-cbind(design, adjustment.table)
	}

	tmpBstar <- (X %*% design %*% solve(t(design)%*%design))
	
	## rescaling the residuals, and addition by Andres
	R <- X-tmpBstar %*% t(design)
	if(rescale.residual){
		R <- t(scale(t(R)))
	}
	d<-EstDimRMT(R)$dim
	
	rnb.info(c("Estimated number of latent components is", d))
	
	test <- RefFreeEwasModel(X, design, d)

	rnb.status("Fitted the RefFreeEWAS model")
	
	if(paired){
		pair.id <- rep(1:n.g1, 2)[order(c(inds.g1,inds.g2))]
		testBoot <- PairsBootRefFreeEwasModel(test, nboot, pair.id)
	}else{
		testBoot <- BootRefFreeEwasModel(test,nboot)
	}
	rnb.status("Pefrormed the bootstrap")
	
	smry<-summary(testBoot)
	rnb.status("Summarized the results")

	tstatBeta<-smry[,2,1,1]/smry[,2,1,2]

	# proper df calculation, added by Andres
	pvals <- pt(-abs(tstatBeta), df=nrow(design)-ncol(design)-nnas)
		
	rnb.logger.completed()
	
	if(ignore.na){
		notna.pvals<-pvals
		pvals<-rep(NA,ncgs)
		pvals[notna.rows]<-notna.pvals
	}
	
	return(pvals)
}

#######################################################################################################################

#' estimatePropotionsCP
#' 
#' Estimates cell type proportions using the constrained projection method from Houseman et al. [1]
#' 
#' @param rnb.set			RnBSet object
#' @param cell.type.column	integer index or character identifier of a column
#' 							in the RnBSet object sample annotation table 
#' 							which gives the mapping to reference cell type
#' 							samples
#' @param n.most.variable	singleton integer specifying how many top variable CpGs should be used for marker selection 
#' @param n.markers			singleton integer specifying how many CpGs should
#' 							be used as markers for fitting the projection model
#' @param constrained		if \code{TRUE} the returned cell type proportion estimates
#' 							are non-negative
#' @param full.output 		if \code{TRUE} not only the estimated proportions 
#' 							but also the intermediate analysis results are returned
#'	
#' @return a matrix of estimated cell type contributions (samples times cell types) or a list with results of the intermetidate steps (see details).
#' 					
#' @details The column specified by \code{cell.type.column} should give 
#' assignment of each reference sample to a cell type and missing values for 
#' all the target samples. First the marker selection model is fit to estimate association of each CpG 
#' with the given reference cell types (first expression in eq. (1) of [1]) . The strength of association is expressed as an F-statistic.
#' Since fitting the marker selection model to all CpGs can take a lot of time, one can limit the marker search only to variable CpG positions  
#' by setting \code{n.most.variable} to non-\code{NA} positive integer. The CpGs will be ranked by decreasing across-sample variance in the 
#' reference data set and \code{n.most.variable} will be taken to fit the marker selection model. 
#' Coefficients of the fit together with the F-statistic value for each CpG are returned in case \code{full.output} is \code{TRUE}. 
#' Thereafter, \code{n.markers} are selected as true quantitative markers and the projection model (eq. [2]) is fit to estimate contributions of each cell type.
#' Depending on the value of \code{constrained} the returned coefficients can be either raw or enforced to attain values between 0 and 1 with within-sample sum 
#' less or equal to 1.
#' 
#' @note Requires the package \pkg{nlme}.
#' 
#' @references 1. Houseman, Eugene and Accomando, William and Koestler, Devin and Christensen, Brock and Marsit, Carmen and Nelson, 
#' 	Heather and Wiencke, John and Kelsey, Karl. DNA methylation arrays as surrogate measures of cell mixture distribution. 
#'  BMC Bioinformatics 2012, 13:86
#' 		
#' @author Pavlo Lutsik 
#' @export
estimateProportionsCP<-function(
		rnb.set, 
		cell.type.column, 
		n.most.variable=NA, 
		n.markers=500L, 
		constrained=TRUE, 
		full.output=FALSE){
	
	if (!suppressPackageStartupMessages(require(nlme))) {
		stop("missing required package nlme")
	}
	
	if(!inherits(rnb.set, "RnBSet")){
		stop("invalid value for rnb.set: object of class RnBSet is expected")
	}
	
	if(!is.integer(cell.type.column) && !is.character(cell.type.column)){
		stop("invalid value for cell.type.column: integer or character singleton is expected")		
	}
	
	if(is.integer(cell.type.column) && 
			(length(cell.type.column)!=1 || cell.type.column<0L || cell.type.column>ncol(pheno(rnb.set)))){
		stop("invalid value for cell.type.column: integer index is out of bounds")		
	}	
	
	if(is.character(cell.type.column) && 
			(length(cell.type.column)!=1 || !cell.type.column %in% colnames(pheno(rnb.set)))){
		stop("invalid value for cell.type.column: integer index is out of bounds")		
	}	
	if(full.output){
		result<-list()
	}
	
	if(full.output){
		result$cell.type.column.name=cell.type.column
		result$cell.type.column=pheno(rnb.set)[,cell.type.column]
	}
	
	cell.types<-as.character(unique(na.omit(pheno(rnb.set)[,cell.type.column])))
	
	if(length(cell.types)<2){
		stop("Found less than two reference cell types")
		if(full.output){
			return(result)
		}else{
			return(NULL)
		}
	}
	
	if(length(cell.types)>length(samples(rnb.set))){
		stop("Found too many reference cell types")
		if(full.output){
			return(result)	
		}else{
			return(NULL)
		}
	}
	
	if(full.output){
		result$all.cell.types<-cell.types	
	}
	
	# Omit the last cell type as a referent, to avoid singularities during the model fitting 
	#referent.ct<-cell.types[length(cell.types)]
	#cell.types<-cell.types[-length(cell.types)]
	
#	if(full.output){
#		result$used.cell.types<-cell.types
#		#result$referent<-referent.ct
#	}
	
	# Load test data
	#load("Example-WBC-Data.RData")
	
	# Define validation model
	cell.types.fla<-paste("ct", 1:length(cell.types), sep="_")
	theModel = paste("y", paste(cell.types.fla, collapse="+"), sep="~")
	# remove the intercept term
	theModel = as.formula(paste(theModel, "1", sep="-"))
	
	
	#sizeModel = length(cell.types)+1
	sizeModel = length(cell.types)
	
	mm<-meth(rnb.set, row.names=TRUE)
	if(is.na(n.most.variable)){
		M <- nrow(mm) # Number of CpGs on array (test data is a subset of 27K)
		validationData<-mm[,!is.na(pheno(rnb.set)[,cell.type.column])]
		targetData<-mm[,is.na(pheno(rnb.set)[,cell.type.column])]
	}else{
		if(n.most.variable<2 || n.most.variable>nrow(mm)){
			stop("invalid value for n.most.variable: either less than 2 or 
							more than CpGs in the RnBSet object")
		}
		M<-n.most.variable
		# get ranking by SD only based on the reference-methylome data 
		sds<-apply(mm, 1, function(mmrow) sd(mmrow[!is.na(pheno(rnb.set)[,cell.type.column])]))
		cands<-order(sds, decreasing=TRUE)[1:n.most.variable]
		validationData<-mm[cands,!is.na(pheno(rnb.set)[,cell.type.column])]
		targetData<-mm[cands,is.na(pheno(rnb.set)[,cell.type.column])]		
	}
	rm(mm)
	# select most variable CpGs 
	
	cell.type.f<-as.character(pheno(rnb.set)[!is.na(pheno(rnb.set)[,cell.type.column]),cell.type.column])
	cell.type.contrasts<-sapply(cell.types, function(ct) as.integer(cell.type.f==ct))
	colnames(cell.type.contrasts)<-cell.types.fla
	plate.col<-grep("Sentrix[| |_]*ID", colnames(pheno(rnb.set)))
	
	if(length(plate.col)==1){
		validationData_Pheno<-data.frame(cell.type.contrasts)
		validationData_Pheno[["PLATE"]]<-pheno(rnb.set)[!is.na(pheno(rnb.set)[,cell.type.column]),plate.col]
	}else{
		validationData_Pheno<-data.frame(cell.type.contrasts)
	}
	
	#################### Detection of the marker set and estimation of the coefficients
	
	# Linear transformation of coefficient vector
	#  representing contrast to test F statistic
	#L.forFstat <- diag(sizeModel)[-1,,drop=F]  #All non-intercept coefficients
	L.forFstat <- diag(sizeModel)[,drop=F]
	
	# Initialize various containers
	sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
	coefEsts <- matrix(NA, M, sizeModel)
	coefVcovs <- list()
	
	for(j in 1:M){ # For each CpG
		
		#Remove missing methylation values
		ii <- !is.na(validationData[j,])
		nObserved[j] <- sum(ii)
		validationData_Pheno[["y"]] <- validationData[j,]
		
		#if(j%%10==0) cat(j,"\n") # Report progress
		
		try({ # Try to fit a mixed model to adjust for plate
					if(length(plate.col)==1){
						fit <- try(lme(theModel, random=~1|PLATE, data=validationData_Pheno[ii,]), silent=TRUE)
					}
					
					if(length(plate.col)!=1 || inherits(fit,"try-error")){ # If LME can't be fit, just use OLS
						fit <- lm(theModel, data=validationData_Pheno[ii,])
						fitCoef <- fit$coef
						sigmaResid[j] <- summary(fit)$sigma
						sigmaIcept[j] <- 0L
						nClusters[j] <- 0L
					}
					else{ 
						fitCoef <- fit$coef$fixed
						sigmaResid[j] <- fit$sigma
						sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
						if(length(plate.col)==1){
							nClusters[j] <- length(fit$coef$random[[1]])
						}else{
							nClusters<-0L
						}
					}
					coefEsts[j,] <- fitCoef
					coefVcovs[[j]] <- vcov(fit)
					
#					useCoef <- L.forFstat[!is.na(fitCoef)[-1],!is.na(fitCoef),drop=F] %*% fitCoef[!is.na(fitCoef)]
					useCoef <- L.forFstat[!is.na(fitCoef),!is.na(fitCoef),drop=F] %*% fitCoef[!is.na(fitCoef)]
#					useV <- L.forFstat[!is.na(fitCoef)[-1],!is.na(fitCoef),drop=F] %*% coefVcovs[[j]] %*% t(L.forFstat[!is.na(fitCoef)[-1],!is.na(fitCoef),drop=F])
					useV <- L.forFstat[!is.na(fitCoef),!is.na(fitCoef),drop=F] %*% coefVcovs[[j]] %*% t(L.forFstat[!is.na(fitCoef),!is.na(fitCoef),drop=F])
					Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
				})
	}
	
	# Name the rows so that they can be easily matched to the target data set
	rownames(coefEsts) <- rownames(validationData)
	#colnames(coefEsts) <- c(names(fitCoef)[1], cell.types)
	colnames(coefEsts) <- c(cell.types)
	
	# Get P values corresponding to F statistics
	Pval <- 1-pf(Fstat, sizeModel, nObserved - nClusters - sizeModel + 1)
	
	if(full.output){
		result$f.stat<-Fstat
		result$f.pval<-Pval
		result$coef.ests<-coefEsts
	}
	
	## order by F statistic since p values become indiscriminatory
	#
	inds<-order(Fstat, decreasing=TRUE)[1:n.markers]
	
	if(full.output){
		result$markers<-inds
		result$n.markers<-n.markers		
	}
	
	# NOTE:  test data consists of CpGs having 500 largest F statistics
	#  in descending order of magnitude
	# For the full 27K array, it is necessary to sort by Pvalue and/or Fstatistic
	#  and choose the top CpGs
	#table(sign(diff(Fstat))) # Fstats decrease as j increases from 1 to M
	
	
	##################### Estimation of the cell type contributions 
	
#	valid.cts<-(colSums(is.na(coefEsts))==0)[-1L]
#	
#	if(full.output){
#		result$valid.cell.types<-valid.cts
#	}
	
	Lwbc <- diag(sizeModel) 
	#Lwbc <- diag(sizeModel)[-1L,,drop=F]
	#Lwbc[,1] <- 1
	rownames(Lwbc) <- colnames(coefEsts)
	colnames(Lwbc) <- colnames(coefEsts)
	
	#Lwbc # View contrast matrix
	
	# CpGSelection = rownames(coefEsts)[inds] 
	# Use the top 100
	# Note:  if the CpGs were scattered throughout the array,
	#    you would want to select them by name as is indicated here.
	#    For this sample version, it would be easier just to use
	#    "[1:100]"
	
	
	####### Projections 
	
	if(full.output){
		
		coefsList<-list()
		
		for(constri in 1:2){
			coefsList[[constri]]<- projectWBC(
					targetData[inds,],
					coefEsts[inds,,drop=FALSE],    
					Lwbc, 
					nonnegative=c(FALSE, TRUE)[constri])
		}
		
		result$contributions<-coefsList[[1]]
		result$contributions.nonneg<-coefsList[[2]]
		
	}else{
		
		coefs <-  projectWBC(
				targetData[inds,],
				coefEsts[inds,,drop=FALSE],    
				Lwbc, nonnegative = constrained)
		
#		if(ncol(coefs)<length(cell.types) && constrained && normalized){
#			
#			pre.coefs<-coefs
#			coefs<-matrix(NA,nrow(pre.coefs),length(cell.types))
#			coefs[,valid.cts]<-pre.coefs/max(rowSums(pre.coefs))
#			coefs[,!valid.cts]<-(1-rowSums(pre.coefs)/(max(rowSums(pre.coefs))))/length(which(!valid.cts))
#			colnames(coefs)<-cell.types
#			rownames(coefs)<-rownames(pre.coefs)
#		}
		
	}
	
	if(full.output){
		return(result)
	}
	
	return(coefs)
}

#######################################################################################################################

#' rnb.plot.marker.fstat
#'
#' Plot the the cell type marker selection based on the reference methylome data  
#' 
#' @param ct.object   Object of class \code{CellTypeInferenceResult} as returned by \link{rnb.execute.ct.estimation}.
#' @param writeToFile If \code{TRUE}, the plot will be written to a file.
#' @param ...         Other arguments to \code{\link{createReportPlot}}.
#' 
#' @details	The F-statistic values from the cell type association model (first part of eqn. (1) in [1]) are plotted in decreasing order 
#' 			for all tested CpG positions. A vertical line gives a cut-off for the number of selected cell type markers.
#' 
#' @return				if \code{writeToFile=TRUE} an object of class \code{\linkS4class{ReportPlot}}, 
#' 						and the plotted reordered F-statistics vector otherwise  
#' 
#' @references  1. Houseman, Eugene and Accomando, William and Koestler, Devin and Christensen, Brock and Marsit, Carmen and Nelson, 
#' 	Heather and Wiencke, John and Kelsey, Karl. DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics 2012, 13:86
#' 
#' @author Pavlo Lutsik
#' @export
rnb.plot.marker.fstat <- function(ct.object, writeToFile=FALSE, ...) {
	
	if(class(ct.object)!="CellTypeInferenceResult"){
		stop("invalid value for ct.object; expected CellTypeInferenceResult")
	}
	if(!isTRUE(ct.object$method == "houseman1")){
		stop("unsupported method in ct.object")
	}
	if (!parameter.is.flag(writeToFile)) {
		stop("invalid value for writeToFile; expected TRUE or FALSE")
	}

	if(writeToFile){
		plot.file<-createReportPlot("CellTypeMarkerFstatPlot", ...)
	}

	dframe <- data.frame(x = 1:length(ct.object$f.stat), y = sort(ct.object$f.stat, decreasing = TRUE))
	print(ggplot2::ggplot(dframe, aes(x = x, y = y)) + ggplot2::geom_point() +
		ggplot2::geom_vline(xintercept = ct.object$n.markers, col = muted("blue")) +
		ggplot2::labs(x = "CpGs", y = "F statistic"))

	if(writeToFile) {
		plot.file <- off(plot.file)
		return(plot.file)
	}
	return(invisible(dframe$y))
}

#######################################################################################################################

#' rnb.plot.ct.heatmap
#'
#' Plot contributions of the cell types  
#' 
#' @param ct.obj		Object of class \code{CellTypeInferenceResult} as returned by \link{rnb.execute.ct.estimation}.
#' @param type			Type of cell type contributions to plot.
#' @param writeToFile	If \code{TRUE}, the plot will be written to a file.
#' @param ...			Other arguments passed to \code{\link{createReportPlot}}.
#' 
#' @details				The cell type contributions are visualized as a heatmap
#' 
#' @return				if \code{writeToFile=TRUE} an object of class \code{\linkS4class{ReportPlot}}, 
#' 						or the protted matrix otherwise   
#' 
#' @author Pavlo Lutsik
#' @export
rnb.plot.ct.heatmap<-function(ct.obj, type="nonnegative", writeToFile=FALSE, ...) {
	
	if(class(ct.obj)!="CellTypeInferenceResult"){
		stop("Invalid value for ct.obj")
	}
	
	if(ct.obj$method!="houseman1"){
		stop("This plot requires the ct.obj obtained usign the houseman1 method")
	}
	
	
	if(type=="raw"){
		tbl<-ct.obj$contributions
	}else if(type=="nonnegative"){
		tbl<-ct.obj$contributions.nonneg
	}
	
	if(writeToFile){
		fun.args <- list(...)
		if (!("fname" %in% names(fun.args))) {
			fun.args$fname <- "CellTypeContributionsPlot"
		}
		if (!("height" %in% names(fun.args))) {
			fun.args$height <- 1.2 + nrow(tbl)* 0.15
		}
		if (!("width" %in% names(fun.args))) {
			fun.args$width <- 3.2 + ncol(tbl) * 0.6
		}
		plot.file<-base::do.call(createReportPlot, fun.args)
	}

	dframe <- data.frame(
		x = factor(rep(colnames(tbl), each = nrow(tbl)), levels = colnames(tbl)),
		y = factor(rep(rownames(tbl), ncol(tbl)), levels = rev(rownames(tbl))),
		v = as.vector(tbl))
	pp <- ggplot(dframe, aes(x, y, fill = v)) +
		geom_tile(colour = "white") +
		scale_fill_gradient(low ="white", high = "steelblue") +
		labs(x = NULL, y = NULL, fill = "Contribution") +
		scale_x_discrete(expand = c(0, 0)) +
		scale_y_discrete(expand = c(0, 0)) +
		theme(panel.background = element_blank(), panel.grid = element_blank()) +
		theme(panel.border = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
		theme(axis.text.x = element_text(size = 8, angle = 330, hjust = 0)) +
		theme(legend.position = c(1, 0.5), legend.justification = c(0, 0.5)) +
		theme(plot.margin = unit(0.1 + c(0, 1, 0, 0), "in"))

	if(writeToFile) {
		pp <- suppressWarnings(ggplot_gtable(ggplot_build(pp)))
		pp$widths[[3]] <- unit(2, "in")
		pp$heights[[length(pp$heights) - 2L]] <- unit(1, "in")
		grid.draw(pp)
		return(off(plot.file))
	}
	return(pp)
}

#######################################################################################################################
 
#' rnb.execute.ct.estimation
#'
#' Perform the estimation of the cell type contributions in each analyzed sample.
#' 
#' @param rnb.set			object of class \code{\linkS4class{RnBSet}}
#' @param cell.type.column	integer index or character identifier of a column in sample annotation table of \code{rnb.set}
#' 							which gives the mapping of samples to reference cell types 
#' @param test.max.markers	maximal amount of CpG positions to use for marker selection
#' @param top.markers		the number of markers to select 
#' @param method 			algorithm used for estmation of the cell type contributions
#' @param verbose			flag specifying whether diagnostic output should be written to the console 
#' 							or to the RnBeads logger in case the latter is initialized
#' 
#' @details The only supported method is the one from Houseman et al BMC Bioinformatics 2012 
#' 
#' @return object of class \code{CellTypeInferenceResult} 
#' 
#' @author Pavlo Lutsik
#' @export
rnb.execute.ct.estimation<-function(rnb.set, 
		cell.type.column=NA,
		test.max.markers=NA,
		top.markers=500,
		method="houseman1",
		verbose=TRUE){
		
	if(method=="houseman1"){
		result<-estimateProportionsCP(rnb.set, 
				cell.type.column, 
				n.most.variable=test.max.markers,
				n.markers=top.markers, full.output=TRUE)
	}else{
		stop("the supplied inference method is not supported yet")
	}
	
	result$method <- method
	class(result)<-"CellTypeInferenceResult"
	return(result)
}

#######################################################################################################################
#
# rnb.section.ct.estimation
#
#
rnb.section.ct.estimation<-function(report, ct.object){
		
	report <- rnb.add.section(report, "Estimation of Cell Type Heterogeneity Effects", NULL)
	
	if(!is.null(ct.object$cell.type.column.name)){
		
		report <- rnb.add.section(report, "Reference methylomes", NULL, level=2)
		
		intro <- c("The dataset contained reference methylomes defined by the sample annotation column <i>",
			ct.object$cell.type.column.name, "</i>.")

		if(length(which(!is.na(ct.object$cell.type.column)))>2){
			intro <- c(intro, " Detected reference methylomes were are summarized in the table below.")
		}
		ref.meth.table <- data.frame(table(na.omit(ct.object$cell.type.column)))
		colnames(ref.meth.table) <- c("Cell type", "Number of samples")

		rnb.add.paragraph(report, intro)

		rnb.add.table(report, ref.meth.table, row.names = FALSE, first.col.header = TRUE)
		
	}
			
	# plot the marker distribution
	if(ct.object$method=="houseman1"){
		
		report <- rnb.add.section(report, "Cell type contributions", NULL, level=2)

		refText <- c("Houseman, E. (2012) DNA methylation arrays as surrogate measures of cell mixture distribution. ",
			"<i>BMC Bioinformatics</i>, <b>13</b>(86)")
		report <- rnb.add.reference(report, refText)

		rnb.add.paragraph(report, c("The contributions of cell types were estimated using the method by Houseman <i>",
			"et al</i> ", rnb.get.reference(report, refText), "."))

		markers.text <- c("In the first step the reference methylomes were used to estimate the association of each ",
			"CpG position to each of the cell types. The stength of association was measured using F-test. ")
		ctMarkers <- rnb.getOption("inference.max.cell.type.markers")
		if(!is.na(ctMarkers)){
			markers.text <- c(markers.text, "To decrease the computation load, only ", ctMarkers, " most variable ",
				"CpGs were considered. ")
		}
		n.markers <- ct.object$n.markers # rnb.getOption("inference.top.cell.type.markers")
		if (n.markers < length(ct.object$f.stat)) {
			borders <- sort(ct.object$f.stat, decreasing = TRUE)[0:1 + n.markers]
			i <- 10^(-1:4)
			cutoff <- floor(borders[1] * i) / i
			i <- which(cutoff > borders[2])
			cutoff <- base::ifelse(length(i) == 0, paste("approximately", cutoff[length(cutoff)]), cutoff[i[1]])
			rm(borders, i)
		} else { # n.markers == length(ct.object$f.stat)
			cutoff <- floor(min(ct.object$f.stat))
		}
		
		markers.text <- c(markers.text, "Finally, only ", n.markers, " CpGs ",
			"with the lowest F-test p-value were used in the contribution estimation. The plot below visualizes the ",
			"distribution of F statistic values for all tested CpGs. Note that selecting the most informative CpGs is ",
			"equivalent to applying an F statistic cut-off of ", cutoff, ".")
		report <- rnb.add.section(report, "Selection of the cell type markers", markers.text, level = 3)

		## Add a plot of F statistic values
		fstat.plot <- rnb.plot.marker.fstat(ct.object, writeToFile = TRUE, report = report)
		txt <- c("Scatter plot visualizing the F statistic of the cell type association model for each CpG position ",
			"from the tested subset. The vertical blue line, if present, reflects the selection of ", n.markers,
			" best markers for the projection.")
		report <- rnb.add.figure(report, txt, fstat.plot)
		rm(n.markers, cutoff)

		## Export the normalized (non-negative) contributions
		fname <- 'contributions_houseman1.csv'
		write.csv(ct.object$contributions.nonneg, file = file.path(rnb.get.directory(report, 'data', TRUE), fname))
		txt <- c("After the marker selection, a projection of the target data onto the space of the marker selection ",
			"model coefficients yields contributions of each reference cell type to each measured DNA methylation ",
			'profile. The resulting cell type contributions are available in a dedicated <a href="',
			rnb.get.directory(report, 'data'), '/', fname, '">comma-separated file</a> accompanying this report. These ',
			'values are also displayed in the heatmap below.', "The contributions are constrained to be greater or equal ",
			"to zero, and the per-sample sums are expected to be close to one, i.e. they are estimates of the cell type ",
			"proportions. Per-sample totals much larger than one may indicate the problems with the procedure, ",
			"e.g. bad correspondence of the target data to the reference methylomes significant batch effects etc.")
	
		report <- rnb.add.section(report, "Cell type contributions via the coefficient projection", txt, level = 3)

		## Add a heatmap of normalized (non-negative) contributions
		contrib.plot <- rnb.plot.ct.heatmap(ct.object, writeToFile = TRUE, report = report)		
		txt <- "Heatmap visualizing estimated cell type contributions, scaled to the range [0, 1]."
		report <- rnb.add.figure(report, txt, contrib.plot)
	}
	
	return(report)
}

#######################################################################################################################
#
# rnb.step.cell.types
#
#
rnb.step.cell.types<-function(rnb.set, report){
	
	logger.start("Estimation of the cell type heterogeneity effects")
	
	logger.start("Performing computations")
		result<-rnb.execute.ct.estimation(rnb.set, 
				cell.type.column=rnb.getOption("inference.reference.methylome.column"),
				test.max.markers=rnb.getOption("inference.max.cell.type.markers"), 
				top.markers=rnb.getOption("inference.top.cell.type.markers"),
				method="houseman1")
	logger.completed()
	
	logger.start("Adding a section to the report")
		new.report<-rnb.section.ct.estimation(report, result)
	logger.completed()
	
	logger.completed()
	
	if(!is.null(result)){
		rnb.set<-set.covariates.ct(rnb.set, result)
		logger.info("Added cell type covariates to the RnBSet object")
	}
	
	return(list(rnb.set=rnb.set, report=new.report))
}

#######################################################################################################################

#' set.covariates.ct
#'
#' Adds the results of cell type estimation to an RnBSet
#' 
#' @param rnb.set The \code{RnBSet} object to which the results should be added
#' @param ct.obj An object of class \code{CellTypeInferenceResult} returned by \code{rnb.execute.ct.estimation}.
#' 
#' @return The modified \code{RnBSet}.
#' 
#' @export set.covariates.ct
set.covariates.ct <- function(rnb.set, ct.obj) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!.hasSlot(rnb.set,"inferred.covariates")) { #.hasSlot ensure backwards compatibility
		rnb.set@inferred.covariates<-list()
	}
	if(is.null(rnb.set@inferred.covariates)){
		rnb.set@inferred.covariates<-list()
	}

	if (!is.null(ct.obj)) {
		## FIXME: Validate ct.obj and its elements
		contributions <- matrix(as.double(NA), nrow = length(samples(rnb.set)), ncol = ncol(ct.obj$contributions),
			dimnames = list(colnames(rnb.set@meth.sites), paste0("ct_", 1:ncol(ct.obj$contributions)))) #colnames(ct.obj$contributions))) #paste0("ct_", ...)
		i.samples <- which(is.na(pheno(rnb.set)[, ct.obj$cell.type.column.name]))
		contributions[i.samples, ] <- ct.obj$contributions
		attr(contributions, "column") <- ct.obj$cell.type.column.name
		rnb.set@inferred.covariates[["cell.types"]] <- contributions
	}

	return(rnb.set)
}

#######################################################################################################################

#' get.covariates.ct
#'
#' Retrieves an NxK matrix of cell type contributions stored in an RnBSet for a given target variable
#' 
#' @param rnb.set \code{RnBSet} object
#' 
#' @return an NxK matrix of K cell types contributions for N samples of the \code{rnb.set}. \code{NULL}
#'		   if the components have not been computed or added to \code{rnb.set}.
#' 
#' @export
get.covariates.ct <- function(rnb.set) {
	if (!.hasSlot(rnb.set,"inferred.covariates")) { #.hasSlot ensure backwards compatibility
		return(NULL)
	}
	rnb.set@inferred.covariates[["cell.types"]]
}

#######################################################################################################################
#' has.covariates.ct
#'
#' Checks whether the given \code{RnBSet} object contains cell type contribution estimates
#' 
#' @param rnb.set \code{RnBSet} object
#' 
#' @return \code{TRUE} if the supplied object contains the cell type covariates information and \code{FALSE} otherwise
#' 
#' @export
has.covariates.ct <- function(rnb.set) {
	if (!.hasSlot(rnb.set,"inferred.covariates")) { #.hasSlot ensure backwards compatibility
		return(FALSE)
	}
	return(!is.null(rnb.set@inferred.covariates[["cell.types"]]))
}
#######################################################################################################################
