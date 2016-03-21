########################################################################################################################
## controlPlotsBiSeq.R
## created: 2012-04-04
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Quality control plots for the reduced-representation and the whole-genome bisulfite sequencing experiments.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.plot.biseq.coverage
#'
#' Plots the sequencing coverage of the RnBiseqSet object across the genomic coordinate 
#'   
#' @param rnbs.set			RnBiseqSet object
#' @param type				\code{character} singleton. If \code{site} the coverage
#' 							information is plotted for each methylation site. Otherwise
#' 							should be one of the regions returned by \code{rnb.region.types}
#' @param sample			unique sample identifier. In case \code{rnb.getOption("identifiers.column")} 
#' 							is not \code{NULL}, \code{sample} should attain values from the corresponding column, 
#' 							or \code{colnames(meth(rnb.set))} otherwise
#' @param writeToFile 		flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param numeric.names if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file name with digits
#' @param covg.lists		if available, the output of \code{\link{rnb.execute.quality}}
#' @param ... 				other arguments to \code{\link{createReportPlot}}
#' 
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class 
#' 			\code{\link{ggplot}} otherwise.
#' TODO add examples
#' @author Pavlo Lutsik
#' @export
rnb.plot.biseq.coverage<-function(
		rnbs.set, 
		sample, 
		type = "sites", 
		writeToFile = FALSE, 
		numeric.names = FALSE, 
		covg.lists=NULL, 
		...){
	
	if(!inherits(rnbs.set, "RnBiseqSet"))
		stop("Invalid value for rnbs.set: object of class RnBiseqSet expected")
	
	if(!(is.character(type) && type %in% c("sites", rnb.region.types(assembly(rnbs.set)))))
		stop("Invalid value for type: \"sites\" or one of the rnb.region.types() expected")
	
	if (!(is.character(sample) && sample %in% samples(rnbs.set))) {
		stop("Invalid value for sample: character expected")
	}
	assembly <- assembly(rnbs.set)	
	covg.rnbs<-covg(rnbs.set, type)
	covg.rnbs[is.na(covg.rnbs)]<-0
	
		
	if(type=="sites"){
		type<-"CpG"
		map<-sites(rnbs.set)
	}else{
		map<-regions(rnbs.set, type)
	}
	
	sample.id<-match(sample, samples(rnbs.set))
	
	if(max(covg.rnbs[,sample.id])==0){
		
		if(writeToFile) plot.file<-createReportPlot(paste("CoveragePlot", ifelse(numeric.names, sample.id, sample), sep="_"), ...)
			print(rnb.message.plot("All CpG have zero coverage in this sample"))
		if(writeToFile) {
			off(plot.file)
			return(plot.file)
		}else{
			return(invisible(NULL))
		}
	}
	
	if(is.null(covg.lists)){
				
		covg.lists<-lapply(names(rnb.get.chromosomes(assembly=assembly)), function(chr) {
					chr.id<-match(chr, names(rnb.get.chromosomes(assembly=assembly)))
					covg.weight<-rep(0, rnb.annotation.size(assembly=assembly)[chr.id])
					covg.weight[map[map[,2]==chr.id,3]]<-covg.rnbs[map[,2]==chr.id,sample.id]
					covg.weight
				})
	}else if(is.list(covg.lists)){
		
		covg.lists<-lapply(1:length(covg.lists), function(chr.id){
					covg.weight<-covg.lists[[chr.id]]
					covg.weight[map[map[,2]==chr.id,3]]<-covg.rnbs[map[,2]==chr.id,sample.id]
					covg.weight
				})
	}
	names(covg.lists)<-names(rnb.get.chromosomes(assembly=assembly))
	
	cv.templ<-coverage(rnb.get.annotation(type,assembly=assembly), weight=covg.lists)
	
	if(writeToFile) plot.file<-createReportPlot(paste("CoveragePlot", ifelse(numeric.names, sample.id, sample), sep="_"), ...)
		
	names(cv.templ)<-NULL
	plot<-ggbio::autoplot(cv.templ, type = "viewMaxs", stat = "slice", lower = 8, ylab="", ylim=2000, labeller=function(variable,value){names(rnb.get.chromosomes(assembly=assembly))[value]})+
					scale_y_continuous(limits=c(0,2000))
		
	if(writeToFile) {
		print(plot)
		off(plot.file)
		return(plot.file)
	}else{
		return(plot)
	}
}

########################################################################################################################

#' rnb.plot.biseq.coverage.hist
#'
#' Plots the histograms of the coverage 
#'   
#' @param rnbs.set			RnBiseqSet object
#' @param sample			unique sample identifier. In case \code{rnb.getOption("identifiers.column")} 
#' 							is not \code{NULL}, \code{sample} should attain values from the corresponding column, 
#' 							or \code{colnames(meth(rnb.set))} otherwise
#' @param type				\code{character} singleton. If \code{site} the coverage
#' 							information is plotted for each methylation site. Otherwise
#' 							should be one of the regions returned by \code{rnb.region.types}
#' @param writeToFile 		a flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param numeric.names if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file name with digits
#' @param covg.max.percentile the maximum percentile of the coverage to be plotted
#' @param ... 				other arguments to \code{\link{createReportPlot}}
#' 
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class 
#' 			\code{\link{ggplot}} otherwise.
#'
#' TODO add examples
#' @author Pavlo Lutsik
#' @export
rnb.plot.biseq.coverage.hist<-function(
		rnbs.set, 
		sample, 
		type = "sites", 
		writeToFile = FALSE, 
		numeric.names = FALSE, 
		covg.max.percentile=1.0, 
		...){
	
	if(!inherits(rnbs.set, "RnBiseqSet"))
		stop("Invalid value for rnbs.set: object of class RnBiseqSet expected")
	
	if(!(is.character(type) && type %in% c("sites", rnb.region.types(assembly(rnbs.set)))))
		stop("Invalid value for type: \"sites\" or one of the rnb.region.types() expected")
	
	if (!(is.character(sample) && sample %in% samples(rnbs.set))) {
		stop("Invalid value for sample: character expected")
	}
	if (!is.double(covg.max.percentile) ||  covg.max.percentile > 1 || covg.max.percentile < 0){
		stop("Invalid value for covg.max.percentile: number between 0 and 1 expected")
	}
	
	sample.id <- match(sample, samples(rnbs.set))
	covg.rnbs <- covg(rnbs.set, type)
	max.covg  <- max(covg.rnbs, na.rm=TRUE)
	if (covg.max.percentile > 0 && covg.max.percentile < 1){
		max.covg <- quantile(covg.rnbs, probs=c(covg.max.percentile), na.rm=TRUE)[1]
	}
	covg.rnbs<-covg.rnbs[,sample.id]
	covg.rnbs[is.na(covg.rnbs)]<-0
	
	
	if(writeToFile) {
		plot.file<-createReportPlot(paste("CoverageHistogram", ifelse(numeric.names, sample.id, sample), sep="_"), ...)
	}
	
	dframe<-data.frame(Coverage=covg.rnbs)
	#TODO: cut after 95 percentile of coverage
	pp <- ggplot(dframe, aes_string(x = "Coverage")) + coord_cartesian(xlim=c(0,max.covg))+ 
		labs(x = "Coverage", y = "# CpGs") + geom_histogram(aes_string(y = "..count.."), binwidth = 10)
	
	if(writeToFile) {
		print(pp)
		off(plot.file)
		return(plot.file)
	}else{
		return(pp)
	}
	
}
########################################################################################################################

#' rnb.plot.biseq.coverage.violin
#'
#' Plots the violin plots of the coverage distribution
#'   
#' @param rnbs.set			RnBiseqSet object
#' @param samples			unique sample identifiers. In case \code{rnb.getOption("identifiers.column")} 
#' 							is not \code{NULL}, \code{samples} should attain values from the corresponding column, 
#' 							or \code{colnames(meth(rnb.set))} otherwise
#' @param fname				base filename for the files to be plotted. If NULL, the plot will not be written to file
#' @param type				\code{character} singleton. If \code{site} the coverage
#' 							information is plotted for each methylation site. Otherwise
#' 							should be one of the regions returned by \code{rnb.region.types}
#' @param covg.range	    Vector of length 2 specifying the range of coverage to be plotted. if \code{NULL} (default) the entire range will be plotted
#' @param ... 				other arguments to \code{\link{createReportPlot}}
#' 
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class 
#' 			\code{\link{ggplot}} otherwise.
#'
#' @author Fabian Mueller
#' TODO add examples
#' @export
rnb.plot.biseq.coverage.violin <- function(rnbs.set, samples, fname=NULL, type = "sites", covg.range=NULL, ...){
	
	if(!inherits(rnbs.set, "RnBiseqSet"))
		stop("Invalid value for rnbs.set: object of class RnBiseqSet expected")
	
	if(!(is.character(type) && type %in% c("sites", rnb.region.types(assembly(rnbs.set)))))
		stop("Invalid value for type: \"sites\" or one of the rnb.region.types() expected")
	
	if (!(is.character(samples) && all(samples %in% samples(rnbs.set)))) {
		stop("Invalid value for samples: character expected")
	}
	if (!is.null(covg.range) && length(covg.range)!=2){
		stop("Invalid value for covg.range. vector of length 2 or NULL expected")
	}
	
	sample.id <- match(samples, samples(rnbs.set))
	covg.rnbs <- covg(rnbs.set, type)
	max.covg  <- max(covg.rnbs, na.rm=TRUE)
	covg.rnbs <- covg.rnbs[,sample.id]
	covg.rnbs[is.na(covg.rnbs)]<-0
	
	dframe <- data.frame(Coverage=as.vector(covg.rnbs),
			Sample=rep(samples,times=rep(nrow(covg.rnbs),length(samples))))
	
	if(!is.null(fname)) plot.file <- createReportPlot(fname, ...)
	
	pp <- ggplot(dframe, aes_string(x = "Sample", y = "Coverage")) 
	if (!is.null(covg.range)) {
		pp <- pp + ylim(covg.range[1],covg.range[2])
	}
	pp <- pp + labs(x = "Sample", y = "Coverage") +
		geom_violin(aes_string(fill="Sample"),scale="count") +
		scale_fill_manual(values=rep(rnb.getOption("colors.category"),length.out=length(samples))) +
		coord_flip() + theme(legend.position="none")
	
	if(!is.null(fname)) {
		print(pp)
		off(plot.file)
		return(plot.file)
	}else{
		return(pp)
	}
}

########################################################################################################################

#' rnb.plot.coverage.thresholds
#'
#' Plots the number of remaining CpGs after applying different thresholds for coverage and support.
#'
#' @param rnb.set       Methylation dataset as an object of type \code{\linkS4class{RnBiseqSet}}.
#' @param min.coverages Non-empty \code{integer} vector storing the unique positive cutoff values to be applied for
#'                      minimal coverage. Names, if present, are interpreted as colors that must be used to denote the
#'                      corresponding values.
#' @param fname         File name to save the generated plot to. See the \emph{Details} section for restrictions.
#' @param ...           Additional named parameters related to saving the plot to files. These can include:
#'                      \code{report}, \code{width}, \code{height}, \code{create.pdf}, \code{low.png} and
#'                      \code{high.png}. These parameters are ignored when \code{fname} is \code{NULL} or \code{NA}.
#' @return If \code{fname} is \code{NULL} or \code{NA} (default), the generated plot as an object of type
#'         \code{ggplot2}; otherwise, the initialized and closed \code{\linkS4class{ReportPlot}} object, invisibly.
#' 
#' @details
#' If \code{fname} is specified, this function calls \code{\link{createReportPlot}} to save the plot to PDF and/or PNG
#' files. See \link[=createReportPlot]{its documentation} for information on acceptable file names. Additional
#' parameters - \code{report}, \code{width}, \code{height}, etc. - can also be given. If image width is not specified,
#' it is set to a value between 4.7 and 9.2 (inches), depending on the number of samples in the dataset. The default
#' image height is fixed to 7.2.
#' 
#' @author Yassen Assenov
#' @export
rnb.plot.coverage.thresholds <- function(rnb.set, min.coverages, fname = NA, ...) {

	## Validate parameter values
	if (!inherits(rnb.set, "RnBiseqSet")) {
		stop("invalid value for rnb.set")
	}
	if (is.double(min.coverages) && isTRUE(all(min.coverages == as.integer(min.coverages)))) {
		cnames <- names(min.coverages)
		min.coverages <- as.integer(min.coverages)
		names(min.coverages) <- cnames
		rm(cnames)
	}
	if (!(is.integer(min.coverages) && length(min.coverages) != 0 && all(!is.na(min.coverages))
		  && all(min.coverages > 0))) {
		stop("invalid value for min.coverages")
	}
	if (anyDuplicated(min.coverages)) {
		stop("invalid value for min.coverage; duplicates found")
	}
	min.coverages <- sort(min.coverages)
	if (is.null(fname)) {
		fname <- NA
	} else if (!is.na(fname)) {
		if (!is.character(fname)) {
			stop("invalid value for fname")
		}
	}

	## Define colors to denote min coverage thresholds
	if (is.null(names(min.coverages))) {
		n <- length(min.coverages)
		colors.covered <- rgb(seq(0, 0.5, length.out = n), seq(0, 0.5, length.out = n), seq(0, 1, length.out = n))
		names(colors.covered) <- min.coverages
	} else {
		colors.covered <- names(min.coverages)
		colors.covered[is.na(colors.covered)] <- "#808080"
		names(colors.covered) <- min.coverages
	}

	## Construct data.frame storing the plot information
	matrix.coverages <- covg(rnb.set)
	n <- ncol(matrix.coverages)
	if (n <= 50) {
		s.values <- n:1
		s.indices <- 1:n
		xlabel <- "Support (number of samples)"
	} else {
		s.values <- seq(0, 1, length.out = 51)[-1]
		s.indices <- as.integer(ceiling(quantile(0:n, probs = s.values)))
		s.indices <- n + 1L - s.indices
		xlabel <- "Support (percentile)"
	}
	dframe <- do.call(rbind, lapply(as.integer(min.coverages), function(mc) {
			interrogated <- as.integer(rowSums(matrix.coverages >= mc))
			cpgs <- cumsum(rev(tabulate(interrogated, ncol(matrix.coverages)))) / 10^6
			data.frame("mc" = mc, "s" = s.values, cpgs = cpgs[s.indices])
		}
	))
	dframe$mc <- factor(as.character(dframe$mc), levels = as.character(min.coverages))

	## Generate the plot
	if (n <= 20) {
		img.width <- 4.7
		xbreaks <- 1:n
		xlabels <- as.character(xbreaks)
		xlimits <- 0.5 + c(n, 0)
	} else if (n <= 50) {
		img.width <- 1.7 + 0.15 * n
		xbreaks <- seq(2L, n, by = 2L)
		xlabels <- as.character(xbreaks)
		xlimits <- 0.5 + c(n, 0)
	} else {
		img.width <- 1.7 + 0.15 * 50
		xbreaks <- seq(0, 1, length.out = 26)
		xlabels <- as.character(xbreaks * 100)
		xlimits <- c(1.01, 0.01)
	}
	pp <- ggplot2::ggplot(dframe, aes_string(x = "s", y = "cpgs", color = "mc", group = "mc")) + ggplot2::geom_line() +
		ggplot2::geom_point() + ggplot2::scale_color_manual(values = colors.covered) +
		ggplot2::scale_x_reverse(breaks = xbreaks, limits = xlimits, labels = xlabels, expand = c(0, 0)) +
		ggplot2::labs(x = xlabel, y = "CpGs (million)", color = "Coverage\nthreshold") +
		ggplot2::theme(legend.position = c(1, 0.5), legend.justification = c(0, 0.5)) +
		ggplot2::theme(plot.margin = unit(0.1 + c(0, 1, 0, 0), "in"))
	if (is.na(fname)) {
		return(pp)
	}

	## Save the plot to file(s)
	additional.params <- list(...)
	if (is.null(additional.params$width)) {
		additional.params$width <- img.width
	}
	if (is.null(additional.params$height)) {
		additional.params$height <- 7.2
	}
	additional.params$fname <- fname
	rplot <- do.call(createReportPlot, additional.params)
	print(pp)
	rplot <- off(rplot)

	## Attach the data as an attribute to the plot
	dframe$mc <- as.integer(as.character(dframe$mc))
	colnames(dframe) <- c("Minimal coverage", xlabel, "CpGs (million)")
	attr(rplot, "data") <- dframe
	invisible(rplot)
}
