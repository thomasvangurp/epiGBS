########################################################################################################################
## qualityControl.R
## created: 2012-05-31
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## 
########################################################################################################################

#' rnb.step.quality
#'
#' Cteates quality control plots section of the quality control report.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report on quality control to contain the generated sections. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return The modified report.
#'
#' @author Pavlo Lutsik
#' @noRd 
rnb.step.quality<-function(rnb.set, report){
	
	if (!inherits(rnb.set, "RnBSet")){
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}

	if(inherits(rnb.set,"RnBeadSet")){
		covg.lists <- NULL
	}else{ # inherits(rnb.set,"RnBiseqSet")
		covg.lists<-rnb.execute.quality(rnb.set)
	}
	report<-rnb.section.quality(report, rnb.set, covg.lists=covg.lists)

	return(report)
}

########################################################################################################################

#' rnb.execute.quality
#'
#' Performs quality control calculations on the loaded DNA methylation data set.
#'
#' @param object            Methylation dataset as an object of class \code{\linkS4class{RnBeadSet}},
#'                          \code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBiseqSet}}.
#' @param type              \code{character} vector of length \code{1} giving the type of genomic regions for which the
#' 				            quality control information is summarized.
#' @param qc.coverage.plots Flag indicating if sequencing coverage information is summarized and returned. This
#'                          parameter is considered only when \code{object} is of type \code{\linkS4class{RnBiseqSet}}.
#' @param verbose	        Flag specifying whether diagnostic output should be written to the console or to the
#'                          RnBeads logger in case the latter is initialized.
#' 
#' @details Currently, summarizing coverage for \code{\linkS4class{RnBiseqSet}} object is the only available function.
#' 
#' @return \code{\linkS4class{RnBeadSet}} object with imputed quality control information
#'
#' @author Pavlo Lutsik
#' @export
rnb.execute.quality<-function(
		object, 
		type="sites", 
		qc.coverage.plots=rnb.getOption("qc.coverage.plots"),
		verbose=TRUE){

  if(verbose){
  	rnb.logger.start("Preparing Quality Control Information")
  }

  if(inherits(object, "RnBiseqSet")){
	  if(qc.coverage.plots){
		  covg.rnbs<-covg(object, type)
		  covg.rnbs[is.na(covg.rnbs)]<-0
		  
		  if(type=="sites"){
			  type<-"CpG"
			  map<-sites(object)
		  }else{
			  map<-regions(object, type)
		  }
		  
		  covg.lists<-lapply(samples(object), function(sample){
		  
		  sample.id<-match(sample, samples(object))
		  assembly <- assembly(object)
		  covg.list<-lapply(names(rnb.get.chromosomes(assembly=assembly)), function(chr) {
					  chr.id<-match(chr, names(rnb.get.chromosomes(assembly=assembly)))
					  covg.weight<-rep(0, rnb.annotation.size(assembly=assembly)[chr.id])
					  covg.weight[map[map[,2]==chr.id,3]]<-covg.rnbs[map[,2]==chr.id,sample.id]
					  covg.weight
				  })
	  	  })
		  names(covg.lists)<-samples(object)  
		  result<-covg.lists
  	  }else{
		  result<-NULL
	  }
  }else{
	  result<-object
  }
  if(verbose){
  	logger.completed()
  }
  return(result)
}

#######################################################################################################################

#' rnb.section.quality
#'
#' Adds quality control probe or coverage section to the quality control report.
#'
#' @param report             Analysis report to contain the newly generated section. This must be an object of type
#'                           \code{\linkS4class{Report}}.
#' @param rnb.set            Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param qc.boxplots        Flag indicating if quality control box plots are to be generated. This option has effect
#'                           only when \code{rnb.set} is a HumanMethylation450 dataset.
#' @param qc.barplots        Flag indicating if quality control bar plots are to be generated. This option has effect
#'                           only when \code{rnb.set} is a HumanMethylation450 dataset.
#' @param qc.negativeboxplot Flag indicating if plots showing the negative control probes for each sample are to be
#'                           generated. This option has effect only when \code{rnb.set} is a
#'                           HumanMethylation450 dataset.
#' @param covg.lists         ...
#' @return The modified report.
#'
#' @author Pavlo Lutsik
#' @noRd 
rnb.section.quality<-function(report, rnb.set, qc.boxplots=rnb.getOption("qc.boxplots"),
	qc.barplots=rnb.getOption("qc.barplots"), qc.negative.boxplot=rnb.getOption("qc.negative.boxplot"),
	qc.coverage.plots=rnb.getOption("qc.coverage.plots"), qc.coverage.histograms=rnb.getOption("qc.coverage.histograms"),
	qc.coverage.violins=rnb.getOption("qc.coverage.violins"),
	qc.coverage.threshold.plot=rnb.getOption("qc.coverage.threshold.plot"),covg.lists=NULL){

	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!parameter.is.flag(qc.boxplots)) {
		stop("invalid value for qc.boxplots; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.barplots)) {
		stop("invalid value for qc.barplots; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.negative.boxplot)) {
		stop("invalid value for qc.negative.boxplots; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.coverage.plots)) {
		stop("invalid value for qc.coverage.plots; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.coverage.histograms)) {
		stop("invalid value for qc.coverage.histograms; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.coverage.violins)) {
		stop("invalid value for qc.coverage.violins; expected TRUE or FALSE")
	}
	## TODO: Add error handling for qc.coverage.threshold.plot
	## TODO: Add error handling for covg.lists

	logger.start("Quality Control Section")
	report <- rnb.add.section(report, "Quality Control", NULL)
	txt <- NULL

	if(inherits(rnb.set, "RnBeadSet")){
		if(is.null(qc(rnb.set))){

			txt <- c("The supplied dataset contains no quality information, therefore, quality control graphics could ",
				"not be generated.")
			rnb.add.paragraph(report, txt)
			logger.info("No quality information present in the dataset")

		}else{
		  
			txt <- "This section contains quality control plots and statistics for the methylation data."
			rnb.add.paragraph(report, txt)

			if(qc.boxplots) {
				txt <- c("Each box plot below shows the signal distribution of quality control probes across all ",
					"samples. The control box plots are separated by control types. Detailed description of the ",
					"control probes is given in the RnBeads vignette.")
				report <- rnb.add.section(report, "Quality Control Box Plots", txt, level = 2)
				report<-add.qc.boxplots(report, rnb.set)
				logger.status("Added quality control box plots")
			}

			if(qc.barplots) {
				txt <- c("The plots below visualize the exact signal levels at each quality control probe. Note that ",
					"the scale is not standardized. Background signal is usualy at the level of 1000 to 2000.")
				report <- rnb.add.section(report, "Quality Control Bar Plots", txt, level = 2)
				report<-add.qc.barplots(report, rnb.set, sample.batch.size=rnb.getOption("qc.sample.batch.size"))
				logger.status("Added quality control bar plots")
			}

			if(qc.negative.boxplot){
				txt <- c("Negative control box plots visualize background intensity distributions of all analyzed ",
					"samples. Samples with skewed distributions and high medians are likely to be of low quality and ",
					"should be discarded from the analysis.")
				report <- rnb.add.section(report, "Negative Control Box Plots", txt, level = 2)
				report<-add.negative.control.boxplot(report, rnb.set, sample.batch.size=rnb.getOption("qc.sample.batch.size"))
				logger.status("Added negative control boxplots")
			}
		}

	}else{ #inherits(rnb.set, "RnBiseqSet")


		if(qc.coverage.plots){
			txt <- c("The sequencing coverage plots visualized effective read coverage at all CpGs in the sequenced ",
					"genome. In case certain samples seem to have significantly decreased coverage, they should be excluded ",
					"from the analysis.")
			report <- rnb.add.section(report, "Sequencing Coverage Plots", txt, level = 2)
			report<-add.seq.coverage.plot(report, rnb.set, covg.lists)
			logger.status("Added sequencing coverage boxplots")
		}

		if(qc.coverage.histograms){
			txt <- c("The sequencing coverage histograms show distribution of coverage across all chromosomes. In case ",
			"certain samples seem to have significantly decreased coverage, they should be excluded from the analysis.")
			report <- rnb.add.section(report, "Sequencing Coverage Histograms", txt, level = 2)
			report<-add.seq.coverage.histograms(report, rnb.set)
			logger.status("Added sequencing coverage histograms")
		}
		
		if(qc.coverage.violins){
			txt <- c("The plots below show an alternative approach to visualizing the coverage distribution.")
			report <- rnb.add.section(report, "Sequencing Coverage Violin Plots", txt, level = 2)
			report<-add.seq.coverage.violins(report, rnb.set)
			logger.status("Added sequencing coverage violin plots")
		}

		if (length(qc.coverage.threshold.plot) != 0) {
			rplot <- rnb.plot.coverage.thresholds(rnb.set, qc.coverage.threshold.plot, fname = "coverage_interrogated",
				report = report)
			dframe.coverages <- attr(rplot, "data")
			fname <- "coverage_interrogated.csv"
			write.csv(dframe.coverages, file = file.path(rnb.get.directory(report, "data", TRUE), fname), row.names = FALSE)
			txt <- sprintf("%1.1f", range(dframe.coverages[dframe.coverages[, 2] == max(dframe.coverages[, 2]), 3]))
			txt <- c('In total, between ', txt[1], ' and ', txt[2], ' million sites are covered in all samples of the ',
				'dataset. The figure below shows the change in supports for different coverage thresholds. The exact ',
				'values are available in a dedicated <a href="', rnb.get.directory(report, 'data'), '/', fname,
				'">comma-separated file</a> accompanying this report.')
			report <- rnb.add.section(report, "Sequencing Coverage Thresholds", txt, level = 2)
			txt <- c("Line plot showing the number of CpG sites with a given support for different thresholds of ",
				"minimal coverage. The support of a CpG site is the minimal number of samples that interrogate it.")
			report <- rnb.add.figure(report, txt, rplot)
			rm(rplot, dframe.coverages, fname)
		}
	}

	if (is.null(txt)) {
		txt <- "No quality control plots are generated because all respective options are disabled."
		rnb.add.paragraph(report, txt)
	}

	logger.completed()
	return(report)
}

#######################################################################################################################

#' rnb.step.mixups
#'
#' Cteates sample mixups section in the quality control report 
#'
#' @param object a \code{\linkS4class{RnBeadSet}} object
#' @param report Report on methylation profiles to contain the dimension reduction section. This must be an object of
#'               type \code{\linkS4class{Report}}.
#'
#' @return the modified report object
#'
#' @author Pavlo Lutsik
#' @noRd 
rnb.step.mixups<-function(object, report){
	
	if(!inherits(object,"RnBSet")){
		stop("Supplied object is not of the class inheriting from RnBSet")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}

	logger.start("Visualizing Sample Mixups")
	if (any(unlist(rnb.options("qc.snp.heatmap", "qc.snp.boxplot", "qc.snp.barplot")))) {
		logger.start("Mixups Visualization Section")
		report<-rnb.section.mixups(report,object)
		logger.completed()
	}
	logger.completed()

	return(report)
}

#######################################################################################################################

#' rnb.section.mixups
#'
#' Adds sample mixups section to quality control report
#'
#' @param report           Report to contain the new section. This must be an object of type \code{\linkS4class{Report}}.
#' @param object           Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param qc.snp.heatmap   ...
#' @param qc.snp.distances Flag indicating if intersample distances based on the beta values of SNP probes are to be
#'                         displayed.
#' @param qc.snpboxplot    Flag indicating if box plots of beta values for each of SNP probes is to be displayed.
#' @param qc.snpbarplot    Flag indicating if bar plots of beta values at SNP probes, one bar plot par sample.
#' @return The (possibly modified) report.
#'
#' @author Pavlo Lutsik
#' @noRd 
rnb.section.mixups<-function(report, object,
	qc.snp.heatmap=rnb.getOption("qc.snp.heatmap"),
	qc.snp.distances=rnb.getOption("qc.snp.distances"),
	qc.snp.boxplot=rnb.getOption("qc.snp.boxplot"),
	qc.snp.barplot=rnb.getOption("qc.snp.barplot")) {

	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (!inherits(object, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!parameter.is.flag(qc.snp.heatmap)) {
		stop("invalid value for qc.snp.heatmap")
	}
	if (!parameter.is.flag(qc.snp.distances)) {
		stop("invalid value for qc.snp.distances")
	}
	if (!parameter.is.flag(qc.snp.boxplot)) {
		stop("invalid value for qc.snp.boxplot")
	}
	if (!parameter.is.flag(qc.snp.barplot)) {
		stop("invalid value for qc.snp.barplot")
	}

	section.title <- "Visualization of Sample Mixups"
	if(object@target=="probes450"){
		snp.probes.present <- any(grepl("^rs", annotation(object, add.names = TRUE)[, "ID"]))
	}else if(object@target=="probes27"){
		snp.probes.present <- !is.null(qc(object))
	}
	if (!snp.probes.present) {
		txt <- c("Overview of SNP-based probes and sample comparison based on them cannot be performed ",
			"because the dataset does not contain such probes.")
		report <- rnb.add.section(report, section.title, txt)
		return(report)
	}
	txt <- c("SNP-based box and bar plots that facilitate identification of the sample mixups.")
	report <- rnb.add.section(report, section.title, txt)

	add.info <- function(stitle, ffunction, txt) {
		result <- rnb.add.section(report, stitle, txt, level = 2)
		result <- ffunction(result, object)
		rnb.status(c("Added", stitle))
		result
	}

	if (qc.snp.heatmap) {
		txt <- "SNP Heatmap allows for identification of sample mixups in particular for genetically matched designs."
		report <- add.info("SNP Heatmap", add.snp.heatmap, txt)
	}
	if (qc.snp.distances) {
		txt <- c("If we inspect the samples in the space defined by the SNP probes only, samples appear close to ",
			"each other are genetically similar. We calculated the Manhattan distance between all vectors of SNP ",
			"probes. The figure below shows the relative distances.")
		report <- add.info("SNP-based Distances", add.snp.distances, txt)
	}
	if (qc.snp.boxplot) {
		txt <- c("SNP box plot allows the identification of sample mixups in case a genetically homogeneous set ",
			"of samples is analyzed.")
		report <- add.info("SNP Box Plot", add.snp.boxplot, txt)
	}
	if (qc.snp.barplot) {
		txt <- "SNP bar plots allow identification of sample mixups, in particular for genetically matched designs."
		report <- add.info("SNP Bar Plots", add.snp.barplot, txt)
	}

	return(report)
}

########################################################################################################################
## UTILS
########################################################################################################################

add.qc.boxplots<-function(report, object){

  descr<-"Quality control box plots."

  if(object@target=="probes450"){
  	ctypes<-rnb.infinium.control.targets(object@target)[c(13,4,3,14,1:2,11:12,6)]
  }else if(object@target=="probes27"){
	ctypes<-rnb.infinium.control.targets(object@target)[c(10,3,2,11,1,9,6)]
  }

  cplots<-lapply(ctypes, rnb.plot.control.boxplot, rnb.set=object, report=report, writeToFile=T, numeric.names=T, width=8, height=6, low.png=100, high.png=300)
  names(cplots)<-1:length(ctypes)

  sn<-list("Control probe type"=ctypes)
  names(sn[[1]])<-1:length(ctypes)

  report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
  report

}

#######################################################################################################################

add.qc.barplots<-function(report, object, sample.batch.size=50){

  descr="Quality control bar plots."
  
  if(object@target=="probes450"){
  	cmd <- rnb.get.annotation("controls450")
  	ctypes<-unique(cmd$Target)[unique(cmd$Target) %in% rnb.infinium.control.targets("probes450")[c(13,4,14,3,1:2,11:12,6)]]
  }else if(object@target=="probes27"){
	cmd <- rnb.get.annotation("controls27")
	ctypes<-unique(cmd$Type)[unique(cmd$Type) %in% rnb.infinium.control.targets("probes27")[c(10,3,2,11,1,9,6)]]
  }
  nsamp<-length(samples(object))
  
  plot.names<-NULL
  
  if(nsamp %% sample.batch.size==1){
	  sample.batch.size<-sample.batch.size-5
  }
  portion.starts<-0:(nsamp %/% sample.batch.size)*sample.batch.size+1
  portion.ends<-portion.starts+sample.batch.size-1
  portion.ends[length(portion.ends)]<-nsamp
  portions<-paste(portion.starts, portion.ends, sep="-")
  
  plots<-lapply(1:length(portions),function(portion.id){
  	 
	  cplots<-lapply(ctypes, function(type){

		if(object@target=="probes450"){
			cmdt <- cmd[cmd[["Target"]] == type, ]
			pn<-paste(type, 1:(dim(cmdt)[1]),  sep=".")
		}else if(object@target=="probes27"){
			cmdt <- cmd[cmd[["Type"]] == type, ]
			pn<-as.character(cmdt$Name)
		}
		
		if(portion.id==1) plot.names<<-c(plot.names, pn)
	    plots<-lapply(pn, rnb.plot.control.barplot, rnb.set=object, 
				sample.subset=portion.starts[portion.id]:portion.ends[portion.id], 
				report=report, writeToFile=T, numeric.names=T, width=8, height=6, low.png=100, high.png=300, verbose=T,
				name.prefix=portions[portion.id])
	
		if(object@target=="probes450"){
			names(plots)<-paste(type, 1:(dim(cmdt)[1]))
		}else if(object@target=="probes27"){
			names(plots)<-as.character(cmdt$Name)
		}
	    plots
		
	  })
	
	  names(cplots)<-NULL
	  
	  cplots<-unlist(cplots)
	  names(cplots)<-1:length(plot.names)
	  cplots
  })

  plots<-unlist(plots)

  sn<-list("Samples #: " = portions, "Control probe ID" = plot.names)
  
    
  names(sn[[1]])<-portions
  if(object@target=="probes450"){
  	names(sn[[2]])<-1:length(plot.names)
  }else if(object@target=="probes27"){
	names(sn[[2]])<-match(plot.names,cmd$Name)
  }
    
  report<-rnb.add.figure(report, description=descr, report.plots=plots, setting.names=sn)

}

#######################################################################################################################

add.negative.control.boxplot<-function(report, object, sample.batch.size=50){
	
	descr<-"Box plots of the negative control probes."
	
#   cplots<-lapply(ctypes, control.boxplot, rnb.set=object, report=report, writeToFile=T)
#   names(cplots)<-gsub(" ", ".", ctypes)
#   
#   sn<-list(tolower(ctypes))
#   names(sn[[1]])<-gsub(" ", ".", ctypes)
	nsamp<-length(samples(object))
	if(nsamp %% sample.batch.size==1){
		sample.batch.size<-sample.batch.size-5
	}
	portion.starts<-0:(nsamp %/% sample.batch.size)*sample.batch.size+1
	portion.ends<-portion.starts+sample.batch.size-1
	portion.ends[length(portion.ends)]<-nsamp
	portions<-paste(portion.starts, portion.ends, sep="-")
	
	cplots<-lapply(1:length(portions), function(portion.id){
				rnb.plot.negative.boxplot(object, sample.subset=portion.starts[portion.id]:portion.ends[portion.id],
						name.prefix=portions[portion.id], writeToFile=T, report=report, width=10, height=6, low.png=75, high.png=100)
	})

	sn<-list("Samples #: " = portions)
	
	names(sn[[1]])<-portions
	
	report<-rnb.add.figure(report, description=descr, report.plots=cplots, sn)
	report
	
}

#######################################################################################################################

add.snp.boxplot<-function(report, object){
  
  descr<-"Box plots of the SNP probes."
  
#   cplots<-lapply(ctypes, control.boxplot, rnb.set=object, report=report, writeToFile=T)
#   names(cplots)<-gsub(" ", ".", ctypes)
#   
#   sn<-list(tolower(ctypes))
#   names(sn[[1]])<-gsub(" ", ".", ctypes)
  cplots<-list(SNPBoxplot=rnb.plot.snp.boxplot(object, writeToFile=T, report=report, width=9, height=6, low.png=100, high.png=300))
  report<-rnb.add.figure(report, description=descr, report.plots=cplots)
  report
  
}

#######################################################################################################################

add.snp.barplot<-function(report, object){
	
	descr<-"Bar plots of the SNP probes beta values. Genetically matched samples should show very similar profiles."

	ids <- samples(object)

	cplots<-lapply(ids, rnb.plot.snp.barplot, rnb.set=object, report=report, writeToFile=T, numeric.names=T, width=9, height=6, low.png=100, high.png=300)
	names(cplots)<-1:length(ids)
	
	sn<-list("Sample labels"=tolower(ids))
	names(sn[[1]])<-1:length(ids)
	
	
	report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
	
	report
	
}

#######################################################################################################################

add.snp.heatmap<-function(report, object){
	
	descr<-"Heat map of the SNP probes. The samples with the same genetic background should cluster together."
	
#   cplots<-lapply(ctypes, control.boxplot, rnb.set=object, report=report, writeToFile=T)
#   names(cplots)<-gsub(" ", ".", ctypes)
#   
#   sn<-list(tolower(ctypes))
#   names(sn[[1]])<-gsub(" ", ".", ctypes)
	cplots<-list(SNPHeatmap=rnb.plot.snp.heatmap(object, writeToFile=T, report=report, width=8, height=9, low.png=100, high.png=300))
	report<-rnb.add.figure(report, description=descr, report.plots=cplots)
	report
	
}

#######################################################################################################################

## add.snp.distances
##
## Adds a section about sample distances based on beta values of SNP probes.
##
## @param report Report to contain the section on SNP probe distances.
## @param object Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
## @return The (possibly modified) report.
## @author Yassen Assenov
add.snp.distances <- function(report, object) {

	## Extract the matrix of beta values on the SNP probes
	if(object@target=="probes450"){
		snp.betas <- meth(object, row.names = TRUE)
		snp.betas <- snp.betas[grep("^rs", rownames(snp.betas)), ]
	}else if(object@target=="probes27"){
		snp.betas <- HM27.snp.betas(qc(object))
	}
	if (!(is.matrix(snp.betas) && nrow(snp.betas) > 1)) {
		return(report)
	}

	## Calculate Manhattan distances between samples based on the SNP probe intensities
	snp.distances <- as.matrix(dist(t(snp.betas), method = "manhattan"))
	rnb.status("Calculated Manhattan distances between samples based on SNP probes")

	report.plots <- list()
	setting.names <- list()
	if (nrow(snp.distances) <= 24) {

		## Create a diagonal heatmap of distances
		tbl.melt <- symmetric.melt(snp.distances)
		colnames(tbl.melt)[3] <- "distance"
		colors.g <- rnb.getOption("colors.gradient")
		pp <- ggplot(tbl.melt, aes_string(x = "x", y = "y", fill = "distance")) + labs(x = NULL, y = NULL) +
			coord_fixed(ratio = 1) + geom_tile(color = "white") +
			scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
			scale_fill_gradient(na.value = "white", low = colors.g[1], high = colors.g[2]) +
			theme(legend.justification = c(0, 1), legend.position = c(1, 1)) +
			theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
			theme(panel.grid.major = element_blank(), panel.background = element_blank()) +
			theme(panel.border = element_blank(), plot.margin = unit(c(0.1, 1.1, 0.1, 0.1), "in"))
		## Fix the areas for x and y axis labels
		i.height <- 2.2 + nrow(snp.distances) * 0.3
		i.width <- 4 + nrow(snp.distances) * 0.3
		report.plots <- createReportPlot("snp_low_dimensional", report, width = i.width, height = i.height)
		pp <- suppressWarnings(ggplot_gtable(ggplot_build(pp)))
		pp$widths[[3]] <- unit(2, "in")
		pp$heights[[length(pp$heights) - 2L]] <- unit(2, "in")
		grid.draw(pp)
		report.plots <- off(report.plots)
		txt <- c("Distances between pairs of samples based on the SNP probes only.")
		rm(tbl.melt, colors.g, i.height, i.width)
		
	} else {
		
		## Create PCA plots
		pr.coords <- prcomp(t(snp.betas), center = TRUE, scale. = FALSE)$x[, 1:2]
		labels.pc <- paste("Principal component", 1:ncol(pr.coords))
		dframe <- data.frame(x = pr.coords[, 1], y = pr.coords[, 2], label = rownames(pr.coords))
		plot.types <- c("point" = "points", "label" = "identifiers")
		for (ptype in names(plot.types)) {
			pp <- ggplot(dframe, aes_string(x = "x", y = "y", label = "label"))
			if (ptype == "point") {
				pp <- pp + geom_point(size = 3)
			} else {
				pp <- pp + geom_text(angle = 45)
			}
			pp <- pp + labs(x = labels.pc[1], y = labels.pc[2]) +
				theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
			fname <- paste0("snp_low_dimensional_", ptype)
			rplot <- createReportPlot(fname, report, width = 6.2, height = 6.2)
			print(pp)
			report.plots <- c(report.plots, off(rplot))
		}
		txt <- c("Scatter plot showing the samples in the first two principal components of the space of their SNP ",
			"probes.")
		setting.names <- list("Display samples as" = plot.types)
		rm(pr.coords, labels.pc, dframe, plot.types, ptype, fname, rplot)
	}
	report <- rnb.add.figure(report, txt, report.plots, setting.names)
	rnb.status("Add plot of SNP distances")
	rm(report.plots, setting.names, pp)

	## Export the distances to file
	fname <- "snp_distances.csv"
	write.csv(snp.distances, file = file.path(rnb.get.directory(report, "data", TRUE), fname))
	txt <- c('The full table of all pairwise distances is stored in a dedicated <a href="',
		rnb.get.directory(report, "data"), '/', fname, '">comma-separated file</a> accompanying this report.')
	rnb.add.paragraph(report, txt)
	rnb.status(c("Exported SNP-based distances to", fname))

	report
}

#######################################################################################################################

add.seq.coverage.plot<-function(report, object, covg.lists=NULL){
	
	descr<-"Sequencing coverage plots visualize effective read coverage over all chromosomes of the genome."
	
	ids <- samples(object)
	
	cplots<-lapply(ids, function(id) rnb.plot.biseq.coverage(rnbs.set=object, sample=id, report=report, writeToFile=T, numeric.names=T,create.pdf=F, width=8, height=16, low.png=100, high.png=200, covg.lists=covg.lists[[id]]))
	
	names(cplots)<-1:length(ids)
	
	sn<-list("Sample labels" = ids)
	names(sn[[1]])<-1:length(ids)
	
	report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
	report
	
}


#######################################################################################################################

add.seq.coverage.histograms<-function(report, object){
	
	descr<-"Sequencing coverage histogram visualize the bulk distribution of read coverage for each sample."
	
	ids <- samples(object)
	
	cplots<-lapply(ids, rnb.plot.biseq.coverage.hist, rnbs.set=object, report=report, writeToFile=T, numeric.names=T, covg.max.percentile=0.99, width=8, height=7, low.png=100, high.png=300)
	
	names(cplots)<-1:length(ids)
	
	sn<-list("Sample labels" = ids)
	names(sn[[1]])<-1:length(ids)
	
	report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
	report
	
}

#######################################################################################################################

add.seq.coverage.violins<-function(report, object, sample.chunk.size=20){
	
	descr<-paste("Sequencing coverage histogram visualized as violin plots. Distributions are based on",
				 nrow(meth(object)),"methylation sites.")
	
	ids <- samples(object)
	n <- length(ids)
	num.intervals <- ceiling(n/sample.chunk.size)
	if (num.intervals>1){
		ssets <- cut(1:n,ceiling(n/sample.chunk.size),labels=FALSE)
	} else {
		ssets <- rep(1,n)
	}
	ssets.n <- length(unique(ssets))
	max.covg.2plot<-quantile(covg(object),probs=c(0.99), na.rm=TRUE)[1]
	
	cplots <- list()
	for (i in 1:ssets.n){
		ssamples <- ids[ssets==i]
		rep.plot <- rnb.plot.biseq.coverage.violin(rnbs.set=object,samples=ssamples,fname=paste("coverageViolin_chunk",i,sep=""),covg.range=c(0,max.covg.2plot),report=report, width=8, height=7, low.png=100, high.png=300)
		cplots<-c(cplots,list(rep.plot))
	}
	
	names(cplots)<-1:ssets.n
	
	sn<-list("Sample chunk" = 1:ssets.n)
	names(sn[[1]])<-paste("chunk",1:ssets.n,sep="")
	
	report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
	report
}
#######################################################################################################################
