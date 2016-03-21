########################################################################################################################
## controlPlots.R
## created: 2012-04-04
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Quality control plots for HumanMethylation450 experiments.
########################################################################################################################

## G L O B A L S #######################################################################################################

.control.plot.theme.samples<-ggplot2::theme(
					plot.margin=grid::unit(c(0.1,1,0.1,0.1), "cm"), 
					axis.title.x = ggplot2::element_blank(),
					axis.text.x=ggplot2::element_text(size=7, angle=-45, hjust=0))

HM450.SAMPLE.INDEPENDENT<-c("STAINING", "EXTENSION", "HYBRIDIZATION", "TARGET REMOVAL")
HM450.SAMPLE.DEPENDENT<-c("BISULFITE CONVERSION I","BISULFITE CONVERSION II",
		"SPECIFICITY I","SPECIFICITY II", "NON-POLYMORPHIC" )

## F U N C T I O N S ###################################################################################################
#'
#' rnb.plot.control.boxplot
#'
#' Box plots of various control probes
#'
#' @param rnb.set \code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBeadSet}} object with valid quality 
#' 				  control information.
#' @param type    type of the control probe; must be one of the \code{"BISULFITE CONVERSION I"},
#'                \code{"BISULFITE CONVERSION II"}, \code{"EXTENSION"}, \code{"HYBRIDIZATION"}, \code{"NEGATIVE"},
#'                \code{"NON-POLYMORPHIC"}, \code{"NORM_A"}, \code{"NORM_C"}, \code{"NORM_G"}, \code{"NORM_T"},
#'                \code{"SPECIFICITY I"}, \code{"SPECIFICITY II"}, \code{"STAINING"}, \code{"TARGET REMOVAL"}.
#' @param writeToFile flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param numeric.names if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file name with digits
#' @param ... other arguments to \code{\link{createReportPlot}}
#' 
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class 
#' 			\code{\link{ggplot}} otherwise.
#'
#' @author Pavlo Lutsik
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.plot.control.boxplot(rnb.set.example)
#' }
rnb.plot.control.boxplot <- function(
		rnb.set, 
		type = rnb.infinium.control.targets(rnb.set@target)[1], 
		writeToFile = FALSE, 
		numeric.names = FALSE, 
		...) {

	if(!type %in% rnb.infinium.control.targets(rnb.set@target)) {
		stop(paste("Unrecognized control type:", type))
	}

	## Extract intensities of the control probes
	if(rnb.set@target=="probes450"){
		meta <- rnb.get.annotation("controls450")
	}else if(rnb.set@target=="probes27"){
		meta <- rnb.get.annotation("controls27")
	}
	
	if(rnb.set@target=="probes450"){
		types<-rnb.infinium.control.targets(rnb.set@target)[c(13,4,3,14,1:2,11:12,6)]
	}else if(rnb.set@target=="probes27"){
		types<-rnb.infinium.control.targets(rnb.set@target)[c(10,3,2,11,1,9,6)]
	}
	
	if(!type %in% types){
		warning("Unoptimized probe type, plotting performance may be decreased")
	}
	
	if(rnb.set@target=="probes450"){
		meta <- meta[type == meta[["Target"]], ]
		ids<-as.character(meta[["ID"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[ids,]
	}else if(rnb.set@target=="probes27"){
		meta <- meta[type == meta[["Type"]], ]
		ids<-as.character(meta[["Address"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[meta$Address %in% ids,]
	}

	if(nrow(meta)==0){
		if(writeToFile) {
			plot.file<-createReportPlot(paste("ControlBoxPlot", ifelse(numeric.names, match(type, types), gsub(" ", ".", type)) , sep="_"), ...)
		}
		plot.obj<-rnb.message.plot(sprintf("Box plot for control type %s not available", type))
		if(writeToFile){
			print(plot.obj)
			off(plot.file)
			return(plot.file)
		}
		return(plot.obj)
	}
	
	intensities <- lapply(c("green" = "Cy3", "red" = "Cy5"), function(channel) {
		result <- qc(rnb.set)[[channel]][ids, ]
		if (is.vector(result)) {
			return(as.matrix(result))
		}
		return(t(result))
	})

	scales<-lapply(qc(rnb.set), get.unified.scale)
	
	if(rnb.set@target=="probes450"){
		## Shorten the words describing probe's expected intensity
		INTENSITIES <- c("Background" = "Bgnd", "High" = "High", "Low" = "Low", "Medium" = "Med")
		levels(meta[, "Expected Intensity"]) <- INTENSITIES[levels(meta[, "Expected Intensity"])]
		meta[, "Expected Intensity"] <- as.character(meta[, "Expected Intensity"])
	
		## Define bar labels in the plot
		plot.names<-sapply(1:length(ids), sprintf, fmt="Probe%2g")
		
		exp.intensity <- lapply(c("green" = "Green", "red" = "Red"), function(channel) {
			ifelse(meta[[paste("Evaluate", channel)]] == "+", meta[["Expected Intensity"]], "Bgnd")
		})
		plot.names <- lapply(exp.intensity, function(exps) { paste(plot.names, exps, sep = ":\n") })
	
		plot.data <- lapply(names(intensities), function(channel) {
			channel.data <- intensities[[channel]]
			data.frame(
				"Intensity" = as.numeric(channel.data),
				"Probe" = as.factor(sapply(plot.names[[channel]], rep, times=nrow(channel.data))))
		})
	}else if(rnb.set@target=="probes27"){
		plot.names<-list()
		plot.names$green <- gsub("conversion |Extension |mismatch ", "", meta$Name)
		plot.names$red <- gsub("conversion |Extension |mismatch ", "", meta$Name)
		plot.data <- lapply(names(intensities), function(channel) {
					channel.data <- intensities[[channel]]
					data.frame(
							"Intensity" = as.numeric(channel.data),
							"Probe" = as.factor(sapply(plot.names[[channel]], rep, times=nrow(channel.data))))
				})
	}
	
	if(writeToFile){
		plot.file<-createReportPlot(paste("ControlBoxPlot", ifelse(numeric.names, match(type, types), gsub(" ", ".", type)) , sep="_"), ...)
	}

	## Define viewports and assign it to grid layout
	#grid.newpage()
	#pushViewport(viewport(layout = grid.layout(2,1)))

	plots<-list()

	for (i in 1:length(intensities)) {
		channel <- names(intensities)[i]
#		print(qplot(Probe, Intensity, data=na.omit(plot.data[[i]]), geom = "boxplot",
		plots[[i]]<-qplot(Probe, Intensity, data=na.omit(plot.data[[i]]), geom = "boxplot",
				main=paste(type, paste(channel, "channel"), sep=": "),
				ylab="Intensity",fill=I(channel)) + scale_y_continuous(limits=scales[[i]])
		#,vp=viewport(layout.pos.row=i, layout.pos.col=1))
	}

	grb<-do.call(arrangeGrob, plots)
	
	if(writeToFile) {
		print(grb)
		off(plot.file)
		return(plot.file)
	}
	return(grb)
}

#######################################################################################################################

#' rnb.plot.negative.boxplot
#'
#' Box plots of negative control probes
#'
#' @param rnb.set 		\code{\linkS4class{RnBeadSet}} object with valid quality control information
#' @param sample.subset an integer vector specifying the subset of samples for which the plotting should be performed
#' @param writeToFile 	flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param name.prefix	in case \code{writeToFile} is \code{TRUE}, a \code{character} singleton specifying a prefix to the variable part of the image file names 
#' @param ... 			other arguments to \code{\link{createReportPlot}}
#' 
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class 
#' 			\code{\link{ggplot}} otherwise.
#' 
#' @author Pavlo Lutsik
#' @export
#' @examples 
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.plot.negative.boxplot(rnb.set.example)
#' }
rnb.plot.negative.boxplot<- function(
		rnb.set,  
		sample.subset=1:length(samples(rnb.set)), 
		writeToFile = FALSE, 
		name.prefix=NULL, 
		...) {
		
	## Extract intensities of the control probes
	if(rnb.set@target=="probes450"){
		meta <- rnb.get.annotation("controls450")
		meta <- meta["NEGATIVE" == meta[["Target"]], ]
		ids<-as.character(meta[["ID"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[ids,]
	}else if(rnb.set@target=="probes27"){
		meta <- rnb.get.annotation("controls27")
		meta <- meta["Negative" == meta[["Type"]], ]
		ids<-as.character(meta[["Address"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[meta$Address %in% ids,]
	}
	if(nrow(meta)==0){
		if(writeToFile) plot.file<-createReportPlot(paste("ControlBoxPlot", ifelse(numeric.names, match(type, types), gsub(" ", ".", type)) , sep="_"), ...)
		print(rnb.message.plot("Box plot for negative control probes is not available"))
		if(writeToFile){
			off(plot.file)
			return(plot.file)
		}
		return(invisible())
	}
	
	intensities <- lapply(c("green" = "Cy3", "red" = "Cy5"), function(channel) {
				result <- qc(rnb.set)[[channel]][ids,sample.subset]
				if (is.vector(result)) {
					return(as.matrix(result))
				}
				return(result)
			})
	
	## Shorten the words describing probe's expected intensity
#	INTENSITIES <- c("Background" = "Bgnd", "High" = "High", "Low" = "Low", "Medium" = "Med")
#	levels(meta[, "Expected Intensity"]) <- INTENSITIES[levels(meta[, "Expected Intensity"])]
#	meta[, "Expected Intensity"] <- as.character(meta[, "Expected Intensity"])
	
	## Define bar labels in the plot
#	plot.names<-sapply(1:length(ids), sprintf, fmt="Probe%2g")
	sample.ids<-samples(rnb.set)[sample.subset]
	
	## FIXME: This might leads to duplicated identifiers and/or reordering of samples 
	sample.labels<-abbreviate.names(sample.ids)
	
#	exp.intensity <- lapply(c("green" = "Green", "red" = "Red"), function(channel) {
#				ifelse(meta[[paste("Evaluate", channel)]] == "+", meta[["Expected Intensity"]], "Bgnd")
#			})
	
#	plot.names <- lapply(exp.intensity, function(exps) { paste(plot.names, exps, sep = ":\n") })
	
	plot.data <- lapply(names(intensities), function(channel) {
				channel.data <- intensities[[channel]]
				data.frame(
						"Intensity" = as.numeric(channel.data),
						"Sample" = factor(sapply(sample.ids, rep, times=nrow(channel.data)), as.character(sample.ids)))
			})
	if(is.null(name.prefix)){
		plot.name<-"NegativeControlBoxPlot"
	}else{
		plot.name<-paste("NegativeControlBoxPlot", name.prefix, sep="_")
	}
	if(writeToFile) {
		plot.file<-createReportPlot(plot.name, ...)
	}
	
	## Define viewports and assign it to grid layout
	#grid.newpage()
	#pushViewport(viewport(layout = grid.layout(2,1)))

	plots<-list()
	
	for (i in 1:length(intensities)) {
		channel <- names(intensities)[i]
		plots[[i]]<-qplot(Sample, Intensity, data=na.omit(plot.data[[i]]), geom = "boxplot",
						main=paste("Negative control probes ", paste(channel, "channel"), sep=": "),
						ylab="Intensity",fill=I(channel))+
						scale_x_discrete(labels=sample.labels)+
						.control.plot.theme.samples
						#,vp=viewport(layout.pos.row=i, layout.pos.col=1))
	}
	
	grb<-do.call(arrangeGrob, plots)
	
	if(writeToFile) {
		print(grb)
		off(plot.file)
		return(plot.file)
	}
	return(grb)
}
#######################################################################################################################
#' rnb.plot.control.barplot
#'
#' Per-sample bar plots of Illumina HumanMethylation control probes
#' 
#' @param rnb.set 		\code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBeadSet}} object with valid 
#' 						quality control information
#' @param probe 		exact id of the control probe consisting of the control probe type (see \code{\link{rnb.plot.control.boxplot}})
#' @param sample.subset an integer vector specifying the subset of samples for which the plotting should be performed
#' @param writeToFile 	flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param numeric.names if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file name with digits
#' @param name.prefix	in case \code{writeToFile} is \code{TRUE}, a \code{character} singleton specifying a prefix to the variable part of the image file names 
#' @param verbose		if \code{TRUE} additional diagnostic output is generated 
#' @param ... 			other arguments to \code{\link{createReportPlot}}
#' 
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class 
#' 			\code{\link{ggplot}} otherwise.
#' 
#' @author Pavlo Lutsik
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' control.meta.data <- rnb.get.annotation("controls450")
#' ctrl.probe<-paste0(unique(control.meta.data[["Target"]])[4], ".5")
#' print(ctrl.probe) # EXTENSION.5
#' rnb.control.barplot(rnb.set.example, ctrl.probe)
#' }
rnb.plot.control.barplot<-function(
		rnb.set, 
		probe, 
		sample.subset=1:length(samples(rnb.set)), 
		writeToFile=FALSE, 
		numeric.names=FALSE, 
		name.prefix=NULL, 
		verbose=F,
		...)
{
	
	if(rnb.set@target=="probes450"){
		control.meta.data <- rnb.get.annotation("controls450")
		type=strsplit(probe,".", fixed=TRUE)[[1]][1]
		index=strsplit(probe,".", fixed=TRUE)[[1]][2]
		if(!(type %in% rnb.infinium.control.targets())){
			rnb.error(c("Unrecognized control probe:", probe))
		}
		id<-subset(control.meta.data, Target==type & Index==index)$ID
		id.col<-"ID"
	}else if(rnb.set@target=="probes27"){
		control.meta.data <- rnb.get.annotation("controls27")		
		if(!probe %in% control.meta.data$Name){
			rnb.error(c("Unrecognized control probe:", probe))
		}
		id<-control.meta.data[probe==control.meta.data[["Name"]],"Address"]
		id.col<-"Address"
	}
	
	#get intensities
	green<-qc(rnb.set)$Cy3[,sample.subset, drop=FALSE]
	red<-qc(rnb.set)$Cy5[,sample.subset, drop=FALSE]
	
	if(! id %in% rownames(green)){
		if(writeToFile){
			if(is.null(name.prefix)){
				plot.name<-paste('ControlBarPlot', ifelse(numeric.names, match(id, control.meta.data[[id.col]] ), probe.name) , sep="_")
			}else{
				plot.name<-paste('ControlBarPlot', name.prefix, ifelse(numeric.names, match(id, control.meta.data[[id.col]] ), probe.name) , sep="_")
			}
			plot.file<-createReportPlot(plot.name, ...)
		}
		plot.obj<-rnb.message.plot(sprintf("Box plot for control type %s not available", type))
		if(writeToFile){
			print(plot.obj)
			off(plot.file)
			return(plot.file)
		}
		return(plot.obj)
	}
		
	green<-green[as.character(id),]
	red<-red[as.character(id),]

#   if(is.null(dim(green))) green<-t(as.matrix(green))
#   if(is.null(dim(red))) red<-t(as.matrix(red))

	## get meta information
	
	if(rnb.set@target=="probes450"){
		meta<-subset(control.meta.data, ID==id)
	}else if(rnb.set@target=="probes27"){
		meta<-subset(control.meta.data, Address==id)
	}

	probe.name<-gsub(" ", ".", probe)
	if(is.null(name.prefix)){
		plot.name<-paste('ControlBarPlot', ifelse(numeric.names, match(id, control.meta.data[[id.col]] ), probe.name) , sep="_")
	}else{
		plot.name<-paste('ControlBarPlot', name.prefix, ifelse(numeric.names, match(id, control.meta.data[[id.col]] ), probe.name) , sep="_")
	}
	if(writeToFile){
		plot.file<-createReportPlot(plot.name, ...)
	}

	Samples<-factor(samples(rnb.set)[sample.subset], as.character(samples(rnb.set)[sample.subset]))
	## shorten sample names 
	
	## FIXME: This might leads to duplicated identifiers and/or reordering of samples 
	sample.labels<-abbreviate.names(samples(rnb.set)[sample.subset])
	
	#grid.newpage()
	# define viewports and assign it to grid layout
	#pushViewport(viewport(layout = grid.layout(2,1)))
	
	### plot green channel
	
	if(rnb.set@target=="probes450"){
		main_txt_grn<-paste(probe, meta[,"Description"], "green channel", if(meta[, "Evaluate Green"]=="+") meta[, "Expected Intensity"] else "Background", sep=": ")
		main_txt_red<-paste(probe, meta[,"Description"], "red channel",if(meta[, "Evaluate Red"]=="+") meta[, "Expected Intensity"] else "Background", sep=": ")
	}else{
		main_txt_grn<-paste(probe, sep="")
		main_txt_red<-paste(probe, sep="")
	}
	
	empty.probe<-sum(sapply(green, is.na))==length(green) 
	
    if(empty.probe) {
		green<-rep(0, length(green))
	}
	green.plot<-qplot(Samples, green, stat="identity", geom="bar",fill=I("green"), 
					main=main_txt_grn)+
					scale_x_discrete(labels=sample.labels)+ylab("Intensity")+
					.control.plot.theme.samples
			
	if(empty.probe) {
		green.plot<-green.plot+annotate("text", label="No intensities found for this quality control probe", x=length(green)%/%2, y=0.5, size=5)
	}
	#print(green.plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
	
	### plot red channel
	if(sum(sapply(red, is.na))==length(red)) empty.probe<-TRUE else empty.probe<-FALSE
	if(empty.probe) {
		red<-rep(0, length(red));
	}
		red.plot<-qplot(Samples, red, stat="identity", geom="bar",fill=I("red"), 
						main=main_txt_red)+ 
					scale_x_discrete(labels=sample.labels)+ylab("Intensity")+
					.control.plot.theme.samples
					
	if(empty.probe){
		red.plot<-red.plot+annotate("text", label="No intensities found for this quality control probe", x=length(red)%/%2, y=0.5, size=5)
	}
	#print(red.plot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
		
  grb<-arrangeGrob(green.plot, red.plot)
  
  if(writeToFile) {
	print(grb)
  	off(plot.file)
  	return(plot.file)
   }
  
  return(grb)
  
}
#######################################################################################################################

## control.probe.PCA
##
## Correlation plots of the principal components of the data set and the quality control information
## 
## @param qc.set Quality control information of an \code{\linkS4class{RnBeadSet}} object.
## @param pca ...
## @param type Desired type of the control probe, one of the c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", 
## "EXTENSION","HYBRIDIZATION", "NEGATIVE", "NON-POLYMORPHIC", "NORM_A", "NORM_C","NORM_G", "NORM_T", "SPECIFICITY I", 
## "SPECIFICITY II", "STAINING","TARGET REMOVAL")
## @param ... other arguments to \code{\link{createReportPlot}}
##
## @export
## @author Pavlo Lutsik
control.probe.PCA <- function(
		qc.set, 
		pca, 
		type, 
		...) {

	if (!(type %in% rnb.infinium.control.targets())){
		rnb.error(c("Unrecognized control probe type:", type))
	}

	control.meta.data <- rnb.get.annotation("controls450")

	green <- qc.set$Cy3
	red <- qc.set$Cy5

	ids <- control.meta.data$ID[grep(paste(type, ".", sep = ""), paste(control.meta.data$Target, ".", sep = ""), fixed = TRUE)]
	green <- green[as.character(ids),]
	red <- red[as.character(ids),]

#	if(is.null(dim(green))) green<-t(as.matrix(green))
#	if(is.null(dim(red))) red<-t(as.matrix(red))

	## Get meta information
#   fd.data<-featureData(cd)@data
#
#   ids<-fd.data[probe, "ProbeID"]
#   meta<-subset(control.meta.data, ID %in% ids)
#

	if (!is.null(red)) {
		if(!is.null(dim(red))) {
			map.red <- apply(pca$x, 2, function(proj) apply(red, 1, cor, y = proj, use = "pairwise.complete.obs"))
			map.green <- apply(pca$x, 2, function(proj) apply(green, 1, cor, y = proj, use = "pairwise.complete.obs"))
		} else {
			map.red <- apply(pca$x, 2, cor, y = red, use = "pairwise.complete.obs")
			map.green <- apply(pca$x, 2, cor, y = green, use = "pairwise.complete.obs")
 		}
	}

	rownames(map.red) <- paste("Probe", 1:dim(map.red)[1])
	rownames(map.green) <- paste("Probe", 1:dim(map.green)[1])

	plot.file1 <- createReportPlot(paste('ControlProbePcaPlot_Red', gsub(" ", ".", type), sep="_"), width=7, height=3, ...)
	corrHeatmap(map.red, 10, -1, margins=c(5,5))
	off(plot.file1)
	plot.file2 <- createReportPlot(paste('ControlProbePcaPlot_Green', gsub(" ", ".", type), sep="_"), width=7, height=3, ...)
    corrHeatmap(map.green, 10, -1, margins=c(5,5))
	off(plot.file2)
	return(list(plot.file1, plot.file2))
}

#######################################################################################################################

#' rnb.plot.snp.boxplot
#'
#' Box plots of beta-values from the genotyping probes
#'
#' @param rnb.set \code{\linkS4class{RnBeadSet}} object
#' @param writeToFile a flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param ... other arguments to \code{\link{createReportPlot}}
#' 
#' @export
#' @author Pavlo Lutsik
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.plot.snp.boxplot(rnb.set.example)
#' }
rnb.plot.snp.boxplot<-function(
		rnb.set, 
		writeToFile=FALSE, 
		...){
	
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!parameter.is.flag(writeToFile)) {
		stop("invalid value for writeToFile; expected TRUE or FALSE")
	}
	
	if(rnb.set@target=="probes450"){
		snp.betas<-meth(rnb.set, row.names=T)
		snp.betas<-t(snp.betas[grep("rs", rownames(snp.betas)),])
	}else if(rnb.set@target=="probes27"){
		snp.betas<-t(HM27.snp.betas(qc(rnb.set)))
	}
	rownames(snp.betas) <- samples(rnb.set)
	
	if (ncol(snp.betas) == 0) {
		stop("invalid value for rnb.set; no SNP probes found")
	}

	if(writeToFile) {
		plot.file<-createReportPlot("SNPBoxplot", ...)
	}
	snp.data<-data.frame(
		Beta.values=as.numeric(snp.betas),
		SNP=as.factor(sapply(colnames(snp.betas),rep,times=dim(snp.betas)[1])))
	plot.obj<-qplot(SNP, Beta.values, data=snp.data, geom = "boxplot", main="",
				  ylab="Beta value")+theme(axis.text.x=element_text(angle=90, vjust=0))
  	#, vp=viewport(layout.pos.row=1, layout.pos.col=1))
  
	if(writeToFile){
		print(plot.obj)
		off(plot.file)
		return(plot.file)
	}
      
	return(plot.obj)
}

#######################################################################################################################

#' rnb.plot.snp.barplot
#'
#' Bar plots of beta-values from the genotyping probes
#'
#' @param rnb.set 	\code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBeadSet}} object
#' @param sample 	unique sample identifier. In case \code{rnb.getOption("identifiers.column")} is not \code{NULL}, 
#' 					\code{sample} should attain values from the corresponding column, or \code{colnames(meth(rnb.set))}
#' 					otherwise.
#' @param writeToFile flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param numeric.names if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file name with digits
#' @param ... other arguments to \code{\link{createReportPlot}}
#'
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class 
#' 			\code{\link{ggplot}} otherwise.
#' 
#' @export
#' @author Pavlo Lutsik
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' samp<-samples(rnb.set.example)[1]
#' rnb.plot.snp.barplot(rnb.set.example, samp)
#' }
rnb.plot.snp.barplot<-function(
		rnb.set, 
		sample, 
		writeToFile=FALSE, 
		numeric.names=FALSE, 
		...){
	
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!parameter.is.flag(writeToFile)) {
		stop("invalid value for writeToFile; expected TRUE or FALSE")
	}
	ids<-samples(rnb.set)
	if (!(sample %in% ids)) {
		stop("invalid value for sample")
	}
	
	if(rnb.set@target=="probes450"){
		snp.betas<-meth(rnb.set, row.names=T)[,match(sample, ids)]
		snp.betas<-snp.betas[grep("rs", names(snp.betas))]
	}else if(rnb.set@target=="probes27"){
		snp.betas<-HM27.snp.betas(qc(rnb.set))[,match(sample, ids)]
	}
	
	if (length(snp.betas) == 0) {
		stop("invalid value for rnb.set; no SNP probes found")
	}

	if(writeToFile) {
		plot.file<-createReportPlot(paste('SNPBarPlot',  ifelse(numeric.names, match(sample, samples(rnb.set)), gsub("[ |_]", ".", sample)) , sep="_"), ...)
	}
	
	Beta.values<-as.numeric(snp.betas)
	SNP<-names(snp.betas)

	plot.obj<-qplot(SNP, Beta.values, geom = "bar", stat="identity", main=sample,
					ylab="Beta value")+scale_y_continuous(limits=c(0,1))+
			.control.plot.theme.samples#, vp=viewport(layout.pos.row=1, layout.pos.col=1))
	
	if(writeToFile) {
		print(plot.obj)
		off(plot.file)
		return(plot.file)
	}
	
	return(plot.obj)
}

#######################################################################################################################

#' rnb.plot.snp.heatmap
#'
#' Heatmap of beta-values from genotyping probes
#'
#' @param rnb.set 		\code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBeadSet}} object
#' @param writeToFile 	flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param ... 			other arguments to \code{\link{createReportPlot}}
#' 
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class 
#' 			\code{\link{ggplot}} otherwise.
#' 
#' @export
#' @author Pavlo Lutsik
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.plot.snp.heatmap(rnb.set.example)
#' }
rnb.plot.snp.heatmap<-function(
		rnb.set, 
		writeToFile=FALSE, 
		...){
	
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!parameter.is.flag(writeToFile)) {
		stop("invalid value for writeToFile; expected TRUE or FALSE")
	}
	
	sample.ids<-abbreviate.names(samples(rnb.set))
	if(rnb.set@target=="probes450"){
		snp.betas<-meth(rnb.set, row.names=T)
		snp.betas<-snp.betas[grep("rs", rownames(snp.betas)),]
	}else if(rnb.set@target=="probes27"){
		snp.betas<-HM27.snp.betas(qc(rnb.set))
	}
		
	if (length(snp.betas) == 0) {
		stop("invalid value for rnb.set; no SNP probes found")
	}
	#threshold for removing rows and columns. If more than this percentage of NAs are present
	#in a row or column, this row or column is removed prior to drawing the heatmap
	thresh.perc.na <- 0.5
	
	rnas<-colMeans(is.na(snp.betas))
	cnas<-rowMeans(is.na(snp.betas))
	
	if(all(rnas>thresh.perc.na)){
		stop("No SNP information present for any of the samples [ERROR_ID: snp_heatmap_rnas]")
	}
	if(any(rnas>thresh.perc.na)){
		snp.betas<-snp.betas[,rnas<=thresh.perc.na]
		sample.ids<-sample.ids[rnas<=thresh.perc.na]
		num.rem <- sum(rnas>thresh.perc.na)
		rnb.warning(c(num.rem,"samples were removed when generating SNP heatmap"))
	}
	if(all(cnas>thresh.perc.na)){
		stop("No SNP information present for any of the samples [ERROR_ID: snp_heatmap_cnas]")
	}
	if(any(cnas>thresh.perc.na)){
		snp.betas<-snp.betas[cnas<=thresh.perc.na,]
		num.rem <- sum(cnas>thresh.perc.na)
		rnb.warning(c(num.rem,"SNP probes were removed when generating SNP heatmap"))
	}
	
	if(writeToFile) {
		plot.file<-createReportPlot('SNPHeatmap', ...)
	}
			
	heatmap.2(snp.betas, 
			scale = "none", trace = "none", margins = c(8,8), #ColSideColors = pheno.colors$colors[,si],
			labRow = rownames(snp.betas), labCol = sample.ids, col = get.methylation.color.panel())
	 
			#lmat=matrix(c(3,4,1,2), ncol=2, byrow=T), lwid=c(0.75, 0.25), lhei=c(0.2,0.8))
	
	if(writeToFile) {
		off(plot.file)
		return(plot.file)
	}
	return(invisible(TRUE))
}

#######################################################################################################################

## point.whisker.ggplot
##
## Creates a ggplot object containing a point-and-whisker plot.
##
## @param dframe Table of values in the form of a \code{data.frame} containig at least the following columns:
##               \code{"Position"}, \code{"Average"}, \code{"Deviation"}.
## @param yrange Requested value range of the y axis.
## @param xlab Title of the x axis.
## @param ylab Title of the y axis.
##
## @author Yassen Assenov
point.whisker.ggplot <- function(dframe, yrange, xlab = "Position on the slide", ylab = "Methylation") {
	ggplot(dframe, aes_string(x = "Position", y = "Average")) +
		coord_cartesian(ylim = yrange) + scale_x_discrete(drop = FALSE) + ggplot2::geom_point(size = 3) +
		geom_errorbar(aes_string(ymin = "Average - Deviation", ymax = "Average + Deviation"), width = 0.4) +
		labs(x = xlab, y = ylab) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

#######################################################################################################################

#' rnb.plot.sentrix.distribution
#'
#' Creates a point-and-whisker plots showing beta value distributions at Sentrix positions for the given slide.
#'
#' @param rnb.set    HumanMethylation450K dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param sentrix.id Slide number (Sentrix ID) as an \code{integer} or \code{character} singleton.
#' @return Generated point-and-whisker plot (an instance of \code{\link{ggplot}}) of mean methylations for the samples
#'         on the specified slide, or \code{FALSE} if the dataset is non-empty but does not contain samples on the
#'         given slide. If the provided dataset does not contain valid Sentrix ID and position information (or is an
#'         empty dataset), this method returns \code{NULL}.
#'
#' 
#' @author Yassen Assenov
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' sid<-as.character(pheno(rnb.set.example)[["Sentrix_ID"]][1])
#' rnb.plot.sentrix.distribution(rnb.set.example,sid)
#' }
rnb.plot.sentrix.distribution <- function(rnb.set, sentrix.id) {
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!((is.character(sentrix.id) || is.integer(sentrix.id)) && length(sentrix.id) == 1 && (!is.na(sentrix.id)))) {
		stop("invalid value for sentrix.id")
	}
	dframe <- rnb.get.sentrix.data(rnb.set)
	if (is.null(dframe) || nrow(dframe) == 0) {
		return(NULL)
	}

	## Find the samples in the dataset on the specified slide
	s.indices <- which(dframe[, "Slide"] == sentrix.id) # this works both for integer and character sentrix.ids
	if (length(s.indices) == 0) {
		return(FALSE)
	}

	## Compute mean methylation and standard deviation
	m.methylations <- apply(as.matrix(meth(rnb.set)), 2, function(x) { c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)) })
	yrange <- range(m.methylations[1, ] - m.methylations[2, ], m.methylations[1, ] + m.methylations[2, ])
	m.methylations <- cbind(dframe[s.indices, ], data.frame(
			"Identifier" = colnames(meth(rnb.set))[s.indices],
			"Average" = m.methylations[1, ],
			"Deviation" = m.methylations[2, ]))
	m.methylations <- m.methylations[order(m.methylations[, "Position"]), ]

	## Create the plot
	point.whisker.ggplot(m.methylations, yrange)
}

#######################################################################################################################

#' rnb.plot.sentrix.distributions
#'
#' Creates one or more point-and-whisker plots showing beta value distributions at Sentrix positions.
#'
#' @param rnb.set HumanMethylation450K dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param fprefix File name prefix to be used in the generated plots. In order to ensure independence of the operating
#'                system, there are strong restrictions on the name of the file. See the documentation of
#'                \code{\link{createReportPlot}} for more information.
#' @param ...     Other arguments passed to \code{\link{createReportPlot}}. These can include the named parameters
#'                \code{report}, \code{width}, \code{height}, and others.
#' @return Point-and-whisker plot (an instance of \code{\linkS4class{ReportPlot}}), or a list of such plots - one per
#'         slide. If the provided dataset does not contain valid Sentrix ID and position information (or is an empty
#'         dataset), this method returns \code{NULL}.
#'
#' @details
#' If no additional parameters are specified, this function creates one PDF and one low-resolution PNG file for every
#' generated plot.
#' 
#' @seealso \code{\link{rnb.plot.sentrix.distribution}} for creating a single plot for a specified slide number
#' 
#' @author Yassen Assenov
#' @export
rnb.plot.sentrix.distributions <- function(rnb.set, fprefix = "sentrix_whisker", ...) {
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!(is.character(fprefix) && length(fprefix) == 1 && (!is.na(fprefix)))) {
		stop("invalid value for fprefix")
	}
	if (!grepl("^[A-Za-z0-9._-]+$", fprefix)) {
		stop("invalid value for fprefix")
	}
	dframe <- rnb.get.sentrix.data(rnb.set)
	if (is.null(dframe) || nrow(dframe) == 0) {
		return(NULL)
	}

	m.methylations <- apply(meth(rnb.set), 2, function(x) { c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)) })
	i.valid <- which(!is.na(m.methylations[2, ]))
	if (length(i.valid) == 0) {
		return(NULL)
	}
	m.methylations <- m.methylations[, i.valid]
	yrange <- range(m.methylations[1, ] - m.methylations[2, ], m.methylations[1, ] + m.methylations[2, ])
	m.methylations <- cbind(dframe[i.valid, ], data.frame(
			"Identifier" = samples(rnb.set)[i.valid],
			"Average" = m.methylations[1, ],
			"Deviation" = m.methylations[2, ]))
	m.methylations <- m.methylations[with(m.methylations, order(Slide, Position)), ]
	slides <- tapply(1:nrow(m.methylations), m.methylations[, "Slide"], identity)

	report.plots <- list()
	for (si in 1:length(slides)) {
		fname <- ifelse(length(slides) == 1, fprefix, paste(fprefix, si, sep = "_"))
		rplot <- createReportPlot(fname = fname, ...)
		print(point.whisker.ggplot(m.methylations[slides[[si]], ], yrange))
		report.plots[[names(slides)[si]]] <- off(rplot)
	}
	if (length(report.plots) == 1) {
		return(report.plots[[1]])
	}
	return(report.plots)
}

#######################################################################################################################

get.unified.scale<-function(intensities){
	
	ir<-range(as.numeric(intensities), na.rm=TRUE)
	
	return(c(0,ir[2]+1000))	
	
}
#######################################################################################################################
HM27.snp.betas<-function(qc.int){
	cmd<-rnb.get.annotation("controls27")
	cmd.snp<-cmd[grep("Genotyping", cmd$Type),]
	snp.ids<-cmd.snp[["Address"]]
	cy3.snp<-qc.int[["Cy3"]][rownames(qc.int[["Cy3"]]) %in% snp.ids,]
	cy5.snp<-qc.int[["Cy5"]][rownames(qc.int[["Cy5"]]) %in% snp.ids,]
	
	rownames(cy3.snp)<-cmd.snp$Name
	rownames(cy5.snp)<-cmd.snp$Name
	
	meth.cy3<-cy3.snp[paste(HM27.CY3.SNP.PROBES, "B", sep="_"),]
	umeth.cy3<-cy3.snp[paste(HM27.CY3.SNP.PROBES, "A", sep="_"),]
	
	meth.cy5<-cy5.snp[paste(HM27.CY5.SNP.PROBES, "B", sep="_"),]
	umeth.cy5<-cy5.snp[paste(HM27.CY5.SNP.PROBES, "A", sep="_"),]
	
	snp.betas.cy3<-meth.cy3/(meth.cy3+umeth.cy3+1)
	snp.betas.cy5<-meth.cy5/(meth.cy5+umeth.cy5+1)
	
	snp.betas<-rbind(snp.betas.cy3, snp.betas.cy5)
	rownames(snp.betas)<-gsub("_B", "", rownames(snp.betas))
	return(snp.betas)
}
#######################################################################################################################