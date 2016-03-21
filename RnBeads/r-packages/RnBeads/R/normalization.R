########################################################################################################################
## normalization.R
## created: 2012-05-31
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the normalization step.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.execute.normalization
#'
#' Performs normalization of the provided HumanMethylation450 data set.
#'
#' @param object Methylation dataset as an object of type \code{\linkS4class{MethyLumiSet}} or
#'               \code{\linkS4class{RnBSet}}.
#' @param method Normalization method, must be one of \code{"none"}, \code{"illumina"}, \code{"swan"},
#'               \code{"minfi.funnorm"}, \code{"bmiq"}, or \code{wm.*} where \code{*} stands for one of the methods
#'               implemented in \pkg{wateRmelon} package.
#'               Note that the execution of methods SWAN and minfi.funnorm requires packages \pkg{minfi} and
#'               \pkg{IlluminaHumanMethylation450kmanifest}. The BMIQ method requires the package \pkg{RPMM}. The
#'               \code{wm.*} methods naturally require \pkg{wateRmelon}.
#' @param bgcorr.method Character singleton specifying which background subtraction should be used. Only methods impemented
#'               in the \pkg{methylumi} package are supported at the moment, namely \code{methylumi.noob}, \code{methylumi.goob}
#'               and \code{methylumi.doob}. See Triche et al. for detailed description of the methods.
#' @param verbose flag specifying whether diagnostic output should be written to the console or to the RnBeads logger 
#' 				 in case the latter is initialized
#'
#' @return Normalized dataset as an object of type \code{\linkS4class{RnBeadSet}}. 
#'
#' @references 1. Triche, Timothy J., Jr., Weisenberger, Daniel J., Van Den Berg, David, Laird, Peter W. and Siegmund, Kimberly D. (2013)
#'             Low-level processing of Illumina Infinium DNA Methylation BeadArrays.
#'             Nucleic Acids Research 41(7):e90-e90.
#'
#' @author Pavlo Lutsik
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.set.norm<-rnb.execute.normalization(rnb.set.example, method="illumina", bgcorr.method="none")
#' }
#' @export
rnb.execute.normalization<-function(
		object, 
		method=rnb.getOption("normalization.method"),
		bgcorr.method=rnb.getOption("normalization.background.method"), 
		verbose=TRUE){
	
	if(!(inherits(object,"MethyLumiSet") || inherits(object,"RnBSet"))) {
		stop("invalid value for object; expected MethyLumiSet, RnBeadSet or RnBiseqSet")
	}
	if (!(is.character(method) && length(method) == 1 && (!is.na(method)) && (method %in% NORMALIZATION.METHODS))) {
		msg <- paste0('"', NORMALIZATION.METHODS, '"', collapse = ", ")
		stop(paste("invalid value for method; expected one of", msg))
	}

	if(inherits(object, "RnBiseqSet")){
		if (method != "none") {
			rnb.options(normalization.method="none")
			rnb.warning("Incompatible values for object and method. Changed the normalization method to \"none\"")
			method <- "none"
		}
		if(bgcorr.method!="none"){
			rnb.options(normalization.background.method="none")
			rnb.warning("Incompatible values for object and bacground correction method. Changed the background correction to \"none\"")
			bgcorr.method <- "none"
		}
	}
	
	if (inherits(object, "RnBeadSet")) {
		if (method == "bmiq") {
			if (object@status$normalized != "none") {
				## Calling BMIQ as a secondary normalization method
				rnb.error("Incompatible values for object and method. Change the normalization method to \"none\"")
			}
		}
		else if (method != "none" && !inherits(object, "RnBeadRawSet")) {
			rnb.warning(c("Incompatible values for object and method: intensity data required to perform normalization with method \"",method,"\". Changed the normalization method to \"none\""))
			rnb.options(normalization.method="none")
			method<-"none"
		}else if (method=="illumina" && is.null(qc(object))){
			warn.txt<-c("Incompatible values for object and method: quality control information is required to to perform normalization with method \"illumina\". Disabled the normalization")
			rnb.warning(warn.txt)
			
			rnb.options(normalization.method="none")
		}
	}
	
	if(inherits(object, "RnBeadSet") && object@target=="probes27" && !method %in% c("illumina")){
		rnb.warning(c("Incompatible values for object and method: \"",method,"\" cannot be applied to HumanMethylation27 data . Changed the normalization method to \"none\""))
		rnb.options(normalization.method="none")
		method<-"none"
	}
	
	if(inherits(object, "MethyLumiSet") && annotation(object)=="IlluminaHumanMethylation27" && !method %in% c("illumina")){
		rnb.warning(c("Incompatible values for object and method: \"",method,"\" cannot be applied to HumanMethylation27 data . Changed the normalization method to \"none\""))
		rnb.options(normalization.method="none")
		method<-"none"
	}
	
	if(inherits(object, "MethyLumiSet") && bgcorr.method!="none"){
		
		if(bgcorr.method=="methylumi.noob" && 
				(!"methylated.OOB" %in% ls(object@assayData) || !"unmethylated.OOB" %in% ls(object@assayData))){
			rnb.warning(c("Incompatible values for object and background correction flag: no out-of-band information found.",
					"Disabled the background correction"))
			bgcorr.method<-"none"
		}
		
		if(bgcorr.method=="methylumi.noob" && annotation(object)=="IlluminaHumanMethylation27"){
			rnb.warning(c("Incompatible methods for object and background correction method: ",
							"no background correction on the Infinium 27k data possible.",
							"Disabled the background correction"))
			bgcorr.method<-"none"
		}
			
		if(grepl("methylumi", bgcorr.method)[1]){
			bgcorr.methylumi<-gsub("methylumi\\.", "", bgcorr.method)
			pheno.columns<-colnames(phenoData(object)@data)
			suppressMessages({
						sinkfile<-ifelse("Windows" %in% Sys.info(),"NUL", "/dev/null")
						sink(sinkfile)
						object<-methylumi.bgcorr(object, method=bgcorr.methylumi)
						sink()			
					})
			#removing the introduced columns
			phenoData(object)<-phenoData(object)[,pheno.columns]
			rnb.status(c("Performed background subtraction with method", bgcorr.method))
		}
		rnb.cleanMem()
	}
	
	
	if(inherits(object, "RnBeadRawSet") && bgcorr.method!="none"){
		
		if(bgcorr.method=="methylumi.noob" && object@target=="probes27"){
			rnb.warning("Incompatible methods for object and background correction method: ]
						no background correction on the Infinium 27k data possible")
			bgcorr.method<-"none"
		}
		
		if(bgcorr.method=="methylumi.noob" && (is.null(M0(object)) || is.null(U0(object)))){
			rnb.warning(c("Incompatible values for object and background correction method: no out-of-band information found.",
					"Disabled the background correction"))
			bgcorr.method<-"none"
		}
		
		if(object@status$normalized=="swan"){
			rnb.warning(c("This RnBeadRawSet object was normalized with method SWAN: no background correction possible"))
			bgcorr.method<-"none"
		}
		
		if(grepl("methylumi", bgcorr.method)[1]){
			bgcorr.methylumi<-gsub("methylumi\\.", "", bgcorr.method)
			pheno.columns<-colnames(pheno(object))
			old.obj<-object
			object<-as(object, "MethyLumiSet")
			if(isTRUE(old.obj@status$discard.ff.matrices)){
				rnb.call.destructor(old.obj)
				rm(old.obj)
			}
			rnb.cleanMem()
			suppressMessages({
						sinkfile<-ifelse("Windows" %in% Sys.info(),"NUL", "/dev/null")
						sink(sinkfile)
						object<-methylumi.bgcorr(object, method=bgcorr.methylumi)
						sink()
					})
			#removing the introduced columns
			phenoData(object)<-phenoData(object)[,pheno.columns]
			object<-as(object, "RnBeadRawSet")
			if(object@status$disk.dump && rnb.getOption("enforce.destroy.disk.dumps")){
				object@status$discard.ff.matrices<-TRUE
			}
			rnb.status(c("Performed background subtraction with method", bgcorr.method))
		}
		object@status$background<-bgcorr.method
		rnb.cleanMem()
	}
	

	if(method=="illumina"){
						
		if(inherits(object, "RnBeadRawSet")){
			object<-as(object,"MethyLumiSet")
		}
		object<-suppressMessages(as(normalizeMethyLumiSet(object), "RnBeadRawSet"))
		rnb.cleanMem()
		
		object@status$normalized<-"illumina"
		object@status$background<-bgcorr.method

	}else if(method=="swan"){

		if(!suppressPackageStartupMessages(require("minfi"))) {
			rnb.error("Missing required package minfi for normalization method \"SWAN\"")
		}
		if(!suppressPackageStartupMessages(require("IlluminaHumanMethylation450kmanifest"))) {
			rnb.error("Missing required package IlluminaHumanMethylation450kmanifest for normalization method \"SWAN\"")
		}
		if(inherits(object,"MethyLumiSet") && (is.null(methylated(object))||is.null(unmethylated(object)))) {
			rnb.error("Invalid value for object; missing intensity information")
		}

		if(inherits(object,"MethyLumiSet")){
			intensities.by.channel<-methylumi.intensities.by.color(object)
		}else if(inherits(object,"RnBeadRawSet")){
			
			intensities.by.channel<-intensities.by.color(object, 
					add.oob=all(!is.null(M0(object)), !is.null(U0(object))),
					add.controls=!is.null(qc(object)))
		}
		
		rg.set<-RGChannelSet(intensities.by.channel$Cy3, intensities.by.channel$Cy5)
		annotation(rg.set)<-c(array="IlluminaHumanMethylation450k")
		suppressMessages({
				sinkfile<-ifelse("Windows" %in% Sys.info(),"%NULL%", "/dev/null")
				sink(sinkfile); methyl.set<-preprocessSWAN(rg.set); sink()	
		})

		meth.minfi<-getMeth(methyl.set)
		umeth.minfi<-getUnmeth(methyl.set)
				
		if(inherits(object, "MethyLumiSet")){		
			methylated(object)<-meth.minfi[match(rownames(meth.minfi), featureNames(object)),]#+methylated(object)[setdiff(featureNames(object), rownames(meth.minfi)),]
			umeth.minfi<-getUnmeth(methyl.set)
			unmethylated(object)<-umeth.minfi[match(rownames(umeth.minfi), featureNames(object)),]#+unmethylated(object)[setdiff(featureNames(object), rownames(umeth.minfi)),]
			rm(rg.set,methyl.set, meth.minfi, umeth.minfi)
			betas(object)<-rbind(methylated(object)/(methylated(object)+unmethylated(object)),
				betas(object)[setdiff(featureNames(object), rownames(methylated(object))),])
			object<-as(object, "RnBeadSet")
		}else if(inherits(object, "RnBeadRawSet")){
			probe.ids<-rownames(annotation(object))
			## rs probes are "lost" during SWAN
			
			meth.minfi<-meth.minfi[rownames(meth.minfi) %in% probe.ids,]
			umeth.minfi<-umeth.minfi[rownames(umeth.minfi) %in% probe.ids,]
			
			object@M[,][match(rownames(meth.minfi),probe.ids),]<-meth.minfi
			object@U[,][match(rownames(umeth.minfi),probe.ids),]<-umeth.minfi
	
			#update.meth(object)
			if(object@status$disk.dump){
				object@meth.sites<-convert.to.ff.matrix.tmp(beta.value(object@M[,], object@U[,]))	
			}else{
				object@meth.sites<-beta.value(object@M, object@U)
			}
			
			
			rm(rg.set,methyl.set, meth.minfi, umeth.minfi)
		}
		
		object@status$normalized<-"swan"
		object@status$background<-bgcorr.method
		rnb.cleanMem()
		
	}else if(method == "bmiq") {

		## Extract methylation value matrix and probe design information
		if (inherits(object, "MethyLumiSet")) {
			m.data <- MethyLumiSet2RnBeadSet(object)
			beta.vals <- m.data$betas
			probe.design <- rnb.annotation2data.frame(rnb.get.annotation("probes450"), add.names = TRUE)
			probe.design <- as.integer(probe.design[rownames(beta.vals), "Design"])
		} else {
			beta.vals <- object@meth.sites[,]
			probe.design <- as.integer(annotation(object)[, "Design"])
		}

		## Perform BMIQ
		samples.skipped <- integer()
		if (parallel.isEnabled()) {
			beta.names <- dimnames(beta.vals)
			beta.vals <- foreach(beta.v = as.data.frame(beta.vals), .combine = cbind, .packages = "RPMM",
					.export = c("BMIQ", "betaEst2", "blc2"),
					.noexport = c("bgcorr.method", "beta.names", "method", "object")) %dopar% {
				i <- which(!is.na(beta.v))
				if (length(i) != 0) {
					p.design <- probe.design[i]
					type1.count <- sum(p.design == 1L)
					if (all(c(type1.count, length(p.design) - type1.count) >= 50000L)) {
						beta.v[i] <- BMIQ(beta.v[i], p.design)$all
					}
				}
				beta.v
			}
			dimnames(beta.vals) <- beta.names
			rm(beta.names)
		} else {
			for (j in 1:ncol(beta.vals)) {
				i <- which(!is.na(beta.vals[, j]))
				if (length(i) != 0) {
					p.design <- probe.design[i]
					type1.count <- sum(p.design == 1L)
					if (any(c(type1.count, length(p.design) - type1.count) < 50000L)) {
						## There are not enough probes of types I and/or II
						samples.skipped <- c(samples.skipped, j)
						rnb.status(c("Skipped sample", j))
						next
					}
					beta.vals[i, j] <- BMIQ(beta.vals[i, j], p.design)$all
				}
				rnb.status(c("Normalized sample", j))
			}
			suppressWarnings(rm(j, i, p.design, type1.count))
		}
		if (length(samples.skipped) != 0) {
			## Some samples were skipped due to not enough observations
			rnb.warning(c("The following samples were not normalized:", paste(samples.skipped, collapse = ", ")))
		}
		rm(probe.design, samples.skipped)

		## Construct the resulting dataset
		if (inherits(object,"MethyLumiSet")) {
			object<-new("RnBeadSet",pheno=m.data$pheno,betas=beta.vals,p.values=m.data$p.values,bead.counts=m.data$bead.counts)
			if ("qc" %in% names(m.data)) {
				qc(object)<-m.data[["qc"]]
			}
		} else {
			if(rnb.getOption("disk.dump.big.matrices")){
				object@meth.sites <- convert.to.ff.matrix.tmp(beta.vals)	
			}else{
				object@meth.sites <- beta.vals
			}

			for (region.type in rnb.region.types.for.analysis(object@assembly)) {
				object <- summarize.regions(object, region.type)
			}
		}
		object@status$normalized<-"bmiq"
		object@status$background<-bgcorr.method

	}else if(grepl("wm\\.",method)[1]){ 
		
		if(!suppressPackageStartupMessages(require("wateRmelon"))){
			rnb.error("Missing required package wateRmelon for normalization method method")
		}
		wm.method<-gsub("wm\\.","", method)
		
		if(inherits(object, "MethyLumiSet")){
			
			object<-do.call(wm.method, list(object))
			
			object<-as(object, "RnBeadSet")	
			
		}else if(inherits(object, "RnBeadRawSet")){
			
			Mmatrix<-M(object, row.names=TRUE)
			Umatrix<-U(object, row.names=TRUE)
			
			ann.full<-rnb.annotation2data.frame(rnb.get.annotation(object@target))
			ann<-annotation(object, add.names=TRUE)
			
			if(nrow(Mmatrix)<nrow(ann.full)){
				filler<-matrix(NA_real_, nrow=nrow(ann.full)-nrow(ann), ncol=length(samples(object)))
				rownames(filler)<-rownames(ann.full[!rownames(ann.full) %in% rownames(ann),])
				ann<-rbind(ann, ann.full[!rownames(ann.full) %in% rownames(ann),colnames(ann.full)[-1]])
				Mmatrix<-rbind(Mmatrix, filler)
				Umatrix<-rbind(Umatrix, filler)
				rm(filler)
			}				
			rm(ann.full)
			
			if(wm.method=="danes"){
				betas.norm<-do.call(wm.method, list(Mmatrix, Umatrix, 
								ann[,"Design"],
								roco=pheno(object)[,grep("Sentrix[ |_]*Position", colnames(pheno(object)))]))
				
				rownames(betas.norm)<-rownames(ann)
				qcl<-qc(object)
				object<-new("RnBeadSet", 
						pheno(object), 
						betas.norm[1:nrow(meth(object)),], 
						dpval(object, row.names=TRUE), 
						covg(object, row.names=TRUE))
				qc(object)<-qcl
				
			}else{
				if(wm.method %in% c("nasen", "naten", "nanet")){
					list.norm<-do.call(wm.method, list(Mmatrix, Umatrix, 
									ann[,"Design"], ret2=TRUE))
				}else{
					list.norm<-do.call(wm.method, list(Mmatrix, Umatrix, 
									ann[,"Design"], ret2=TRUE,roco=pheno(object)[,grep("Sentrix[ |_]Position", colnames(pheno(object)))]))
				}
				
				object@M[,]<-list.norm[[1]][1:nrow(meth(object)),]
				object@U[,]<-list.norm[[2]][1:nrow(meth(object)),]
				
				rm(list.norm)
				
				#update.meth(object)
				if(object@status$disk.dump){
					object@meth.sites<-convert.to.ff.matrix.tmp(beta.value(object@M[,], object@U[,]))	
				}else{
					object@meth.sites<-beta.value(object@M, object@U)
				}
			}
		}
		object@status$normalized<-method
		object@status$background<-bgcorr.method
		rnb.cleanMem()
			
	}else if(method == "minfi.funnorm"){ 
		
		if(!suppressPackageStartupMessages(require("minfi"))) {
			rnb.error("Missing required package minfi for normalization method \"SWAN\"")
		}
		if(!suppressPackageStartupMessages(require("IlluminaHumanMethylation450kmanifest"))) {
			rnb.error("Missing required package IlluminaHumanMethylation450kmanifest for normalization method \"SWAN\"")
		}
		if(inherits(object,"MethyLumiSet") && (is.null(methylated(object))||is.null(unmethylated(object)))) {
			rnb.error("Invalid value for object; missing intensity information")
		}
		
		if(inherits(object,"MethyLumiSet")){
			
			if(nrow(methylated(object))<sum(elementLengths(rnb.get.annotation("probes450"))) || 
					sum(is.na(methylated(object)))>0 || sum(is.na(unmethylated(object)))>0){
				rnb.warning("Funtional normalization is only supported for unfiltered data sets where intensity values are present for all probes. Skipping normalization")
				object<-as(object, "RnBeadRawSet")
				object@status$normalized<-"none"
				object@status$background<-bgcorr.method
				rnb.cleanMem()
				return(object)
			}
			intensities.by.channel<-methylumi.intensities.by.color(object)

		}else if(inherits(object,"RnBeadRawSet")){

			if(nrow(sites(object))<sum(elementLengths(rnb.get.annotation("probes450"))) ||
					sum(is.na(M(object)))>0 || sum(is.na(U(object)))>0){
				rnb.warning("Funtional normalization is only supported for unfiltered data sets where intensity values are present for all probes. Skipping normalization")
				object@status$normalized<-"none"
				object@status$background<-bgcorr.method
				rnb.cleanMem()
				return(object)
			}
			intensities.by.channel<-intensities.by.color(object, 
					add.oob=all(!is.null(M0(object)), !is.null(U0(object))),
					add.controls=!is.null(qc(object)))
		}
		
		rg.set<-RGChannelSet(intensities.by.channel$Cy3, intensities.by.channel$Cy5)
		annotation(rg.set)<-c(array="IlluminaHumanMethylation450k", annotation=minfi:::.default.450k.annotation)
		suppressMessages({
					#sinkfile<-ifelse("Windows" %in% Sys.info(),"%NULL%", "/dev/null")
					#sink(sinkfile); 
					#methyl.set<-preprocessFunnorm(rg.set);
					rg.set <- updateObject(rg.set)
					gmSet <- mapToGenome(rg.set)
					extractedData <- minfi:::.extractFromRGSet450k(rg.set)
					
					gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
					sex <- rep(1L, length(gmSet$predictedSex))
					sex[gmSet$predictedSex == "F"] <- 2L
					
					rm(rg.set)
					CN <- getCN(gmSet)
					methyl.set <- minfi:::.normalizeFunnorm450k(object = gmSet, 
							extractedData = extractedData, 
							sex = NULL, nPCs = 2, verbose = 0)
					
					#sink()	
				})
		
		meth.minfi<-getMeth(methyl.set)
		umeth.minfi<-getUnmeth(methyl.set)
		
		meth.minfi<-pmax(meth.minfi, 0)
		umeth.minfi<-pmax(umeth.minfi, 0)
		
		if(inherits(object, "MethyLumiSet")){		
			
			methylated(object)<-meth.minfi[match(rownames(meth.minfi), featureNames(object)),]#+methylated(object)[setdiff(featureNames(object), rownames(meth.minfi)),]
			umeth.minfi<-getUnmeth(methyl.set)
			unmethylated(object)<-umeth.minfi[match(rownames(umeth.minfi), featureNames(object)),]#+unmethylated(object)[setdiff(featureNames(object), rownames(umeth.minfi)),]
			rm(methyl.set, meth.minfi, umeth.minfi)
			betas(object)<-rbind(methylated(object)/(methylated(object)+unmethylated(object)),
					betas(object)[setdiff(featureNames(object), rownames(methylated(object))),])
			object<-as(object, "RnBeadSet")
			
		}else if(inherits(object, "RnBeadRawSet")){
			
			probe.ids<-rownames(annotation(object, add.names=TRUE))
			## rs probes are "lost" during SWAN
			
			meth.minfi<-meth.minfi[rownames(meth.minfi) %in% probe.ids,]
			umeth.minfi<-umeth.minfi[rownames(umeth.minfi) %in% probe.ids,]
			
			object@M[,][match(rownames(meth.minfi),probe.ids),]<-meth.minfi
			object@U[,][match(rownames(umeth.minfi),probe.ids),]<-umeth.minfi
			
			#update.meth(object)
			if(object@status$disk.dump){
				object@meth.sites<-convert.to.ff.matrix.tmp(beta.value(object@M[,], object@U[,]))	
			}else{
				object@meth.sites<-beta.value(object@M, object@U)
			}
			rm(methyl.set, meth.minfi, umeth.minfi)
		}
		
		object@status$normalized<-method
		object@status$background<-bgcorr.method
		rnb.cleanMem()
	
	}else { # method == "none"
		if(inherits(object,"MethyLumiSet")){	
			object<-as(object, "RnBeadSet")
			object@status$normalized<-"none"
			object@status$background<-bgcorr.method
			
		}
		if(is.null(object@status$normalized)) {
			object@status$normalized<-"none"
			object@status$background<-bgcorr.method
		}
		
	}
	
	rnb.cleanMem()
	object
}

########################################################################################################################

#' rnb.section.normalization.shifts
#'
#' Generates a report sub-section on methylation beta values before and after correction.
#'
#' @param report       The report to contain the generated section.
#' @param betas.before Vector of raw beta values (before correction).
#' @param shifts       Vector of methylation value shifts (corrections). The \code{i}-th element in this vector must
#'                     correspond to the shift of the \code{i}-th element in \code{betas.before}.
#' @param bgcorr       Flag indicating if background subtraction was performed.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.normalization.shifts <- function(report, betas.before, shifts, bgcorr) {

	## Calculate frequencies of observed shifts and beta values
	shift.break.width <- 0.01
	shift.max <- ceiling(max(abs(shifts)) * 10) / 10
	shift.break.max <- shift.max + shift.break.width / 2
	shift.breaks <- seq(-shift.break.max, shift.break.max, by = shift.break.width)
	shift.mids <- seq(-shift.max, shift.max, length.out = length(shift.breaks) - 1)
	beta.break.width <- 0.025
	beta.breaks <- seq(0, 1, by = beta.break.width)
	beta.mids <- (beta.breaks[-length(beta.breaks)] + beta.breaks[-1]) / 2
	i.sb <- ceiling(betas.before / beta.break.width)
	i.sb[i.sb == 0] <- 1L
	i.sb <- as.integer(ceiling((shifts - shift.breaks[1]) / shift.break.width) - 1) * length(beta.mids) + i.sb
	i.sb <- tabulate(i.sb, length(shift.mids) * length(beta.mids))
	i.sb <- data.frame(
		xmin = rep(beta.breaks[-length(beta.breaks)], length(shift.mids)),
		xmax = rep(beta.breaks[-1], length(shift.mids)),
		ymin = rep(shift.breaks[-length(shift.breaks)], each = length(beta.mids)),
		ymax = rep(shift.breaks[-1], each = length(beta.mids)),
		freq = i.sb)

	## Plot histogram of shifts (differences)
	dframe <- data.frame(x = shift.mids, y = tapply(i.sb$freq, i.sb$ymin, sum))
	pp <- ggplot(dframe, aes_string(x = "x", y = "y")) + ggplot2::geom_bar(stat = "identity", width = 0.01) +
		scale_x_continuous(limits = c(-shift.break.max, shift.break.max), expand = c(0, 0)) +
		scale_y_continuous(expand = c(0, 0)) + labs(x = 'Shift', y = 'Frequency')
	rplot <- createReportPlot("correction_shifts", report, width = 7, height = 5.2)
	print(pp)
	off(rplot)
	txt <- c("The next figure gives an idea of the magnitude of the correction by showing the distribution of shifts, ",
		"i.e. degrees of modification of the raw methylation values.")
	rnb.add.paragraph(report, txt)
	txt <- "Histogram of observed magnitude of &beta; value correction."
	report <- rnb.add.figure(report, txt, rplot)
	logger.status("Added histogram of observed beta shifts (magnitude of correction)")
	rm(shift.break.width, shift.max, shift.breaks, shift.mids, beta.break.width, beta.breaks, beta.mids)
	rm(dframe, pp)

	## Plot 2D histogram of beta value (before normalization) and shift
	pp <- ggplot(i.sb, aes_string(xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax", fill = "freq")) +
		labs(x = expression(beta), y = "Shift", fill = "frequency") + geom_rect() +
		scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
		scale_fill_gradient(low = "#FFFFFF", high = "#000000") +
#		scale_fill_gradient(low = rnb.getOption("colors.gradient")[1], high = rnb.getOption("colors.gradient")[2]) +
		scale_y_continuous(limits = c(-shift.break.max, shift.break.max), expand = c(0, 0)) +
		theme(legend.justification = c(0, 0.5), legend.position = c(1, 0.5)) +
		theme(panel.border = element_blank(), panel.background = element_blank(), panel.grid = element_blank()) +
		theme(plot.margin = unit(c(0.1, 1.1, 0.1, 0.1), "in"))
	rplot <- createReportPlot("betas_shifts", report, width = 7.2, height = 6.2)
	print(pp)
	off(rplot)
	txt <- c("The figure below gives a more detailed view. This color-coded 2D histogram shows the uncorrected &beta; ",
		"values and their respective shifts after performing the normalization procedure.")
	rnb.add.paragraph(report, txt)
	txt <- "2D histogram showing the raw &beta; values and the magnitude of the corrections."
	report <- rnb.add.figure(report, txt, rplot)
	logger.status("Added 2D histogram of observed beta values and shifts")

	return(report)
}

########################################################################################################################

#' rnb.section.normalize.regions
#'
#' Adds a section summarizing the requested regions.
#'
#' @param report  Report on loading to contain the newly constructed section.
#' @param rnb.set Methylation dataset to be analyzed.
#' @param regions Non-empty \code{character} vector of region names.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.normalize.regions <- function(report, rnb.set, regions) {
	msg <- function(txt, e.class = "disabled") {
		paste0('<span class="', e.class, '">', txt, '</span>')
	}
	reg.counts <- sapply(regions, function(region) {
			tryCatch(as.character(nrow(meth(rnb.set, region))), error = function(e) { msg("not supported") })
		}
	)
	get.r.description <- function(reg) {
		result <- attr(rnb.get.annotation(reg, assembly = assembly(rnb.set)), "description")
		if (is.null(result) || identical("", result)) {
			result <- msg("n.a.")
		} else {
			result <- paste0('<p style="font-weight:normal;text-align:left">', result, '</p>')
		}
		result
	}
	reg.descriptions <- sapply(regions, function(region) {
			tryCatch(get.r.description(region), error = function(e) { msg("not imported", "outdated") })
		}
	)
	table.statistics <- data.frame(
		"Annotation" = regions, "Description" = reg.descriptions, "Regions in the Dataset" = reg.counts,
		check.names = FALSE, stringsAsFactors = FALSE)
	txt <- c(ifelse(rnb.getOption("analyze.sites"), "In addition to CpG sites, there", "There"),
		ifelse(length(regions) == 1, " is one set", paste(" are", length(regions), "sets")),
		" of genomic regions to be covered in the analysis. The table below gives a summary of these annotations.")
	report <- rnb.add.section(report, "Region Annotations", txt, level = 2)
	table.header <- c("<colgroup>", paste0('\t<col width="', c(210, 420, 150), 'px" />'), "</colgroup>")
	rnb.add.table(report, table.statistics, row.names = FALSE, thead = table.header)
	
	return(report)
}

#######################################################################################################################

#' rnb.section.normalization
#'
#' Performs normalization of the loaded methylation dataset.
#'
#' @param report    Report to contain the normalization section. This must be an object of type
#'                  \code{\linkS4class{Report}}.
#' @param rnb.set   Normalized methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param betas.raw Matrix of methylation beta values for CpG sites, as they would appear without (prior to)
#'                  normalization. If this parameter is \code{NULL}, the effect of normalization on beta values is
#'                  not examined and visualized.
#' @return The modified report.
#'
#' @details
#' If specified \code{betas.raw}, must be a matrix of methylation values with dimensions identical to the one
#' encapsulated in \code{rnb.set}. The former matrix is expected to store beta values before normalization,
#' whereas the latter one - after normalization. If the provided matrices contain sufficient amounts of non-missing
#' values, this function creates figures that compare the two distributions of values and examines the magnitute of
#' modifications.   
#' 
#' @author Pavlo Lutsik
#' @noRd
rnb.section.normalization <- function(report, rnb.set, betas.raw = NULL) {
	if (!inherits(report, "Report")) {
		stop("Invalid value for report")
	}
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!is.null(betas.raw)) {
		if (!(is.matrix(betas.raw) && (is.double(betas.raw) || is.integer(betas.raw)))) {
			stop("invalid value for betas.raw")
		}
		betas.after <- meth(rnb.set, row.names = TRUE)
		if (!(nrow(betas.raw) == nrow(betas.after) && ncol(betas.raw) == ncol(betas.after) &&
			setequal(rownames(betas.raw), rownames(betas.after)) &&
			setequal(colnames(betas.raw), colnames(betas.after)))) {
			stop("incompatible values for rnb.set and betas.raw")
		}
		betas.after <- betas.after[rownames(betas.raw), colnames(betas.raw)]
	}

	txt.methylumi <- "<a href=\"http://www.bioconductor.org/packages/release/bioc/html/methylumi.html\">methylumi</a>"

	bgcorr.method <- rnb.set@status$background
	bgcorr <- (bgcorr.method != "none")
	bgcorr.method <- gsub("methylumi\\.","", bgcorr.method)
	if (bgcorr) {
		refText <- c("Triche, T.J. Jr, Weisenberger, D.J., Van Den Berg, D., Laird, P.W., and Siegmund, K.D. (2013)",
			"Low-level processing of Illumina Infinium DNA Methylation BeadArrays. ",
			"<i>Nucleic Acids Research</i> <b>41</b>(7), e90")
		report <- rnb.add.reference(report, refText)
		txt <- c("The background was subtracted using the ", txt.methylumi, " package (method \"",
			bgcorr.method, "\") ", rnb.get.reference(report, refText), ".")
		txt.methylumi <- "methylumi"
	} else {
		txt <- ""
	}
	method <- rnb.set@status$normalized
	if (method == "illumina") {
		txt <- c(txt, "The signal intensity values were normalized using Illumina's recommended normalization method, ",
			"as implemented in the ", txt.methylumi, " package.")
	} else if (method == "swan") {
		txt <- c(txt, "The signal intensity values were normalized using the SWAN normalization method, as ",
			"implemented in the ",
			"<a href=\"http://www.bioconductor.org/packages/release/bioc/html/minfi.html\">minfi</a> package.")
	} else if (method == "bmiq") {
		refText <- c("Teschendorff, A. E., Marabita, F., Lechner, M., Bartlett, T., Tegner, J., Gomez-Cabrero, D. and ",
			"Beck, S. (2013) A beta-mixture quantile normalization method for correcting probe design bias in ",
			"Illumina Infinium 450 k DNA methylation data. <i>Bioinformatics</i> <b>29</b>(2), 189-196")
		report <- rnb.add.reference(report, refText)
		txt <- c(txt, "The methylation &beta; values were normalized using the BMIQ normalization method ",
			rnb.get.reference(report, refText), ".")
	} else if (grepl("wm\\.", method)){
		refText <- c("Pidsley, R., Wong, C.,  Volta, M., Lunnon, K., Mill, J., and Schalkwyk, L. (2013)",
				"A data-driven approach to preprocessing Illumina 450K methylation array data."," <i>BMC Genomics</i>, <b>14<b>(1), 293")
		report <- rnb.add.reference(report, refText)
		wm.method<-gsub("wm\\.", "", method)
		txt<-c(txt, sprintf("The data was normalized using method %s from %s.", wm.method, rnb.get.reference(report,refText)))
	}else if (method == "minfi.funnorm"){ 
		refText<-c("Fortin, J., Labbe, A., Lemire, M., Zanke, B. W., Hudson, T. J., Fertig, E. J., Greenwood, C.M.T., ",
			"Hansen, K. D. (2014) Functional normalization of 450k methylation array data improves replication in ",
			"large cancer studies. <i>BioRxiv</i>.")
		report <- rnb.add.reference(report, refText)
		txt<-c(txt, sprintf("The data was normalized using the functional normalization method from %s.", rnb.get.reference(report,refText)))
		
	}else{ # method == "none"
		txt <- c(txt, "The measurements in this dataset were not normalized after ",
			ifelse(txt == "", "loading", "background subtraction"), ".")
		if(rnb.getOption("normalization.method")=="minfi.funnorm"){
			if(nrow(sites(rnb.set))<sum(elementLengths(rnb.get.annotation("probes450")))){
				txt<-c(txt, "The desired normalization method \"minfi.funnorm\" was not applied since the data set did not 
					contain information for all probes. This is most likely because of the preceding quality filtering steps.", 
					"In order to apply minfi.funnorm in the full RnBeads analysis you should disable SNP filtering and Greedycut.",
					"Otherwise, you apply rnb.execute.normalization() function to an unfiltered RnBead(Raw)Set after loading and start
					the full pipeline with the returned object as input.")
			}else if(sum(is.na(M(rnb.set)))>0 || sum(is.na(U(rnb.set)))>0){
				txt<-c(txt, "The desired normalization method \"minfi.funnorm\" was not applied because intensity information contains 
						missing values.", "Try to apply rnb.execute.normalization() function to an unfiltered RnBead(Raw)Set after loading and start
								the full pipeline with the returned object as input.")
			}
		}
	}
	report <- rnb.add.section(report, "Normalization", txt)

	## Add comparison of betas before and after correction
	if ((!is.null(betas.raw)) && (bgcorr || method != "none")) {
		betas.raw <- as.vector(betas.raw)
		betas.after <- as.vector(betas.after)
		min.observations <- 501L

		i.before <- !is.na(betas.raw); i.after <- !is.na(betas.after)
		i.both <- which(i.before & i.after)
		i.before <- which(i.before); i.after <- which(i.after)

		section.title <- "Effect of Correction"
		if (length(i.before) >= min.observations && length(i.after) >= min.observations) {
			txt <- c("This section shows the influence of the applied normalization procedure on CpG methylation ",
				"values. The following figure compares the distributions of the &beta; values before and after ",
				"performing normalization.")
			report <- rnb.add.section(report, section.title, txt, level = 2)

			## Distributions before and after normalization
			beta.values <- list("Before correction" = betas.raw[i.before], "After correction" = betas.after[i.after])
			report.plots <- rnb.plot.beta.comparison(beta.values, "correction_comparison", report, min.observations)
			txt <- "Comparison of &beta; values before and after correction."
			txt <- c(txt, add.text.subsampling(attr(report.plots, "subsampled"), paste("betas", names(beta.values))))
			setting.names <- list("Plot type" =
				c("density" = "density estimation", "histogram" = "histograms", "qq" = "quantile-quantile plot"))
			report <- rnb.add.figure(report, txt, report.plots, setting.names)
			rm(i.before, i.after, beta.values, report.plots, setting.names)
			rnb.cleanMem()
			logger.status("Added comparison between non-normalized and normalized beta values")

			if (rnb.getOption("normalization.plot.shifts")) {
				if (length(i.both) != 0) {
					betas.raw <- betas.raw[i.both]
					betas.after <- betas.after[i.both] - betas.raw
					report <- rnb.section.normalization.shifts(report, betas.raw, betas.after, bgcorr)
				} else {
					txt <- c("No plots on &beta; value shifts were created because there are not enough data to ",
						"visualize.")
					rnb.add.paragraph(report, txt)
				}
			}
		} else {
			txt <- c("The effects of the applied procedure on CpG methylation values is not summarized because there ",
				"are not enough &beta; values to accurately present the observed distributions before and after ",
				"normalization.")
			report <- rnb.add.section(report, section.title, txt)
		}
		rm(betas.raw, betas.after)
	}

	if (!inherits(rnb.set, "RnBeadSet")) {
		return(report)
	}

	## Add point-and-whisker plots of mean methylation for every slide
	report.plots <- rnb.plot.sentrix.distributions(rnb.set, report = report, width = 6, height = 6)
	if (is.null(report.plots)) {
		txt <- "Sample average methylation cannot be visualized because "
		if (length(samples(rnb.set)) == 0) {
			txt <- c(txt, "the normalized dataset is empty.")
		} else if (all(is.na(meth(rnb.set)))) {
			txt <- c(txt, "the dataset contains no valid methylation beta values.")
		} else {
			txt <- c(txt, "no valid Sentrix ID and Sentrix Position information could be extracted from the ",
				"sample annotation table.")
		}
	} else {
		txt <- "The following figure visualizes the average methylation per sample. Samples are grouped by slide."
	}
	report <- rnb.add.section(report, "Sample Mean Methylations", txt, level = 2)
	if (!is.null(report.plots)) {
		description <- "Point-and-whisker plot showing mean and standard deviation among all beta values in a sample."
		if (inherits(report.plots, "ReportPlot")) {
			setting.names <- list()
		} else { # is.list(report.plots)
			setting.names <- list("Slide number" = names(report.plots))
			names(setting.names[[1]]) <- 1:length(report.plots)
		}
		report <- rnb.add.figure(report, description, report.plots, setting.names)
	}

	r.types <- rnb.getOption("region.types")
	if (is.null(r.types)) {
		r.types <- summarized.regions(rnb.set)
	}
	if (length(r.types) != 0) {
		report <- rnb.section.normalize.regions(report, rnb.set, r.types)
	}

	return(report)
}

#######################################################################################################################

#' rnb.step.normalization
#'
#' Performs normalization of the loaded HumanMethylation450 dataset and adds a corresponding section to the given
#' report.
#'
#' @param object Methylation dataset as an object of type inheriting \code{\linkS4class{MethyLumiSet}} or
#'               \code{\linkS4class{RnBSet}}.
#' @param report Report to contain the normalization section. This must be an object of type
#'               \code{\linkS4class{Report}}.
#' @param method Normalization method, must be one of \code{"illumina"}, \code{"swan"} or \code{"none"}. Note that
#'               the execution of method SWAN requires packages \pkg{minfi} and
#'               \pkg{IlluminaHumanMethylation450kmanifest}.
#'
#' @return List with two elements:
#'         \describe{
#'           \item{\code{dataset}}{Normalized dataset as an object of type \code{\linkS4class{RnBeadSet}}.}
#'           \item{\code{report}}{the modified report.}
#'         }
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.step.normalization<-function(object, report, method = rnb.getOption("normalization.method")){
	
	if(!(inherits(object,"MethyLumiSet") || inherits(object,"RnBSet"))) {
		stop("Invalid value for object; expected MethyLumiSet, RnBeadSet or RnBiseqSet")
	}
	if (!inherits(report, "Report")) {
		stop("Invalid value for report")
	}
	if (!(is.character(method) && length(method) == 1 && (!is.na(method)) &&
			(method %in% NORMALIZATION.METHODS))) {
		msg <- paste0('"', NORMALIZATION.METHODS, '"', collapse = ", ")
		stop(paste("invalid value for method; expected one of", msg))
	}

	logger.start("Normalization Procedure") 


	if (method == "none") {
		betas.raw <- NULL
	} else if (inherits(object, "MethyLumiSet")) {
		betas.raw <- betas(object)
		pheno.table <- phenoData(object)@data
		id.column <- rnb.getOption("identifiers.column")
		if (!(is.null(pheno.table) || is.null(id.column))) {
			ids <- NA
			if (is.character(id.column) && (id.column %in% colnames(pheno.table))) {
				ids <- pheno.table[, id.column]
			} else if (1 <= id.column && id.column <= ncol(pheno.table)) {
				ids <- pheno.table[, id.column]
			}
			if (any(is.na(ids)) == FALSE && anyDuplicated(ids) == 0) {
				colnames(betas.raw) <- as.character(ids)
			}
			rm(ids)
		}
		rm(pheno.table, id.column)
		rnb.cleanMem()
	} else { # inherits(object, "RnBSet")
		betas.raw <- meth(object, row.names = TRUE)
	}
	#only destroy an object if the new object is not the same as the old object,
	#i.e. if normalization or background subtraction took place
	#if ((rnb.set@status$normalized != "none") || (rnb.set@status$background != "none")){
	#	rnb.call.destructor(object)
	#}

	if(object@status$disk.dump && rnb.getOption("enforce.destroy.disk.dumps")){
		object@status$discard.ff.matrices<-TRUE
	}
	object <- rnb.execute.normalization(object, method)
	rnb.cleanMem()
	if(isTRUE(object@status$discard.ff.matrices)){
		object@status$discard.ff.matrices<-NULL
	}
#	if(object@status$background!="none"){
#		logger.status(sprintf("Performed background correction with method \"%s\"", object@status$background))
#	}
	if(object@status$normalized!="none"){
		logger.status(sprintf("Performed normalization with method \"%s\"", object@status$normalized))
	}
	
	report<-rnb.section.normalization(report, object, betas.raw)
	logger.status("Added normalization section")

	logger.completed()
	return(list(dataset=object, report=report))
}
