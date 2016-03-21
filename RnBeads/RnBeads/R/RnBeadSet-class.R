########################################################################################################################
## RnBeadSet-class.R
## created: 2012-04-06
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## RnBeadSet class definition.
########################################################################################################################
## GLOBALS

RNBEADSET.SLOTNAMES<-c("pval.sites")

## ---------------------------------------------------------------------------------------------------------------------
## CLASS DEFINITIONS
## ---------------------------------------------------------------------------------------------------------------------

#' RnBeadSet Class
#'
#' Stores the preprocessed information from HumanMethylation experiments
#'
#' @details
#' There are multiple ways to create an object of type \code{RnBeadSet}:
#' \describe{
#'  \item{Loading from files}{Dataset can be loaded from text or binary files. See the function
#'       \code{\link{rnb.execute.import}} for more details.}
#'  \item{Downloading from GEO}{See the function \code{\link{read.geo}} for details.}
#'  \item{Converting from \code{MethyLumiSet}}{...}
#' } 
#'
#' @section Slots:
#' \describe{
#'   \item{\code{pval.sites}}{\code{matrix} of detection p-values with the same dimensions as \code{betas}, or
#'        \code{NULL} if the detection p-values are not available.}
#'   \item{\code{pval.regions}}{\code{list} of methylation \code{matrix} objects, one per available region type. Every row in a 
#' 		matrix corresponds to a methylation site, and every column - to a sample.}
#'   \item{\code{covg.sites}}{\code{matrix} of bead counts per probe with the same dimensions as \code{betas}, or
#'        \code{NULL} if this data are not available.}
#'   \item{\code{qc}}{Quality control probe information in the form of a \code{list} of two elements - \code{"Cy3"} and
#'        \code{"Cy5"}, storing intensities of probes on the green and red channels, respectively. This slot's value is
#'        \code{NULL} if no control probe information is available.}
#' }
#'
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{\link[=samples,RnBSet-method]{samples}}}{Gets the identifiers of all samples in the dataset.}
#'   \item{\code{\link[=pheno,RnBSet-method]{pheno}}}{Gets the phenotypic and processing data of the dataset.}
#'   \item{\code{\link[=meth,RnBSet-method]{meth}}}{Gets the \code{matrix} of methylation beta-values of the dataset.}
#'   \item{\code{\link[=dpval,RnBeadSet-method]{dpval}}}{Gets the \code{matrix} of detection p-values of the dataset.}
#'   \item{\code{\link[=covg,RnBSet-method]{covg}}}{Gets the \code{matrix} of bead counts of the dataset.}
#'   \item{\code{\link[=qc,RnBeadSet-method]{qc}}}{Gets the intensities of the quality control probes.}
#'   \item{\code{\link[=remove.sites,RnBSet-method]{remove.sites}}}{Removes probes from the dataset.}
#'   \item{\code{\link[=remove.samples,RnBeadSet-method]{remove.samples}}}{Removes samples from the dataset.}
#'   \item{\code{\link[BiocGenerics]{combine}}}{Combines two datasets.}
#' }
#'
#' @name RnBeadSet-class
#' @rdname RnBeadSet-class
#' @author Pavlo Lutsik
#' @exportClass RnBeadSet
setClass("RnBeadSet",
         representation(pval.sites="matrixOrffOrNULL",
				 		pval.regions="listOrNULL",
				 		qc="listOrNULL",
						status="listOrNULL"
		#				bead.counts="listOrNULL"
		),
		 contains="RnBSet",
         prototype(pheno=data.frame(),
				 meth.sites=matrix(),
				 pval.sites=NULL,
				 pval.regions=NULL,
				 qc=NULL, 
				 status=NULL
#				 bead.counts=NULL
 		 ), package = "RnBeads"
)

## ---------------------------------------------------------------------------------------------------------------------
## VALIDITY
## ---------------------------------------------------------------------------------------------------------------------

##
## validRnBeadSetObject
##
## Valididty check for RnBeadSet-class
##
validRnBeadSetObject<-function(object){

  if(is.null(object@pheno) || is.null(object@meth.sites)) return("Pheno and betas slots can not be NULL")
  if(dim(object@pheno)[1L]!=dim(object@meth.sites)[2L]) return("The number of samples is ambiguous")

  if(!is.null(object@pval.sites)){
    if(!setequal(rownames(object@pval.sites), rownames(object@meth.sites))) return("Rows of the beta and p-value tables do not match")
    if(!setequal(colnames(object@pval.sites), colnames(object@meth.sites))) return("Columns of the beta and p-value tables do not match")

  }

  if(!is.null(object@covg.sites)){
    if(!setequal(rownames(object@covg.sites), rownames(object@meth.sites))) return("Rows of the beta and bead count tables do not match")
    if(!setequal(colnames(object@covg.sites), colnames(object@meth.sites))) return("Columns of the beta and bead count tables do not match")
  }

  return(TRUE)
}

setValidity("RnBeadSet", method=validRnBeadSetObject)

## ---------------------------------------------------------------------------------------------------------------------
## CONSTRUCTORS
## ---------------------------------------------------------------------------------------------------------------------

#' initialize.RnBeadSet
#'
#' Initializes \code{RnBeadSet} object with data from supplied data matrices.
#'
#' @param pheno       phenotypic data.
#' @param betas       matrix of beta values. Should contain Infinium probe identifiers as row names.
#' @param p.values    matrix of detection p-values.
#' @param bead.counts matrix of bead counts per probe.
#' @param platform	  \code{character} singleton specifying the microarray platform: \code{"450k"} corresponds to HumanMethylation450 microarray, and \code{"27k"} stands for HumanMethylation27.
#' @param region.types	a \code{character} vector specifying the region types, for which the methylation infromation will be summarized.
#' @param useff		  if \code{TRUE} the data matrices will be stored as \code{ff} objects 
#'
#' @details Row names of the \code{betas} matrix are used to match the data to annotation. All matrices 
#' should have identical dimensions.
#' 
#' @export
#' 
#' @docType methods
#' @rdname RnBeadSet-class
setMethod("initialize", "RnBeadSet",
	function(.Object, 
			pheno, 
			betas, 
			p.values=NULL, 
			bead.counts = NULL, 
			platform = "450k",
			region.types = rnb.region.types.for.analysis("hg19"),
			useff=rnb.getOption("disk.dump.big.matrices")) {
		
		## A hack necessary for correct processing of classes 
		## that inherit from RnBeadSet
		##
		if(missing(pheno) || missing(betas)){
			#warning("Valid RnBeadSet object was not created: pheno and/or betas missing")
			return(.Object)
		}
		
		if (platform == "450k") {
			.Object@target <- "probes450"
		} else if(platform == "27k"){ 
			.Object@target <- "probes27"
		}else{
			stop("Invalid value for platform")
		}

		## Load probe annotation table
		probe.annotation <- rnb.get.annotation(.Object@target, .Object@assembly)

		## Construct the matrix of site indices
		#x.data<-GenomicRanges::as.data.frame(probe.annotation)
		
		chrs<-rnb.get.chromosomes(.Object@assembly)
		
		annotated.probes<-do.call("c", lapply(probe.annotation, names))
		if(length(which(rownames(betas) %in% annotated.probes))<2){
			err<-"Annotations could be found for less than two rows from the supplied beta value table"
			rnb.error(err)		
		}
						
		if(!all(rownames(betas) %in% annotated.probes)){
			warn<-"Some of the supplied probes are missing annotation and will be discarded"
			rnb.warning(warn)			
		}
		
		x.data<-rnb.annotation2data.frame(probe.annotation)
		rownames(x.data)<-annotated.probes
		rm(annotated.probes)
		
		site.ids<-character()
		p.indices <- lapply(unique(x.data[["Chromosome"]]), function(chr) {
					chr.ids<-rownames(x.data)[x.data[["Chromosome"]]==chr]
					table<-match(rownames(betas), chr.ids)
					table<-sort(na.omit(table))
					site.ids<<-c(site.ids, chr.ids[table])
					table
				})
		
		p.infos<-matrix(c(as.integer(x.data[site.ids, "Context"]), match(x.data[site.ids,"Chromosome"], names(chrs)), unlist(p.indices)),
				ncol = 3, dimnames = list(site.ids, c("context", "chr", "index")))
		
		.Object@sites <- p.infos
		
		## Assign tables of methylation and other values
		.Object@pheno <- pheno
		if(useff){
			.Object@meth.sites <- convert.to.ff.matrix.tmp(betas[site.ids, , drop = FALSE])
			.Object@status<-list(disk.dump=TRUE)
		}else{
			.Object@meth.sites <- betas[site.ids, , drop = FALSE]
			.Object@status<-list(disk.dump=FALSE)
		}
		if (!is.null(p.values)) {
			if(useff){
				.Object@pval.sites <- convert.to.ff.matrix.tmp(p.values[site.ids, , drop = FALSE])
			}else{
				.Object@pval.sites <- p.values[site.ids, , drop = FALSE]
			}
		}
		if (!is.null(bead.counts)) {
			if(useff){
				.Object@covg.sites <- convert.to.ff.matrix.tmp(bead.counts[site.ids, , drop = FALSE])
			}else{
				.Object@covg.sites <- bead.counts[site.ids, , drop = FALSE]
			}
		}

		for (region.type in region.types) {
			if (region.type %in% rnb.region.types("hg19")) {
				.Object <- summarize.regions(.Object, region.type)
			}
		}

		.Object@status[["normalized"]]<-"none"
		.Object@status[["background"]]<-"none"

		.Object@inferred.covariates <- list()
		
		.Object
	}
)

########################################################################################################################

## Prints a summary of a RnBeadSet dataset. This is used in the \code{show} methods.
## 
## @param object The RnBeadSet object to be presented.
##
rnb.show.rnbeadset <- function(object) {
	cat("Object of class ", class(object), "\n", sep = "")
	cat(sprintf("%8d samples\n", nrow(object@pheno)))
	cat(sprintf("%8d probes\n", nrow(object@sites)))
	probe.types <- rownames(object@sites)
	if (!is.null(probe.types)) {
		probe.types <- sapply(c("^cg", "^ch", "^rs"), function(type) { sum(grepl(type, probe.types)) })
		cat(sprintf("\tof which: %g CpG, %g CpH, and %g rs\n", probe.types[1], probe.types[2], probe.types[3]))
	}
	if (!is.null(object@regions)) {
		cat("Region types:\n")
		for (rn in names(object@regions)) {
			cat(sprintf("\t%8d regions of type %s\n", nrow(object@regions[[rn]]), rn))
		}
	}
	cat(sprintf("Intensity information is %s\n", ifelse(inherits(object, "RnBeadRawSet"), "present", "absent")))
	cat(sprintf("Detection p-values are %s\n", ifelse(is.null(object@pval.sites), "absent", "present")))
	cat(sprintf("Bead counts are %s\n", ifelse(is.null(object@covg.sites), "absent", "present")))
	cat(sprintf("Quality control information is %s\n", ifelse(is.null(qc(object)), "absent", "present")))
	cat(sprintf("Summary of normalization procedures:\n"))
	
	if(!is.null(object@status) && !is.null(object@status$normalized)){
		if(object@status$normalized=="none"){
			cat(sprintf("\tThe methylation data was not normalized.\n"))
		}else{
			cat(sprintf("\tThe methylation data was normalized with method %s.\n", object@status$normalized))			
		}
	}
	
	if(!is.null(object@status) && !is.null(object@status$background)){
		if(object@status$background=="none"){
			cat(sprintf("\tNo background correction was performed.\n"))
		}else{
			cat(sprintf("\tBackground correction was performed with method %s.\n", object@status$background))			
		}
	}
}

setMethod("show", "RnBeadSet", rnb.show.rnbeadset)

########################################################################################################################
## MethyLumi2RnBeadSet
##
## Extracts all methylation data from a MethyLumiSet object, necessary to instantiate RnBeadSet
##
## @param methySet \code{MethylumiSet} instance to be converted.
## @return \code{list} with up to 5 elements: \code{"pheno"}, \code{"betas"}, \code{"bead.counts"}, \code{"p.values"},
##         \code{"qc"}
##
## @author Pavlo Lutsik
MethyLumiSet2RnBeadSet <- function(methySet) {
	result <- list()
	result[["pheno"]] <- as(phenoData(methySet), "data.frame")
	beta.vals <- betas(methySet)

	if(all(c("methylated.N", "unmethylated.N") %in% assayDataElementNames(methySet))){

		## Set beta values with low coverage (bead.counts) to NA
		result[["bead.counts"]]<-methylated.N(methySet)+unmethylated.N(methySet)
		beta.vals[methylated.N(methySet)<2]<-NA
		beta.vals[unmethylated.N(methySet)<2]<-NA
		
	}
	result[["betas"]] <- beta.vals
	result[["p.values"]] <- pvals(methySet)

	if(!is.null(controlData(methySet))){

		## Extract quality control data
		cd<-controlData(methySet)
		green<-intensitiesByChannel(cd,"Cy3")
		red<-intensitiesByChannel(cd,"Cy5")

		if(length(grep('[A-z]', rownames(green)))==dim(green)[1L]){

			if("ProbeID" %in% varLabels(featureData(cd)))
				probe.id.col<-"ProbeID" else if ("Address" %in% varLabels(featureData(cd))) probe.id.col<-"Address"

			probe.mapping<-as(featureData(cd)[,probe.id.col], "data.frame")
			rownames(green)<-probe.mapping[rownames(green), probe.id.col]
			rownames(red)<-probe.mapping[rownames(red), probe.id.col]
			
		}

		result[["qc"]]<-list(Cy3=green, Cy5=red)
	}
	result
}

#' as("RnBeadSet", "MethyLumiSet")
#'
#' Convert a \code{\linkS4class{RnBeadSet}} object to \code{\linkS4class{MethyLumiSet}}
#' 
#' @name coercion-methods
#' 
setAs("MethyLumiSet", "RnBeadSet",

	function(from, to){

		if(!inherits(from,"MethyLumiSet")){
			stop("not a MethyLumiSet object:", deparse(substitute(methylumi.set)))
		}

		m.data <- MethyLumiSet2RnBeadSet(from)
		object<-new(to,pheno=m.data$pheno,betas=m.data$betas,p.values=m.data$p.values,bead.counts=m.data$bead.counts)
		if ("qc" %in% names(m.data)) {
			qc(object)<-m.data[["qc"]]
		}

		object
})

## ---------------------------------------------------------------------------------------------------------------------
## GETTERS
## ---------------------------------------------------------------------------------------------------------------------

if(!isGeneric("dpval")) setGeneric('dpval',
           function(object, ...) standardGeneric('dpval'))

#' dpval-methods
#'
#' Extract detection p-values from an object of \code{\linkS4class{RnBeadSet}} class.
#'
#' @param object 		\code{\linkS4class{RnBeadSet}} or \code{\linkS4class{RnBeadRawSet}} object
#' @param type 			\code{character} singleton. If \code{sites} detection p-values per each available 
#' 						site is returned. Otherwise should be one of region types for for which the summarized 
#' 						p-values are available
#' @param row.names	    Flag indicating of row names are to be generated in the result.

#' 
#' @return detection p-values available for the dataset in the form of a \code{matrix}.
#'
#' @rdname dpval-methods
#' @docType methods
#' @export
#' @aliases dpval
#' @aliases dpval,RnBeadSet-method
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' dp<-dpval(rnb.set.example, row.names=TRUE)
#' head(dp)
#' }
setMethod("dpval", signature(object="RnBeadSet"),
          function(object, type="sites", row.names=FALSE){
			  get.dataset.matrix(object, type, row.names, object@pval.sites, object@meth.regions)
          })

########################################################################################################################
  
setGeneric("qc", function(object) standardGeneric("qc"))

#' qc-methods
#'
#' Extracts HumanMethylation quality control information
#'
#' @param object Dataset of interest.
#' @return Quality control information available for the dataset in the form of a \code{list} with two elements:
#' \code{Cy3} and \code{Cy5}.
#'
#' @rdname qc-methods
#' @aliases qc
#' @aliases qc,RnBeadSet-method
#' @docType methods
#' @export 
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' qcinf<-dpval(rnb.set.example, row.names=TRUE)
#' head(qcinf$Cy3)
#' head(qcinf$Cy5)
#' }
setMethod("qc", signature(object="RnBeadSet"), function(object){ object@qc })

setGeneric("qc<-", function(object, value) standardGeneric("qc<-"))
setMethod("qc<-", signature(object = "RnBeadSet", value = "listOrNULL"),
	function(object, value) {
		## TODO: Add validation and export, or simply use @qc in our code
		object@qc <- value
		object
	}
)

## ---------------------------------------------------------------------------------------------------------------------
## MODIFIERS
## ---------------------------------------------------------------------------------------------------------------------

if (!isGeneric("remove.sites")) {
	setGeneric("remove.sites", function(object, probelist, verbose = TRUE) standardGeneric("remove.sites"))
}

#' @rdname remove.sites-methods
#' @aliases remove.sites,RnBeadSet-method
#' @docType methods
#' @export
setMethod("remove.sites", signature(object = "RnBeadSet"),
	function(object, probelist, verbose = TRUE) {
		inds <- get.i.vector(probelist, rownames(object@meth.sites))
		if (length(inds) != 0) {
			if (!is.null(object@pval.sites)) {
				if(object@status$disk.dump){
					new.matrix<-object@pval.sites[-inds,]
					if(isTRUE(object@status$discard.ff.matrices)){
						delete(object@pval.sites)
					}
					object@pval.sites <- convert.to.ff.matrix.tmp(new.matrix)
					rm(new.matrix); rnb.cleanMem()
				}else{
					object@pval.sites <- object@pval.sites[-inds, ]
				}
			}
		}

		callNextMethod()
	}
)

########################################################################################################################

if (!isGeneric("remove.samples")) {
	setGeneric("remove.samples", function(object, samplelist) standardGeneric("remove.samples"))
}

#' @rdname remove.samples-methods
#' @aliases remove.samples,RnBeadSet-method
#' @docType methods
#' @export
setMethod("remove.samples", signature(object = "RnBeadSet"),
	function(object, samplelist) {
		inds <- get.i.vector(samplelist, samples(object))
		if (length(inds) != 0) {
			if (!is.null(object@pval.sites)) {
				if(object@status$disk.dump){
					new.matrix<-object@pval.sites[,-inds, drop=F]
					if(isTRUE(object@status$discard.ff.matrices)){
						delete(object@pval.sites)
					}
					object@pval.sites <- convert.to.ff.matrix.tmp(new.matrix)
					rm(new.matrix); rnb.cleanMem()
				}else{
					object@pval.sites <- object@pval.sites[,-inds, drop=F]
				}
			}
			if (!is.null(object@qc)) {
				object@qc$Cy3 <- object@qc$Cy3[,-inds, drop=F]
				object@qc$Cy5 <- object@qc$Cy5[,-inds, drop=F]
			}
		}
		callNextMethod()
	}
)

########################################################################################################################

setMethod("save.matrices", signature(object="RnBeadSet", path="character"),
		function(object, path){
			
			if(!file.exists(path)){
				dir.create(path)
			}
			if(!is.null(object@pval.sites)){
				
				if(!is.null(object@status) && object@status$disk.dump){
					
					if("ff" %in% class(object@pval.sites)){
						ffmatrix<-object@pval.sites
						ffsave(ffmatrix, file=file.path(path, "rnb.pvals"),rootpath=getOption('fftempdir'))
						rm(ffmatrix)
						
					}
				}
				
				if(!is.null(object@regions)){
					
					for(rgn in names(object@regions)){
						
						rgnpath<-file.path(path,rgn)
						if(!file.exists(rgnpath)){
							dir.create(rgnpath)	
						}
						
						if("ff" %in% class(object@pval.regions[[rgn]])){
							
							ffmatrix<-object@meth.regions[[rgn]]
							ffsave(ffmatrix, file=file.path(path, rgn, "rnb.pval"),rootpath=getOption('fftempdir'))
							rm(ffmatrix)
							
						}			
					}
				}
			}
			callNextMethod(object, path)
			
		})

########################################################################################################################
		
setMethod("load.matrices", signature(object="RnBeadSet", path="character"),
		
		function(object, path, temp.dir=tempdir()){
			
		if(!is.null(object@pval.sites)){
			
			if(sum(grepl("rnb.pvals",list.files(path)))==2){
				load_env<-new.env()
				suppressMessages(ffload(file=file.path(path, "rnb.pvals"), envir=load_env,rootpath=getOption("fftempdir")))
				object@pval.sites<-get("ffmatrix", envir=load_env)
				rm(load_env)
			}
			
			if(!is.null(object@regions)){
				
				for(rgn in names(object@regions)){
					
					if(sum(grepl("rnb.pvals",list.files(file.path(path, rgn))))==2){
						load_env<-new.env()
						suppressMessages(ffload(file=file.path(path, rgn, "rnb.pvals"), envir=load_env,rootpath=getOption("fftempdir")))
						object@pval.regions[[rgn]]<-get("ffmatrix", envir=load_env)
						rm(load_env)
					}		
				}
			}
		}
			
			callNextMethod(object=object, path=path, temp.dir=temp.dir)
			
		})

########################################################################################################################

#' @rdname destroy-methods
#' @aliases destroy,RnBeadSet-method
#' @docType methods
#' @export
setMethod("destroy", signature(object="RnBeadSet"),
		function(object){
			
			if(object@status$disk.dump){
				
				if(!is.null(object@pval.sites)){
					delete(object@pval.sites)
				}
			}
			callNextMethod()
			
		}
)

## ---------------------------------------------------------------------------------------------------------------------
## HELPER ROUTINES
## ---------------------------------------------------------------------------------------------------------------------
get.relative.covg<-function(m, design.vector){
	
	des1.quant<-sort(m[design.vector=="I",])[round(length(which(design.vector=="I"))*ncol(m)*0.001)]
	des2.quant<-sort(m[design.vector=="II",])[round(length(which(design.vector=="II"))*ncol(m)*0.001)]
	m[design.vector=="I",]<-m[design.vector=="I",]/as.double(des1.quant)
	m[design.vector=="II",]<-m[design.vector=="II",]/as.double(des2.quant)
	m
	
}
########################################################################################################################
summarize.bead.counts<-function(bead.counts.M, bead.counts.U, method="min"){
	
	nrow.M<-nrow(bead.counts.M); ncol.M<-ncol(bead.counts.M)
	nrow.U<-nrow(bead.counts.U); ncol.U<-ncol(bead.counts.U)
	
	if(nrow.M!=nrow.U || ncol.M!=ncol.U){
		stop("Dimensions of bead count matrices differ")
	}else{
		bead.counts<-matrix(NA,nrow.M,ncol.M)
		rownames(bead.counts)<-rownames(bead.counts.M)
	}
	
	if(method=="min"){
		index1<-bead.counts.M<=bead.counts.U
		index1[is.na(index1)]<-FALSE
		bead.counts[index1]<-bead.counts.M[index1]
		index2<-bead.counts.U<bead.counts.M
		index2[is.na(index2)]<-FALSE
		bead.counts[index2]<-bead.counts.U[index2]
		
	}
	
	bead.counts
	
}
########################################################################################################################
combine.sites<-function(sites1, sites2){
	
	common.chr<-intersect(unique(sites1[,2]), unique(sites2[,2]))
	
	common.sites<-lapply(common.chr, function(chr){
				
				sts<-intersect(sites1[sites1[,2]==chr,3],sites2[sites2[,2]==chr,3])				
				
				rbind(rep(1,length(sts)), rep(chr,length(sts)), sts)
				
			})
	
		
	if("ff_matrix" %in% c(class(sites1), class(sites2))){
	
		new.sites<-ff(vmode="integer", dim=c(sum(sapply(common.sites, nrow)),3))

		ixx<-1
		for(sts in common.sites){
			new.sites[,][ixx:(ixx+nrow(sts)),]<-sts
			ixx<-ixx+nrow(sts)+1
		}
		
	}else{
		new.sites<-do.call("rbind", common.sites)
	}
	
	
	new.sites
		
}
########################################################################################################################
