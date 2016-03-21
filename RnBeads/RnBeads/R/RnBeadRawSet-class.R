########################################################################################################################
## RnBeadRawSet-class.R
## created: 2013-xx-xx
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## RnBeadRawSet class definition.
########################################################################################################################

## ---------------------------------------------------------------------------------------------------------------------
## CLASS DEFINITIONS
## ---------------------------------------------------------------------------------------------------------------------

#' RnBeadRawSet-class 
#' 
#' Main class for storing HumanMethylation micorarray data which includes intensity information 
#' 
#' @section Slots:
#' \describe{
#'	\item{\code{pheno}}{Phenotypic data.}
#'	\item{\code{M}}{\code{matrix} of intensities for the probes measuring the abundance of methylated molecules.} 
#'	\item{\code{U}}{\code{matrix} of intensities for the probes measuring the abundance of unmethylated molecules.} 
#' 	\item{\code{M0}}{\code{matrix} of "out-of-band" intensities for the probes measuring the abundance of methylated molecules.} 
#' 	\item{\code{U0}}{\code{matrix} of "out-of-band" intensities for the probes measuring the abundance of unmethylated molecules.}
#' 	\item{\code{bead.counts.M}}{\code{matrix} of bead counts per probe.}
#' 	\item{\code{bead.counts.U}}{\code{matrix} of bead counts per probe.}
#' }
#' 
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{samples}}{Gets the identifiers of all samples in the dataset.}
#'   \item{\code{\link[=M,RnBeadRawSet-method]{M}}}{Get the matrix of intensities for the probes measuring the abundance of methylated molecules.}
#'   \item{\code{\link[=U,RnBeadRawSet-method]{U}}}{Get the matrix of intensities for the probes measuring the abundance of unmethylated molecules.}
#'   \item{\code{\link{intensities.by.color}}}{Get probe intensities in each color channel.}
#' } 
#' 
#' @name RnBeadRawSet-class
#' @rdname RnBeadRawSet-class
#' @author Pavlo Lutsik
#' @exportClass RnBeadRawSet
#' @include RnBeadSet-class.R
setClass("RnBeadRawSet",
		representation(M="matrixOrffOrNULL",
				U="matrixOrffOrNULL",
				M0="matrixOrffOrNULL",
				U0="matrixOrffOrNULL",
				bead.counts.M="matrixOrffOrNULL",
				bead.counts.U="matrixOrffOrNULL"
		),
		contains="RnBeadSet",
		prototype(pheno=data.frame(),
				betas=matrix(),
				meth.sites=matrix(),
				pval.sites=NULL,
				pval.regions=NULL,
				qc=NULL, 
				status=NULL,
				M=matrix(),
				U=matrix(),
				M0=NULL,
				U0=NULL,
				bead.counts.M=NULL,
				bead.counts.U=NULL
		),
		package = "RnBeads"
)

RNBRAWSET.SLOTNAMES<-c("M","U","M0","U0","bead.counts.M", "bead.counts.U")

########################################################################################################################

#' initialize.RnBeadRawSet
#'
#' Direct slot filling.
#'
#' @param pheno       		Phenotypic data.
#' @param M       	  		Matrix of intensities for the probes measuring the abundance of methylated molecules 
#' @param U       	  		Matrix of intensities for the probes measuring the abundance of unmethylated molecules 
#' @param M0       	  		Matrix of "out-of-band" intensities for the probes measuring the abundance of methylated molecules 
#' @param U0       	  		Matrix of "out-of-band" intensities for the probes measuring the abundance of unmethylated molecules
#' @param bead.counts.M 	Matrix of bead counts per probe.
#' @param bead.counts.U 	Matrix of bead counts per probe.
#' @param p.values    		Matrix of detection p-values.
#' @param platform	   		\code{character} singleton specifying the microarray platform: \code{"450k"} corresponds to HumanMethylation450 microarray, and \code{"27k"} stands for HumanMethylation27.
#' @param region.types		A \code{character} vector specifying the region types, for which the methylation infromation will be summarized.
#' @param useff		  		If \code{TRUE} the data matrices will be stored as \code{ff} objects
#' @param beta.offset		A regularization constant which is added to the denominator at beta-value calculation
#' @param summarize.bead.counts	If \code{TRUE} the coverage slot is filled by summarizing the \code{bead.counts.M} and \code{bead.counts.U} matrices. For type I probes the summarization is done using \code{min} operation, while for type II probes the bead counts should be identical in both supplied matrices
#'
#' 
#' @details Row names of the \code{M} matrix are used to match the data to annotation. All matrices should have identical dimensions.
#' 
#' @export
#' @docType methods
#' @rdname RnBeadRawSet-class
setMethod("initialize", "RnBeadRawSet",
		function(.Object,
				pheno,
				M,
				U,
				M0=NULL,
				U0=NULL,
				bead.counts.M = NULL, 
				bead.counts.U = NULL,
				p.values=NULL,
				platform = "450k",
				region.types = rnb.region.types.for.analysis("hg19"),
				useff=rnb.getOption("disk.dump.big.matrices"),
				beta.offset=100,
				summarize.bead.counts=TRUE) {
			
						
			if(is.null(rownames(M)) || is.null(rownames(U))){
				err<-"M and U matrices do not have rownames"
				rnb.error(err)
			}
			
			if(!setequal(rownames(M), rownames(U))){
				err<-"Inconsistent M and U matrices, the rownames do not match"
				rnb.error(err)
			}
			
			rn<-rownames(M)
			
			if (platform =="450k") {
				.Object@target <- "probes450"
			} else if(platform == "27k"){
				.Object@target <- "probes27"
			}else{
				rnb.error("Invalid value for platform")
			}
			
			probe.annotation <- rnb.get.annotation(.Object@target, .Object@assembly)
			
			## Construct the matrix of site indices
			#x.data<-GenomicRanges::as.data.frame(probe.annotation)
			
			chrs<-rnb.get.chromosomes(.Object@assembly)
			
			annotated.probes<-do.call("c", lapply(probe.annotation, names))
			if(length(which(rn %in% annotated.probes))<2){
				err<-"Annotations could be found for less than two rows from the supplied beta value table"
				rnb.error(err)		
			}
			
			if(!all(rn %in% annotated.probes)){
				warn<-"Some of the supplied probes are missing annotation and will be discarded"
				rnb.warning(warn)		
			}
			
			x.data<-rnb.annotation2data.frame(probe.annotation)
			rownames(x.data)<-annotated.probes
			rm(annotated.probes)
			
			#site.ids<-character()
			#site.ids<-integer()
			site.ids<- lapply(unique(x.data[["Chromosome"]]), function(chr) {
						chr.ids<-rownames(x.data)[x.data[["Chromosome"]]==chr]
						#table<-match(rn, chr.ids)
						table<-match(chr.ids,rn)
						#table<-sort(na.omit(table))
						table<-na.omit(table)
						#site.ids<<-c(site.ids, chr.ids[table])
						#site.ids<<-c(site.ids,table)
						table
					})
			site.ids<-do.call("c", site.ids)
			
			if(!is.null(bead.counts.M) && !is.null(bead.counts.U) && summarize.bead.counts){
				bead.counts<-summarize.bead.counts(bead.counts.M,bead.counts.U)
			}else{
				bead.counts<-NULL
			}
			
			betas<-beta.value(M[site.ids,,drop=FALSE],U[site.ids,,drop=FALSE], beta.offset)
		
			if(useff){
				
				.Object@M<-convert.to.ff.matrix.tmp(M[site.ids,,drop=FALSE])
				.Object@U<-convert.to.ff.matrix.tmp(U[site.ids,,drop=FALSE])
				if(!is.null(M0)){
					.Object@M0<-convert.to.ff.matrix.tmp(M0[site.ids,,drop=FALSE])
				}
				if(!is.null(U0)){
					.Object@U0<-convert.to.ff.matrix.tmp(U0[site.ids,,drop=FALSE])
				}
				if(!is.null(bead.counts.M)){
					.Object@bead.counts.M<-convert.to.ff.matrix.tmp(bead.counts.M[site.ids,,drop=FALSE])
				}
				if(!is.null(bead.counts.U)){
					.Object@bead.counts.U<-convert.to.ff.matrix.tmp(bead.counts.U[site.ids,,drop=FALSE])
				}
				
			}else{
				
				.Object@M<-M[site.ids,,drop=FALSE]
				.Object@U<-U[site.ids,,drop=FALSE]
				if(!is.null(M0)){
					.Object@M0<-M0[site.ids,,drop=FALSE]
				}
				if(!is.null(U0)){
					.Object@U0<-U0[site.ids,,drop=FALSE]
				}
				if(!is.null(bead.counts.M)){
					.Object@bead.counts.M<-bead.counts.M[site.ids,,drop=FALSE]
				}
				if(!is.null(bead.counts.U)){
					.Object@bead.counts.U<-bead.counts.U[site.ids,,drop=FALSE]
				}
				
			}
			
			rownames(.Object@M)<-NULL			
			rownames(.Object@U)<-NULL
			
			if(!is.null(.Object@M0)){
				rownames(.Object@M0)<-NULL
			}
			if(!is.null(.Object@U0)){
				rownames(.Object@U0)<-NULL
			}
			if(!is.null(.Object@bead.counts.M)){
				rownames(.Object@bead.counts.M)<-NULL
			}
			if(!is.null(.Object@bead.counts.U)){
				rownames(.Object@bead.counts.U)<-NULL
			}
			
			callNextMethod(.Object,
					pheno=pheno,
					betas=betas, 
					p.values=p.values, 
					bead.counts=bead.counts,
					platform=platform,
					region.types=region.types,
					useff=useff)
						
		})

########################################################################################################################
		
## FIXME: dummy validity method
validRnBeadRawSetObject<-function(object){
	return(TRUE)
}

setValidity("RnBeadRawSet", method=validRnBeadRawSetObject)

########################################################################################################################

setMethod("show", "RnBeadRawSet", rnb.show.rnbeadset)

########################################################################################################################

#' as("MethyLumiSet", "RnBeadRawSet")
#' 
#' Convert a \code{\linkS4class{MethyLumiSet}} object to \code{\linkS4class{RnBeadRawSet}}
#' 
#' @name as.RnBeadRawSet
setAs("MethyLumiSet", "RnBeadRawSet",
		
		function(from, to){
			
			if(!inherits(from,"MethyLumiSet")){
				stop("not a MethyLumiSet object:", deparse(substitute(methylumi.set)))
			}
			
			m.data <- MethyLumiSet2RnBeadSet(from)
			
			if("methylated.N" %in% ls(from@assayData) && "unmethylated.N" %in% ls(from@assayData) ){
				meth.N.element<-"methylated.N"
				umeth.N.element<-"unmethylated.N"
			}else if("Avg_NBEADS_A" %in% ls(from@assayData) && "Avg_NBEADS_B" %in% ls(from@assayData)){
				meth.N.element<-"Avg_NBEADS_B"	
				umeth.N.element<-"Avg_NBEADS_A"
			}else{
				meth.N.element<-NULL
				umeth.N.element<-NULL
			}
			
			if("methylated.OOB" %in% ls(from@assayData) && "unmethylated.OOB" %in% ls(from@assayData)){
				meth.oob.element<-"methylated.OOB"
				umeth.oob.element<-"unmethylated.OOB"
			}else{
				meth.oob.element<-NULL
				umeth.oob.element<-NULL
			}
			
			if(annotation(from)=="IlluminaHumanMethylation450k"){
				platform="450k"
			}else if(annotation(from)=="IlluminaHumanMethylation27k"){
				platform="27k"
			}
			
			object<-new(to,
					pheno=m.data$pheno,
					M=methylated(from),
					U=unmethylated(from),
					p.values=m.data$p.values,
					M0=if(!is.null(meth.oob.element)) get(meth.oob.element,from@assayData) else NULL,
					U0=if(!is.null(umeth.oob.element)) get(umeth.oob.element,from@assayData) else NULL,
					bead.counts.M=if(!is.null(meth.N.element)) get(meth.N.element,from@assayData) else NULL,
					bead.counts.U=if(!is.null(umeth.N.element)) get(umeth.N.element,from@assayData) else NULL,
					platform=platform
					)
			
			if ("qc" %in% names(m.data)) {
				qc(object)<-m.data[["qc"]]
			}
						
			object
			
		})

########################################################################################################################
		
#' as("RnBeadRawSet", "MethyLumiSet")
#'
#' Convert a \code{\linkS4class{RnBeadRawSet}} object to \code{\linkS4class{MethyLumiSet}}
#' 
#' @name as.RnBeadRawSet 
setAs("RnBeadRawSet","MethyLumiSet",
		
		function(from, to){
			
			if(!inherits(from,"RnBeadRawSet")){
				stop("not a RnBeadRawSet object:", deparse(substitute(methylumi.set)))
			}
				
			assd<-new.env()
			
			assign("betas", meth(from),  envir=assd)
			assign("pvals", dpval(from),  envir=assd)
			assign("methylated", M(from),  envir=assd)
			assign("unmethylated", U(from),  envir=assd)
			
			if(!is.null(M0(from))){
				assign("methylated.OOB", M0(from),  envir=assd)
			}
			if(!is.null(U0(from))){
				assign("unmethylated.OOB", U0(from),  envir=assd)
			}
			if(!is.null(bead.counts.M(from))){
				assign("methylated.N", bead.counts.M(from),  envir=assd)
			}
			if(!is.null(bead.counts.U(from))){
				assign("unmethylated.N", bead.counts.U(from),  envir=assd)
			}
			
			pd<-pheno(from)
			rownames(pd)<-colnames(meth(from))
			
			mset<-new(to, assd, as(pd, "AnnotatedDataFrame"))		
			
			rm(assd)
			
			ann<-annotation(from, add.names=TRUE)
			featureData(mset)<-as(data.frame(COLOR_CHANNEL=ann$Color, rownames=rownames(ann)), "AnnotatedDataFrame")
			featureNames(mset)<-ann[["ID"]]
			
			if (!is.null(qc(from))){
				
				assd<-new.env()
				
				assign("methylated", qc(from)$Cy3,  envir=assd)
				assign("unmethylated", qc(from)$Cy5,  envir=assd)
				
				mset@QC<-new("MethyLumiQC", assd)
						
				if(from@target == "probes450"){
					probeIDs<-rnb.get.annotation("controls450")[,"Target"]
					probeIDs<-paste(probeIDs, unlist(sapply(table(probeIDs)[unique(probeIDs)], seq, from=1 )), sep=".")
				}else if(from@target == "probes27"){
					probeIDs<-rnb.get.annotation("controls27")[,"Name"]
				}
				
				featureData(mset@QC)<-as(data.frame(Address=rownames(qc(from)$Cy3), rownames=probeIDs), "AnnotatedDataFrame")
				featureNames(mset@QC)<-probeIDs
								
				if(from@target == "probes450"){
					annotation(mset@QC) <- "IlluminaHumanMethylation450k"
				}else if(from@target == "probes27"){
					annotation(mset) <- "IlluminaHumanMethylation27k"
				}
				
			}
			if(from@target == "probes450"){
				annotation(mset) <- "IlluminaHumanMethylation450k"
			}else if(from@target == "probes27"){
				annotation(mset) <- "IlluminaHumanMethylation27k"
			}
			
			mset
			
		})


########################################################################################################################
		
## ---------------------------------------------------------------------------------------------------------------------
## ACCESSORS
## ---------------------------------------------------------------------------------------------------------------------

if(!isGeneric("M")) setGeneric('M',
			function(object, ...) standardGeneric('M'))

#' M-methods
#'
#' Extract raw methylated probe intensity from an object of \code{RnBeadRawSet} class.
#'
#' @param object 		Dataset of interest.
#' @param row.names		Flag indicating whether the resulting matrix will be assigned row names
#'  
#' @return \code{matrix} of the methylated probe intensities
#'
#' @rdname M-methods
#' @docType methods
#' @export
#' @aliases M
#' @aliases M,RnBeadRawSet-method
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' M.intensity<-M(rnb.set.example)
#' head(M.intensity)
#' } 
#' 
setMethod("M", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@M, object@meth.regions)
		})


if(!isGeneric("U")) setGeneric('U',
			function(object, ...) standardGeneric('U'))

########################################################################################################################

#' U-methods
#'
#' Extract raw unmethylated probe intensity from an object of \code{RnBeadRawSet} class.
#'
#' @param object 		Dataset of interest.
#' @param row.names		Flag indicating whether the resulting matrix will be assigned row names
#'  
#' @return \code{matrix} of the unmethylated probe intensities
#'
#' @rdname U-methods
#' @docType methods
#' @export
#' @aliases U
#' @aliases U,RnBeadRawSet-method
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' U.intensity<-U(rnb.set.example)
#' head(U.intensity)
#' } 
setMethod("U", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@U, object@meth.regions)
		})

########################################################################################################################

setGeneric('M0',
			function(object, ...) standardGeneric('M0'))

setMethod("M0", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@M0, object@meth.regions)
		})

########################################################################################################################

setGeneric('U0',
			function(object, ...) standardGeneric('U0'))

setMethod("U0", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@U0, object@meth.regions)
		})
########################################################################################################################

setGeneric('bead.counts.M',
			function(object, ...) standardGeneric('bead.counts.M'))


setMethod("bead.counts.M", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@bead.counts.M, object@meth.regions)
		})

########################################################################################################################

setGeneric('bead.counts.U',
			function(object, ...) standardGeneric('bead.counts.U'))


setMethod("bead.counts.U", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@bead.counts.U, object@meth.regions)
		})

## ---------------------------------------------------------------------------------------------------------------------
## MODIFIERS
## ---------------------------------------------------------------------------------------------------------------------

setGeneric('M<-',
			function(object, value) standardGeneric('M<-'))

setMethod("M<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@M)
				object@M<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@M<-value
			}
			
		})
########################################################################################################################

setGeneric('U<-',
			function(object, value) standardGeneric('U<-'))

setMethod("U<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@U)
				object@U<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@U<-value
			}
		})

########################################################################################################################

setGeneric('M0<-',
			function(object, value) standardGeneric('M0<-'))

setMethod("M0<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@M0)
				object@M0<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@M0<-value
			}
		})
########################################################################################################################

setGeneric('U0<-',
			function(object, value) standardGeneric('U0<-'))

setMethod("U0<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@U0)
				object@U0<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@U0<-value
			}
		})
########################################################################################################################

setGeneric('bead.counts.M<-',
			function(object, value) standardGeneric('bead.counts.M<-'))

setMethod("bead.counts.M<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@bead.counts.M)
				object@bead.counts.M<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@bead.counts.M<-value
			}
		})
########################################################################################################################

setGeneric('bead.counts.U<-',
			function(object, value) standardGeneric('bead.counts.U<-'))

setMethod("bead.counts.U<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@bead.counts.U)
				object@bead.counts.U<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@bead.counts.U<-value
			}
		})
########################################################################################################################

if (!isGeneric("remove.sites")) {
	setGeneric("remove.sites", function(object, probelist, verbose = TRUE) standardGeneric("remove.sites"))
}

#' @rdname remove.sites-methods
#' @aliases remove.sites,RnBeadRawSet-method
#' @docType methods
#' @export
setMethod("remove.sites", signature(object = "RnBeadRawSet"),
		function(object, probelist, verbose = TRUE) {
			inds <- get.i.vector(probelist, rownames(object@meth.sites))
			if (length(inds) != 0) {
				for(sl in RNBRAWSET.SLOTNAMES){
					if(!is.null(slot(object,sl))){
						if(!is.null(object@status) && object@status$disk.dump){
							new.matrix<-slot(object,sl)[-inds,]
							if(isTRUE(object@status$discard.ff.matrices)){
								delete(slot(object,sl))
							}
							slot(object,sl)<-convert.to.ff.matrix.tmp(new.matrix)
							rm(new.matrix); rnb.cleanMem()
						}else{
							slot(object,sl)<-slot(object,sl)[-inds,]
						}
					
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
#' @aliases remove.samples,RnBeadRawSet-method
#' @docType methods
#' @export
setMethod("remove.samples", signature(object = "RnBeadRawSet"),
		function(object, samplelist) {
			inds <- get.i.vector(samplelist, samples(object))
			if (length(inds) != 0) {
				for(sl in RNBRAWSET.SLOTNAMES){
					if(!is.null(slot(object,sl))){
						if(!is.null(object@status) && object@status$disk.dump){
							new.matrix<-slot(object,sl)[,-inds, drop=F]
							if(isTRUE(object@status$discard.ff.matrices)){
								delete(slot(object,sl))
							}
							slot(object,sl)<-convert.to.ff.matrix.tmp(new.matrix)
							rm(new.matrix); rnb.cleanMem()
						}else{
							slot(object,sl)<-slot(object,sl)[,-inds, drop=F]
						}
					}
				}
			}
			callNextMethod()
		}
)

#######################################################################################################################

#if (!isGeneric("update.meth")) {
	setGeneric("update.meth", function(object) standardGeneric("update.meth"))
#}

##  
## update.meth
## 
## Update the methylation calls, after the change of intensity values
##
## param object 		RnBeadRawSet object
##
## return Updated RnBeadRawSet object
##
setMethod("update.meth", signature(object="RnBeadRawSet"),
		function(object){
			
			if(object@status$disk.dump){
				object@meth.sites<-convert.to.ff.matrix.tmp(beta.value(object@M[,], object@U[,]))	
			}else{
				object@meth.sites<-beta.value(object@M, object@U)
			}
			return(object)
		})

#######################################################################################################################
## save, load and destroy 

setMethod("save.matrices", signature(object="RnBeadRawSet", path="character"),
		function(object, path){
			
			if(!file.exists(path)){
				dir.create(path)
			}
			
			for(sl in RNBRAWSET.SLOTNAMES){
				if(!is.null(slot(object,sl))){
					
					if(!is.null(object@status) && object@status$disk.dump){
						
						if("ff" %in% class(slot(object,sl))){
							ffmatrix<-slot(object,sl)
							ffsave(ffmatrix, file=file.path(path, paste("rnb", sl, sep=".")),
									rootpath=getOption('fftempdir'))
							rm(ffmatrix)
							
						}
					}
				}
			}
			callNextMethod(object, path)
			
		})

#######################################################################################################################

setMethod("load.matrices", signature(object="RnBeadRawSet", path="character"),
		
		function(object, path, temp.dir=tempdir()){
			slot.names <- RNBRAWSET.SLOTNAMES
			for(sl in slot.names){
				if(!is.null(slot(object, sl))){
					
					if(paste("rnb",sl,"RData", sep=".") %in% list.files(path) &&
							paste("rnb",sl,"ffData", sep=".") %in% list.files(path)){
						load_env<-new.env()
						suppressMessages(ffload(file=file.path(path, paste("rnb", sl, sep=".")), 
										envir=load_env,rootpath=getOption("fftempdir")))
						slot(object, sl)<-get("ffmatrix", envir=load_env)
						rm(load_env)
					}
					
				}

			}
			
			callNextMethod(object=object, path=path, temp.dir=temp.dir)
			
		})

#######################################################################################################################

#' @rdname destroy-methods
#' @aliases destroy,RnBeadRawSet-method
#' @docType methods
#' @export
setMethod("destroy", signature(object="RnBeadRawSet"),
		function(object){
			
			if(object@status$disk.dump){
				for(sl in RNBRAWSET.SLOTNAMES){
					if(!is.null(slot(object,sl))){
						delete(slot(object, sl))	
					}
				}
			}
			callNextMethod()
			
		}
)

## ---------------------------------------------------------------------------------------------------------------------
## HELPER ROUTINES
## ---------------------------------------------------------------------------------------------------------------------

beta.value<-function(M,U,offset=100){
	
	betas<-M/(M+U+offset)
	return(betas)
	
}
#######################################################################################################################

m.value<-function(M,U,offset=100){
	
	mvals<-log2((M+offset)/(U+offset))
	return(mvals)
	
}

#######################################################################################################################

#' intensities.by.color
#' 
#' Rearranges information from "M" and "U" slots of a RnBeadsRawSet object by color channer.
#'
#' @param raw.set 			RnBeadRawSet object
#' @param address.rownames  if \code{TRUE} the rows of the returned matrices are named with the with the correspoding Illumina probe addresses 
#' @param add.oob			if \code{TRUE} the "out-of-band" intensities are included
#' @param add.controls		if \code{TRUE} the control probe intensities are included
#' @param add.missing		if \code{TRUE} the rows for the probes missing in \code{raw.set} is imputed with \code{NA} values
#' 
#' 
#' @author Pavlo Lutsik
intensities.by.color<-function(raw.set, 
		address.rownames=TRUE, 
		add.oob=TRUE, 
		add.controls=TRUE, 
		add.missing=TRUE
		){
	
	if(!require("IlluminaHumanMethylation450kmanifest")){
		rnb.error("IlluminaHumanMethylation450kmanifest should be installed")
	}		
	
	Mmatrix<-M(raw.set, row.names=TRUE)
	Umatrix<-U(raw.set, row.names=TRUE)
	if(add.oob){
		M0matrix<-M0(raw.set, row.names=TRUE)
		U0matrix<-U0(raw.set, row.names=TRUE)
	}
	
	pinfos <- annotation(raw.set, add.names=TRUE)

	if(add.missing){
		full.ann<-rnb.annotation2data.frame(rnb.get.annotation(raw.set@target))
		ann.missing<-full.ann[!rownames(full.ann)%in%rownames(pinfos),]
		pinfos<-rbind(pinfos, ann.missing[,colnames(full.ann)])
		filler<-matrix(NA_real_, nrow=nrow(ann.missing), ncol=length(samples(raw.set)))
		rownames(filler)<-rownames(ann.missing)

		Mmatrix<-rbind(Mmatrix, filler)
		Umatrix<-rbind(Umatrix, filler)
		
		if(add.oob){
			M0matrix<-rbind(M0matrix, filler)
			U0matrix<-rbind(U0matrix, filler)
		}
		rm(ann.missing, filler, full.ann)
	}
	
	rnb.set.probe.ids<-pinfos[["ID"]]
	
	dII.probes <- rnb.set.probe.ids[pinfos[,"Design"] == "II"]
	
	#dII.probes <- dII.probes[!grepl("rs", dII.probes)]
	
	if(address.rownames){
	
		tII<-rbind(as.data.frame(IlluminaHumanMethylation450kmanifest@data$TypeII[,c("Name", "AddressA")]),
			as.data.frame(IlluminaHumanMethylation450kmanifest@data$TypeSnpII[,c("Name", "AddressA")]))
	
		tII<-tII[match(dII.probes, tII$Name),]
	}
	dII.grn<-Mmatrix[pinfos[,"Design"] == "II",]
	if(address.rownames) rownames(dII.grn)<-tII$AddressA
	
	dII.red<-Umatrix[pinfos[,"Design"] == "II",]
	if(address.rownames) rownames(dII.red)<-tII$AddressA
	
	dI.red.probes <- rnb.set.probe.ids[pinfos[, "Color"] == "Red"]
	#dI.red.probes <- dI.red.probes[!grepl("rs", dI.red.probes)]
	dI.green.probes <- rnb.set.probe.ids[pinfos[, "Color"] == "Grn"]
	#dI.green.probes <- dI.green.probes[!grepl("rs", dI.green.probes)]
	
	if(address.rownames){
	
		tI<-rbind(as.data.frame(IlluminaHumanMethylation450kmanifest@data$TypeI[,c("Name","Color", "AddressA", "AddressB")]),
				as.data.frame(IlluminaHumanMethylation450kmanifest@data$TypeSnpI[,c("Name","Color", "AddressA", "AddressB")]))
	
	
		tI.red<-tI[tI$Color=="Red",]
		tI.red<-tI.red[match(dI.red.probes, tI.red$Name),]
	
		tI.grn<-tI[tI$Color=="Grn",]
		tI.grn<-tI.grn[match(dI.green.probes, tI.grn$Name),]
	}
	
	dI.red.meth<-Mmatrix[pinfos[, "Color"] == "Red",]
	
	if(address.rownames) rownames(dI.red.meth)<-tI.red[,"AddressB"]
	
	dI.red.umeth<-Umatrix[pinfos[, "Color"] == "Red",]
	if(address.rownames) rownames(dI.red.umeth)<-tI.red[,"AddressA"]
	
	if(add.oob){
		dI.red.meth.oob<-M0matrix[pinfos[, "Color"] == "Red",]
		if(address.rownames) rownames(dI.red.meth.oob)<-tI.red[,"AddressB"]
		
		dI.red.umeth.oob<-U0matrix[pinfos[, "Color"] == "Red",]
		if(address.rownames) rownames(dI.red.umeth.oob)<-tI.red[,"AddressA"]
	}
	
	dI.grn.meth<-Mmatrix[pinfos[, "Color"] == "Grn",]
	if(address.rownames) rownames(dI.grn.meth)<-tI.grn[,"AddressB"]
	
	dI.grn.umeth<-Umatrix[pinfos[, "Color"] == "Grn",]
	if(address.rownames) rownames(dI.grn.umeth)<-tI.grn[,"AddressA"]
	
	if(add.oob){
		dI.grn.meth.oob<-M0matrix[pinfos[, "Color"] == "Grn",]
		if(address.rownames) rownames(dI.grn.meth.oob)<-tI.grn[,"AddressB"]
		
		dI.grn.umeth.oob<-U0matrix[pinfos[, "Color"] == "Grn",]
		if(address.rownames) rownames(dI.grn.umeth.oob)<-tI.grn[,"AddressA"]
	}
	
	intensities.by.channel <- list(
			Cy3=rbind(dII.grn, dI.grn.meth,dI.grn.umeth, 
					if(add.oob) dI.red.meth.oob else NULL, if(add.oob) dI.red.umeth.oob else NULL),
			Cy5=rbind(dII.red, dI.red.meth, dI.red.umeth, 
					if(add.oob) dI.grn.meth.oob else NULL, if(add.oob) dI.grn.umeth.oob else NULL))
	
	rm(dII.grn, dI.grn.meth, dI.grn.umeth, dI.red.meth.oob, dI.red.umeth.oob, 
			dII.red, dI.red.meth, dI.red.umeth, dI.grn.meth.oob, dI.grn.umeth.oob)

	gc()
	
	if(address.rownames) intensities.by.channel$Cy5<-intensities.by.channel$Cy5[rownames(intensities.by.channel$Cy3),]
	if(add.controls){
		ncd<-rnb.get.annotation("controls450")
		#ncd<-ncd[ncd[["Target"]] == "NEGATIVE", ]
		ncd$Target<-tolower(ncd$Target)
		controls.by.channel<-qc(raw.set)
			
		controls.by.channel$Cy3<-controls.by.channel$Cy3[as.character(ncd$ID),]
		controls.by.channel$Cy5<-controls.by.channel$Cy5[as.character(ncd$ID),]
		
		intensities.by.channel$Cy3<-rbind(intensities.by.channel$Cy3, controls.by.channel$Cy3)
		intensities.by.channel$Cy5<-rbind(intensities.by.channel$Cy5, controls.by.channel$Cy5)
	}
	return(intensities.by.channel)
}
########################################################################################################################
