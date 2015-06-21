########################################################################################################################
## RnBiseqSet-class.R
## created: 2012-10-29
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## RnBiseqSet class definition.
########################################################################################################################

#' RnBiseqSet Class
#'
#' A class for storing the DNA methylation and quality information from bisulfite sequencing experiments 
#'
#' @details TBA
#'
#' @section Slots:
#' \describe{
#'   \item{\code{status}}{Normalization status.}
#' }
#'
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{\link[BiocGenerics]{combine}}}{Combines two datasets.}
#' }
#'
#' @name RnBiseqSet-class
#' @rdname RnBiseqSet-class
#' @author Pavlo Lutsik
#' @exportClass RnBiseqSet
setClass("RnBiseqSet",
		representation(status="listOrNULL"),
		contains="RnBSet",
		prototype(pheno=data.frame(), meth.sites=matrix(), meth=matrix(), coverage=NULL, status=NULL),
		package = "RnBeads")

## ---------------------------------------------------------------------------------------------------------------------
## CONSTRUCTORS
## ---------------------------------------------------------------------------------------------------------------------

#' initialize.RnBiseqSet
#'
#' Direct slot filling.
#'
#' @param pheno       	phenotypic data.
#' @param sites       	matrix of CpG site mappings.
#' @param meth.sites  	matrix of summarized methylation calls 
#' @param covg.sites 	matrix with read coverage information 
#' @param assembly		the genome assembly
#' @param target		target DNA methylation features (CpG sites)
#' @param regions		list of region mapping matrices
#' @param meth.regions	list of matrices with summarized region-level methylation calls  
#' @param covg.regions  list of matrices with summarized region-level coverage data
#' @param region.types	region annotations for which the methylation data should be summarized 
#' @param useff			a flag specifying whether the ff functionality should be used
#' @param verbose		a flag specifying whether the output 
#' 
#' @export
#' @docType methods
#' @rdname RnBiseqSet-class
setMethod("initialize", "RnBiseqSet",
		function(.Object,
				pheno,
				sites,
				meth.sites,
				covg.sites,
				assembly="hg19",
				target="CpG",
				regions=list(),
				meth.regions=list(),
				covg.regions=NULL,
				region.types=rnb.region.types.for.analysis(assembly),
				useff=rnb.getOption("disk.dump.big.matrices"),
				verbose=FALSE){

			if (nrow(sites)!= nrow(meth.sites)){
				msg<-c("Inconsistent values for sites and meth.sites")
				rnb.error(msg)
			}
			
			if(useff){
				
				#subsampling for large datasets (ff currently supports only objects of size .Machine$integer.max)
				if (prod(dim(meth.sites))>.Machine$integer.max){
					sites.allowed <- as.integer(.Machine$integer.max/ncol(meth.sites))
					sample.site.inds <- sort(sample.int(nrow(meth.sites),sites.allowed))
					msg<-c("Full dataset is too large to be supported by ff. --> downsampling to",sites.allowed,"( of",nrow(meth.sites),") sites")
					rnb.warning(msg)
					meth.sites <- meth.sites[sample.site.inds,]
					if(!is.null(covg.sites)){
						covg.sites <- covg.sites[sample.site.inds,]
					}
					sites <- sites[sample.site.inds,]
				}
					
				if("ff_matrix" %in% class(meth.sites)){
					.Object@meth.sites<-meth.sites
				}else if(class(meth.sites)=="matrix"){
					.Object@meth.sites<-convert.to.ff.matrix.tmp(meth.sites)
				}else{
					stop("Invalid value for meth.sites supplied")
				}
				
				if(!is.null(covg.sites)){
					if("ff_matrix" %in% class(covg.sites)){
						.Object@covg.sites<-covg.sites
					}else if(class(covg.sites)=="matrix"){
						.Object@covg.sites<-convert.to.ff.matrix.tmp(covg.sites)
					}else {
						stop("Invalid value for covg.sites supplied")
					}
				}
				
				.Object@status<-list(disk.dump=TRUE)
			}else{
				.Object@meth.sites<-meth.sites
				.Object@covg.sites<-covg.sites
				.Object@status<-list(disk.dump=FALSE)
			}
			.Object@pheno<-pheno	
			.Object@sites<-sites

			.Object@regions<-regions
			.Object@meth.regions<-meth.regions
			.Object@covg.regions<-covg.regions
			.Object@assembly<-assembly
			.Object@target<-target

			.Object@inferred.covariates <- list()

			if(!rnb.getOption("strand.specific")){
				msg<-c("Summarizing","strand","methylation")
				if(verbose){
					rnb.status(msg)
				}
				if(is.null(covg.sites)){
					aggr<-"mean"
				}else{
					aggr<-"coverage.weighted"
				}
				.Object <- summarize.regions(.Object, "strands", aggregation=aggr)
			}
			
			for (region.type in region.types) {
				if (region.type %in% rnb.region.types(assembly)) {
					msg<-c("Summarizing",region.type,"methylation")
					if(verbose){
						rnb.status(msg)
					}
					.Object <- summarize.regions(.Object, region.type)
				}
			}
			.Object
		})

## ---------------------------------------------------------------------------------------------------------------------
## VALIDITY
## ---------------------------------------------------------------------------------------------------------------------
#
# check.rnb.biseq.set
#
# A routine that checks the validity of a freshly loaded RnBiseqSet object
# 
# #param object RnBiseqSet object to test
# #param verbose A flag specifying whether the logger is initialized
#
# #return TRUE if all checks are successful
# 
# Pavlo Lutsik
#
check.rnb.biseq.set<-function(object, verbose=TRUE){
	
	
	if(!inherits(object, "RnBiseqSet")){
		if(verbose){ 
			rnb.info("The supplied object is not of class RnBiseqSet. Breaking the check...")
		}
		return(FALSE)
	}else{
		if(verbose){
			rnb.info("Checking the supplied RnBiseqSet object")
		}
		valid<-TRUE
	}
	
	#check the methylation data
	
	mm<-meth(object)
	
	if(nrow(mm)<1){
		if(verbose){
			rnb.info("The object contains information for 0 methylation sites")
		}
		valid<-FALSE
	}else{
		if(verbose){
			rnb.info(sprintf("The object contains information for %d methylation sites",nrow(mm)))
		}
	}
	
	if(ncol(mm)<1){
		if(verbose){
			rnb.info("The object contains information for 0 samples")
		}
		valid<-FALSE
	}else{
		if(verbose){
			rnb.info(sprintf("The object contains information for %d samples",ncol(mm)))
		}
	}
	
	if(sum(!is.na(mm))<1){
		if(verbose){
			rnb.info("All methylation values are missing")
		}
		valid<-FALSE
	}else{
		if(verbose){
			rnb.info(sprintf("The object contains %d missing methylation values",sum(is.na(mm))))
		}
	}
	
	if(sum(mm[!is.na(mm)]<0 & mm[!is.na(mm)]>1)>0){
		if(verbose){
			rnb.info("The object contains incorrect methylation values (less than 0 or greater than 1)")
		}
		valid<-FALSE
	}else{
		if(verbose){
			rnb.info("Methylation values are within the expected range")
		}
	}
	
	#check the coverage data
	
	if(is.null(object@covg.sites)){
		if(verbose){
			rnb.info("No coverage information found")
		}
		valid<-FALSE
		return(valid)
	}else{
		if(verbose){
			rnb.info("The object contains coverage information")
		}
		cvg<-covg(object)
	}
	
	if(sum(cvg[!is.na(cvg)]<0)>0){
		if(verbose){
			rnb.info("The object contains incorrect coverage values (negative values)")
		}
		valid<-FALSE
	}else if(sum(cvg[!is.na(cvg)]>1000)>0.5*length(as.numeric(cvg))){
		if(verbose){
			rnb.info("The object contains incorrect coverage values (half of the sites have coverage > 1000 )")
		}
	}else{
		if(verbose){
			rnb.info("Coverage values are within the expected range")
		}
	}
	
	return(valid)
	
}
########################################################################################################################

setMethod("show", "RnBiseqSet",
	function(object) {
		cat("Object of class RnBiseqSet\n")
		cat(sprintf("%8d samples\n", nrow(object@pheno)))
		cat(sprintf("%8d methylation sites\n", nrow(object@sites)))
		if (!is.null(object@regions)){
			cat("Region types:\n")
			for (rn in names(object@regions)) {
				cat(sprintf("\t%8d regions of type %s\n", nrow(object@regions[[rn]]), rn))
			}
		}
		cat(sprintf("Coverage information is %s\n", ifelse(is.null(object@covg.sites), "absent", "present")))
	}
)

########################################################################################################################