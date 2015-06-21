########################################################################################################################
## RnBDiffMeth-class.R
## created: 2012-08-30
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## RnBDiffMeth class definition.
########################################################################################################################

setClassUnion("characterOrNULL", c("character", "NULL"))

#' RnBDiffMeth Class
#'
#' A class for storing differential methylation data.
#' 
#' @details
#' Contains differential methylation tables (DMT) for multiple comparisons and region types. DMTs can be stored in memory
#' as R objects or on disk
#'
#' @section Slots:
#' \describe{
#'   \item{\code{sites}}{List of differential methylation tables on site level (see \code{computeDiffMeth.bin.site} for details).
#' 					Indexed by comparison.}
#'   \item{\code{regions}}{List of lists of differential methylation tables on region levels (see \code{computeDiffMeth.bin.region} for details).
#' 					Indexed by region type on the top level and comparison on the lower level.} 
#'   \item{\code{comparisons}}{character vector of all comparisons stored in the objects. Vector indices correspond to indices in the \code{sites} and
#' 					\code{regions} list slots.}
#'   \item{\code{region.types}}{character vector of all region types stored in the objects. Vector indices correspond to indices in
#' 					the \code{regions} list slot.}
#'   \item{\code{comparison.grouplabels}}{A character matrix with 2 columns containing group labels of all comparisons in the object}
#'   \item{\code{comparison.info}}{A list containing comparison information for each comparison. See \code{\link{get.comparison.info}} for details.}
#'   \item{\code{site.test.method}}{method which was applied to obtain the site-level p-values.}
#'   \item{\code{covg.thres}}{coverage threshold. Important for certain columns of the differential methylation tables.
#' 					See \code{computeDiffMeth.bin.site} and \code{computeDiffMeth.bin.region} for details.}
#'   \item{\code{disk.dump}}{Flag indicating whether the tables should be stored on disk rather than in the main memory}
#'   \item{\code{disk.path}}{path on the disk for DMTs.Only meaningful if \code{disk.dump} is \code{TRUE}}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{destroy,RnBDiffMeth-method}}}{remove tables stored to disk from the file system}
#'   \item{\code{\link{get.region.types,RnBDiffMeth-method}}}{Gets all region types represented in the object as character vector}
#' 	 \item{\code{\link{get.comparisons,RnBDiffMeth-method}}}{Gets all comparisons represented in the object as character vector}
#'   \item{\code{\link{get.comparison.grouplabels,RnBDiffMeth-method}}}{Gets all comparison group names as a matrix}
#'   \item{\code{\link{get.covg.thres,RnBDiffMeth-method}}}{Gets the coverage threshold employed for obtaining statistics in the differential methylation tables}
#'   \item{\code{\link{get.table,RnBDiffMeth-method}}}{Gets a differential methylation table}
#'   \item{\code{\link{addDiffMethTable,RnBDiffMeth-method}}}{Adds a differential methylation table}
#' 	 \item{\code{\link{reload,RnBDiffMeth-method}}}{relink disk dumped tables. Useful if the files are manually copied or if the object is loaded again}
#'   \item{\code{\link{save.tables,RnBDiffMeth-method}}}{save disk dumped tables as binaries and zip them.
#' 						Useful if the files are copied or shared.}
#'   \item{\code{\link{join.diffMeth}}}{Merges two disjoint RnBDiffMeth objects into one}
#' }
#'
#' @name RnBDiffMeth-class
#' @rdname RnBDiffMeth-class
#' @author Fabian Mueller
#' @exportClass RnBDiffMeth
setClass("RnBDiffMeth",
	slots = list(
		sites="list",
		regions="list",
		comparisons="character",
		region.types="character",
		comparison.grouplabels="matrix",
		comparison.info="list",
		site.test.method="characterOrNULL",
		covg.thres="integer",
		disk.dump="logical",
		disk.path="characterOrNULL"
	),
	prototype = prototype(
		sites=list(),
		regions=list(),
		comparisons=character(),
		region.types=character(),
		comparison.grouplabels=matrix(ncol=2,nrow=0),
		comparison.info=list(),
		site.test.method=NULL,
		covg.thres=-1L,
		disk.dump=FALSE,
		disk.path=NULL
	),
	package = "RnBeads")

#' initialize.RnBDiffMeth
#'
#' Initialize an RnBDiffMeth object

#' @param site.test.method method which was applied to obtain the site-level p-values.
#' @param covg.thres  coverage threshold. Important for certain columns of the differential methylation tables.
#' 					  See \code{computeDiffMeth.bin.site} and \code{computeDiffMeth.bin.region} for details.
#' @param disk.dump   Flag indicating whether the tables should be stored on disk rather than in the main memory
#' @param disk.path   Path on the disk for DMTs.Only meaningful if \code{disk.dump} is \code{TRUE}
#'
#' @export
#' @author Fabian Mueller
#' @docType methods
setMethod("initialize", "RnBDiffMeth",
	function(.Object,
		site.test.method=rnb.getOption("differential.site.test.method"),
		covg.thres=rnb.getOption("filtering.coverage.threshold"),
		disk.dump=FALSE,disk.path=NULL){
			#create directory in which to dump the matrices to file. must be a non-existing directory
			if (disk.dump){
				if (!is.null(disk.path)){
					if (!create.path(disk.path, FALSE, showWarnings = FALSE)) {
						stop(paste("Could not create directory (RnBDiffMeth):",disk.path))
					}
				}
			}
			
			.Object@sites <- list()
			.Object@regions <- list()
			.Object@comparisons <- character()
			.Object@region.types <- character()
			.Object@comparison.grouplabels <- matrix(ncol=2,nrow=0)
			.Object@site.test.method <- site.test.method
			.Object@covg.thres <- covg.thres
			.Object@comparison.info <- list()
			.Object@disk.dump <- disk.dump
			.Object@disk.path <- disk.path
			.Object
	}
)

#remove the directory containing the dumped matrices
if (!isGeneric("destroy")) setGeneric("destroy", function(object) standardGeneric("destroy"))
#' destroy-methods
#'
#' remove tables stored to disk from the file system. Useful for cleaning up disk dumped objects.
#' CAUTION: currently only works with reloaded objects
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @return Nothing of particular interest
#'
#' @rdname destroy-RnBDiffMeth-methods
#' @docType methods
#' @aliases destroy,RnBDiffMeth-method
#' @author Fabian Mueller
#' @export
setMethod("destroy", signature(object="RnBDiffMeth"),
	function(object){
		n.comps <- length(object@comparisons)
		n.region.types <- length(object@region.types)
		if (object@disk.dump){			
			logger.start("Deleting RnBDiffMeth disk dump files from disk")
			for (cci in 1:n.comps) {
				cc <- object@comparisons[cci]
				if (!is.null(object@sites[[cci]])){
					delete(object@sites[[cci]])
				}
				for (rri in 1:n.region.types){
					rr <- object@region.types[rri]
					if (!is.null(object@regions[[rri]][[cci]])){
						delete(object@regions[[rri]][[cci]])
					}
				}
			}
			unlink(object@disk.path,recursive=TRUE)
			logger.completed()
		}
	}
)

if (!isGeneric("get.region.types")) setGeneric("get.region.types", function(object) standardGeneric("get.region.types"))
#' get.region.types-methods
#'
#' Gets all region types represented in the object as character vector
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @return character vector containing region types
#'
#' @rdname get.region.types-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases get.region.types
#' @aliases get.region.types,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' get.region.types(dm)
#' }
setMethod("get.region.types", signature(object="RnBDiffMeth"),
	function(object){
		return(object@region.types)
	}
)

if (!isGeneric("get.comparisons")) setGeneric("get.comparisons", function(object) standardGeneric("get.comparisons"))
#' get.comparisons-methods
#'
#' Gets all comparisons represented in the object as character vector
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @return character vector containing comparisons
#'
#' @rdname get.comparisons-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases get.comparisons
#' @aliases get.comparisons,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' get.comparisons(dm)
#' }
setMethod("get.comparisons", signature(object="RnBDiffMeth"),
	function(object){
		return(object@comparisons)
	}
)

if (!isGeneric("get.comparison.grouplabels")) setGeneric("get.comparison.grouplabels", function(object) standardGeneric("get.comparison.grouplabels"))
#' get.comparison.grouplabels-methods
#'
#' Gets all comparison grouplabels represented in the object as character matrix of dimension n.comparisons x 2
#' where the columns specify group names 1 and 2 respectively
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @return character matrix containing comparison group names
#'
#' @rdname get.comparison.grouplabels-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases get.comparison.grouplabels
#' @aliases get.comparison.grouplabels,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' get.comparison.grouplabels(dm)
#' }
setMethod("get.comparison.grouplabels", signature(object="RnBDiffMeth"),
		function(object){
			return(object@comparison.grouplabels)
		}
)

if (!isGeneric("get.site.test.method")) setGeneric("get.site.test.method", function(object) standardGeneric("get.site.test.method"))
#' get.site.test.method-methods
#'
#' Gets the site testing method used to obtain the p-values in the differential methylation tables
#'
#' @param object RnBDiffMeth object
#' @return character describing the site test method
#'
#' @rdname get.site.test.method-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases get.site.test.method
#' @aliases get.site.test.method,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' get.site.test.method(dm)
#' }
setMethod("get.site.test.method", signature(object="RnBDiffMeth"),
		function(object){
			if (.hasSlot(object,"site.test.method")) { #.hasSlot ensure backwards compatibility
				return(object@site.test.method)
			} else {
				return(rnb.getOption("differential.site.test.method"))
			}
		}
)

if (!isGeneric("get.covg.thres")) setGeneric("get.covg.thres", function(object) standardGeneric("get.covg.thres"))
#' get.covg.thres-methods
#'
#' Gets the coverage threshold employed for obtaining statistics in the differential methylation tables
#'
#' @param object RnBDiffMeth object
#' @return integer coverage threshold
#'
#' @rdname get.covg.thres-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases get.covg.thres
#' @aliases get.covg.thres,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' get.covg.thres(dm)
#' }
setMethod("get.covg.thres", signature(object="RnBDiffMeth"),
		function(object){
			return(object@covg.thres)
		}
)

if (!isGeneric("get.table")) setGeneric("get.table", function(object,...) standardGeneric("get.table"))
#' get.table-methods
#'
#' Gets a differential methylation table
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @param comparison character or index of the comparison of the table to retrieve
#' @param region.type character or index of the region type of the table to retrieve
#' @param undump Flag indicating whether to convert the table into a matrix instead of using the file descriptor.
#' 		         Only meaningful if the if the objects's \code{disk.dump} slot is true.
#' @param return.data.frame should a data.frame be returned instead of a matrix?
#' @return differential methylation table. See \code{computeDiffMeth.bin.site} and \code{computeDiffMeth.bin.region} for details.
#'
#' @rdname get.table-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases get.table
#' @aliases get.table,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' dm.promoters <- get.table(dm,get.comparisons(dm)[1],"promoters",return.data.frame=TRUE)
#' summary(dm.promoters)
#' }
setMethod("get.table", signature(object="RnBDiffMeth"),
	function(object,comparison,region.type,undump=TRUE,return.data.frame=FALSE){
		if (!is.element(comparison,object@comparisons)) {
			stop(paste("invalid comparison:",comparison))
		}
		if (!is.element(region.type,c("sites",object@region.types))) {
			stop(paste("invalid region.type:",region.type))
		}
		if (region.type == "sites"){
			res <- object@sites[[comparison]]
		} else {
			res <- object@regions[[region.type]][[comparison]]
		}
		if (object@disk.dump && undump && !is.null(res)){
			res <- res[,]
		}
		if (return.data.frame){
			res <- data.frame(res)
		}
		return(res)
	}
)

if (!isGeneric("addDiffMethTable")) setGeneric("addDiffMethTable", function(object,...) standardGeneric("addDiffMethTable"))
#' addDiffMethTable-methods
#'
#' Adds a differential methylation table
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @param dmt    Differential methylation table to add
#' @param comparison character or index of the comparison of the table to retrieve
#' @param grp.labs character vector of length 2 specifying the names of the groups being compared
#' @param region.type character or index of the region type of the table to retrieve
#' @return the updated RnBDiffMeth object
#'
#' @note Caveat: if disk dumping is enabled the resulting object tables will be stored in the initial location of the object.
#' @rdname addDiffMethTable-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases addDiffMethTable
#' @aliases addDiffMethTable,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols="Sample_Group",region.types=c("genes","tiling"))
#' sample.groups <- rnb.sample.groups(rnb.set.example,"Sample_Group")[[1]]
#' dmt.sites <- computeDiffTab.extended.site(meth(rnb.set.example),sample.groups[[1]],sample.groups[[2]])
#' map.regions.to.sites <- regionMapping(rnb.set.example,"promoters")
#' dmt.promoters <- computeDiffTab.default.region(dmt.sites,map.regions.to.sites)
#' cmp.name <- get.comparisons(dm)[1]
#' grp.labs <- get.comparison.grouplabels(dm)[1,]
#' #add the promoter level differential methylation table
#' dm.add <- addDiffMethTable(dm,dmt.promoters,cmp.name,"promoters",grp.labs)
#' get.region.types(dm.add)
#' }
setMethod("addDiffMethTable", signature(object="RnBDiffMeth"),
	function(object,dmt,comparison,region.type,grp.labs=c("group1","group2")){
		if(!(is.character(comparison) && length(comparison)==1)) stop("Invalid argument: comparison")
		if(!(is.character(grp.labs) && length(grp.labs)==2)) stop("Invalid argument: grp.labs")
		if(!(is.character(region.type) && length(region.type)==1)) stop("Invalid argument: region.type")
		if(!(is.data.frame(dmt) || is.matrix(dmt))) stop("Invalid argument: dmt")
		
		disk.dump <- object@disk.dump
		disk.path <- object@disk.path
		
		#check if there is already an entry for the diff meth table
		if (region.type == "sites"){
			if (!is.null(object@sites[[comparison]])) {
				warning(paste("DiffMethTable already exists for combination:",comparison," , sites . Returning unmodified object"))
				return(object)
			}
		} else {
			if (!is.null(object@regions[[region.type]][[comparison]])) {
				warning(paste("DiffMethTable already exists for combination:",comparison,",",reg.type,". Returning unmodified object"))
				return(object)
			}
		}
		
		#check if the comparison is already there
		if(!is.element(comparison,object@comparisons)){
			object@comparisons <- c(object@comparisons,comparison)
			n.comps <- length(object@comparisons)
			names(object@comparisons)[n.comps] <- paste0("cmp",n.comps)
			object@comparison.grouplabels <- rbind(object@comparison.grouplabels,grp.labs)
			rownames(object@comparison.grouplabels) <- object@comparisons
			cmp.i <- length(object@comparisons)
			#append comparison info
			object@comparison.info <- c(object@comparison.info,list(NULL))
			names(object@comparison.info)[cmp.i] <- comparison
			#append the comparison to all sites and regions objects
			object@sites <- c(object@sites,list(NULL))
			names(object@sites)[cmp.i] <- comparison
			for (rr in object@region.types){
				object@regions[[rr]] <- c(object@regions[[rr]],list(NULL))
				names(object@regions[[rr]])[cmp.i] <- comparison
			}
		}
		cmp.i <- which(object@comparisons == comparison) #index of the comparison
		cmp.fname <- paste0("cmp",cmp.i)
		
		if (region.type == "sites"){
			reg.i <- 0
			table.obj <- as.matrix(dmt)
			rownames(table.obj) <- NULL #make certain that rownames do not exist. they are inefficient when using ff
			if (disk.dump){
				#create the ff matrix object
				fileN <- paste0(paste("sites",cmp.fname,sep="_"),".ff")
				table.obj <- ff(table.obj,dim=dim(table.obj),dimnames=dimnames(table.obj),filename=file.path(disk.path,fileN))
			}
			object@sites[[comparison]] <- table.obj

		} else {
			#check if the region type is already there
			if(!is.element(region.type,object@region.types)){
				object@region.types <- c(object@region.types,region.type)
				n.reg.types <- length(object@region.types)
				names(object@region.types)[n.reg.types] <- paste0("reg",n.reg.types)
				empty.cmp.list.4.regs <- rep(list(NULL),length(object@comparisons))#empty list of comparisons
				names(empty.cmp.list.4.regs) <- object@comparisons
				object@regions <- c(object@regions,list(empty.cmp.list.4.regs))
				names(object@regions)[length(object@regions)] <- region.type
			}
			reg.i <- which(object@region.types == region.type)
			reg.dir.name <- paste0("reg",reg.i)
			table.obj <- as.matrix(dmt)
			rownames(table.obj) <- NULL #make certain that rownames do not exist. they are inefficient when using ff
			if (disk.dump){
				#create the ff matrix object
				fileN <- paste0(paste("regions",reg.dir.name,cmp.fname,sep="_"),".ff")
				table.obj <- ff(table.obj,dim=dim(table.obj),dimnames=dimnames(table.obj),filename=file.path(disk.path,fileN))
			}
			object@regions[[region.type]][[comparison]] <- table.obj
		}
		
		return(object)
	}
)

if (!isGeneric("save.tables")) setGeneric("save.tables", function(object,...) standardGeneric("save.tables"))
#' save.tables-methods
#'
#' save the disk dumped tables to an ff archive for later reloading
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @param file path on the disk to save to.
#' @return success
#'
#' @rdname save.tables-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases save.tables
#' @aliases save.tables,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"),disk.dump=TRUE,disk.dump.dir=tempfile())
#' save.tables(dm,tempfile())
#' }
setMethod("save.tables", signature(object="RnBDiffMeth"),
	function(object,file){
		n.comps <- length(object@comparisons)
		n.region.types <- length(object@region.types)
		if (object@disk.dump){
			ee <- new.env()
			for (cci in 1:n.comps) {
				cc <- object@comparisons[cci]
				ccn <- paste0("cmp",cci)
				if (!is.null(object@sites[[cci]])){
					ee[[paste("sites",ccn,sep=".")]] <- object@sites[[cci]]
				}
				for (rri in 1:n.region.types){
					rr <- object@region.types[rri]
					rrn <- paste0("reg",rri)
					if (!is.null(object@regions[[rri]][[cci]])){
						ee[[paste("regions",rrn,ccn,sep=".")]] <- object@regions[[rri]][[cci]]
					}
				}
			}
			ffsave(list=ls(ee),envir=ee,rootpath=object@disk.path,file=file)
		}
	}
)
			
if (!isGeneric("reload")) setGeneric("reload", function(object,...) standardGeneric("reload"))
#' reload-methods
#'
#' reload disk dumped tables. Useful if the table files are manually copied or if the object is loaded again.
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @param save.file location of the ff data saved to disk (i.e. save in save.RData and save.ffData)
#' @param disk.path path on the disk for DMTs. can be new or be the same as in the original object
#' @return the updated RnBDiffMeth object
#'
#' @rdname reload-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases reload
#' @aliases reload,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' #compute differential methylation
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"),disk.dump=TRUE,disk.dump.dir=tempfile(pattern="working"))
#' #get temporary file names
#' fn.save.tabs <- tempfile(pattern="saveTables")
#' fn.save.obj  <- tempfile(pattern="saveObject")
#' #save the object and the tables to disk
#' save(dm,file=fn.save.obj)
#' save.tables(dm,fn.save.tabs)
#' #delete the object from the workspace
#' destroy(dm)
#' rm(dm)
#' #reload the object and tables
#' load(fn.save.obj)
#' dm.new <- reload(dm,fn.save.tabs)
#' }
setMethod("reload", signature(object="RnBDiffMeth"),
	function(object,save.file,disk.path=tempfile(pattern="diffmeth_", tmpdir=getOption('fftempdir'))){
		if (!object@disk.dump){
			warning("RnBDiffMeth object is not dumped to disk. Returning unmodified object")
			return(object)
		}
#		if (!file.exists(disk.path)) {
#			stop(paste("RnBDiffMeth dump directory does not exist:",disk.path))
#		}
		n.comps <- length(object@comparisons)
		n.region.types <- length(object@region.types)
		
		logger.start("Relinking RnBDiffMeth object")
		object@disk.path <- disk.path
		
		ee <- new.env(parent=emptyenv())
		ffload(save.file,rootpath=disk.path,envir=ee)
		
		for (cci in 1:n.comps) {
			cc <- object@comparisons[cci]
			ccn <- paste0("cmp",cci)
			site.obj.name <- paste("sites",ccn,sep=".")
			if (exists(site.obj.name,ee)){
				object@sites[[cci]] <- get(site.obj.name,ee)
			} else {
				logger.warning(c("Could not relink:","sites","--",cc))
				object@sites[cci] <- list(NULL)
			}
			for (rri in 1:n.region.types){
				rr <- object@region.types[rri]
				rrn <- paste0("reg",rri)
				reg.obj.name <- paste("regions",rrn,ccn,sep=".")
				if (exists(reg.obj.name,ee)){
					object@regions[[rri]][[cci]] <- get(reg.obj.name,ee)
				} else {
					logger.warning(c("Could not relink:",rr,"--",cc))
					object@regions[[rri]][cci] <- list(NULL)
				}
			}
		}
		logger.completed()
		
		return(object)
	}
)

#' save.rnb.diffmeth
#'
#' save an \code{\linkS4class{RnBDiffMeth}} object to disk
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @param path path on the disk to save to.
#' @author Fabian Mueller
#' @aliases save.rnb.diffmeth
#' @export
save.rnb.diffmeth <- function(object,path){
	if (file.exists(path)){
		stop("Saving unsuccessfull. Target directory already exists")
	}

	dir.create(path, showWarnings=FALSE, recursive=FALSE)
	save(object,file=file.path(path,"rnbDiffMeth.RData"))
	if (object@disk.dump){
		if(.Platform$OS == "windows" && Sys.getenv("R_ZIPCMD")==""){
			rnb.warning(c("Zip not found on this Windows system, this RnBDiffMeth object will not be saved.",
					"See the instructions for installing ZIP on Windows in the FAQ section of the RnBeads website."))
			return(invisible(FALSE))
		}
		save.tables(object,file.path(path,"rnbDiffMeth_tables"))
	}
}
#' load.rnb.diffmeth
#'
#' load a saved \code{\linkS4class{RnBDiffMeth}} object from disk
#'
#' @param path path of the saved object (a directory containing a corresponding \code{rnbDiffMeth.RData} file and possibly \code{rnbDiffMeth_tables} files)
#' @return the loaded \code{\linkS4class{RnBDiffMeth}} object
#' @author Fabian Mueller
#' @aliases load.rnb.diffmeth
#' @export
load.rnb.diffmeth <- function(path){
	if (!file.exists(path)){
		stop("Loading unsuccessfull. Path does not exist")
	}
	load.env <- new.env(parent=emptyenv())
	load(file.path(path,"rnbDiffMeth.RData"),envir=load.env)
	object <- get("object", load.env)
	if (object@disk.dump){
		object <- reload(object,file.path(path,"rnbDiffMeth_tables"))
	}
	return(object)
}

setGeneric("addComparisonInfo", function(object,...) standardGeneric("addComparisonInfo"))

#' addComparisonInfo-methods
#'
#' Adds a differential methylation comparison information.
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @param cmp.info	Comparison information as returned by \code{\link{get.comparison.info}}
#' @return the updated RnBDiffMeth object
#'
#' @rdname addComparisonInfo-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases addComparisonInfo,RnBDiffMeth-method
#' @noRd
setMethod("addComparisonInfo", signature(object="RnBDiffMeth"),
	function(object,cmp.info){
		if (!is.list(cmp.info)) {
			stop("invalid value for cmp.info. expected list")
		}
		if (length(cmp.info) < 1 || length(names(cmp.info)) < 1) {
			stop("invalid value for cmp.info: list is too short or is unnamed")
		}
		if (!all(names(cmp.info) %in% names(object@comparison.info))) {
			stop(paste("cmp.info does not match comparisons. The following comparison names are not",
					   "present in the RnBDiffMeth object: ",
					   paste(setdiff(names(cmp.info),names(object@comparison.info)),collapse=","))
			)
		}
		for (cn in names(cmp.info)){
			if (!is.null(cmp.info[[cn]])) {
				if (is.null(object@comparison.info[[cn]])) {
					object@comparison.info[[cn]] <- cmp.info[[cn]]
				} else {
					warning(paste("Did not overwrite comparison information for comparison:",cn,""))
				}
			}
		}
		return(object)
	}
)


if (!isGeneric("join.diffMeth")) setGeneric("join.diffMeth", function(obj1,obj2,...) standardGeneric("join.diffMeth"))
#' join.diffMeth-methods
#'
#' Merges two disjoint RnBDiffMeth objects into one. Disjoint here means, that no differential methylation table is specified in both
#' objects.
#'
#' @param obj1 \code{\linkS4class{RnBDiffMeth}} object. Its base properties will be used to create the joint object
#' 			   this is particularly imported for disk dumped objects as its path will be used and tables from the second
#' 			   object will be copied there
#' @param obj2 \code{\linkS4class{RnBDiffMeth}} object
#' @return the merged \code{\linkS4class{RnBDiffMeth}} object
#'
#' @note Caveat: if disk dumping is enabled the resulting object tables will be stored in the initial location of the first object to be joined
#' I.e. deleting the first object will lead to a broken joined object and deleting the joined object will lead to an broken first object.
#' @rdname join.diffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases join.diffMeth
#' @aliases join.diffMeth,RnBDiffMeth,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm1 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group"),region.types=c("genes","tiling"))
#' dm2 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"),region.types=c("promoters"))
#' dm.join1 <- join.diffMeth(dm1,dm2)
#' #the following joint object is invalid, because some region type - comparison combinations are missing
#' is.valid(dm.join1)
#' dm3 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Treatment"),region.types=c("genes","tiling"))
#' dm.join2 <- join.diffMeth(dm.join1,dm3)
#' #after joining the missing information, the new object is valid
#' is.valid(dm.join2)
#' }
setMethod("join.diffMeth", signature(obj1="RnBDiffMeth",obj2="RnBDiffMeth"),
	function(obj1,obj2){
		is.compatible <- (obj1@site.test.method == obj2@site.test.method) &&
						 (obj1@covg.thres == obj2@covg.thres) && 
						 (obj1@disk.dump == obj2@disk.dump)
		if (!is.compatible){
			stop("incompatible RnBDiffMeth objects")
		}
#		logger.start("Combining RnBDiffMeth objects")
		res <- obj1
		new.regs <- setdiff(obj2@region.types,obj1@region.types)
		n.new.regs <- length(new.regs)
		res@region.types <- c(res@region.types,new.regs)
		names(res@region.types) <- paste0("reg",1:length(res@region.types))
		
		new.comps <- setdiff(obj2@comparisons,obj1@comparisons)
		n.new.comps <- length(new.comps)
		res@comparisons <- c(res@comparisons,new.comps)
		names(res@comparisons) <- paste0("cmp",1:length(res@comparisons))
		new.comps.gns <- obj2@comparison.grouplabels[new.comps,]
		res@comparison.grouplabels <- rbind(res@comparison.grouplabels,new.comps.gns)
		rownames(res@comparison.grouplabels) <- res@comparisons
		
		#check for conflicting definitions, i.e. if anything is defined in both sets
		comps.intersect <- intersect(obj1@comparisons,obj2@comparisons)
		region.types.intersect <- intersect(obj1@region.types,obj2@region.types)
		#sites
		is.def.obj1 <- !vapply(obj1@sites[comps.intersect],is.null,FALSE)
		is.def.obj2 <- !vapply(obj2@sites[comps.intersect],is.null,FALSE)
		#comparison information: don't overwrite obj1
		cmp.info.names.obj1 <- names(obj1@comparison.info)
		cmp.info.names.obj2 <- names(obj2@comparison.info)
		cmp.info.names.intersect <- intersect(cmp.info.names.obj1,cmp.info.names.obj2)
		if (length(cmp.info.names.intersect)>0) {
			warning(paste("Join RnBDiffMeth: Comparison info defined in both objects:",cmp.info.names.intersect,
						  "--> keeping obj1"))
		}
		res@comparison.info <- rep(list(NULL),length(res@comparisons))
		names(res@comparison.info) <- res@comparisons
		res <- suppressWarnings(addComparisonInfo(res,obj1@comparison.info))
		res <- suppressWarnings(addComparisonInfo(res,obj2@comparison.info))

		#is any comparison defined on both objects: site level
		is.def.both <- (is.def.obj1+is.def.obj2)>1
		if (any(is.def.both)) {
			warning(paste("Join RnBDiffMeth: Comparison on site level defined in both objects:",comps.intersect[is.def.both],
						  "--> keeping obj1"))
		}
		#regions
		for (rr in region.types.intersect){
			is.def.obj1 <- !vapply(obj1@regions[[rr]][comps.intersect],is.null,FALSE)
			is.def.obj2 <- !vapply(obj2@regions[[rr]][comps.intersect],is.null,FALSE)
			#is any comparison defined on both objects
			is.def.both <- (is.def.obj1+is.def.obj2)>1
			if (any(is.def.both)) {
				warning(paste("Join RnBDiffMeth: Comparison on region (",rr,") level defined in both objects:",comps.intersect[is.def.both],
							  "--> keeping obj1"))
			}
		}
		#add empty lists for sites and regions
		new.comp.list <- rep(list(NULL),n.new.comps)
		names(new.comp.list) <- new.comps
		res@sites <- c(res@sites,new.comp.list)
		for (rr in obj1@region.types){
			res@regions[[rr]] <- c(res@regions[[rr]],new.comp.list) 
		}
		n.all.comps <- length(res@comparisons)
		new.reg.list <- rep(list(NULL),n.all.comps)
		names(new.reg.list) <- res@comparisons
		new.reg.list.list <- rep(list(new.reg.list),n.new.regs)
		names(new.reg.list.list) <- new.regs
		res@regions <- c(res@regions,new.reg.list.list)
		
		iis <- match(obj2@comparisons,res@comparisons)
		jjs <- match(obj2@region.types,res@region.types)
		for (i.old in 1:length(obj2@comparisons)){
			i.new <- iis[i.old]
			cc <- obj2@comparisons[i.old]
			
			cmp.fname <- paste0("cmp",i.new)
			if (!is.null(obj2@sites[[cc]]) && is.null(res@sites[[cc]])) {
				#if dumped, copy the matrix of obj2 to the location of obj1
				if (obj2@disk.dump) {
					base.path <- res@disk.path
					fileN <- paste0(paste("sites",cmp.fname,sep="_"),".ff")
					res@sites[[cc]] <- clone(obj2@sites[[cc]],vmode=vmode(obj2@sites[[cc]]),filename=file.path(base.path,fileN))
				} else {
					res@sites[[cc]] <- obj2@sites[[cc]]
				}
			}
			for (j.old in 1:length(obj2@region.types)){
				j.new <- jjs[j.old]
				rr <- obj2@region.types[j.old]
				reg.dir.name <- paste0("reg",j.new)
				if (!is.null(obj2@regions[[rr]][[cc]]) && is.null(res@regions[[rr]][[cc]])){
					#if dumped, copy the matrix of obj2 to the location of obj1
					if (obj2@disk.dump) {
						base.path <- res@disk.path
						fileN <- paste0(paste("regions",reg.dir.name,cmp.fname,sep="_"),".ff")
						res@regions[[rr]][[cc]] <- clone(obj2@regions[[rr]][[cc]],vmode=vmode(obj2@regions[[rr]][[cc]]),filename=file.path(base.path,fileN))
					} else {
						res@regions[[rr]][[cc]] <- obj2@regions[[rr]][[cc]]
					}
				}
			}
		}
#		logger.completed()
		return(res)
	}
)

if (!isGeneric("is.valid")) setGeneric("is.valid", function(object,...) standardGeneric("is.valid"))
#' is.valid-methods
#'
#' Validate an RnBDiffMeth object, ie. verify that all differential methylation tables are specified
#' and accounted for
#'
#' @param object \code{\linkS4class{RnBDiffMeth}} object
#' @param verbose print more info to the logger
#' @return TRUE iff all differential methylation tables are present and accounted for
#'
#' @rdname is.valid-RnBDiffMeth-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases is.valid
#' @aliases is.valid,RnBDiffMeth-method
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm1 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group"),region.types=c("genes","tiling"))
#' dm2 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"),region.types=c("promoters"))
#' dm.join1 <- join.diffMeth(dm1,dm2)
#' #the following joint object is invalid, because some region type - comparison combinations are missing
#' is.valid(dm.join1)
#' dm3 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Treatment"),region.types=c("genes","tiling"))
#' dm.join2 <- join.diffMeth(dm.join1,dm3)
#' #after joining the missing information, the new object is valid
#' is.valid(dm.join2)
#' }
setMethod("is.valid", signature(object="RnBDiffMeth"),
	function(object,verbose=FALSE){
		n.comps <- length(object@comparisons)
		n.region.types <- length(object@region.types)
		
		for (cci in 1:n.comps) {
			cc <- object@comparisons[cci]
			ccn <- paste0("cmp",cci)
			if (is.null(object@sites[[cc]])){
				if (verbose) logger.info(paste0("No table found for comparison '",cc,"' (sites)"))
				return(FALSE)
			}
			if (object@disk.dump){
				fileN <- file.path(object@disk.path,paste0(paste("sites",ccn,sep="_"),".ff"))
				if (!file.exists(fileN)){
					if (verbose) logger.info(paste0("Disk dump file ['",fileN,"'] not found for comparison '",cc,"' (sites)"))
					return(FALSE)
				}
			}
			for (rri in 1:n.region.types){
				rr <- object@region.types[rri]
				rrn <- paste0("reg",rri)
				if (is.null(object@regions[[rr]][[cc]])){
					if (verbose) logger.info(paste0("No table found for comparison '",cc,"' (region: '",rr,"')"))
					return(FALSE)
				}
				if (object@disk.dump){
					fileN <- file.path(object@disk.path,paste0(paste("regions",rrn,ccn,sep="_"),".ff"))
					if (!file.exists(fileN)){
						if (verbose) logger.info(paste0("Disk dump file ['",fileN,"'] not found for comparison '",cc,"' (region: '",rr,"')"))
						return(FALSE)
					}
				}
			}
		}
		return(TRUE)
	}
)
