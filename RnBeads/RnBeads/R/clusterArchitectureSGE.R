################################################################################
# Cluster Architecture Descriptions
################################################################################
################################################################################
# Concrete implementations for the SGE environment
# have to implement:
# Method: getSubCmdTokens
################################################################################

#' ClusterArchitectureSGE Class
#'
#' A child class of \code{\linkS4class{ClusterArchitecture}} implementing specifications of Sun Grid Engine (SGE) architectures.
#'
#' @details
#' Follow this template if you want to create your own ClusterArchitecture class.
#'
#' @section Slots:
#' see \code{\linkS4class{ClusterArchitecture}}
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getSubCmdTokens,ClusterArchitectureSGE-method}}}{Returns a vector of command line tokens corresponding to submitting
#'   a job with the given command to the cluster}
#' }
#'
#' @name ClusterArchitectureSGE-class
#' @rdname ClusterArchitectureSGE-class
#' @author Fabian Mueller
#' @exportClass ClusterArchitecture
setClass("ClusterArchitectureSGE",
	contains = "ClusterArchitecture"
)

#' initialize.ClusterArchitectureSGE
#'
#' Initialize an ClusterArchitecture object for a Sun Grid Engine (SGE)
#' 
#' @param name A name or identifier
#' @param ... arguments passed on to the constructor of \code{\linkS4class{ClusterArchitecture}} (the parent class)
#'
#' @export
#' @author Fabian Mueller
#' @docType methods
setMethod("initialize","ClusterArchitectureSGE",
	function(
		.Object,
		name="ClusterArchitectureSGE",
		...
	) {
		.Object <- callNextMethod(.Object=.Object, name=name, ...)
		.Object <- setExecutable(.Object,"R","R")
		.Object <- setExecutable(.Object,"Rscript","Rscript")
		.Object <- setExecutable(.Object,"python","python")
		.Object@getSubCmdTokens.optional.args <- c("sub.binary","quote.cmd")
		.Object
	}
)

#' getSubCmdTokens-methods
#'
#' Returns a string for the of command line corresponding to submitting
#' a job with the given command to the cluster.
#' @details
#' For a concrete child class implementation for a sun grid architecture specification see \code{\linkS4class{ClusterArchitectureSGE}}
#'
#' @param object \code{\linkS4class{ClusterArchitectureSGE}} object
#' @param cmd.tokens a character vector specifying the executable command that should be wrapped in the cluster submission command
#' @param log file name and path of the log file that the submitted job writes to
#' @param job.name name of the submitted job
#' @param res.req character vector specifying required resources. The resource requirements should be the values of the vector,
#'                the names should specify the resource name
#' @param depend.jobs character vector containg names or ids of jobs the submitted job will depend on.
#' @param sub.binary treat the command as binary (see \code{-b} flag of \code{qsub} of the SGE documentation)
#' @param quote.cmd Flag indicating whether the submitted cammed should also be wrapped in quotes
#' @return A character vector containing the submission command tokens
#'
#' @rdname getSubCmdTokens-ClusterArchitectureSGE-methods
#' @docType methods
#' @aliases getSubCmdTokens,ClusterArchitectureSGE-method
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' arch <- new("ClusterArchitectureSGE",
#' 	name="my_sge_architecture"
#' )
#' getSubCmdTokens(arch,c("Rscript","my_great_script.R"),"my_logfile.log")
#' }
setMethod("getSubCmdTokens",
	signature(
		object="ClusterArchitectureSGE"
	),
	function(
		object,
		cmd.tokens,
		log,
		job.name = "",
		res.req = character(0),
		depend.jobs = character(0),
		sub.binary = TRUE,
		quote.cmd = TRUE
	) {
		res.req.tokens <- NULL
		if (length(res.req) > 0){
			if (any(is.na(names(res.req))) || is.null(names(res.req))){
				stop("Invalid resource requirement specification. Need names for all requirements")
			}
			rr <- as.vector(rbind(rep("-l",length(res.req)),paste0(names(res.req),"=",res.req)))
			res.req.tokens <- c(res.req.tokens,rr)
		}
		log.token <- NULL
		if (nchar(log)>0) {
			log.token <- c("-o",log)
		}
		job.name.token <- NULL
		if (nchar(job.name)>0) {
			job.name.token <- c(job.name.token,"-N",job.name)
		}
		dependency.token <- NULL
		if (length(depend.jobs)>0){
			dependency.token <- c(dependency.token, "-hold_jid", paste0(paste(depend.jobs,collapse=",")))
		}
		if (quote.cmd){
			cmd.tokens <- paste0("'",paste(cmd.tokens,collapse=" "),"'")
		}
		bin.token <- ifelse(sub.binary,"y","n")
		res <- c(
			"qsub",
			res.req.tokens,
			log.token,
			"-j","y",
			job.name.token,
			dependency.token,
			"-b",bin.token,
			cmd.tokens
		)
		return(res)
	}
)
