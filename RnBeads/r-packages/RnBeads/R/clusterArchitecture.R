################################################################################
# Cluster Architecture Descriptions
################################################################################
# General Virtual class template
################################################################################
#' ClusterArchitecture Class
#'
#' A virtual class for storing specifications of architectures for different compute clusters. It is designed to
#' let other classes inherit from it
#' 
#' @details
#' For a concrete child class for a sun grid architecture specification see \code{\linkS4class{ClusterArchitectureSGE}}
#' If you want to implement your own child class be sure to at least implement the following functions:
#' \code{\link{getSubCmdTokens,ClusterArchitecture-method}}.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{name}}{A name or identifier}
#'   \item{\code{executables}}{A NAMED character vector of executables that can be used by the cluster. For instance, the \code{R} executable is important} 
#'   \item{\code{getSubCmdTokens.optional.args}}{character vector containing the valid optional arguments to the \code{\link{getSubCmdTokens,ClusterArchitecture-method}} function.}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getSubCmdTokens,ClusterArchitecture-method}}}{Returns a vector of command line tokens corresponding to submitting
#'   a job with the given command to the cluster}
#'   \item{\code{\link{getSubCmdStr,ClusterArchitecture-method}}}{Returns a string for the of command line corresponding to submitting
#'   a job with the given command to the cluster}
#'   \item{\code{\link{setExecutable,ClusterArchitecture,character,character-method}}}{Tells the cluster architecture about an executable that can be submitted as job}
#'   \item{\code{\link{getExecutable,ClusterArchitecture,character-method}}}{Gets the location of an executable associated with a name}
#' }
#'
#' @name ClusterArchitecture-class
#' @rdname ClusterArchitecture-class
#' @author Fabian Mueller
#' @exportClass ClusterArchitecture
setClass("ClusterArchitecture",
	slots = list(
		name = "character",
		executables = "character",
		getSubCmdTokens.optional.args = "character"
	),
	contains = "VIRTUAL"
)

#' initialize.ClusterArchitecture
#'
#' Initialize an ClusterArchitecture object
#' 
#' @param name A name or identifier
#'
#' @author Fabian Mueller
#' @docType methods
setMethod("initialize","ClusterArchitecture",
	function(
		.Object,
		name="ClusterArchitecture"
	) {
		.Object@name           <- name
		.Object@executables    <- character(0)
		.Object@getSubCmdTokens.optional.args <- character(0)
		.Object
	}
)

#implement for each inheriting class
if (!isGeneric("getSubCmdTokens")) {
	setGeneric(
		"getSubCmdTokens",
		function(object, cmd.tokens, log, job.name="", res.req=character(0), depend.jobs=character(0), ...) standardGeneric("getSubCmdTokens"),
		signature=c("object","cmd.tokens","log","job.name","res.req","depend.jobs")
	)
}
#' getSubCmdTokens-methods
#'
#' Returns a string for the of command line corresponding to submitting
#' a job with the given command to the cluster.
#' @details
#' For a concrete child class implementation for a sun grid architecture specification see \code{\link{getSubCmdTokens,ClusterArchitectureSGE-method}}
#'
#' @param object \code{\linkS4class{ClusterArchitecture}} object
#' @param cmd.tokens a character vector specifying the executable command that should be wrapped in the cluster submission command
#' @param log file name and path of the log file that the submitted job writes to
#' @param job.name name of the submitted job
#' @param res.req character vector specifying required resources. The resource requirements should be the values of the vector,
#'                the names should specify the resource name
#' @param depend.jobs character vector containg names or ids of jobs the submitted job will depend on.
#' @return A character vector containing the submission command tokens
#'
#' @rdname getSubCmdTokens-ClusterArchitecture-methods
#' @docType methods
#' @aliases getSubCmdTokens
#' @aliases getSubCmdTokens,ClusterArchitecture-method
#' @author Fabian Mueller
#' @export
setMethod("getSubCmdTokens",
	signature(
		object="ClusterArchitecture"
	),
	function(
		object,
		cmd.tokens,
		log,
		job.name = "",
		res.req = character(0),
		depend.jobs = character(0)
	) {
		return(character(0))
	}
)

if (!isGeneric("getSubCmdStr")) setGeneric("getSubCmdStr", function(object,...) standardGeneric("getSubCmdStr"))
#' getSubCmdStr-methods
#'
#' Returns a string for the of command line corresponding to submitting
#' a job with the given command to the cluster.
#'
#' @param object \code{\linkS4class{ClusterArchitecture}} object
#' @param ... arguments passed on to \code{\link{getSubCmdTokens,ClusterArchitecture-method}}
#' @return A string containing the submission command
#'
#' @rdname getSubCmdStr-ClusterArchitecture-methods
#' @docType methods
#' @aliases getSubCmdStr
#' @aliases getSubCmdStr,ClusterArchitecture-method
#' @author Fabian Mueller
#' @export
setMethod("getSubCmdStr","ClusterArchitecture",
	function(
		object,
		...
	) {
		res <- getSubCmdTokens(object,...)
		res <- paste(res, collapse=" ")
		return(res)
	}
)

if (!isGeneric("setExecutable")) setGeneric("setExecutable", function(object,exec.name,exec.loc) standardGeneric("setExecutable"))
#' setExecutable-methods
#'
#' Tells the cluster architecture about an executable that can be submitted as job
#'
#' @param object \code{\linkS4class{ClusterArchitecture}} object
#' @param exec.name A name/identifier that will be associated with the given executable
#' @param exec.loc The executable's location
#' @return The modified object
#'
#' @rdname setExecutable-ClusterArchitecture-methods
#' @docType methods
#' @aliases setExecutable
#' @aliases setExecutable,ClusterArchitecture-method
#' @author Fabian Mueller
#' @export
setMethod("setExecutable",
	signature(
		object="ClusterArchitecture",
		exec.name="character",
		exec.loc="character"
	),
	function(
		object,
		exec.name,
		exec.loc
	) {
		object@executables[exec.name] <- exec.loc
		return(object)
	}
)

if (!isGeneric("getExecutable")) setGeneric("getExecutable", function(object,exec.name) standardGeneric("getExecutable"))
#' getExecutable-methods
#'
#' Retrieves the executable associated with a name/identifier
#'
#' @param object \code{\linkS4class{ClusterArchitecture}} object
#' @param exec.name The executable's name/identifier
#' @return The executable. If the name is not associated with any executable, the names will be returned and a warning will be raised
#'
#' @rdname getExecutable-ClusterArchitecture-methods
#' @docType methods
#' @aliases getExecutable
#' @aliases getExecutable,ClusterArchitecture-method
#' @author Fabian Mueller
#' @export
setMethod("getExecutable",
	signature(
		object="ClusterArchitecture",
		exec.name="character"
	),
	function(
		object,
		exec.name
	) {
		res <- object@executables[exec.name]
		if (is.na(res)){
			warning(paste0("Executable '",exec.name,"' not found. Returning the name as value"))
			res <- exec.name
		}
		return(res)
	}
)
