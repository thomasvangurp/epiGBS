########################################################################################################################
## parallelProcessing.R
## created: 2012-05-04
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for parallel processing.
########################################################################################################################

## G L O B A L S #######################################################################################################

.parallel <- new.env()
.parallel[["do.par"]] <- FALSE
.parallel[["num.cores"]] <- -1

## F U N C T I O N S ###################################################################################################

#' parallel.setup
#'
#' Sets up parallel processing. Requires the \pkg{foreach} and \pkg{doParallel} packages
#' @param ... Parameters for \code{registerDoParallel} from the \pkg{doParallel} package.
#' 			  This allows, for instance, for specificying the number of workers.
#' @return \code{TRUE} (invisible) to indicate that parallelization is set up.
#' @note Requires the packages \pkg{foreach} and \pkg{doParallel}.
#'
#' @author Fabian Mueller
#' @export parallel.setup
#' @examples
#' \dontrun{
#' parallel.setup(2)
#' }
parallel.setup <- function(...){
	logger.start("Setting up Multicore")
	require(foreach)
	require(doParallel)
	registerDoParallel(...)
	.parallel[["num.cores"]] <- getDoParWorkers()
	.parallel[["do.par"]] <- TRUE
	logger.info(c("Using",.parallel[["num.cores"]],"cores"))
	logger.completed()
	invisible(.parallel[["do.par"]])
}

########################################################################################################################

#' parallel.disable
#'
#' Disables parallel processing.
#'
#' @author Fabian Mueller
#' @export parallel.disable
#' @examples
#' \dontrun{
#' parallel.getNumWorkers()
#' parallel.setup(2)
#' parallel.getNumWorkers()
#' parallel.disable()
#' parallel.getNumWorkers()
#' }
parallel.disable <- function(){
	.parallel[["do.par"]] <- FALSE
	.parallel[["num.cores"]] <- -1
	invisible(TRUE)
}

########################################################################################################################

#' parallel.getNumWorkers
#'
#' Return the number of workers used for parallel processing.
#'
#' @author Fabian Mueller
#' @export parallel.getNumWorkers
#' @examples
#' \dontrun{
#' parallel.getNumWorkers()
#' parallel.setup(2)
#' parallel.getNumWorkers()
#' parallel.disable()
#' parallel.getNumWorkers()
#' }
parallel.getNumWorkers <- function(){
	return(.parallel[["num.cores"]])
}

########################################################################################################################

#' parallel.isEnabled
#'
#' Checks if whether parallel processing is enabled.
#'
#' @return \code{TRUE} if multicore processing is enabled, \code{FALSE} otherwise.
#'
#' @author Fabian Mueller
#' @export parallel.isEnabled
#' @examples
#' \dontrun{
#' parallel.isEnabled()
#' parallel.setup(2)
#' parallel.isEnabled()
#' parallel.disable()
#' parallel.isEnabled()
#' }
parallel.isEnabled <- function(){
	return(.parallel[["do.par"]])
}
