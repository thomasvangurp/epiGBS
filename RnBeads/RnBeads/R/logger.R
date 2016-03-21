########################################################################################################################
## logger.R
## created: 2012-03-28
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Logging functionality of the RnBeads modules.
########################################################################################################################

## G L O B A L S #######################################################################################################

.LOG.INFO <- new.env()

.LOG.INDENT <- "    "
.LOGGER <- "RnBeads"
.LOG.STAT <- c("INFO" = "INFO", "STATUS" = "STATUS", "WARNING" = "WARNING", "ERROR" = "ERROR")

## F U N C T I O N S ###################################################################################################

## get.memory.usage
##
## Gets the memory used by this process.
##
## @return Memory, in Gigabytes, used in this R session.
## @details
## In Windows, the returned value measures only the memory allocated by this session. It does not include the memory
## used by the R system itself (and many of the loaded libraries). In Linux, the returned value is the total memory used
## by the R process that runs this script.
## @author Yassen Assenov
get.memory.usage <- function() {
	if (.Platform$OS == "windows") {
		return(memory.size() / 1024)
	}
	processinfo <- paste("/proc", Sys.getpid(), "status", sep = "/")
	processinfo <- scan(processinfo, what = character(), sep = "\n", quiet = TRUE)
	memused <- grep("^VmSize\\:(.+)kB", processinfo, value = TRUE)
	memused <- as.double(gsub("\\s+", "", substr(memused, 8, nchar(memused) - 2))) / 1048576
	return(memused)
}

########################################################################################################################

## get.disk.usage
##
## Gets the space used by all files and subdirectories of the given path.
##
## @param path Base directory to scan.
## @return Combined size, in Gigabytes, used by the files in the given path and in its subdirectories. The disk
##         space used by the directory entries themselves is not included.
## @author Yassen Assenov
get.disk.usage <- function(path = getOption("fftempdir")) {
	if (!isTRUE(file.info(path)[, "isdir"])) {
		stop("invalid value for path; expected existing directory")
	}
	sum(file.info(dir(path, full.names = TRUE, recursive = TRUE))[, "size"] / 1048576) / 1024
}

########################################################################################################################

format.usage <- function(usage.extractor) {
	tryCatch(format(usage.extractor(), digits = 1L, nsmall = 1L, justify = "right", width = 6),
		error = function(e) { "      " })
}

########################################################################################################################

## logger.format
##
## Memory usage formatter for the RnBeads text logger.
##
## @param record Log message in the form of a record. This message is expected to start with a status word, e.g. with
##               \code{WARNING}.
## @return       Character containing the formatted message.
logger.format <- function(record) {

	memory.used <- ifelse(.LOG.INFO[["memory"]], paste0(format.usage(get.memory.usage), " "), "")
	disk.used <- ifelse(.LOG.INFO[["disk"]], paste0(format.usage(get.disk.usage), " "), "")

	## Right-align the status word
	status.word <- sub("^([A-Z]+) .+", "\\1", record$msg)
	if (status.word %in% .LOG.STAT) {
		prefix <- paste(rep(" ", max(nchar(.LOG.STAT)) - nchar(status.word)), collapse = "")
	} else {
		prefix <- ""
	}
	## Prepend date, time, memory used and temporary disk space used
	paste0(record$timestamp, " ", memory.used, disk.used, prefix, record$msg)
}

########################################################################################################################

## Transforms a vector of message elements to a character.
##
## @param txt    Character vector with message elements.
## @param indent Flag indicating if the message must be prepended by \code{.LOG.INDENT} to indicate it belongs to
##               a specific section.
## @return       Single-element character vector concatenating all elements in \code{txt}, possibly with indentation.
logger.transform <- function(txt, indent = TRUE) {
	if (!(is.character(txt) && is.vector(txt))) {
		stop("invalid value for parameter txt")
	}
	if (length(txt) > 1) {
		txt <- paste(txt, collapse = " ")
	}
	if (indent && length(.LOG.INFO[["titles"]]) != 0) {
		txt <- paste(paste(rep(.LOG.INDENT, length(.LOG.INFO[["titles"]])), collapse = ""), txt, sep = "")
	}
	return(txt)
}

########################################################################################################################

logger.paste <- function(status.word, txt, logger = .LOGGER) {
	record <- list(
		timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S "),
		msg = paste(.LOG.STAT[status.word], logger.transform(txt)))
	txt <- paste(.LOG.INFO[["formatter"]](record), "\n", sep = "")
	for (fname in .LOG.INFO[["logfiles"]]) {
		cat(txt, file = ifelse(is.na(fname), "", fname), append = TRUE)
	}
}

########################################################################################################################

#' logger.getfiles
#'
#' Gets the files currently used by the logger.
#'
#' @return Vector storing the full names of the files that are being used by the logger. This vector contains \code{NA}
#'         as an element if the logger is (also) using the console for its output. If logging functionality is disabled
#'         (see \code{\link{rnb.options}}) or the logger is not initialized, this function returns \code{NULL}.
#'
#' @examples
#' \dontrun{
#' if (NA %in% logger.getfiles())
#'   cat("Console logger is enabled\n")
#' }
#' @seealso \code{\link{logger.isinitialized}} to check if logging is activated;
#'   \code{\link{logger.start}} for initializing a logger or starting a section
#' @author Yassen Assenov
#' @export
logger.getfiles <- function() {
	if (rnb.getOption("logging") && exists("logfiles", envir = .LOG.INFO, inherits = FALSE)) {
		return(get("logfiles", envir = .LOG.INFO, inherits = FALSE))
	}
	return(NULL)
}

########################################################################################################################

#' logger.isinitialized
#'
#' Checks if the logger is initialized.
#'
#' @return \code{TRUE} if the logger was initialized and is in use; \code{FALSE} otherwise.
#'
#' @examples
#' \dontrun{
#' if (!logger.isinitialized())
#'   logger.start(fname = NA)
#' }
#' @seealso \code{\link{logger.start}} for initializing a logger or starting a section
#' @author Yassen Assenov
#' @export
logger.isinitialized <- function() {
	if (rnb.getOption("logging")) {
		return(exists("logfiles", envir = .LOG.INFO, inherits = FALSE))
	}
	return(FALSE)
}

########################################################################################################################

#' Writing text messages to the log file.
#'
#' Appends a single-line status message to the log text file. The message is prepended by its type, which is one of
#' \code{STATUS}, \code{INFO}, \code{WARNING} or \code{ERROR}.
#'
#' @rdname loggerMessages
#' @aliases logger.status
#' @aliases logger.info
#' @aliases logger.warning
#' @aliases logger.error
#'
#' @param txt       Text to add to the log file. This must be a \code{character} vector; its elements are concatenated
#'                  using a single space (\code{" "}) as a separator.
#' @param terminate Flag indicating if the execution is to be terminated after this error message is added to the log.
#'
#' @examples
#' \dontrun{
#' if (!logger.isinitialized())
#'   logger.start(fname = NA)
#' logger.status(c("Reached step", 2))
#' logger.info(c("Provided email:", rnb.getOption("email")))
#' }
#' @seealso \code{\link{logger.isinitialized}}to check if logging is activated;
#'   \code{\link{logger.start}} for initializing a logger or starting a section
#' @author Yassen Assenov
#' @export
logger.status <- function(txt) {
	if (rnb.getOption("logging")) {
		if (!logger.isinitialized()) {
			stop("logger is not initialized")
		}
		logger.paste("STATUS", txt)
	}
}

########################################################################################################################

#' @rdname loggerMessages
#' @export
logger.info <- function(txt) {
	if (rnb.getOption("logging")) {
		if (!logger.isinitialized()) {
			stop("logger is not initialized")
		}
		logger.paste("INFO", txt)
	}
}

########################################################################################################################

#' @rdname loggerMessages
#' @export
logger.warning <- function(txt) {
	if (rnb.getOption("logging")) {
		if (!logger.isinitialized()) {
			stop("logger is not initialized")
		}
		logger.paste("WARNING", txt)
	}
}

########################################################################################################################

#' @rdname loggerMessages
#' @export
logger.error <- function(txt, terminate = rnb.getOption("logging.exit.on.error")) {
	if (rnb.getOption("logging")) {
		if (!logger.isinitialized()) {
			stop("logger is not initialized")
		}
		logger.paste("ERROR", txt)
		if (terminate) {
			for (logfile in get("logfiles", envir = .LOG.INFO, inherits = FALSE)) {
				if (is.na(logfile)) {
					cat("\n")
				} else {
					cat("\n", file = logfile, append = TRUE)
				}
			}
			quit(save = "no", status = 1L)
		}
		stop(logger.transform(txt, FALSE))
	}
}

########################################################################################################################

#' Log File Management
#'
#' Functions for logger management.
#'
#' @rdname loggerManagement
#' @aliases logger.start
#' @aliases logger.completed
#' @aliases logger.close
#'
#' @param txt   Description to add to the log file. The words \code{STARTED} and \code{COMPLETED} are prepended to the
#'              message upon initialization and completion of the section, respectively.
#' @param fname Name of the log file and/or console. Note that at most one file name can be specified. The function
#'              \code{logger.start} normalizes the given name, that is, it converts it to an absolute name. If this
#'              parameter is \code{NA}, logger messages are printed to the console. If it is a two-element vector
#'              containing one file name and \code{NA}, the logger is (re)initialized to print messages both to the
#'              given file name and the console. A value of \code{NULL} (default) indicates the logger should continue
#'              using the previously specified file.
#'
#' @examples
#' \dontrun{
#' if (!logger.isinitialized())
#'   logger.start(fname = NA)
#' logger.start("Tests for Significance")
#' logger.completed()
#' logger.close()
#' }
#' @section Details:
#' \code{logger.start} initializes the logger and/or starts a new section. \code{logger.completed} completes the last
#' (innermost) open section in the log. \code{logger.close} deinitializes the logger. Note that after reinitialization
#' or deinitialization, the information about the current output file, as well as any open sections, is deleted.
#'
#' @seealso logger.isinitialized
#' @author Yassen Assenov
#' @export
logger.start <- function(txt = character(0), fname = NULL) {
	if (rnb.getOption("logging")) {
		if (!logger.isinitialized()) {
			if (is.null(fname)) {
				stop("logger is not initialized")
			}
		}
		if (!is.null(fname)) {
			if (!((length(fname) == 1 && (is.character(fname) || is.na(fname))) ||
				(length(fname) == 2 && (is.character(fname) && sum(is.na(fname)) == 1 &&
						sum(fname == "", na.rm = TRUE) == 0)))) {
				stop("invalid value for fname")
			}
			if (exists("logfiles", envir = .LOG.INFO, inherits = FALSE)) {
				logger.close()
			}
			.LOG.INFO[["memory"]] <- rnb.getOption("logging.memory")
			.LOG.INFO[["disk"]] <- (rnb.getOption("logging.disk") && rnb.getOption("logging.disk"))
			.LOG.INFO[["formatter"]] <- logger.format
			for (i in 1:length(fname)) {
				if (!is.na(fname[i])) {
					fname[i] <- normalizePath(fname[i], mustWork = FALSE)
				}
			}
			.LOG.INFO[["logfiles"]] <- as.character(fname)
			.LOG.INFO[["titles"]] <- character(0)
		}
		txt <- logger.transform(txt, indent = FALSE)
		if (length(txt) != 0) {
			logger.status(paste("STARTED", txt))
			.LOG.INFO[["titles"]] <- c(.LOG.INFO[["titles"]], txt)
		}
	}
}

########################################################################################################################

#' @rdname loggerManagement
#' @export
logger.completed <- function() {
	if (rnb.getOption("logging")) {
		if (!logger.isinitialized()) {
			stop("logger is not initialized")
		}
		N <- length(.LOG.INFO[["titles"]])
		if (length(.LOG.INFO[["titles"]]) == 0) {
			logger.error("No section to complete")
		}
		txt <- paste("COMPLETED ", .LOG.INFO[["titles"]][N], ifelse(N == 1, "\n", ""), sep = "")
		.LOG.INFO[["titles"]] <- .LOG.INFO[["titles"]][-N]
		logger.status(txt)
	}
}

########################################################################################################################

## logger.addfile
##
## Adds a new file (or console) to contain the messages of the logger.
##
## @param fname Name of the log file. This function normalizes the given file name, that is, it converts it to an
##              absolute name. Set this to \code{NA} in order to print log messages to the console.
##
## @seealso \code{\link{logger.start}} for re-initializing the logger
##
## @author Yassen Assenov
logger.addfile <- function(fname) {
	if (rnb.getOption("logging")) {
		if (!(length(fname) == 1 && (is.na(fname) || is.character(fname)))) {
			stop("invalid value for fname")
		}
		if (!logger.isinitialized()) {
			stop("logger is not initialized")
		}
		logfiles <- get("logfiles", envir = .LOG.INFO, inherits = FALSE)
		if (!(fname %in% logfiles)) {
			if (any(sapply(logfiles, function(x) { is.character(x) && (!is.na(x)) }))) {
				## Adding a file when there is already a file
				stop("logger is already initialized to file")
			}
			.LOG.INFO[["logfiles"]] <- c(logfiles, as.character(fname))
		}
	}
}

########################################################################################################################

#' @rdname loggerManagement
#' @export
logger.close <- function() {
	if (rnb.getOption("logging")) {
		rm(list = ls(envir = .LOG.INFO), envir = .LOG.INFO)
	}
}

########################################################################################################################

#' logger.argument
#'
#' Reads a command-line argument supplied to a script.
#'
#' @param arg.names       \code{character} vector of acceptable argument names. This function scans the provided
#'                        arguments and performs a case insensitive match.
#' @param full.name       One-element \code{character} vector giving the argument's full name or description. This is
#'                        used in a log message in case of an error.
#' @param arg.type        Variable type of the argument. Must be one of \code{"character"}, \code{"logical"},
#'                        \code{"integer"}, \code{"double"}, \code{"numeric"} or \code{"real"}. The last three types are
#'                        all synonyms.
#' @param accepted.values Vector of accepted values for the argument. This must be of the type given in \code{arg.type}.
#'                        Set this to \code{NULL} if there are no restrictions on the argument values.
#' @param default         Default value for the argument in case it is not specified. Setting this to \code{NULL} makes
#'                        the argument required, that is, an error is generated if the argument is not specified. Set
#'                        this to \code{NA} if is not a required argument and it shouldn't default to a specific value.
#'                        Otherwise, if \code{accepted.values} is provided, this must be one of its elements.
#' @param arg.list        Vector of arguments provided at the execution of the script. The arguments should be provided
#'                        as \emph{name=value} pairs.
#' @return                Argument's value, or \code{NULL} if such is not provided.
#'
#' @details This is convenience function for reading parameters supplied to the script in the form \emph{name = value}.
#'          It expects that logging is enabled (see \code{\link{rnb.options}}). The function fails if this condition is
#'          not met.
#' @examples
#' \dontrun{
#' n.iterations <- logger.argument("iterations", "number of iterations", "integer",
#'   accepted.values = 1:100, default = 1L)
#' logger.close()
#' }
#' @author Yassen Assenov
#' @export
logger.argument <- function(arg.names, full.name, arg.type = "character", accepted.values = NULL, default = NULL,
	arg.list = commandArgs()) {
	if (rnb.getOption("logging") == FALSE) {
		stop("logging is disabled")
	}
	if (!(is.character(arg.names) && length(arg.names) != 0 && (!any(is.na(arg.names))))) {
		stop("invalid value for arg.names")
	}
	if (!(is.character(full.name) && length(full.name) == 1 && (!is.na(full.name)))) {
		stop("invalid value for full.name")
	}
	if (!(is.character(arg.type) && length(arg.type) == 1 && (!is.na(arg.type)))) {
		stop("invalid value for arg.type")
	}
	if (!(arg.type %in% c("character", "logical", "integer", "double", "numeric", "real"))) {
		stop("invalid value for arg.type")
	}
	if (!(is.null(default) || length(default) == 1)) {
		stop("invalid value for default; expected single value or NULL")
	}
	if (!(is.null(accepted.values) || is.null(default))) {
		accepted.values <- c(NA, accepted.values)
		if (!(default %in% accepted.values)) {
			stop("invalid value for default; expected one of the accepted values")
		}
	}
	if (!(is.character(arg.list) && length(arg.list) != 0 && (!any(is.na(arg.list))))) {
		stop("invalid value for arg.list")
	}

	if (!logger.isinitialized()) {
		logger.start(fname = NA)
	}
	value <- NULL
	wrongType <- function(e) {
		logger.error(c("Invalid value for argument", full.name, ":", value))
	}
	arg.names <- tolower(arg.names)
	for (arg in strsplit(arg.list, split = "=", fixed = TRUE)) {
		if (length(arg) != 2) {
			next
		}
		arg[1] <- tolower(arg[1])
		if (arg[1] %in% arg.names) {
			if (!is.null(value)) {
				logger.error(c("Redefinition of", full.name))
			}
			## Perform type check
			value <- arg[2]
			if (arg.type == "character") {
				## No type conversion necessary
			} else if (arg.type == "logical") {
				value <- as.logical(toupper(arg[2]))
				if (is.na(value)) {
					value <- arg[2]
					wrongType("")
				}
			} else if (arg.type == "integer") {
				value <- tryCatch(as.integer(arg[2]), warning = wrongType, error = wrongType)
			} else { # arg.type %in% c("double", "numeric", "real")
				value <- tryCatch(as.numeric(arg[2]), warning = wrongType, error = wrongType)
			}
		}
	}
	if (is.null(value)) {
		if (is.null(default)) {
			logger.error(c("Missing required argument", full.name))
		}
		if (is.na(default)) {
			if (arg.type == "character") {
				value <- as.character(NA)
			} else if (arg.type == "logical") {
				value <- as.logical(NA)
			} else if (arg.type == "integer") {
				value <- as.integer(NA)
			} else {
				value <- as.numeric(NA)
			}
		} else {
			value <- default
		}
	}
	if (!(is.null(accepted.values) || (value %in% accepted.values))) {
		wrongType("")
	}
	return(value)
}

########################################################################################################################

#' logger.validate.file
#'
#' Validates the specified file or directory exists. Prints an error or a warning message to the log if it does not
#' exist, it is not of the accepted type or is not accessible.
#'
#' @param file      Name of file or directory to validate.
#' @param is.file   Flag indicating if the given name must denote an existing file. If this is \code{FALSE}, the given
#'                  name must denote a directory. Set this to \code{NA} if both types are an acceptable scenario.
#' @param terminate Flag indicating if the execution is to be terminated in case the validation fails. This parameter
#'                  determines if an error message (\code{terminate} is \code{TRUE}) or a warning message
#'                  (\code{terminate} is \code{FALSE}) is to be sent to the log when the specified file or directory
#'                  does not exist, is not of the accepted type or is not accessible.
#' @return Whether the validation succeeded or not, invisibly. Note that when \code{terminate} is \code{TRUE} and the
#'         validation fails, the R session is closed and thus no value is returned.
#'
#' @examples
#' \dontrun{
#' if (!logger.isinitialized())
#'   logger.start(fname = NA)
#' # Validate the current working directory exists
#' logger.validate.file(getwd(), FALSE)
#' }
#' @author Yassen Assenov
#' @export
logger.validate.file <- function(file, is.file = TRUE, terminate = TRUE) {
	if (!(is.character(file) && length(file) == 1 && (!is.na(file)))) {
		stop("invalid value for file; expected single character")
	}
	if (!parameter.is.flag(is.file)) {
		stop("invalid value for is.file; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(terminate)) {
		stop("invalid value for terminate; expected TRUE or FALSE")
	}
	file <- file[1]
	is.file <- is.file[1]
	terminate <- terminate[1]
	if (rnb.getOption("logging")) {
		if (!file.exists(file)) {
			msg <- ifelse(is.na(is.file), "File / directory", ifelse(is.file, "File", "Directory"))
			msg <- c(msg, "not found:", file)
			if (terminate) {
				logger.error(msg)
			}
			logger.warning(msg)
			return(invisible(FALSE))
		}
		if (!is.na(is.file)) {
			is.dir <- file.info(file)[1, "isdir"]
			if (is.file == is.dir) {
				msg <- c(file, "is a", ifelse(is.dir, "directory", "file"))
				if (terminate) {
					logger.error(msg)
				}
				logger.warning(msg)
				return(invisible(FALSE))
			}
		}
		return(invisible(TRUE))
	}
}
