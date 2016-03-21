########################################################################################################################
## utilities.R
## created: 2012-05-10
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Collection of (mostly internal) helper constants and functions.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.beta2mval
#'
#' Transforms beta values to M values, adjusting for +infinity and -infinity.
#' 
#' @param betas   \code{numeric} vector or matrix of beta values to be transformed.
#' @param epsilon Single \code{numeric} in the range [0, 0.5], giving the threshold of beta values to use when
#'                adjusting for potential M values close to +infinity or -infinity. Setting this parameter to 0 (zero)
#'                disables stabilization; in which case M values of -infinity or +infinity could be returned.
#' @return The calculated and adjusted M values.
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' mvals <- rnb.beta2mval(meth(rnb.set.example))
#' summary(mvals)
#' }
rnb.beta2mval <- function(betas, epsilon = 0.00001) {
	if (!is.numeric(betas)) {
		stop("invalid value for betas")
	}
	if (!(is.numeric(epsilon) && length(epsilon) == 1 && (!is.na(epsilon)))) {
		stop("invalid value for epsilon")
	}
	if (epsilon < 0 || epsilon > 0.5) {
		stop("invalid value for epsilon; expected 0 <= epsilon <= 0.5")
	}
	betas[betas < epsilon] <- epsilon
	betas[betas > (1 - epsilon)] <- 1 - epsilon
	return(log2(betas / (1 - betas)))
}

########################################################################################################################

#' rnb.mval2beta
#'
#' Transforms M values to beta values.
#' 
#' @param mvals \code{numeric} vector or matrix of M values to be transformed.
#' 
#' @return The calculated beta values.
#'
#' @author Pavlo Lutsik
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' mvals <- rnb.beta2mval(meth(rnb.set.example))
#' bvals <- rnb.mval2beta(mvals)
#' all((bvals-meth(rnb.set.example))<1e-10)
#' }
rnb.mval2beta <- function(mvals){
	if (!is.numeric(mvals)) {
		stop("invalid value for mvals")
	}
	intmf<-2^(mvals)
	return(intmf/(intmf+1))
}

########################################################################################################################

#' rnb.require
#'
#' Loads the specified required package.
#'
#' @param pkg Package name to be loaded.
#'
#' @author Yassen Assenov
#' @noRd
rnb.require <- function(pkg) {
	if (!suppressPackageStartupMessages(do.call(require, list(package = pkg)))) {
		rnb.error(paste("Missing required package", pkg, "or its dependency"))
	}
}

########################################################################################################################

## rnb.cleanMem
##
## perform memory cleanup, i.e. invoke R's garbage collector
## does nothing, if \code{rnb.getOption("enforce.memory.management")==FALSE}
## 
## @param iter.gc number of times to invoke the garbage
## @return invisible \code{TRUE}
##
## @author Fabian Mueller
rnb.cleanMem <- function(iter.gc=10) {
	if (rnb.getOption("enforce.memory.management")){
		#it might be helpful to invoke gc() multiple times
		#see http://stackoverflow.com/questions/1467201/forcing-garbage-collection-to-run-in-r-with-the-gc-command
	#	logger.start("cleaning memory")
		for (i in 1:iter.gc){gc()}
	#	logger.completed()
	}
	invisible(TRUE)
}

########################################################################################################################

#' rnb.call.destructor
#'
#' calls the destructor of an RnBSet, RnBeadSet or RnBeadRawSet object
#' conditionally on whether the \code{enforce.destroy.disk.dumps} option is enabled.
#' 
#' @param object object to be destroyed
#' @param ...	 further arguments to the method \code{\link[=destroy,RnBSet-method]{destroy}}
#' @return invisible \code{TRUE}
#' 
#' @author Fabian Mueller
rnb.call.destructor <- function(object,...) {
	if (rnb.getOption("enforce.destroy.disk.dumps")){
		destroy(object,...)
	}
	invisible(TRUE)
}

########################################################################################################################

#' rnb.sample.groups
#'
#' Identifies sample subgroups defined in the given annotation information.
#'
#' @param annotations     Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}},
#'                        or its sample annotations in the form of a \code{data.frame}. If this parameter is
#'                        a dataset, the annotation information is extracted using the method
#'                        \code{\link[=pheno,RnBSet-method]{pheno}}.
#' @param columns         Optional; predefined column names (in the form of a \code{character} vector) or indices (an
#'                        \code{integer} vector) to consider. All other columns in the annotation table will be ignored.
#' @param columns.pairs	  Optional; a NAMED vector containing for each column name for which paired comparisons should
#'                        be performed (say columnA) the name or index of another column (say columnB) in which same
#'                        values indicate the same pairing. columnA should be the name of the value columnB in this
#'                        vector.
#' @param min.group.size  Minimum number of samples in each subgroup. This must be a positive \code{integer}.
#' @param max.group.count Maxumum number of subgroups defined by a trait. This must be an \code{integer} greater than 1.
#' @return List of traits that define subgroups in the dataset. For each trait, the defined subgroups are represented by
#'         a list of \code{integer} vectors storing the corresponding sample indices.
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' str(rnb.sample.groups(rnb.set.example))
#' }
#'
#' @author Yassen Assenov
#' @export
rnb.sample.groups <- function(annotations, columns = NULL, columns.pairs = NULL, min.group.size = rnb.getOption("min.group.size"),
		max.group.count = rnb.getOption("max.group.count")) {
	if (!is.null(annotations)) {
		if (inherits(annotations, "RnBSet")) {
			annotations <- pheno(annotations)
		} else if (!is.data.frame(annotations)) {
			stop("invalid value for annotations; expected object of type RnBSet or a data.frame")
		}
	}
	if (!(is.null(columns) || is.character(columns) || is.integer(columns))) {
		stop("invalid value for columns; expected character or integer")
	}
	if (!(is.null(columns.pairs) || is.character(columns.pairs) || is.integer(columns.pairs))) {
		stop("invalid value for columns.pairs; expected character or integer")
	}
	if (!is.null(columns.pairs)) {
		if (is.null(names(columns.pairs))) {
			stop("invalid value for columns.pairs; expected NAMED vector")
		}
		if (any(!(names(columns.pairs) %in% colnames(annotations)))) {
			stop("invalid value for columns.pairs; not all names are present in the annotation")
		}
	}
	if (is.numeric(min.group.size) && all(min.group.size == as.integer(min.group.size), na.rm = TRUE)) {
		min.group.size <- as.integer(min.group.size)
	}
	if (!(is.integer(min.group.size) && length(min.group.size) == 1 && (!is.na(min.group.size)))) {
		stop("invalid value for min.group.size; expected single integer")
	}
	min.group.size <- min.group.size[1]
	if (min.group.size < 1) {
		stop("invalid value for min.group.size; expected positive integer")
	}
	if (is.numeric(max.group.count) && all(max.group.count == as.integer(max.group.count), na.rm = TRUE)) {
		max.group.count <- as.integer(max.group.count)
	}
	#set the maximum number of groups to be N-1 for NULL values
	if (is.null(max.group.count)) {
		max.group.count <- nrow(annotations) - 1
		max.group.count <- as.integer(max(c(2,max.group.count)))
	}
	if (!(is.integer(max.group.count) && length(max.group.count) == 1 && (!is.na(max.group.count)))) {
		stop("invalid value for max.group.count; expected single integer")
	}
	max.group.count <- max.group.count[1]
	if (max.group.count < 2) {
		stop("invalid value for max.group.count; expected at least 2")
	}
	
	if (is.null(annotations)) {
		result <- list()
		names(result) <- character(0)
	} else {
		annotations.all <- annotations
		if (!is.null(columns)) {
			if (is.character(columns)) {
				columns <- intersect(columns, colnames(annotations))
			} else { # is.integer(columns)
				columns <- intersect(columns, 1:ncol(annotations))
			}
			pp <- annotations
			annotations <- annotations[, columns]
			if (length(columns) == 1){
				annotations <- data.frame(annotations)
				rownames(annotations) <- rownames(pp)
				colnames(annotations) <- columns
			}
		}
		result <- lapply(annotations, function(trait) {
					res <- tapply(1:nrow(annotations), trait, identity)
					if (length(res) < 2) return(NULL)
					res <- res[sapply(res, length) != 0] # ignore levels that are missing in the dataset
					res <- res[sapply(res, function(x){!any(is.na(x))})] # ignore levels that are missing in the dataset
					if (length(res) < 2 || max.group.count < length(res) || any(sapply(res, length) < min.group.size)) {
						return(NULL)
					}
					res
				}
		)
		result <- result[!sapply(result, is.null)]
		has.pairing <- rep(FALSE,length(result))
		#find pairing information
		if (!is.null(columns.pairs)){
			has.pairing <- names(result) %in% names(columns.pairs)
			pair.list <- lapply(names(columns.pairs),FUN=function(cc){
				annotations.all[,columns.pairs[cc]]
			})
			names(pair.list) <- names(columns.pairs)
			for (i in 1:length(result)){
				cc <- names(result)[i]
				if (has.pairing[i]){
					ggs <- result[[i]]
					ggs.pair.ids <- lapply(ggs,FUN=function(iis){
						pair.list[[cc]][iis]
					}) #pairing identifiers for the individual groups
					all.pair.ids <- factor(unlist(ggs.pair.ids))
					pair.id.order <- levels(all.pair.ids)
					is.valid.pairing <- TRUE
					ggs.sorted <- list()
					for (j in 1:length(ggs)) {
						ggn <- names(ggs)[j]
						if (any(duplicated(ggs.pair.ids[[ggn]])) || !all(pair.id.order %in% ggs.pair.ids[[ggn]])){
							is.valid.pairing <- FALSE
							break
						}
						ggs.cur <- ggs[[ggn]]
						names(ggs.cur) <- ggs.pair.ids[[ggn]]
						ggs.sorted[[j]] <- ggs.cur[pair.id.order]
					}
					if (is.valid.pairing){
						names(ggs.sorted) <- names(ggs)
						dd <- data.frame(ggs.sorted)
						rownames(dd) <- pair.id.order
						result[[i]] <- dd
					} else {
						logger.warning(paste("Invalid pairing information for column",cc,". --> Treating as unpaired."))
						has.pairing[i] <- FALSE
					}
				}
			}
		}
		attr(result,"paired") <- has.pairing
	}
	return(result)
}

########################################################################################################################

#' rnb.sample.replicates
#'
#' Identifies sample replicates defined in the given sample annotation table.
#'
#' @param rnb.set          Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param replicate.id.col Trait (column name in the sample annotation table) that indicates sample replicates.
#' 						   Replicates should have the same value for this trait, while samples without replicates are
#'                         expected to have unique values or missing values.
#' @return List of length of the number of replicates in the dataset. Each element is an \code{integer} vector storing
#'         the corresponding sample indices.
#'
#' @author Fabian Mueller
#' @export
rnb.sample.replicates <- function(rnb.set, replicate.id.col) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	result <- list()
	names(result) <- character(0)
	if (is.null(replicate.id.col)) {
		return(result)
	}
	if (!(is.character(replicate.id.col) && length(replicate.id.col) == 1 && (!is.na(replicate.id.col)))) {
		stop("invalid value for replicate.id.col")
	}
	if (!(replicate.id.col %in% colnames(pheno(rnb.set)))) {
		warning("specified replicate.id.column not found in rnb.set")
		return(result)
	}

	ids <- pheno(rnb.set)[, replicate.id.col]
	return(tapply(1:length(ids), ids, identity))
}

########################################################################################################################

#' rnb.show.report
#' 
#' Opens the given HMTL report file in the browser.
#' 
#' @param report \code{\linkS4class{Report}} object to open.
#' @author Pavlo Lutsik
#' @export
rnb.show.report <- function(report) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	
	tryCatch(browseURL(sprintf("file://%s", report@fname)),
		error = function(err) {
			rnb.warning(c("Could not display the report due to the following error message:", err$message))
			invisible(NULL)
		}
	)
}

########################################################################################################################

#' rnb.load.sitelist
#'
#' Loads a list of probe or site identifiers. This function is used in the preprocessing module for loading a whitelist
#' and/or a blacklist of identifiers.
#'
#' @param fname   File listing the identifiers, one per line.
#' @param verbose Flag indicating if messages are to be printed. If the values is \code{TRUE} and a logger is
#'                initialized, this function adds a message to the log.
#' @return The loaded list of identifiers, or \code{NULL} if \code{fname} could not be open.
#'
#' @seealso \code{\link{logger.start}} for initializing a logger
#' @author Yassen Assenov
#' @export
rnb.load.sitelist <- function(fname, verbose = FALSE) {
	result <- tryCatch(scan(fname, "", sep = "\n", quiet = TRUE), error = function(er) { NULL })
	if (verbose) {
		if (is.null(result)) {
			msg <- paste("Could not load site list from", fname)
			if (logger.isinitialized()) {
				logger.error(msg, terminate = FALSE)
			} else {
				message(msg)
			}
		} else {
			msg <- paste("Loaded", length(result), "sites from", fname)
			if (logger.isinitialized()) {
				rnb.status(msg)
			} else {
				message(msg)
			}
		}
	}
	result
}

########################################################################################################################

## rnb.process.sitelist
##
## Loads and processes a list of probe or site identifiers.
##
## @param fname      File containing the list of identifiers to be processed.
## @param anno.table Probe or site annotation table.
## @return \code{integer} vector containing the indices in the annotation table that are targeted by the loaded list.
##
## @author Yassen Assenov
rnb.process.sitelist <- function(fname, anno.table) {
	if (length(fname) == 0) {
		return(integer())
	}
	id.list <- rnb.load.sitelist(fname, TRUE)
	if (is.null(id.list)) {
		return(integer())
	}
	which(rownames(anno.table) %in% id.list)
}

########################################################################################################################

#' rnb.write.table
#'
#' Writes a table to a file. Different formats and compression options are available.
#' 
#' @param tt     Table to be written to file, usually in the form of a \code{matrix} or \code{data.frame}.
#' @param fname  Target file name. If this file already exists, it will be overwritten.
#' @param fpath  Target file path. If "" (default value), \code{fname} is assumed to contain the absolute path.
#' @param format Target format; one of \code{"csv"}, \code{"tab"} or \code{"txt"}, denoting comma-separated,
#'               tab-separated and default text format, respectively. The last format allows for a user-specified
#'               delimiter through an additional parameter \code{sep}. See the documentation of
#'               \code{\link{write.table}} for more details.
#' @param gz     Flag indicating whether the file should be zipped in \code{gz} format.
#' @param ...    Any additional arguments to be passed on to \code{write.table} or \code{utils::write.csv}.
#' @return The (possibly updated) target file name, invisibly. If \code{gz} is \code{TRUE}, the string \code{".gz"} will
#'         be appended to \code{fname}.
#'
#' @seealso \code{\link{write.table}}
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' data(mtcars)
#' rnb.write.table(mtcars,tempfile(pattern="cars",fileext=".csv"))
#' }
rnb.write.table <- function(tt, fname, fpath = "", format = "csv", gz = FALSE, ...) {
	if (!(is.character(fname) && length(fname) == 1 && (!is.na(fname)))) {
		stop("invalid value for fname; expected single character")
	}
	if (!(is.character(fname) && length(fname) == 1)) {
		stop("invalid value for fpath; expected single character")
	} else if (is.na(fpath)) {
		fpath <- ""
	}
	if (!(is.character(format) && length(format) == 1)) {
		stop("invalid value for format")
	}
	if (!parameter.is.flag(gz)) {
		stop("invalid value for gz; expected TRUE or FALSE")
	}

	if (gz && grepl("\\.gz$", fname) == FALSE) {
		fname <- paste(fname, "gz", sep = ".")
	}
	fname.full <- ifelse(fpath != "", file.path(fpath, fname), fname)
	if (gz) {
		fname.full <- gzfile(fname.full, open = "w")
	}
	if (format == "csv") {
		utils::write.csv(tt, file = fname.full, ...)
	} else if (format == "tab") {
		write.table(tt, file = fname.full, sep = "\t", ...)
	} else {
		write.table(tt, file = fname.full, ...)
	}
	if (gz) {
		close(fname.full)
	}
	invisible(fname)
}

########################################################################################################################

#' rnb.region.types.for.analysis
#'
#' Identifies the region types that are summarized by the given dataset and pointed to for analysis.
#' 
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @return List of all region types to be analyzed in the current dataset in the form of a \code{character} vector.
#'
#' @details
#' This function intersects the value of the analysis option \code{"region.types"} with the region types that are
#' summarized in the provided dataset. In case the option's value is \code{NULL}, this function returns all summarized
#' region types in \code{rnb.set}.
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' "promoters" %in% rnb.region.types.for.analysis(rnb.set.example)
#' }
#' @seealso \code{\link{rnb.getOption}} for checking the value of the \code{"region.types"} option;
#'   \code{\link[=summarized.regions,RnBSet-method]{summarized.regions}} for obtaining the region types summarized in a
#'   dataset
#' @author Yassen Assenov
#' @export
rnb.region.types.for.analysis <- function(rnb.set) {
	region.types <- rnb.getOption("region.types")
	if (is.character(rnb.set) && length(rnb.set) == 1 && (rnb.set %in% rnb.get.assemblies())) {
		## Handle the special (undocumented) case when an assembly is specified instead of a dataset
		if (is.null(region.types)) {
			region.types <- rnb.region.types(rnb.set)
		} else {
			region.types <- intersect(region.types, rnb.region.types(rnb.set))
		}
		return(region.types)
	}
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	region.types.dataset <- summarized.regions(rnb.set)
	if (is.null(region.types)) {
		region.types <- region.types.dataset
	} else {
		r.lost <- setdiff(region.types, region.types.dataset)
		if (logger.isinitialized() && length(r.lost) != 0) {
			txt <- c("The following", length(r.lost), "region types will not be included in the analysis:")
			logger.warning(c(txt, paste(r.lost, collapse = ", ")))
		}
		region.types <- intersect(region.types, region.types.dataset)
	}
	region.types
}

########################################################################################################################

## rnb.regions2sites
##
## Creates a mapping from regions to sites.
##
## @param regions GRangesList object storing genomic regions per chromosome. The number of elements and their names must
##                match the names of \code{sites}.
## @param sites   GRangesList object storing genomic loci (sites) per chromosome.
## @return Mapping from the specified regions to the sites contained in them, in the form of a list of mappings per
##         chromosome. Each entry in this list is an object of type \code{IRanges}, storing indices of \code{sites} on
##         the respective chromosome that are contained in the corresponding region. Regions that do not contain sites
##         are left out of the result.
##
## @author Fabian Mueller
rnb.regions2sites <- function(regions, sites) {
	if (length(setdiff(names(regions), names(sites))) != 0) {
		stop("sites is not defined for all chromosomes in regions")
	}
	ovl <- list()
	for (chrom in names(regions)) {
		ov <- GenomicRanges::findOverlaps(regions[[chrom]], sites[[chrom]], ignore.strand = TRUE)
		if (length(ov) != 0) {
			indexMap <- tapply(subjectHits(ov), queryHits(ov), range)
			reg.indices <- names(indexMap)
			indexMap <- matrix(unlist(indexMap, use.names = FALSE), nrow = 2)
			indexMap <- IRanges(start = indexMap[1, ], end = indexMap[2, ], names = reg.indices)
			ovl[[chrom]] <- indexMap
		}
	}
	return(ovl)
}

########################################################################################################################

## rnb.get.sentrix.data
##
## Extracts information on Sentrix ID and Sentrix position from the given sample information table.
##
## @param dframe Sample annotation in the form of a \code{data.frame}; or a dataset as an object of type inheriting
##               \code{\linkS4class{RnBeadSet}}.
## @return Table (\code{data.frame}) with two columns - \code{"Slide"} and \code{"Position"}, listing the slide number
##         location on the slide for each sample. If the provided annotation contains no valid sentrix ID and position
##         information, the returned value is \code{NULL}.
##
## @author Yassen Assenov
rnb.get.sentrix.data <- function(dframe) {
	if (inherits(dframe, "RnBeadSet")) {
		dframe <- pheno(dframe)
	}
	if (all(c("Sentrix_ID", "Sentrix_Position") %in% colnames(dframe))) {
		sentrix.ids <- dframe[, "Sentrix_ID"]
		sentrix.pos <- as.character(dframe[, "Sentrix_Position"])
	} else if (all(c("Sentrix ID", "Sentrix Position") %in% colnames(dframe))) {
		sentrix.ids <- dframe[, "Sentrix ID"]
		sentrix.pos <- as.character(dframe[, "Sentrix Position"])
	} else if ("barcode" %in% colnames(dframe)) {
		barcodes <- strsplit(as.character(dframe[, "barcode"]), "_", fixed = TRUE)
		if (all(sapply(barcodes, length) == 2L)) {
			sentrix.ids <- sapply(barcodes, '[', 1L)
			sentrix.pos <- sapply(barcodes, '[', 2L)
		}
	}
	if (exists("sentrix.ids", inherits = FALSE)) {
		positions <- paste0("R0", 1:6, "C0", rep(1:2, each = 6))
		if ((!(any(c(is.na(sentrix.ids), is.na(sentrix.pos)))) || any(sentrix.ids == "")) &&
			all(sentrix.pos %in% positions)) {
			return(data.frame("Slide" = sentrix.ids, "Position" = factor(sentrix.pos, levels = positions),
					stringsAsFactors = TRUE, check.names = FALSE))
		}
	}
	return(NULL)
}

########################################################################################################################

#' rnb.step.analyze.targets
#'
#' Gets the list of all site and region types that are both requested and supported.
#'
#' @param rnb.set       Dataset of interest.
#' @param report        Report to contain the section once the analysis is complete.
#' @param section.title Title of the report section in human-readable form.
#' @param analyze.sites Flag indicating if the sites are requested to be analyzed.
#' @return List of supported site and region types in the form of a \code{character} vector. If there are region types
#'         that are requested but not summarized in the dataset, this list contains a sentence as an attribute named
#'         \code{"warning"}.
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.analyze.targets <- function(rnb.set, report, section.title, analyze.sites = rnb.getOption("analyze.sites")) {
	if (analyze.sites) {
		targets <- "sites"
	} else {
		targets <- character()
	}
	regions.requested <- rnb.region.types.for.analysis(rnb.set)
	regions.supported <- intersect(regions.requested, names(rnb.set@regions))
	regions.unsupported <- setdiff(regions.requested, regions.supported)
	targets <- c(targets, regions.supported)

	if (length(targets) == 0) {
		if (length(regions.unsupported) != 0) {
			txt <- c("This procedure was skipped because none of the requested region types are summarized in the ",
				"dataset.")
		} else {
			txt <- "This procedure was skipped because neither site nor region types are specified for analysis."
		}
		rnb.add.section(report, section.title, txt)
	} else if (length(regions.unsupported) != 0) {
		txt <- paste("<b>", regions.unsupported, "</b>", sep = "", collapse = ", ")
		txt <- paste("The following region types are not analyzed because they are not summarized in the dataset:", txt)
		attr(targets, "warning") <- txt
	}
	return(targets)
}

########################################################################################################################

## capitalize
##
## Makes the first letter of each of the given strings capital.
##
## @param x          \code{character} vector containing strings to be capitalized.
## @param every.word Flag indicating if the first letter of every word must be capitalized. Note that this method
##                   splits to words using the space character only, therefore, other whitespaces (tabulation or newline
##                   characters) are not recognized as word separators. If this parameter is \code{FALSE} (default),
##                   only the first letter of every string is capitalized. 
## @return \code{character} vector of the same length as \code{x}, in which every element is capitalized.
##
## @details
## If \code{every.word} is \code{TRUE}, this method removes a trailing space (if such exists) from every element of
## \code{x}. This is a side effect of splitting into words.
##
## @author Yassen Assenov
capitalize <- function(x, every.word = FALSE) {
	if (every.word) {
		result <- sapply(strsplit(x, split = " "), capitalize, USE.NAMES = !is.null(names(x)))
		result <- sapply(result, paste, collapse = " ")
	} else {
		result <- paste(toupper(substr(x, 1, 1)), substring(x, 2), sep = "")
		result[is.na(x)] <- NA
	}
	return(result)
}

########################################################################################################################

## symmetric.melt
##
## Constructs a \code{data.frame} for plotting the lower triangular values in the given symmetric matrix.
## 
## @param tbl.symmetric    Symmetric \code{matrix} of values with at least 2 rows.
## @param include.diagonal Flag indicating if the diagonal of \code{tbl.symmetric} should also be displayed.
## @return Table, in the form of a \code{data.frame} with three columns, to be passed to \pkg{ggplot}. Missing values
##         (\code{NA}) in the original matrix are omitted from the resulting table.
## @author Yassen Assenov
symmetric.melt <- function(tbl.symmetric, include.diagonal = TRUE) {
	tbl <- tbl.symmetric
	tbl[upper.tri(tbl)] <- NA
	if (!include.diagonal) {
		tbl <- tbl[-1, -ncol(tbl)]
		if (nrow(tbl.symmetric) == 2) {
			tbl <- matrix(tbl, nrow = 1, ncol = 1,
				dimnames = list(rownames(tbl.symmetric)[2], colnames(tbl.symmetric)[1]))
		}
	}
	tbl.melt <- melt(tbl, varnames = c("x", "y"))
	tbl.melt[[1]] <- factor(as.character(tbl.melt[[1]]), levels = rev(rownames(tbl)))
	tbl.melt[[2]] <- factor(as.character(tbl.melt[[2]]), levels = colnames(tbl))
	tbl.melt[!is.na(tbl.melt[[3]]), ]
}

########################################################################################################################

## get.i.vector
##
## Transforms the given sublisting of items to indices (if necessary) and validates its scope.
##
## @param i.list Sublisting of items in the form of a \code{logical}, \code{integer} or \code{character} vector.
## @param items  Full list of items to be indexed.
## @return       Sorted \code{integer} vector representing indices in \code{items}.
## @author Yassen Assenov
get.i.vector <- function(i.list, items) {
	if (is.double(i.list) && all(i.list == as.integer(i.list), na.rm = TRUE)) {
		i.list <- as.integer(i.list)
	}
	if (is.logical(i.list)) {
		if (!is.null(items) && length(i.list) != length(items)) {
			stop("specified logical vector is of unexpected length")
		}
		inds <- which(i.list)
	} else if (is.character(i.list)) {
		inds <- which(items %in% i.list)
		if (length(i.list) != length(inds)) {
			stop("unknown names specified")
		}
	} else if (is.integer(i.list)) {
		inds <- sort(i.list)
		if (length(inds) != 0 && (min(inds) < 1 || (!is.null(items) && length(items) < max(inds)))) {
			stop("invalid indices specified")
		}
	} else {
		stop("invalid item list specified")
	}
	inds
}

########################################################################################################################

## parameter.is.flag
##
## Checks if the provided parameter value is a flag.
## 
## @param value Value to be tested.
## @return \code{TRUE} if \code{values} is \code{TRUE} or \code{FALSE}, \code{FALSE} otherwise.
## @author Yassen Assenov
parameter.is.flag <- function(value) {
	is.logical(value) && length(value) == 1 && (!is.na(value))
}

########################################################################################################################

## Validates the given vector or list contains a single element. This function is used in validating function or method
## arguments.
##
## @param x          Value vector or list to validate.
## @param param.name Name of parameter or slot that is validated. This is used in the generation of failing message.
## @return Short message that encodes the result of the validation in the form of a \code{character}. It is either the
##         string \code{ok}, or a short phrase describing the divergence from the "single value assumption".
## @author Yassen Assenov
validate.single <- function(x, param.name = "x") {
	if (is.null(x) || length(x) == 0) {
		result <- paste("missing value for", param.name)
	} else if (length(x) > 1) {
		result <- paste("multiple values for", param.name)
	} else if (is.na(x)) {
		result <- paste("missing value for", param.name)
	} else {
		result <- "ok"
	}
	return(result)
}

########################################################################################################################

## If there is a logger initialized, validates that the given directory exists.
##
## @param dname Name of directory to be validated.
## @author Yassen Assenov
validate.dir <- function(dname) {
	if (logger.isinitialized()) {
		logger.validate.file(dname, is.file = FALSE)
	}
}

########################################################################################################################

## write.line
##
## Writes a line of text to the specified text file. This function is used in the generation of HTML reports.
##
## @param txt    Character vector storing the text to be written. The elements of this vector are concatenated without
##               a separator.
## @param fname  Name of the file to write the text to.
## @param indent Indentation of the text, given as number of \code{TAB} characters.
## @param append Flag indicating if the line is to be appended to the text file. If this is \code{FALSE}, the file's
##               contents are overwritten.
## @author Yassen Assenov
write.line <- function(txt, fname, indent = 0, append = TRUE) {
	strprefix <- paste(rep("\t", times = indent), collapse = "")
	cat(strprefix, paste0(txt, collapse = ""), "\n", file = fname, sep = "", append = append)
}

########################################################################################################################

## write.table.options
##
## Generates HTML code for a table of module option values in the specified report.
##
## @param optionlist List of \pkg{RnBeads} module options. The attribute \code{"enabled"}, if present, is used to mark 
##                   which of these options were actually used.
## @param report     Report to which the table is to be added.
## @param indent     Default indentation, in number of tabulation characters, to apply to the \code{<table>} tag.
## @author Yassen Assenov
write.table.options <- function(optionlist, report, indent = 0L) {
	if (is.null(attr(optionlist, "enabled"))) {
		attr(optionlist, "enabled") <- rep.int(TRUE, length(optionlist))
	}

	wline <- function(txt, indentation = indent) {
		write.line(txt, report@fname, indent = indentation)
	}
	wline("<table class=\"tindex\">")
	wline("<thead>")
	wline("<tr>")
	wline("<th>Option</th>", indent + 1L)
	wline("<th>Value</th>", indent + 1L)
	wline("</tr>")
	wline("</thead>")
	wline("<tbody>")
	for (i in 1:length(optionlist)) {
		wline(c("<tr", ifelse(i %% 2 == 0, " class=\"darker\"", ""), ">"))
		tdopen <- c("<td", ifelse(attr(optionlist, "enabled")[i], "", " class=\"disabled\""), ">")
		wline(c(tdopen, names(optionlist)[i], "</td>"), indent + 1L)
		optionvalue <- optionlist[[i]]
		if (is.null(optionvalue)) {
			optionvalue <- "<span class=\"missing\">default</span>"
		} else if (length(optionvalue) != 1) {
			optionvalue <- paste(optionvalue, collapse = ", ")
		} else if (is.na(optionvalue)) {
			optionvalue <- "<span class=\"disabled\">n/a</span>"
		} else if (is.logical(optionvalue)) {
			optionvalue <- ifelse(optionvalue, "yes", "no")
		}
		wline(c(tdopen, optionvalue, "</td>"), indent + 1L)
		wline("</tr>")
	}
	wline("</tbody>")
	wline("</table>")
}

########################################################################################################################

#' rnb.add.exported.tables
#'
#' Adds a paragraph informing the readers that data are exported as CSV files, following by a table of links.
#'
#' @param report Report to be written to.
#' @param tbl    Table containing links to the exported CSV files.
#'
#' @author Yassen Assenov
#' @noRd
rnb.add.exported.tables <- function(report, tbl) {
	txt <- c("The values presented in the figure above are avaialable in CSV (comma-separated value) files ",
		"accompanying this report.")
	rnb.add.paragraph(report, txt)
	rnb.add.table(report, tbl, row.names = FALSE)
}

########################################################################################################################

#' rnb.get.reliability.matrix
#' 
#' Gets a matrix of reliability indications for every measurement in the given dataset.
#' 
#' @param rnb.set   Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param row.names	Flag indicating of row names are to be generated in the result.
#' @return \code{logical} matrix in which every row corresponds to a CpG site or probe and every column - to a patient.
#'         If the dataset does not contain coverage or detection p-value information, the returned value is \code{NULL}.
#' 
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.options(identifiers.column = "Sample_ID")
#' str(rnb.get.reliability.matrix(rnb.set.example))
#' }
#' @author Yassen Assenov
#' @export
rnb.get.reliability.matrix <- function(rnb.set, row.names = FALSE) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (inherits(rnb.set, "RnBeadSet")) {
		result <- dpval(rnb.set, row.names = row.names)
		if (!is.null(result)) {
			result <- (result < rnb.getOption("filtering.greedycut.pvalue.threshold"))
		}
	} else { # inherits(rnb.set, "RnBiseqSet")
		result <- covg(rnb.set, row.names = row.names)
		if (!is.null(result)) {
			result <- (result >= rnb.getOption("filtering.coverage.threshold"))
		}
	}
	result
}

########################################################################################################################

## methylumi.intensities.by.color
## 
## Rearranges information from "methylated" and "unmethylated" slots of a MethyLumiSet object by color channer.
##
## @param mset MethyLumiSet object
## @param address.rownames If \code{TRUE} the returned intensity matrices have Infinium probe address as row names   
## 
## @author Pavlo Lutsik
methylumi.intensities.by.color<-function(mset,address.rownames=TRUE){

	if(!require("IlluminaHumanMethylation450kmanifest"))
		stop("IlluminaHumanMethylation450kmanifest should be installed")

	pinfos <- rnb.annotation2data.frame(rnb.get.annotation("probes450"), add.names=TRUE)[featureNames(mset), ]
	dII.probes <- featureNames(mset)[pinfos[,"Design"] == "II"]
		
	dII.probes <- dII.probes[!grepl("rs", dII.probes)]
			
	if(address.rownames) tII<-IlluminaHumanMethylation450kmanifest@data$TypeII 
	
	if(address.rownames) tII<-tII[match(dII.probes, tII$Name),]
	
	dII.grn<-methylated(mset)[dII.probes,]
	if(address.rownames) rownames(dII.grn)<-tII$Address
	
	dII.red<-unmethylated(mset)[dII.probes,]
	if(address.rownames) rownames(dII.red)<-tII$Address

	dI.red.probes <- featureNames(mset)[pinfos[, "Color"] == "Red"]
	dI.red.probes <- dI.red.probes[!grepl("rs", dI.red.probes)]
	dI.green.probes <- featureNames(mset)[pinfos[, "Color"] == "Grn"]
	dI.green.probes <- dI.green.probes[!grepl("rs", dI.green.probes)]
	if(address.rownames) tI<-IlluminaHumanMethylation450kmanifest@data$TypeI
		
	if(address.rownames) tI.red<-tI[tI$Color=="Red",]
	if(address.rownames) tI.red<-tI.red[match(dI.red.probes, tI.red$Name),]
	
	if(address.rownames) tI.grn<-tI[tI$Color=="Grn",]
	if(address.rownames) tI.grn<-tI.grn[match(dI.green.probes, tI.grn$Name),]
	
	dI.red.meth<-methylated(mset)[dI.red.probes,]
	
	if(address.rownames) rownames(dI.red.meth)<-tI.red[,"AddressB"]
	
	dI.red.umeth<-unmethylated(mset)[dI.red.probes,]
	if(address.rownames) rownames(dI.red.umeth)<-tI.red[,"AddressA"]
	
	
	dI.red.meth.oob<-assayDataElement(mset, "methylated.OOB")[dI.red.probes,]
	if(address.rownames) rownames(dI.red.meth.oob)<-tI.red[,"AddressB"]
	
	dI.red.umeth.oob<-assayDataElement(mset, "unmethylated.OOB")(mset)[dI.red.probes,]
	if(address.rownames) rownames(dI.red.umeth.oob)<-tI.red[,"AddressA"]
	
	
	dI.grn.meth<-methylated(mset)[dI.green.probes,]
	if(address.rownames) rownames(dI.grn.meth)<-tI.grn[,"AddressB"]
	
	dI.grn.umeth<-unmethylated(mset)[dI.green.probes,]
	if(address.rownames) rownames(dI.grn.umeth)<-tI.grn[,"AddressA"]
	
	dI.grn.meth.oob<-assayDataElement(mset, "methylated.OOB")(mset)[dI.green.probes,]
	if(address.rownames) rownames(dI.grn.meth.oob)<-tI.grn[,"AddressB"]
	
	dI.grn.umeth.oob<-assayDataElement(mset, "unmethylated.OOB")[dI.green.probes,]
	if(address.rownames) rownames(dI.grn.umeth.oob)<-tI.grn[,"AddressA"]

	
	intensities.by.channel <- list(
			Cy3=rbind(dII.grn, dI.grn.meth,dI.grn.umeth, dI.red.meth.oob, dI.red.umeth.oob),
			Cy5=rbind(dII.red, dI.red.meth, dI.red.umeth, dI.grn.meth.oob, dI.grn.umeth.oob))

#	intensities.by.channel$Cy3<-rbind(intensities.by.channel$Cy3, 
#			matrix(0, nrow=length(setdiff(rownames(intensities.by.channel$Cy5),rownames(intensities.by.channel$Cy3))),
#					ncol=ncol(intensities.by.channel$Cy3),
#					dimnames=list(rownames=setdiff(rownames(intensities.by.channel$Cy5),rownames(intensities.by.channel$Cy3)))))
#	intensities.by.channel$Cy5<-rbind(intensities.by.channel$Cy5, 
#			matrix(0, nrow=length(setdiff(rownames(intensities.by.channel$Cy3),rownames(intensities.by.channel$Cy5))), 
#					ncol=ncol(intensities.by.channel$Cy5), 
#					dimnames=list(rownames=setdiff(rownames(intensities.by.channel$Cy3),rownames(intensities.by.channel$Cy5)))))

	intensities.by.channel$Cy5<-intensities.by.channel$Cy5[rownames(intensities.by.channel$Cy3),]

	ncd<-rnb.get.annotation("controls450")
	ncd<-ncd[ncd[["Target"]] == "NEGATIVE", ]
	ncd$Target<-tolower(ncd$Target)
	controls.by.channel<-intensitiesByChannel(controlData(mset))
	rownames(controls.by.channel$Cy3)<-tolower(rownames(controls.by.channel$Cy3))
	rownames(controls.by.channel$Cy5)<-tolower(rownames(controls.by.channel$Cy5))

	ncd<-ncd[1:length(grep("negative", rownames(controls.by.channel$Cy3))),]
	controls.by.channel$Cy3<-controls.by.channel$Cy3[sub("\\.0", "" ,paste(ncd$Target, 
							ncd$Index-ifelse("negative" %in% rownames(controls.by.channel$Cy3),1,0), sep=".")),]
	controls.by.channel$Cy5<-controls.by.channel$Cy5[sub("\\.0", "" ,paste(ncd$Target, 
							ncd$Index-ifelse("negative" %in% rownames(controls.by.channel$Cy3),1,0), sep=".")),]

	rownames(controls.by.channel$Cy3)<-ncd$ID
	rownames(controls.by.channel$Cy5)<-ncd$ID

	intensities.by.channel$Cy3<-rbind(intensities.by.channel$Cy3, controls.by.channel$Cy3)
	intensities.by.channel$Cy5<-rbind(intensities.by.channel$Cy5, controls.by.channel$Cy5)

	return(intensities.by.channel)
}

########################################################################################################################

## abbreviate.names
## 
## Creates short versions, if necessary, of the supplied vector labels.
##
## @param labels ...
## @param max.length ...
## @return ...
## @author Pavlo Lutsik
abbreviate.names <- function(labels, max.length = 15L) {
	N <- (max.length - 3) %/% 2
	ifelse(nchar(labels) <= max.length, labels,
		paste0(substring(labels, 1, N), "...", substring(labels, nchar(labels) - N + 1)))
}

########################################################################################################################

#' rnb.logger.start
#'
#' Starts a new section in the log (if it is initialized).
#' @param txt \code{character} vector storing the section title.
#' @details If the logger is not initialized, calling this method has no effect.
#'
#' @author Yassen Assenov
#' @noRd
rnb.logger.start <- function(txt) {
	if (logger.isinitialized()) {
		logger.start(txt)
	}
}

########################################################################################################################

#' rnb.logger.completed
#'
#' Completes the last open section in the log (if it is initialized).
#' @details If the logger is not initialized, calling this method has no effect.
#'
#' @author Yassen Assenov
#' @noRd
rnb.logger.completed <- function() {
	if (logger.isinitialized()) {
		logger.completed()
	}
}

########################################################################################################################

#' rnb.info
#'
#' Prints an information message to the log (if it is initialized).
#' @param txt \code{character} vector storing the status message.
#' @details If the logger is not initialized, calling this function has no effect.
#'
#' @author Yassen Assenov
#' @noRd
rnb.info <- function(txt) {
	if (logger.isinitialized()) {
		logger.info(txt)
	} else {
		message(paste(txt, collapse = " "))
	}
}

########################################################################################################################

#' rnb.status
#'
#' Prints a status message to the log (if it is initialized).
#' @param txt \code{character} vector storing the status message.
#' @details If the logger is not initialized, calling this function has no effect.
#'
#' @author Yassen Assenov
#' @noRd
rnb.status <- function(txt) {
	if (logger.isinitialized()) {
		logger.status(txt)
	} else {
		message(paste(txt, collapse = " "))
	}
}

########################################################################################################################

#' rnb.warning
#'
#' Prints a warning message to the log (if it is initialized) or the console (otherwise).
#' @param txt \code{character} vector storing the warning message.
#'
#' @author Yassen Assenov
#' @noRd
rnb.warning <- function(txt) {
	if (logger.isinitialized()) {
		logger.warning(txt)
	} else {
		warning(paste(txt, collapse = " "), call.=FALSE, immediate.=TRUE)
	}
}

########################################################################################################################

#' rnb.error
#'
#' Prints an error message to the log (if it is initialized) or executes an error action (otherwise).
#' @param txt \code{character} vector storing the error message.
#'
#' @author Yassen Assenov
#' @noRd
rnb.error <- function(txt) {
	if (logger.isinitialized()) {
		logger.error(txt)
	} else {
		stop(paste(txt, collapse = " "), call.=FALSE)
	}
}

########################################################################################################################

## convert.to.ff.tmp.file
##
## converts a table to an ff object
##
## @author Fabian Mueller
convert.to.ff.matrix.tmp <- function(x,...) {
	ff(x,dim=dim(x),dimnames=dimnames(x),finalizer="delete",...)
}

## additions with respect to a pached ff
convert.to.ff.matrix.tmp2 <- function(x,...) {
	if(is.integer(x)){
		vm<-'integer'
	}else{
		vm<-'double'
	}
	
	ffobj<-ff(vmode=vm, dim=dim(x),dimnames=dimnames(x),finalizer="delete",...)
	for(i in 1:ncol(x)){
		ff[,i]<-x[,i]
	}
	ffobj
}

## creating an empty ff matrix
create.empty.ff.matrix.tmp <- function(vm,dim, dimnames=NULL,...) {
	ff(vmode=vm,dim=dim,dimnames=dimnames,finalizer="delete",...)
}


get.matrix.from.ff<-function(x){

	nm<-matrix(NA,ncol=ncol(x), nrow=nrow(x))
	for(i in 1:ncol(x)){
		nm[,i]<-x[,i]
	}
	nm
}

########################################################################################################################

## prepare.idat.dir
##
## Creates a temporary directory with symbolic links to all \code{idat} files located in the given directory or its
## subdirectories.
##
## @param base.dir Base directory to scan for \code{idat} files.
## @param temp.dir Temporary directory to be created for storing the symbolic links.
## @return \code{NULL} if the temporary directory could not be created or if no symbolic links could be created at all;
##         otherwise, the value of \code{temp.dir}. Additionally, an attribute named \code{"failed"} contains the number
##         of idat files in \code{rootdir} that could not be linked to.
## @author Pavlo Lutsik
prepare.idat.dir <- function(base.dir, temp.dir = tempfile(pattern = "idat")) {
	if (!dir.create(temp.dir)) {
		return(NULL)
	}
	idat.files <- list.files(base.dir, pattern = "^.+\\.idat$", recursive = TRUE)
	failed <- !sapply(idat.files, function(fname) {
		suppressWarnings(file.symlink(file.path(base.dir, fname), file.path(temp.dir, basename(fname))))
	})
	if (all(failed)) {
		return(NULL)
	}
	result <- temp.dir
	attr(result, "failed") <- as.integer(sum(failed))
	return(result)
}

########################################################################################################################

## check.idat.subdirs
##
## Checks whether the given directory contains subdirectories with IDAT files.
##
## @param rootdir Base directory to scan.
## @return \code{TRUE} if \code{base.dir} contains a subdirectory that contains one or more \code{idat} files;
##         \code{FALSE} otherwise.
## @author Pavlo Lutsik
check.idat.subdirs <- function(base.dir) {
	subdirectories <- dir(base.dir, full.names = TRUE)
	subdirectories <- subdirectories[file.info(subdirectories)[, "isdir"] == FALSE]
	for (subdir in subdirectories) {
		if (length(dir(subdir, pattern = "^.+\\.idat$", recursive = TRUE, include.dirs = FALSE)) != 0) {
			return(TRUE)
		}
	}
	return(FALSE)
}

########################################################################################################################

## mergeColumns
##
## Given a matrix X and a list containing column indices for each group in the matrix apply mergeFun to each submatrix
## defined by the column indices in the list
## @param X the input matrix
## @param mergeList a list containing column indices defining subgroups in matrix X
## @param mergeFun the function to be applied to each submatrix
## @return the combined matrix
## @note Requires the packages \pkg{foreach} and \pkg{doParallel}.
##
## @author Fabian Mueller
mergeColumns <- function(X,mergeList,mergeFun=function(X.sub){rowMeans(X.sub,na.rm=TRUE)}){
	res <- sapply(mergeList,FUN=function(iis){
		do.call(mergeFun,list(X.sub=as.matrix(X[,iis])))
	})
	return(res)
}
