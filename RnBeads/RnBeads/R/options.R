########################################################################################################################
## options.R
## created: 2012-05-10
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## RnBeads options and their accessor and modifier functions.
########################################################################################################################

## G L O B A L S #######################################################################################################

.rnb.options <- new.env()

## Supported methods to aggregate several methylation values of CpGs to a single value (region methylation)
AGGREGATION.METHODS <- c("min", "max", "mean", "median", "sum", "coverage.weighted")

## Supported normalization methods
NORMALIZATION.METHODS <- c("none", "bmiq", "illumina", "swan", "minfi.funnorm",
	"wm.dasen", "wm.nasen", "wm.betaqn", "wm.naten", "wm.nanet", "wm.nanes", "wm.danes",
	"wm.danet", "wm.danen", "wm.daten1", "wm.daten2", "wm.tost", "wm.fuks", "wm.swan")

## SNP-based filtering criteria
FILTERING.SNPS <- c("no", "3", "5", "any", "yes")

PROFILE.FILES <- c("450k_minimal", "450k_moderate", "450k_full", "bs_minimal", "bs_moderate", "bs_full")

INFERENCE.SVA.NUM.METHODS <- c("leek","be")

## F U N C T I O N S ###################################################################################################

#' rnb.load.options
#'
#' Loads the option names, types and default values from the \code{"options.txt"} file in the package. This function is
#' called upon initialization only.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.load.options <- function() {
	fname <- system.file(file.path("extdata", "options.txt"), package = "RnBeads")
	table.options <- read.delim(fname, quote = "", stringsAsFactors = FALSE)
	rownames(table.options) <- table.options[, "Name"]

	parse.value <- function(v.type, txt) {
		if (identical(txt, "NULL")) {
			return(NULL)
		}
		txt.names <- NULL
		regex.vector <- "^(.+)\\.vector$"
		if (grepl(regex.vector, v.type)) {
			regex.brackets <- "^\\((.*)\\)$"
			if (!grepl(regex.brackets, txt)) {
				stop("invalid value for txt")
			}
			txt <- gsub(regex.brackets, "\\1", txt)
			txt <- strsplit(txt, ",", fixed = TRUE)[[1]]
			regex.names <- "^(.+)=(.+)$"
			if (length(txt) != 0) {
				if (all(grepl(regex.names, txt))) {
					txt.names <- gsub(regex.names, "\\1", txt)
					txt <- gsub(regex.names, "\\2", txt)
				}
			}
			v.type <- gsub(regex.vector, "\\1", v.type)
		}
		
		if (v.type == "character" && identical(txt, "")) {
			result <- character()
		} else {
			result <- tryCatch(do.call(paste0("as.", v.type), list(txt)),
				warning = function(x) { stop("invalid value for txt") })
			names(result) <- txt.names
		}
		result
	}

	OPTION.TYPES <- table.options[, "Type"]
	names(OPTION.TYPES) <- rownames(table.options)
	assign("OPTION.TYPES", OPTION.TYPES, .rnb.options)
	assign("OPTIONS.NULL.ACCEPTED", rownames(table.options)[table.options[, "Null"] == "yes"], .rnb.options)
	option.values <- mapply(parse.value, as.list(OPTION.TYPES), as.list(table.options[, "Default"]),
			SIMPLIFY = FALSE)
	assign("current", option.values, .rnb.options)
	assign("previous", list(), .rnb.options)
}
rnb.load.options()
rm(rnb.load.options)

########################################################################################################################

## rnb.validate.option
##
## Validates the provided values for an option is acceptable, and converts it if necessary.
##
## @return The (possibly modified) option value.
## @author Yassen Assenov
rnb.validate.option <- function(oname, ovalue) {
	if (!(oname %in% names(.rnb.options[["OPTION.TYPES"]]))) {
		stop(paste(oname, "is invalid option"))
	}
	if (is.null(ovalue)) {
		if (!(oname %in% .rnb.options[["OPTIONS.NULL.ACCEPTED"]])) {
			stop(paste("invalid value for option", oname))
		}
		## Option value is NULL and this is considered valid
		if (oname == "qc.coverage.threshold.plot") {
			return(integer())
		}
		if (oname == "filtering.context.removal") {
			return(character())
		}
		return(NULL)
	}
	otype <- .rnb.options[["OPTION.TYPES"]][oname]
	## Infer type for options that accept multiple types of values
	if (oname == "identifiers.column") {
		if (is.double(ovalue) || is.integer(ovalue)) {
			otype <- "integer"
		}
	} else if (oname == "import.bed.columns") {
		if (is.character(ovalue) && length(ovalue) >= 1) {
			v.names <- as.character(ovalue)
			ovalue <- 1:length(ovalue)
			names(ovalue) <- v.names
			rm(v.names)
			supported.columns <- c("chr", "start", "end", "strand", "meth", "coverage", "c", "t")
			ovalue <- ovalue[v.names %in% supported.columns]
		}
	} else if (oname == "filtering.snp") {
		if (is.integer(ovalue) && length(ovalue) == 1 && (as.character(ovalue) %in% FILTERING.SNPS)) {
			ovalue <- as.character(ovalue)
		}
	} else if (oname %in% c("columns.pairing", "exploratory.columns", "differential.comparison.columns",
			"differential.comparison.columns.all.pairwise", "covariate.adjustment.columns")) {
		if (is.character(ovalue)) {
			otype <- "character.vector"
		}
	}

	## Validate type
	ovnames <- names(ovalue)
	if (otype == "character") {
		valid <- (is.character(ovalue) && length(ovalue) == 1 && (!is.na(ovalue)))
		if (valid) { ovalue <- ovalue[1] }
	} else if (otype == "double") {
		if (is.integer(ovalue)) {
			ovalue <- as.double(ovalue)
			names(ovalue) <- ovnames
		}
		valid <- (is.double(ovalue) && length(ovalue) == 1 && (!is.na(ovalue)))
		if (valid) { ovalue <- ovalue[1] }
	} else if (otype == "integer") {
		if (is.double(ovalue) && all(ovalue == as.integer(ovalue), na.rm = TRUE)) {
			ovalue <- as.integer(ovalue)
			names(ovalue) <- ovnames
		}
		valid <- (is.integer(ovalue) && length(ovalue) == 1 && (!is.na(ovalue)))
		if (valid) { ovalue <- ovalue[1] }
	} else if (otype == "logical") {
		valid <- parameter.is.flag(ovalue)
		if (valid) { ovalue <- ovalue[1] }
	} else if (otype == "character.vector") {
		valid <- (is.character(ovalue) && all(!is.na(ovalue)))
	} else { # otype == "integer.vector"
		if (is.double(ovalue) && all(ovalue == as.integer(ovalue), na.rm = TRUE)) {
			ovalue <- as.integer(ovalue)
			names(ovalue) <- ovnames
		}
		valid <- TRUE
		if (oname != "import.bed.columns") {
			valid <- (is.integer(ovalue) && all(!is.na(ovalue)))
		}
	}
	if (oname == "columns.pairing"){
		if (is.null(names(ovalue))){
			valid <- FALSE
		}
	}
	if (!valid) {
		stop(paste("invalid value for option", oname))
	}
	## Validate value(s)
	if (oname == "assembly") {
		valid <- ovalue %in% rnb.get.assemblies()
	} else if (oname == "region.aggregation") {
		valid <- ovalue %in% AGGREGATION.METHODS
	} else if (oname %in% c("points.category", "colors.category", "colors.gradient")) {
		valid <- (length(ovalue) > 1)
	} else if (oname == "min.group.size") {
		valid <- (ovalue > 0)
	} else if (oname == "max.group.count") {
		valid <- (ovalue >= 2)
	} else if (oname == "import.default.data.type") {
		valid <- (ovalue %in% c("data.dir", "idat.dir", "GS.report", "GEO", "bed.dir"))
	} else if (oname == "import.bed.style") {
		valid <- (ovalue %in% c("BisSNP", "Encode", "EPP", "bismarkCytosine", "bismarkCov", "custom"))
	} else if (oname == "import.bed.columns") {
		v.names <- names(ovalue)
		supported.columns <- c("chr", "start", "end", "strand", "meth", "coverage", "c", "t")
		if (is.null(v.names) || (!all(v.names %in% supported.columns))) {
			valid <- FALSE
		} else {
			ovalue <- ovalue[supported.columns]
			names(ovalue) <- supported.columns
			valid <- !(any(is.na(ovalue[c("chr", "start")])) ||
							(is.na(ovalue["meth"]) && sum(is.na(ovalue[c("coverage", "c", "t")])) > 1))
		}
		valid <- valid && all(ovalue > 0, na.rm = TRUE) && anyDuplicated(ovalue,incomparables=NA) == 0
	} else if (oname == "import.idat.chunk.size") {
		valid <- is.null(ovalue) || ovalue > 1
	} else if (oname == "qc.coverage.threshold.plot") {
		ovalue <- sort(unique(as.integer(ovalue)))
		valid <- isTRUE(all(ovalue > 0))
	} else if (oname == "normalization.method") {
		ovalue <- tolower(ovalue)
		valid <- (ovalue %in% NORMALIZATION.METHODS)  
	} else if (oname == "normalization.background.method") {
		valid <- (ovalue %in% c("none", "methylumi.noob", "methylumi.goob", "methylumi.lumi"))
	} else if (oname == "filtering.context.removal") {
		valid <- all(ovalue %in% c("CG", "CC", "CAG", "CAH", "CTG", "CTH", "Other"))
		ovalue <- intersect(c("CG", "CC", "CAG", "CAH", "CTG", "CTH", "Other"), ovalue)
	} else if (oname == "filtering.snp") {
		valid <- (ovalue %in% FILTERING.SNPS)
	} else if (oname %in% c("filtering.coverage.threshold", "exploratory.correlation.permutations",
		"differential.permutations")) {
		valid <- (ovalue >= 0)
	} else if (oname %in% c("filtering.greedycut.pvalue.threshold", "filtering.deviation.threshold",
		"filtering.missing.value.quantile", "exploratory.correlation.pvalue.threshold")) {
		valid <- (0 <= ovalue && ovalue <= 1)
	} else if (oname == "filtering.greedycut.rc.ties") {
		valid <- (ovalue %in% c("row", "column", "any"))
	} else if (oname == "inference.sva.num.method") {
		ovalue <- tolower(ovalue)
		valid <- (ovalue %in% INFERENCE.SVA.NUM.METHODS)  
	} else if (oname == "exploratory.principal.components") {
		valid <- (0 <= ovalue && ovalue <= 10)
	} else if (oname == "exploratory.clustering") {
		valid <- (ovalue %in% c("all", "top", "none"))
	} else if (oname == "exploratory.clustering.top.sites") {
		valid <- (length(ovalue) > 0 && all(ovalue > 0))
		if (valid) { ovalue <- unique(sort(ovalue)) }
	} else if (oname == "differential.site.test.method") {
		valid <- (ovalue %in% c("ttest","limma","refFreeEWAS"))
	} else if (oname == "colors.3.gradient") {
		valid <- (length(ovalue) == 3)
	} else if (oname == "colors.meth") {
		valid <- (length(ovalue) %in% c(2, 3))
	}
	if (!valid) {
		stop(paste("invalid value for option", oname))
	}
	## Modify values of linked options
	if (oname == "import.bed.style") {
		.rnb.options[["previous"]][["import.bed.columns"]] <- .rnb.options[["current"]][["import.bed.columns"]]
		if (ovalue == "Encode") {
			.rnb.options[["current"]][["import.bed.columns"]] <-
				c(chr = 1L, start = 2L, end = 3L, strand = 6L, meth = 4L, coverage = 5L, c = NA, t = NA)
		} else if (ovalue == "BisSNP"){ # ovalue == "BisSNP"
			.rnb.options[["current"]][["import.bed.columns"]] <-
				c(chr = 1L, start = 2L, end = 3L, strand = 6L, meth = 4L, coverage = 5L, c = NA, t = NA)
		} else if (ovalue == "bismarkCytosine"){
			.rnb.options[["current"]][["import.bed.columns"]] <-
				c(chr = 1L, start = 2L, end = NA, strand = 3L, meth = NA, coverage = NA, c = 4L, t = 5L)
		} else if (ovalue == "bismarkCov"){
			.rnb.options[["current"]][["import.bed.columns"]] <-
				c(chr = 1L, start = 2L, end = NA, strand = NA, meth = NA, coverage = NA, c = 5L, t = 6L)
		} #else assume custom style or EPP style which is handled in read.bissnp.bed() 
	}
	
	return(ovalue)
}

########################################################################################################################

## rnb.get.option
##
## Gets or sets the value of a single RnBeads option.
##
## @param oname    Name of the option to get or set.
## @param ovalue   New value for the option. This is used only if \code{setvalue} is \code{TRUE}.
## @param setvalue Flag indicating if the value of the option is to be changed.
## @return Empty \code{character} string if the operation was successful; the text of an error message otherwise.
## @author Yassen Assenov
rnb.get.option <- function(oname, ovalue = NULL, setvalue = FALSE) {
	if (!(oname %in% names(.rnb.options[["current"]]))) {
		return(paste(oname, "is invalid option"))
	}
	result <- ""
	cvalue <- .rnb.options[["current"]][[oname]]
	if (is.null(cvalue)) {
		.rnb.options[["previous"]][oname] <- list(NULL) 
	} else {
		.rnb.options[["previous"]][[oname]] <- cvalue
	}
	if (setvalue) {
		ovalue <- tryCatch(rnb.validate.option(oname, ovalue), error = function(e) { e })
		if (inherits(ovalue, "error")) {
			result <- ovalue$message
		} else if (is.null(ovalue)) {
			.rnb.options[["current"]][oname] <- list(NULL)
		} else {
			.rnb.options[["current"]][[oname]] <- ovalue
		}
	}
	return(result)
}

########################################################################################################################

#' rnb.is.option
#'
#' Checks if the specified text is an option name.
#' 
#' @param txt Potential option name. This should be a one-element \code{character} vector.
#' @return \code{TRUE} if the specified parameter is a valid analysis option name; \code{FALSE} otherwise.
#' 
#' @seealso \code{\link{rnb.options}} for getting and setting option values
#' @examples
#' \dontrun{
#' rnb.is.option("logging") # TRUE
#' rnb.is.option("Logging") # FALSE
#' }
#' @author Yassen Assenov
#' @export
rnb.is.option <- function(txt) {
	return(is.character(txt) && length(txt) == 1 && txt %in% names(.rnb.options[["current"]]))
}

########################################################################################################################

#' RnBeads Options
#'
#' Allows the user to set and examine a variety of \pkg{RnBeads} global options. They affect the way in which the
#' package computes and displays its results.
#'
#' @rdname rnb.options
#' @aliases rnb.getOption
#'
#' @param ... Option names as \code{character}s, or new option values given in the form \code{name = value}.
#' @param x   Option name in the form of a \code{character} vector of length 1.
#' @return For \code{rnb.getOption}, the current value for \code{x}. For \code{rnb.options()}, a list of all
#'         \pkg{RnBeads} options and their current values. If option names are given, a list of all requested options
#'         and their values. If option values are set, \code{rnb.options} returns the previous values of the modified
#'         options, invisibly.
#'
#' @details
#' Invoking \code{rnb.options()} with no arguments returns a list with the current values of the options. To access the
#' value of a single option, one should use, e.g., \code{rnb.getOption("filtering.greedycut")}, rather than
#' \code{rnb.options("filtering.greedycut")} which is a \emph{list} of length one. Also, only a limited set of options
#' is available (see below). Attempting to get or set the value of a non-existing option results in an error.
#'
#' @section Options used in RnBeads:
#' \describe{
#'   \item{\bold{\code{analysis.name}}\code{ = NULL}}{
#'        One-element \code{character} vector storing a short title of the analysis. If specified, this name appears at
#'        the page title of every report.}
#'   \item{\bold{\code{logging}}\code{ = TRUE}}{
#'        Flag indicating if logging functionality is enabled in the automatic runs of the pipeline.}
#'   \item{\bold{\code{email}}\code{ = NULL}}{
#'        Email address associated with the analyses.}
#'   \item{\bold{\code{assembly}}\code{ = "hg19"}}{
#'        Genome assembly to be used. Currently only important for bisulfite mode. The supported genomes returned by the
#'        function \code{\link{rnb.get.assemblies}}.}
#'   \item{\bold{\code{analyze.sites}}\code{ = TRUE}}{
#'        Flag indicating if analysis on site or probe level is to be conducted. Note that the preprocessing module
#'        always operates on the site level (only), regardless of the value of this option.}
#'   \item{\bold{\code{region.types}}\code{ = NULL}}{
#'        Region types to carry out analysis on, in the form of a \code{character} vector. \code{NULL} (default value)
#'        signifies that all available region annotations (as returned by \code{\link{rnb.region.types}}) are summarized 
#'        upon loading and normalization, and the other modules analyze all regions summarized in the dataset. If this
#'        option is set to an empty vector, analysis on the region level is skipped.}
#'   \item{\bold{\code{region.aggregation}}\code{ = "mean"}}{
#'        Aggregation function to apply when calculating the methylation value for a region based on the values of the
#'        CpGs associated with that region. Accepted values for this function are \code{"min"}, \code{"max"},
#'        \code{"mean"} (default), \code{"median"}, \code{"sum"}, \code{"coverage.weighted"}. The last method is
#'        applicable only for sequencing-based methylation datasets. It computes the weighted average of the values of
#'        the associated CpGs, whereby weights are calculated based on the coverages of the respective sites.}
#'   \item{\bold{\code{region.subsegments}}\code{ = 0}}{
#'        If a number larger than 1 is specified, \pkg{RnBeads} will subdivide each region specified in the
#'        \code{region.types} option into subsegments containing on average \code{region.subsegments} sites per 
#'        subsegment. This is done by clustering the sites within each regions according to their genomic coordinates.
#'        These subsegments are then used for subsequent analysis. 
#'        Use cautiously as this will significantly increase the runtime of the pipeline.}
#'   \item{\bold{\code{region.subsegments.types}}\code{ = NULL}}{
#'        The region types to which subsegmentation will be applied. Defaults to \code{region.types} when set to
#'        \code{NULL}.}
#'   \item{\bold{\code{identifiers.column}}\code{ = NULL}}{
#'        Column name or index in the table of phenotypic information to be used when plotting sample identifiers. If
#'        this option is \code{NULL}, it points to a non-existing column or a column that does not list IDs, the default
#'        identifiers are used. These are the row names of the sample phenotype table (and the column names of the beta
#'        value matrix).}
#'   \item{\bold{\code{colors.category}}\code{ = c("#1B9E77","#D95F02",...)}}{
#'        \code{character} vector of length 2 or more giving the color scheme for displaying categorical trait values in
#'        plots. RnBeads denotes missing values (\code{NA}) by grey, therefore, it is not recommended to include shades
#'        of grey in this vector. The default value of this option is the result of the \code{"Dark2"} palette of
#'        \emph{RColorBrewer} with 8 values.}
#'   \item{\bold{\code{colors.gradient}}\code{ = c("#132B43","#56B1F7")}}{
#'        \code{character} vector of length 2 or more giving the color scheme for displaying continuous (gradient) trait
#'        values in plots. \pkg{RnBeads} interpolates between the color values.}
#'   \item{\bold{\code{min.group.size}}\code{ = 2}}{
#'        Minimum number of samples each subgroup defined by a trait, in order for this trait to be considered in the
#'        methylation profiles and in the differential methylation modules. This must be a positive \code{integer}.}
#'   \item{\bold{\code{max.group.count}}\code{ = NULL}}{
#'        Maximum number of subgroups defined by a trait, in order for this trait to be considered in the methylation
#'        profiles and in the differential methylation modules. This must be an \code{integer} of value \code{2} or
#'        more. As a special case, a value of \code{NULL} (default) indicates that the maximum number of subgroups is
#'		  the number of samples in an analysis minus \code{1}, i.e. traits with all unique values will be ignored.}
#'   \item{\bold{\code{replicate.id.column}}\code{ = NULL}}{
#'        Column name in the sample annotation table that indicates sample replicates. Replicates are expected to
#'        contain the same value. Samples without replicates should contain unique or missing values. If this option is
#'        \code{NULL} (default), replicate handling is disabled.}
#'   \item{\bold{\code{gz.large.files}}\code{ = FALSE}}{
#'        Flag indicating whether large output files should be compressed (in \code{.gz} format).}
#'   \item{\bold{\code{import}}\code{ = TRUE}}{
#'        Flag controlling whether data import report should be generated. This option be set to \code{FALSE} only when
#'        the provided data source is an object of type \linkS4class{RnBSet}, i.e. the data has been previously loaded
#'        by \pkg{RnBeads}.}
#'   \item{\bold{\code{import.default.data.type}}\code{ = "infinium.idat.dir"}}{
#'        Type of data assumed to be supplied by default (Infinium 450k microarray). 
#' 		  For sequencing data set this to \code{bs.bed.dir} and save the options. 
#' 		  See \code{\link{rnb.execute.import}} for further details.}
#'   \item{\bold{\code{import.table.separator}}\code{ = ","}}{
#'        Separator used in the plain text data tables. See \code{\link{rnb.execute.import}} for details.}
#'   \item{\bold{\code{import.bed.style}}\code{ = "BisSNP"}}{
#' 		  Preset for bed-like formats. \code{"BisSNP", "Encode","EPP", "bismarkCytosine", "bismarkCov"} are currently
#' 		  supported. See the \pkg{RnBeads} vignette and the FAQ section on the website for more details.}
#'   \item{\bold{\code{import.bed.columns}}}{Column indices in the supplied BED file with DNA methylation information.
#'        These are represented by a named \code{integer} vector, in which the names are: \code{"chr"}, \code{"start"},
#'        \code{"end"}, \code{"strand"}, \code{"meth"}, \code{"coverage"}, \code{"c"} and \code{"t"}. These names
#'        correspond the columns for chromosome, start position, end position, strand, methylation degree, read
#'        coverage, number of reads with C and number of reads with T, respectively. Methylation degree and/or read
#'        coverage, if not specified, are inferred from the values in the columns \code{"c"} and \code{"t"}.
#' 		  Further details and examples of BED files can be found in Section 4.1 of the RnBeads vignette.}
#'   \item{\bold{\code{import.bed.frame.shift}}\code{ = 1}}{Singleton of type \code{integer} specifying the frame shift between
#'        the coordinates in the input BED file and the corresponding genomic reference. This (\code{integer}) value 
#'        is added to the coordinates from the BED file before matching the methylation sites to the annotated ones.}
#'   \item{\bold{\code{import.bed.test}}\code{ = TRUE}}{
#' 		  Perform a small loading test, by reading 1000 rows from each BED file, after which normal loading is performed.
#'        See \pkg{RnBeads} vignette and the FAQ section on the website for more details.}
#'   \item{\bold{\code{import.bed.test.only}}\code{ = FALSE}}{
#' 		  Perform only the small loading test, and skip loading all the data.}
#'   \item{\bold{\code{preprocessing}}\code{ = TRUE}}{Flag controlling whether the data should be preprocessed 
#' 		  (whether quality filtering and in case of Infinium microarray data normalization should be applied).}
#'   \item{\bold{\code{normalization}}\code{ = NULL}}{
#'        Flag controlling whether the data should be normalized and normalization report generated. Setting this to
#'        \code{NULL} (default) enables this step for analysis on Infinium datasets, but disables it in case of
#'        sequencing-based datasets. Note that normalization is never applied in sequencing datasets; if this flag is
#'        enabled, it will lead to a warning message.}
#'   \item{\bold{\code{normalization.method}}\code{ = "swan"}}{
#'        Normalization method to be applied, or \code{"none"}. Multiple normalization methods are supported:
#'        \code{"illumina"} -
#'        \href{http://www.bioconductor.org/packages/devel/bioc/html/methylumi.html}{methylumi}-implemented
#'        Illumina scaling normalization; \code{"swan"} (default) - SWAN-normalization by Gordon et al., as implemented
#'        in \href{http://www.bioconductor.org/packages/release/bioc/html/minfi.html}{minfi}; \code{"bmiq"} -
#'        beta-mixture quantile normalization method by Teschendorff et al; as well as \code{"wm.dasen"},
#'        \code{"wm.nasen"}, \code{"wm.betaqn"}, \code{"wm.naten"}, \code{"wm.nanet"}, \code{"wm.nanes"},
#'        \code{"wm.danes"}, \code{"wm.danet"}, \code{"wm.danen"}, \code{"wm.daten1"}, \code{"wm.daten2"},
#'        \code{"wm.tost"}, \code{"wm.fuks"} and \code{"wm.swan"} - all normalization methods implemented in the
#'        \href{http://www.bioconductor.org/packages/release/bioc/html/wateRmelon.html}{wateRmelon} package. When
#'        setting this option to a specific algorithm, make sure its dedicated package is installed.}
#'   \item{\bold{\code{normalization.background.method}}\code{ = "methylumi.noob"}}{
#'        A character singleton specifying which background subtraction is to be performed during normalization.
#'        The \href{http://www.bioconductor.org/packages/devel/bioc/html/methylumi.html}{methylumi} background
#'        correction methods are supported. The following values are accepted: \code{"none"}, \code{"methylumi.noob"},
#'        \code{"methylumi.goob"} and \code{"methylumi.lumi"}.}
#'   \item{\bold{\code{normalization.plot.shifts}}\code{ = TRUE}}{
#'        Flag indicating if the report on normalization should include plots of shifts (degrees of beta value
#'        correction).}
#'   \item{\bold{\code{qc}}\code{ = TRUE}}{
#'        Flag indicating if the quality control module is to be executed.}
#'   \item{\bold{\code{qc.boxplots}}\code{ = TRUE}}{
#'        [Infinium 450k] Add boxplots for all types of quality control probes to the quality control report. The boxplots give signal
#'        distribution across samples.}
#'   \item{\bold{\code{qc.barplots}}\code{ = TRUE}}{
#'        [Infinium 450k] Add barplots for each quality control probes to the quality control report.}
#'   \item{\bold{\code{qc.negative.boxplot}}\code{ = TRUE}}{
#'        [Infinium 450k] Add boxplot of negative control probe intensities for all samples.}
#'   \item{\bold{\code{qc.snp.distances}}\code{ = TRUE}}{
#'        [Infinium 450k] Flag indicating if intersample distances based on the beta values of SNP probes are to be displayed. This
#'        can help identify or validate genetically similar or identical samples.}
#'   \item{\bold{\code{qc.snp.boxplot}}\code{ = FALSE}}{
#'        [Infinium 450k] Add boxplot of beta-values for the SNP-calling probes. Can be useful for detection of sample mix-ups.}
#'   \item{\bold{\code{qc.snp.barplot}}\code{ = FALSE}}{
#'        [Infinium 450k] Add bar plots of beta-values for the SNP-calling probes in each profiled sample.}
#'   \item{\bold{\code{qc.sample.batch.size}}\code{ = 50}}{
#'        [Infinium 450k] Maximal number of samples included in a single quality control barplot and negative control boxplot.}
#'   \item{\bold{\code{qc.coverage.plots}}\code{ = FALSE}}{
#'        [Bisulfite sequencing] Add genome-wide sequencing coverage plot for each sample.}
#'   \item{\bold{\code{qc.coverage.threshold.plot}}\code{ = 1:10}}{
#'        [Bisulfite sequencing] Values for coverage cutoffs to be shown in a coverage thresholds plot. This must be an \code{integer}
#'        vector of positive values. Setting this to an empty vector disables the coverage thresholds plot.}
#'   \item{\bold{\code{qc.coverage.histograms}}\code{ = FALSE}}{
#'        [Bisulfite sequencing] Add sequencing coverage histogram for each sample.}
#'   \item{\bold{\code{qc.coverage.violins}}\code{ = FALSE}}{
#'        [Bisulfite sequencing] Add sequencing coverage violin plot for each sample.}
#'   \item{\bold{\code{filtering.whitelist}}\code{ = NULL}}{Name of a file specifying site or probe identifiers to be
#'        whitelisted. Every line in this file must contain exactly one identifier. The whitelisted sites are always
#'        retained in the analysed datasets, even if filtering criteria or blacklisting requires their removal.
#'        For Infinium studies, the file must contain Infinium probe identifiers. For bisulfite sequencing studies,
#'        the file must contain CpG positions in the form "chromosome:coordinate" (1-based coordinate of the cytosine),
#'        e.g. \code{chr2:48607772}. Unknown identifiers are silently ignored.}
#'   \item{\bold{\code{filtering.blacklist}}\code{ = NULL}}{Name of a file specifying site or probe identifiers to be
#'        blacklisted. Every line in this file must contain exactly one identifier. The blacklisted sites are removed
#'        from the analysed datasets as a first step in the preprocessing module. For Infinium studies, the file must
#'        contain Infinium probe identifiers. For bisulfite sequencing studies, the file must contain CpG positions in
#'        the form "chromosome:coordinate" (1-based coordinate of the cytosine), e.g. \code{chr2:48607772}. 
#'        Unknown identifiers are silently ignored.}
#'   \item{\bold{\code{filtering.context.removal}}\code{ = c("CC","CAG",...)}}{
#'        \code{character} vector giving the list of probe context types to be removed as a filtering step. Possible
#'        context values are \code{"CC"}, \code{"CG"}, \code{"CAG"}, \code{"CAH"}, \code{"CTG"}, \code{"CTH"} and
#'        \code{"Other"}. Probes in the second context measure CpG methylation; the last context denotes probes
#'        dedicated to SNP detection. Setting this option to \code{NULL} or an empty vector effectively disables the
#'        step of context-specific probe removal.}
#'   \item{\bold{\code{filtering.snp}}\code{ = "3"}}{
#'        Removal of sites or probes based on overlap with SNPs. The accepted values for this option are:
#'        \describe{
#'           \item{\code{"no"}}{no SNP-based filtering;}
#'           \item{\code{"3"}}{filter out a probe when the last 3 bases in its target sequence overlap with SNP;}
#'           \item{\code{"5"}}{filter out a probe when the last 5 bases in its target sequence overlap with SNP;}
#'           \item{\code{"any"} or \code{"yes"}}{filter out a CpG site or probe when any base in its target sequence
#'                overlaps with SNP.}}
#'        Bisulfite sequencing datasets operate on sites instead of probes, therefore, the values \code{"3"} and
#'        \code{"5"} are treated as \code{"yes"}.}
#'   \item{\bold{\code{filtering.sex.chromosomes.removal}}\code{ = FALSE}}{
#'        Flag indicating if the removal of probes located on sex chromosomes should be performed as a filtering step.}
#'   \item{\bold{\code{filtering.missing.value.quantile}}\code{ = 1}}{
#'        Number between 0 and 1, indicating the fraction of allowed missing values per site. A site is filtered out
#'        when its methylation beta values are \code{NA}s in a larger fraction of samples than this threshold. Setting
#'        this option to 1 (default) retains all sites, and thus effectively disables the missing value filtering step
#'        in the preprocessing module. If this is set to 0, all sites that contain missing values are filtered out.}
#'   \item{\bold{\code{filtering.coverage.threshold}}\code{ = 5}}{
#'        Threshold for minimal acceptable coverage. This must be a non-negative value. Setting this option to 0 (zero)
#'        effectively considers any known or unknown read coverage for sufficiently deep.}
#'   \item{\bold{\code{filtering.low.coverage.masking}}\code{ = FALSE}}{
#'        Flag indicating whether methylation values for low coverage sites should be set to missing. In combination
#'        with \code{filtering.missing.value.quantile} this can lead to the removal of sites.}
#'   \item{\bold{\code{filtering.high.coverage.outliers}}\code{ = FALSE}}{
#'        (Bisulfite sequencing mode) Flag indicating whether methylation sites with a coverage of more than 10 times
#'        the 95-percentile of coverage should be removed.}
#'   \item{\bold{\code{filtering.greedycut}}\code{ = TRUE}}{
#'        Flag indicating if the Greedycut procedure should be run as part of the preprocessing module.}
#'   \item{\bold{\code{filtering.greedycut.pvalue.threshold}}\code{ = 0.05}}{
#'        Threshold for the detection p-value to be used in Greedycut. This is a value between 0 and 1. This option has
#'        effect only when \code{filtering.greedycut} is \code{TRUE}.}
#'   \item{\bold{\code{filtering.greedycut.rc.ties}}\code{ = "row"}}{
#'        Indicator of what the behaviour of Greedycut should be in case of ties between the scores of rows (probes) and
#'        columns (samples). The value of this option must be one of \code{"row"}, \code{"column"} or \code{"any"}; the
#'        last one indicating random choice. This option has effect only when \code{filtering.greedycut} is
#'        \code{TRUE}.}
#'   \item{\bold{\code{filtering.deviation.threshold}}\code{ = 0}}{
#'        Threshold used to filter probes based on the variability of their assigned beta values. This must be a real
#'        value between 0 and 1, denoting minimum standard deviation of the beta values in one site across all samples.
#'        Any sites that have standard deviation lower than this threshold are filtered out. Note that sites with
#'        undetermined varibility, that is, sites for which there are no measurements (all beta values are \code{NA}s),
#'        are retained. Setting this option to 0 (default) disables filtering based on methylation variability.}
#'   \item{\bold{\code{inference}}\code{ = FALSE}}{
#'        Flag indicating if the covariate inference analysis module is to be executed.}
#'   \item{\bold{\code{inference.targets.sva}}\code{ = character()}}{
#'        Column names in the sample annotation table for which surrogate variable analysis (SVA) should be conducted.
#'		  An empty vector (default) means that SVA is skipped.}
#'   \item{\bold{\code{inference.reference.methylome.column}}\code{ = character()}}{
#'        Column name in the sample annotation table giving the assignment of samples to reference methylomes.
#' 		  The target samples should have \code{NA} values in this column.}
#'   \item{\bold{\code{inference.max.cell.type.markers}}\code{ = 10000}}{
#'        A number of most variable CpGs which are tested for association with the reference cell types.}
#'   \item{\bold{\code{inference.top.cell.type.markers}}\code{ = 500}}{
#'        The number of top cell type markers used for determining cell type contributions to the target 
#' 		  DNA methylation profiles using the projection method of Houseman et al.}
#'   \item{\bold{\code{inference.sva.num.method}}\code{ = "leek"}}{
#'        Name of the method to be used for estimating the number of surrogate variables. 
#'		  must be either 'leek' or 'be', See \code{sva} function for details.}
#'   \item{\bold{\code{exploratory}}\code{ = TRUE}}{
#'        Flag indicating if the exploratory analysis module is to be executed.}
#'   \item{\bold{\code{exploratory.columns}}\code{ = NULL}}{
#'        Traits, given as column names or indices in the sample annotation table, to be used in the exploratory
#'        analysis. These traits are used in multiple steps in the module: they are visualized using point types and
#'        colors in the dimension reduction plots; tested for strong correlations and associations with principal
#'        components in a methylation space; used to define groups when plotting beta distributions and/or inter-sample
#'        methylation variability. The default value of this parameter - \code{NULL} - indicates that columns should be
#'        automatically selected; see \code{\link{rnb.sample.groups}} for how this is done.}
#'   \item{\bold{\code{exploratory.top.dimensions}}\code{ = 0}}{
#'        Number of most variable probes, sites or regions to select prior to performing dimension reduction techniques
#'        and tests for associations. Preselection can significantly reduce the running time and memory usage in the
#'        exploratory analysis module. Setting this number to zero (default) disables preselection.}
#'   \item{\bold{\code{exploratory.principal.components}}\code{ = 8}}{
#'        Maximum number of principal components to be tested for associations with other factors, such as control probe
#'        states and sample traits. This must be an \code{integer} value between \code{0} and \code{10}. Setting this
#'        option to \code{0} disables such tests.}
#'   \item{\bold{\code{exploratory.correlation.pvalue.threshold}}\code{ = 0.01}}{
#'        Significance threshold for a p-value resulting from applying a test for association. This is a value between
#'        0 and 1.}
#'   \item{\bold{\code{exploratory.correlation.permutations}}\code{ = 10000}}{
#'        Number of permutations in tests performed to check for associations between traits, and between control probe
#'        intensities and coordinates in the prinicipal component space. This must be a non-negative \code{integer}.
#'        Setting this option to \code{0} disables permutation tests.}
#'   \item{\bold{\code{exploratory.correlation.qc}}\code{ = TRUE}}{
#'        [Infinium 450k] Flag indicating if quality-associated batch effects should be studied. This amounts to testing for
#'        associations between intensities of quality control probes and principal components. This option has effect
#'        only when \code{exploratory.principal.components} is non-zero.}
#'   \item{\bold{\code{exploratory.beta.distribution}}\code{ = TRUE}}{
#'        Flag indicating whether beta value distributions for sample groups and probe or site categories should be
#'        computed.}
#'   \item{\bold{\code{exploratory.intersample}}\code{ = TRUE}}{
#'        Flag indicating if methylation variability in sample groups should be computed as part of the exploratory
#'        analysis module.}
#'   \item{\bold{\code{exploratory.deviation.plots}}\code{ = NULL}}{
#'        Flag indicating if the inter-sample methylation variability step in the exploratory analysis module should
#'        include deviation plots. Deviation plots show intra-group methylation variability at the covered sites and
#'        regions. Setting this option to \code{NULL} (default) enables deviation plots on Infinium datasets, but
#'        disables them in case of sequencing-based datasets, because their generation can be very computationally
#'        intensive. This option has effect only when \code{exploratory.intersample} is \code{TRUE}.}
#'   \item{\bold{\code{exploratory.clustering}}\code{ = "all"}}{
#'        Which sites should be used by clustering algorithms in the exploraroty analysis module.
#'        \pkg{RnBeads} performs several algorithms that cluster the samples in the dataset. If this option is set to
#'        \code{"all"} (default), clustering is performed using all sites; a value of \code{"top"} indicates that only
#'        the most variable sites are used (see the option \code{exploratory.clustering.top.sites}); and \code{"none"}
#'        disables clustering.}
#'   \item{\bold{\code{exploratory.clustering.top.sites}}\code{ = 1000}}{
#'        Number of most variable sites to use when visualizing heatmaps. This must be a non-empty \code{integer} vector
#'        containing positive values. This option is ignored when \code{exploratory.clustering} is \code{"none"}.}
#'   \item{\bold{\code{exploratory.clustering.heatmaps.pdf}}\code{ = FALSE}}{
#'        Flag indicating if the generated methylation value heatmaps in the clustering section of the exploratory
#'        analysis module should be saved as PDF files. Enabling this option is not recommended for large values of
#'        \code{exploratory.clustering.top.sites} (more than 200), because heatmaps might generate very large PDF files.}
#'   \item{\bold{\code{exploratory.region.profiles}}\code{ = NULL}}{
#'        Region types for generating regional methylation profiles. If \code{NULL} (default), regional methylation
#'        profiles are created only for the region types that are available for the targeted assembly and summarized in
#'        the dataset of interest. Setting this option to an empty vector disables the region profiles step in the
#'        exploratory analysis module.}
#'   \item{\bold{\code{exploratory.gene.symbols}}\code{ = NULL}}{
#'        A list of gene symbols to be used for custom locus profiling. Locus views will be generated for these genes.
#'        }
#'   \item{\bold{\code{exploratory.custom.loci.bed}}\code{ = NULL}}{
#'        Path to a bed file containing custom genomic regions. Locus views will be generated for these regions.
#'        }
#'   \item{\bold{\code{differential}}\code{ = TRUE}}{
#'        Flag indicating if the differential methylation module is to be executed.}
#'   \item{\bold{\code{differential.site.test.method}}\code{ = "limma"}}{
#'        Method to be used for calculating p-values on the site level. Currently supported options are "ttest" for a (paired)
#'        t-test and "limma" for a linear modeling approach implemented in the \code{limma} package for differential expression
#'        in microarrays.} 
#'   \item{\bold{\code{differential.permutations}}\code{ = 0}}{
#'        Number of permutation tests performed to compute the p-value of rank permutation tests in the differential
#'        methylation analysis. This must be a non-negative \code{integer}. Setting this option to \code{0} (default)
#'        disables permutation tests for rank permutations. Note that p-values for differential methylation are
#'        computed and also considered for the ranking in any case.}
#'   \item{\bold{\code{differential.comparison.columns}}\code{ = NULL}}{
#'        Column names or indices in the table of the sample annotation table to be used for group definition in the
#'        differential methylation analysis. The default value - \code{NULL} - indicates that columns should be
#'        automatically selected. See\code{\link{rnb.sample.groups}} for how this is done. By default, 
#'		  the comparisons are done in a one vs. all manner if there are multiple
#'		  groups defined in a column. }
#'   \item{\bold{\code{differential.comparison.columns.all.pairwise}}\code{ = NULL}}{
#'        Column names or indices in the table of sample annotation table to be used for group definition in the
#'        differential methylation analysis in which all pairwise comparisons between groups should be conducted (the default
#'		  is one vs all if multiple groups are specified in a column).
#'        Caution: for large numbers of sample groups this can lead to combinatorial explosion and thus to huge runtimes.
#'        A value of \code{NULL} (default) indicates that no column is selected for all pairwise comparisons explicitely.
#'        If specified, the selected columns must be a subset of the columns that will be selected according to the 
#'        \code{differential.comparison.columns} option.}
#'   \item{\bold{\code{covariate.adjustment.columns}}\code{ = NULL}}{
#'        Column names or indices in the table of phenotypic information to be used for confounder adjustment in the
#'        differential methylation analysis. Currently this is only supported for \code{differential.site.test.method=="limma"}.
#'        }
#'	 \item{\bold{\code{columns.pairing}}\code{ = NULL}}{
#' 		  A NAMED vector containing for each column name for which paired analysis
#' 		  should be performed (say columnA) the name or index of another column (say columnB) in which same values indicate
#' 		  the same pairing. columnA should be the name of the value columnB in this vector.
#' 		  For more details see \code{\link{rnb.sample.groups}}}
#'   \item{\bold{\code{differential.adjustment.sva}}\code{ = TRUE}}{
#'        Flag indicating if the differential methylation analysis should account for Surrogate Variables. If
#'        \code{TRUE}, \pkg{RnBeads} looks for overlaps between the \code{differential.comparison.columns} and
#'        \code{inference.targets.sva} options and include the surrogate variables as confounding factors only for these
#'        columns. In other words, it will only have an effect if the corresponding inference option
#'		  (see \code{inference.targets.sva} option for details) is enabled.
#'		  Currently this is only supported for \code{differential.site.test.method=="limma"}.}
#'   \item{\bold{\code{differential.adjustment.celltype}}\code{ = TRUE}}{
#'        Should the differential methylation analysis account for celltype using the reference based Houseman method.
#'        It will only have an effect if the corresponding inference option is enabled (see \code{inference.reference.methylome.column}
#'		  option for details). Currently this is only supported for \code{differential.site.test.method=="limma"}.
#'        }
#'   \item{\bold{\code{differential.enrichment}}\code{ = FALSE}}{
#'        Flag indicating whether \href{http://www.geneontology.org/}{Gene Ontology} (GO)-enrichment analysis is to be
#'        conducted on the identified differentially methylated regions.}
#'   \item{\bold{\code{export.to.bed}}\code{ = TRUE}}{
#'        Flag indicating whether the data should be exported to bed files.}
#'   \item{\bold{\code{export.to.trackhub}}\code{ = c("bigBed","bigWig")}}{
#'        \code{character} vector specifying which data types should be exported to
#'        \href{http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html}{Track hub directories}. Possible values
#'        in the vector are \code{"bigBed"} and \code{"bigWig"}. When this options is set to \code{NULL}, track hub
#'        export is disabled. Note that if \code{"bigBed"} is contained in this option, bed files are created
#'        automatically.}
#'   \item{\bold{\code{export.to.csv}}\code{ = FALSE}}{
#'        Flag indicating whether methylation value matrices are to be exported to comma-separated value (CSV) files.}
#'   \item{\bold{\code{export.to.ewasher}}\code{ = FALSE}}{
#'        Flag indicating whether methylation values and differential methylation analysis settings should be exported to
#'		  a format compatible with FaST-LMM-EWASher, a tool for adjusting for cell-type compositions.
#'		  See \href{http://www.nature.com/nmeth/journal/v11/n3/full/nmeth.2815.html}{Zou, J., et al., Nature Methods, 2014} for further details on the tool.}
#'   \item{\bold{\code{export.types}}\code{ = "sites"}}{
#'        \code{character} vector of sites and region names to be exported. If \code{NULL}, no region methylation values
#'        are exported.}
#'   \item{\bold{\code{disk.dump.big.matrices}}\code{ = FALSE}}{
#'        Flag indicating whether big tables should be stored on disk rather than in main memory in order to keep memory
#'        requirements down. May slow down analysis!}
#'   \item{\bold{\code{logging.exit.on.error}}\code{ = FALSE}}{
#'        Flag indicating if the active R session should be terminated when an error is encountered during execution.}
#'   \item{\bold{\code{distribution.subsample}}\code{ = 1000000}}{
#'        When plotting methylation value distributions, this threshold specifies the number of observations drawn per
#'        group. Distributions are estimated and plotted based on these random subsamples. This approach can
#'        significantly reduce the memory requirements of the preprocessing and exploratory analysis modules, where
#'        methylation value distributions are plotted. Setting this to \code{0} disables subsampling. More information
#'        is presented the Details section of \code{\link{rnb.step.betadistribution}}}.
#'   \item{\bold{\code{enforce.memory.management}}\code{ = FALSE}}{
#'        Flag indicating whether in some places of the code memory management should actively being enforced in order to
#'        achieve a better memory profile. I.e. garbage collection, variable removal is conducted actively.
#'        May slow down analysis.}
#'   \item{\bold{\code{enforce.destroy.disk.dumps}}\code{ = FALSE}}{
#'        Flag indicating whether disked dumped big matrices (see \code{disk.dump.big.matrices} option) should actively
#'        be deleted when RnBSets are modified. You should switch it to \code{TRUE} when \code{disk.dump.big.matrices}
#' 		  is \code{TRUE} and the amount of hard drive space is also limited.}
#' }
#'
#' @examples
#' \dontrun{
#' str(rnb.options())
#' rnb.getOption("filtering.greedycut")
#' }
#'
#' @author Yassen Assenov
#' @export
rnb.options <- function(...) {
	optlist <- list(...)
	if (length(optlist) == 0) {
		## Return all options and their current values
		return(.rnb.options[["current"]])
	}
	if (is.null(names(optlist))) {
		## Extract the value of one or more options
		for (oname in optlist) {
			if (!(is.character(oname) && length(oname) == 1 && (!is.na(oname)))) {
				stop("invalid option name specified")
			}
			oresult <- rnb.get.option(oname[1])
			if (oresult != "") {
				.rnb.options[["previous"]] <- list()
				stop(oresult)
			}
		}
		return.invisible <- FALSE
		names(optlist) <- unlist(optlist, use.names = FALSE)
	} else {
		return.invisible <- TRUE
		## Set the values of one or more options
		for (i in 1:length(optlist)) {
			oname <- names(optlist)[i]
			if (oname == "") {
				return.invisible <- FALSE
				oname <- optlist[[i]]
				if (!(is.character(oname) && length(oname) == 1 && (!is.na(oname)))) {
					stop("invalid option name specified")
				}
				names(optlist)[i] <- oname
				oresult <- rnb.get.option(oname)
			} else {
				oresult <- rnb.get.option(oname, optlist[[oname]], setvalue = TRUE)
			}
			if (oresult != "") {
				## An error was reached, rollback
				if (length(.rnb.options[["previous"]]) != 0) {
					.rnb.options[["current"]][names(.rnb.options[["previous"]])] <- .rnb.options[["previous"]]
					.rnb.options[["previous"]] <- list()
				}
				stop(oresult)
			}
		}
	}
	result <- .rnb.options[["previous"]]
	.rnb.options[["previous"]] <- list()
	if (return.invisible) {
		return(invisible(result))
	}
	return(result)
}

########################################################################################################################

#' @rdname rnb.options
#' @export
rnb.getOption <- function(x) {
	if (!(is.character(x) && length(x) == 1 && (!is.na(x)))) {
		stop("invalid option name specified")
	}
	if (!(x[1] %in% names(.rnb.options[["current"]]))) {
		stop(paste(x[1], "is invalid option"))
	}
	return(.rnb.options[["current"]][[x[1]]])
}

########################################################################################################################

#' rnb.options2xml
#'
#' Exports all option values to an XML document.
#'
#' @param pretty Flag indicating if the document should be formatted to be easily readable. For example, if this is set
#'               to \code{TRUE} (default), every element is located on separate line. Formatting does not affect the
#'               validity of the generated XML tree.
#' @return XML document in the form of a \code{character} that encodes all options and their current values.
#'
#' @examples
#' \dontrun{
#' cat(rnb.options2xml(), file = "rnbeads_options.xml")
#' }
#' @author Yassen Assenov
#' @export
rnb.options2xml<- function(pretty = TRUE) {
	if (!parameter.is.flag(pretty)) {
		stop("invalid value for pretty; expected TRUE or FALSE")
	}
	o.value2character <- function(o.value) {
		if (is.null(o.value)) {
			return(' null="null">')
		}
		result <- as.character(o.value)
		if (is.logical(o.value)) {
			result <- tolower(result)
		}
		if (length(result) > 1) {
			result <- paste0(result, collapse = ",")
		}
		result <- gsub("&", "&amp;", result, fixed = TRUE)
		result <- gsub("<", "&lt;", result, fixed = TRUE)
		result <- gsub(">", "&gt;", result, fixed = TRUE)
		nns <- names(o.value)
		names.attr <- ""
		if (!is.null(nns)){
			if (!any(is.na(nns))){
				names.attr <- TRUE
				names.attr <- paste0(' names="',paste(nns,collapse=","),'"')
				names.attr <- gsub("&", "&amp;", names.attr, fixed = TRUE)
				names.attr <- gsub("<", "&lt;", names.attr, fixed = TRUE)
				names.attr <- gsub(">", "&gt;", names.attr, fixed = TRUE)
			}
		}
		if (nchar(names.attr)>0){
			paste0(names.attr, ">", result)
		} else {
			paste0(">", result)
		}
	}
	NL <- ifelse(pretty, "\n", "")
	TAB <- ifelse(pretty, "\t", "")
	o.values <- sapply(rnb.options(), o.value2character)
	o.values <- paste0(TAB, "<", names(o.values), o.values, "</", names(o.values), ">", collapse = NL)
	paste0("<rnb.xml>", NL, o.values, NL, "</rnb.xml>", NL)
}

#' rnb.performance.profile
#'
#' Enables one of the pre-installed option profiles
#'
#' @param profile A character singleton specifying the option profile. Can attain values
#' 				  "450k_minimal", "450k_moderate", "450k_full", "bs_minimal", "bs_moderate",
#' 				  "bs_full".
#' 
#' @author Pavlo Lutsik
#' @export
rnb.performance.profile<-function(profile){

	if(missing(profile)){
		stop("One of the pre-installed profiles has to be specified")
	}
	if(!profile %in% PROFILE.FILES){
		stop("Unknown profile supplied")
	}

	fname <- system.file(paste0("extdata/optionProfiles/", profile, ".xml"), package="RnBeads")
	invisible(rnb.xml2options(fname))
}
