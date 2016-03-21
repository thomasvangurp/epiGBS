########################################################################################################################
## filtering.R
## created: 2012-06-04
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the SNP removal and NA removal steps of the filtering steps of the preprocessing module.
########################################################################################################################

## G L O B A L S #######################################################################################################

HIGH.COVER.OUTLIER.QUANTILE <- 0.95 #used for high coverage outlier filtering
HIGH.COVER.OUTLIER.FACTOR <- 50 #used for high coverage outlier filtering

## F U N C T I O N S ###################################################################################################

## rnb.filtering.results
##
## Constructs the list with common elements returned by execute functions in the filtering steps.
##
## @param rnb.set  Methylation dataset before filtering.
## @param filtered \code{logical} vector signifying which of the sites in \code{rnb.set} will be removed.
## @return List of three or four elements:
##         \describe{
##           \item{\code{"dataset.before"}}{Copy of \code{rnb.set}.}
##           \item{\code{"dataset"}}{The (possibly modified) dataset after performing site removal.}
##           \item{\code{"filtered"}}{Indices of sites in \code{rnb.set} that were removed.}
##           \item{\code{"betas"}}{\code{matrix} of the beta values for the sites that were removed from the dataset.
##                This element is added to the list only when \code{filtered} is non-empty.}
##         }
## @author Yassen Assenov
rnb.filtering.results <- function(rnb.set, filtered = NULL) {
	if (is.null(filtered)) {
		filtered <- integer(0)
	} else {
		filtered <- which(filtered)
	}
	if (length(filtered) != 0) {
		return(list(dataset.before = rnb.set, dataset = remove.sites(rnb.set, filtered, verbose = TRUE),
				filtered = filtered, betas = meth(rnb.set)[filtered, ]))
	}
	rnb.cleanMem()
	return(list(dataset.before = rnb.set, dataset = rnb.set, filtered = filtered))
}

########################################################################################################################

## validate.stats
##
## Validates that the given statistics are (partially) produced by \code{\link{rnb.filtering.results}}.
##
## @param stats Filtering result statistics.
## @ author Yassen Assenov
validate.stats <- function(stats) {
	if (!is.list(stats) && all(c("dataset.before", "dataset", "filtered") %in% names(stats))) {
		stop("invalid value for stats")
	}
	if (!inherits(stats$dataset.before, "RnBSet")) {
		stop("invalid value for stats$dataset.before; expected methylation dataset")
	}
	if (!inherits(stats$dataset, "RnBSet")) {
		stop("invalid value for stats$dataset; expected methylation dataset")
	}
	samples.before <- colnames(meth(stats$dataset.before))
	samples.after <- colnames(meth(stats$dataset))
	if (!identical(samples.after, samples.before)) {
		stop("invalid value for stats; incompatible dataset.before and dataset")
	}
	site.count.before <- nrow(meth(stats$dataset.before))
	site.count.after <- nrow(meth(stats$dataset))
	if (!(site.count.after <= site.count.before)) {
		stop("invalid value for stats; incompatible dataset.before and dataset")
	}
	filtered <- stats$filtered
	if (!(is.integer(filtered) && all(!is.na(filtered)) && all(1L <= filtered) && anyDuplicated(filtered) == 0)) {
		stop("invalid value for stats$filtered")
	}
	if (!(all(filtered <= site.count.before) && site.count.before == length(filtered) + site.count.after)) {
		stop("invalid value for stats; incompatible dataset.before, dataset and filtered")
	}
}

########################################################################################################################

## rnb.save.removed.sites
##
## Creates a table with information on the specified removed sites or probes and saves it to a file.
##
## @param anno.table     Table of site or probe annotation.
## @param report         Report to link to the generated file with an annotation table.
## @param fname          Name of the file to store the created annotation table.
## @param p.columns      Column names from the complete annotation table to include in the created one.
## @param p.columns.file Column names to use for the created table.
## @param filtered       Indices to include in the created annotation table.
## @param ...            Any additional annotation columns, specified as named parameters.
## @return Name of the generated file, as a local path with respect to the given report.
##
## @author Yassen Assenov
rnb.save.removed.sites <- function(anno.table, report, fname, p.columns = c("ID", "Chromosome", "Start", "End"),
	p.columns.file = p.columns, ...) {

	p.columns <- intersect(p.columns, colnames(anno.table))
	p.infos <- anno.table[, p.columns]

	additional.values <- list(...)
	for (i in names(additional.values)) {
		p.infos[, i] <- additional.values[[i]]
	}
	fname.full <- file.path(rnb.get.directory(report, "data", absolute = TRUE), fname)
	utils::write.csv(p.infos, file = fname.full, row.names = FALSE)
	rnb.status(c("Saved removed sites to", fname.full))
	return(paste(rnb.get.directory(report, "data"), fname, sep = "/"))
}

########################################################################################################################

#' rnb.execute.context.removal
#'
#' Removes all probes that belong to specific context from the given dataset.
#'
#' @param rnb.set  Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param contexts Probe contexts to be filtered out.
#' @return List of three or four elements:
#'         \describe{
#'           \item{\code{"dataset.before"}}{Copy of \code{rnb.set}.}
#'           \item{\code{"dataset"}}{The (possibly modified) \code{RnBeadSet} object after performing the missing
#'                value removal.}
#'           \item{\code{"filtered"}}{\code{integer} vector storing the indices of all removed probes in
#'                \code{dataset.before}.}
#'           \item{\code{"contexts"}}{The value of the parameter \code{contexts}.}
#'         }
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' contexts.to.ignore <- c("CC", "CAG", "CAH")
#' rnb.set.filtered <- rnb.execute.context.removal(rnb.set.example, contexts.to.ignore)$dataset
#' identical(rnb.set.example, rnb.set.filtered) # FALSE
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.context.removal <- function(rnb.set, contexts = rnb.getOption("filtering.context.removal")) {
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	## TODO: Validate contexts

	filtered <- rnb.execute.context.removal.internal(integer(), contexts, annotation(rnb.set, add.names = TRUE))
	return(list(dataset.before = rnb.set, dataset = remove.sites(rnb.set, filtered), filtered = filtered,
			contexts = contexts))
}

rnb.execute.context.removal.internal <- function(sites2ignore, contexts, anno.table) {
	setdiff(which(anno.table[, "Context"] %in% contexts), sites2ignore)
}

########################################################################################################################

#' rnb.section.context.removal
#'
#' Adds a section on context-specific probe removal to the specified report.
#'
#' @param report Report to summarize the outcome of the procedure. This must be an object of type
#'               \code{\linkS4class{Report}}.
#' @param stats  Statistics on context-specific probe filtering, as returned by
#'               \code{\link{rnb.execute.context.removal}}. See the documentation of the function for more details.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.context.removal <- function(report, stats) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	validate.stats(stats)
	## TODO: Validate stats$contexts
	rnb.section.context.removal.internal(report, stats, annotation(stats$dataset.before, add.names = TRUE))
}

rnb.section.context.removal.internal <- function(report, contexts, filtered, anno.table) {

	if (length(contexts) == 0) {
		return(report)
	}
	f.count <- length(filtered)
	txt.title <- "Context-specific Probe Removal"
	txt <- paste("The studied dataset contains",
		ifelse(f.count == 0, "no probes", ifelse(f.count == 1, "1 (one) probe", paste("in total", f.count, "probes"))),
		"of the specified", ifelse(length(contexts) == 1, "context.", "contexts."))
	if (f.count == 0) {
		report <- rnb.add.section(report, txt.title, txt)
		return(report)
	}

	## Construct a table of selected information for all removed probes
	p.columns <- c("ID", "Chromosome", "Start", "End", "Context")
	fname <- "removed_sites_context.csv"
	fname <- rnb.save.removed.sites(anno.table[filtered, ], report, fname, p.columns)
	txt <- c(txt, ifelse(f.count == 1, " This (removed) probe is", " All these (removed) probes are "),
		"available in a <a href=\"", fname, "\">dedicated table</a> accompanying this report.")
	rm(p.columns, fname)

	## Construct a summary table of number of removed probes per context
	probe.infos <- as.character(anno.table[filtered, "Context"])
	removed.summary <- data.frame(
		"Context" = contexts,
		"Probes" = sapply(contexts, function(x) { sum(probe.infos == x) }))

	if (nrow(removed.summary) == 1) {
		report <- rnb.add.section(report, txt.title, txt)
		return(report)
	}
	txt <- c(txt, " The table below summarizes the number of removed probes per context.")
	report <- rnb.add.section(report, txt.title, txt)
	rnb.add.table(report, removed.summary, row.names = FALSE)
	return(report)
}

########################################################################################################################

#' rnb.step.context.removal
#'
#' Performs the procedure for context-specific probe removal (if applicable) to filter the given methylation dataset
#' and adds a corresponding section to the specified report.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return List of up to three elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The (possibly modified) \code{RnBeadSet} object after performing context-specific
#'                probe removal.}
#'           \item{\code{"report"}}{The modified report.}
#'           \item{\code{"filtered"}}{}
#'           \item{\code{"betas"}}{\code{matrix} of the beta values for probes that have undesired contexts. These
#'                values were essentially removed from the dataset. This element is added to the list only when
#'                \code{rnb.set} contains such probes.}
#'         }
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.context.removal <- function(rnb.set, report) {
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	result <- rnb.step.context.removal.internal(integer(), report, annotation(rnb.set, add.names = TRUE))
	return(list(dataset = remove.sites(rnb.set, result$filtered), report = result$report))
}

rnb.step.context.removal.internal <- function(sites2ignore, report, anno.table) {
	logger.start("Probe Context Removal")
	contexts <- rnb.getOption("filtering.context.removal")
	filtered <- rnb.execute.context.removal.internal(sites2ignore, contexts, anno.table)
	logger.status(c("Removed", length(filtered), "probe(s) having not acceptable context"))
	report <- rnb.section.context.removal.internal(report, contexts, filtered, anno.table)
	logger.status("Added a corresponding section to the report")
	logger.completed()
	list(report = report, filtered = filtered)
}

########################################################################################################################

#' rnb.execute.snp.removal
#'
#' Removes all probes overlapping with single nucleotide polymorphisms (SNPs) from the given dataset.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param snp     Criterion for the removal of sites or probes based on overlap with SNPs. Possible values are
#'                \code{"no"}, \code{"3"}, \code{"5"}, \code{"any"} or \code{"yes"}. See the documentation of
#'                \code{\link{rnb.options}} for a detailed explanation of the procedures these values encode.
#' @return \code{list} of four elements:
#'         \describe{
#'           \item{\code{"dataset.before"}}{Copy of \code{rnb.set}.}
#'           \item{\code{"dataset"}}{The (possibly) modified dataset object after removing probes that overlap
#'                with SNPs.}
#'           \item{\code{"filtered"}}{\code{integer} vector storing the indices (in beta matrix of the unfiltered
#'                dataset) of all removed sites or probes.}
#'           \item{\code{"snp"}}{The value of the \code{snp} parameter.}
#'         }
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.set.filtered <- rnb.execute.snp.removal(rnb.set.example, "any")$dataset
#' identical(meth(rnb.set.example), meth(rnb.set.filtered)) # FALSE
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.snp.removal <- function(rnb.set, snp = rnb.getOption("filtering.snp")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if ((is.double(snp) || is.integer(snp)) && length(snp) == 1 && (as.character(snp) %in% FILTERING.SNPS)) {
		snp <- as.character(snp)
	}
	if (!(is.character(snp) && length(snp) == 1 && (snp %in% FILTERING.SNPS))) {
		msg <- paste0('"', FILTERING.SNPS, '"', collapse = ", ")
		stop(paste("invalid value for snp; expected one of", msg))
	}
	filtered <- rnb.execute.snp.removal.internal(integer(), snp[1], annotation(rnb.set))
	
	list(dataset.before = rnb.set, dataset = remove.sites(rnb.set, filtered), filtered = filtered, snp = snp)
}

rnb.execute.snp.removal.internal <- function(sites2ignore, snp, anno.table) {
	if (snp == "no") {
		filtered <- NULL
	} else {
		snp.overlap.column <- paste("SNPs", ifelse(snp %in% c("3", "5"), snp, "Full"))
		if (snp.overlap.column %in% colnames(anno.table)) {
			filtered <- setdiff(which(anno.table[, snp.overlap.column] > 0), sites2ignore)
		} else {
			stop("no SNP annotation")
		}
	}
	filtered
}

########################################################################################################################

#' rnb.section.snp.removal
#'
#' Adds a section on removing SNP-overlapping sites or probes to the specified report.
#'
#' @param report  Report to contain the new section. This must be an object of type \code{\linkS4class{Report}}.
#' @param rnb.set Methylation dataset before filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param stats   Statistics on SNP-overlapping probe filtering, as returned by
#'                \code{\link{rnb.execute.snp.removal}}. See the documentation of the function for more details.
#' @return The possibly modified report.
#'
#' @seealso \code{\link{rnb.execute.snp.removal}}, \code{\link{rnb.step.snp.removal}}, \code{\link{rnb.run.filtering}}
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.snp.removal <- function(report, rnb.set, stats) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	validate.stats(stats)
	if (!inherits(stats$dataset, "RnBSet")) {
		stop("invalid value for stats$dataset; expected RnBset")
	}
	anno.table <- annotation(stats$dataset.before, add.names = inherits(stats$dataset, "RnBeadSet"))
	filtered <- stats$filtered
	if (!(is.integer(filtered) && anyDuplicated(filtered) == 0 && 1 <= min(filtered) &&
		  	max(filtered) <= nrow(anno.table))) {
		stop("invalid value for stats$filtered")
	}
	snp <- stats$snp
	if ((is.double(snp) || is.integer(snp)) && length(snp) == 1 && (as.character(snp) %in% FILTERING.SNPS)) {
		snp <- as.character(snp)
	}
	if (!(is.character(snp) && length(snp) == 1 && (snp %in% FILTERING.SNPS))) {
		msg <- paste0('"', FILTERING.SNPS, '"', collapse = ", ")
		stop(paste("invalid value for stats$snp; expected one of", msg))
	}
	rnb.section.snp.removal.internal(report, class(stats$dataset), filtered, anno.table, snp)
}

rnb.section.snp.removal.internal <- function(report, dataset.class, filtered, anno.table, snp) {
	txt.site <- rnb.get.row.token(dataset.class)
	txt.sites <- rnb.get.row.token(dataset.class, plural = TRUE)
	txt.title <- paste("Removal of SNP-enriched", capitalize(txt.sites))
	if (snp == "no") {
		return(report)
	}
	if (is.null(filtered)) {
		txt <- "SNP-based filtering was not performed because the site annotation does not include SNP information."
		report <- rnb.add.section(report, txt.title, txt)
		return(report)
	}
	N <- length(filtered)

	if (dataset.class == "RnBeadSet") {
		if (N == 0) {
			txt <- "No probes were found for which"
		} else if (N == 1) {
			txt <- "<b>One</b> probe was removed because"
		} else {
			txt <- paste0("<b>", N, "</b> probes were removed because")
		}
		if (snp %in% c("3", "5")) {
			txt <- paste(txt, "the last", snp, "bases of")
		}
		txt <- paste(txt, "their sequences overlap with SNPs.")
	} else { # dataset.class == "RnBiseqSet"
		if (N == 0) {
			txt <- "No sites were found that overlap"
		} else if (N == 1) {
			txt <- "<b>One</b> site was removed because it overlaps"
		} else {
			txt <- paste0("<b>", N, "</b> sites were removed because they overlap")
		}
		txt <- paste(txt, "with SNPs.")
	}
	
	if (N != 0) {
		snp.overlap.column <- paste("SNPs", ifelse(snp %in% c("3", "5"), snp, "Full"))
		p.columns <- c("ID", "Chromosome", "Start", "End", snp.overlap.column)
		names(p.columns) <- c(p.columns[1:(length(p.columns) - 1)], "SNPs")
		fname <- "removed_sites_snp.csv"
		fname <- rnb.save.removed.sites(anno.table[filtered, ], report, fname, p.columns, names(p.columns))
		txt <- paste(txt, "The", ifelse(N == 1, paste("removed", txt.site), paste("list of removed", txt.sites)),
					 ' is available in a <a href="', fname, '">dedicated table</a> accompanying this report.')
	}
	report <- rnb.add.section(report, txt.title, txt)
	return(report)
}

########################################################################################################################

#' rnb.step.snp.removal
#'
#' Performs the procedure for removal of sites or probes (if any) that overlap with many SNPs to filter the given
#' methylation dataset and adds a corresponding section to the specified report.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return List of two elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The (possibly modified) dataset after removing sites that overlap with SNPs.}
#'           \item{\code{"report"}}{The modified report.}
#'         }
#'
#' @seealso \code{\link{rnb.execute.snp.removal}}, \code{\link{rnb.section.snp.removal}}, \code{\link{rnb.run.filtering}}
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.snp.removal <- function(rnb.set, report) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	result <- rnb.step.snp.removal.internal(class(rnb.set), integer(), report,
		annotation(rnb.set, add.names = TRUE))
	return(list(dataset = remove.sites(rnb.set, result$filtered), report = result$report))
}

rnb.step.snp.removal.internal <- function(dataset.class, sites2ignore, report, anno.table) {
	snp <- rnb.getOption("filtering.snp")
	logger.start("Removal of SNP-enriched Sites")
	filtered <- tryCatch(
		rnb.execute.snp.removal.internal(sites2ignore, snp, anno.table), error = function(x) { NULL })
	if (is.null(filtered)) {
		filtered <- integer()
		report <- rnb.section.snp.removal.internal(report, dataset.class, NULL, anno.table, snp)
	} else {
		msg <- paste("Removed", length(filtered), ifelse(length(filtered) == 1, "site", "sites"))
		logger.status(paste0(msg, ' using SNP criterion "', snp, '"'))
		report <- rnb.section.snp.removal.internal(report, dataset.class, filtered, anno.table, snp)
	}
	logger.status("Added a corresponding section to the report")
	logger.completed()
	list(report = report, filtered = filtered)
}

########################################################################################################################

#' rnb.execute.sex.removal
#'
#' Removes all sites in sex chromosomes from the given dataset.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @return List of three elements:
#'         \describe{
#'           \item{\code{"dataset.before"}}{Copy of \code{rnb.set}.}
#'           \item{\code{"dataset"}}{The (possibly) modified dataset after retaining sites on autosomes only.}
#'           \item{\code{"filtered"}}{\code{integer} vector storing the indices (in beta matrix of the unfiltered
#'                dataset) of all removed probes.}
#'         }
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.example)$dataset
#' identical(meth(rnb.set.example), meth(rnb.set.filtered)) # FALSE
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.sex.removal <- function(rnb.set) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	filtered <- rnb.execute.sex.removal.internal(integer(),
		annotation(rnb.set, add.names = inherits(rnb.set, "RnBeadSet")))

	if (length(filtered) != 0) {
		dataset <- remove.sites(rnb.set, filtered)
	} else {
		dataset <- rnb.set
	}
	return(list(dataset.before = rnb.set, dataset = dataset, filtered = filtered))
}

rnb.execute.sex.removal.internal <- function(sites2ignore, anno.table) {
	setdiff(which(anno.table[, "Chromosome"] %in% c("chrX", "chrY")), sites2ignore)
}

########################################################################################################################

#' rnb.section.sex.removal
#'
#' Adds a section on removing sex chromosome sites to the specified report.
#'
#' @param report Report to contain the new section. This must be an object of type \code{\linkS4class{Report}}.
#' @param stats  Statistics on sex chromosome site filtering, as returned by \code{\link{rnb.execute.sex.removal}}.
#'               See the documentation of the function for more details.
#' @return The modified report.
#'
#' @seealso \code{\link{rnb.execute.sex.removal}}, \code{\link{rnb.step.sex.removal}}, \code{\link{rnb.run.filtering}}
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.sex.removal <- function(report, stats) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	validate.stats(stats)
	rnb.section.sex.removal.internal(report, class(stats$dataset), stats$filtered,
		annotation(stats$dataset.before, add.names = inherits(stats$dataset, "RnBeadSet")))
}

rnb.section.sex.removal.internal <- function(report, dataset.class, filtered, anno.table) {
	txt.site <- rnb.get.row.token(dataset.class)
	txt.sites <- rnb.get.row.token(dataset.class, plural = TRUE)
	pcount <- length(filtered)
	if (pcount == 0) {
		txt <- c("No ", txt.sites, " located on sex chromosomes were found.")
	} else {
		fname <- "removed_sites_sex.csv"
		fname <- rnb.save.removed.sites(anno.table[filtered, ], report, fname)
		txt <- c(pcount, " ", ifelse(pcount == 1, txt.site, txt.sites), " on sex chromosomes ",
			ifelse(pcount == 1, "was", "were"), " removed at this step. The ",
			ifelse(pcount == 1, paste("removed", txt.site), paste("list of removed", txt.sites)),
			" is available in a <a href=\"", fname, "\">dedicated table</a> accompanying this report.")
	}
	report <- rnb.add.section(report, paste("Removal of", capitalize(txt.sites), "on Sex Chromosomes"), txt)
	return(report)
}

########################################################################################################################

#' rnb.step.sex.removal
#'
#' Performs the procedure for removal of sites or probes on sex chromosomes (if any) to filter the given methylation
#' dataset and adds a corresponding section to the specified report.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return List of two elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The (possibly modified) dataset after removing sites on sex chromosomes.}
#'           \item{\code{"report"}}{The modified report.}
#'         }
#'
#' @seealso \code{\link{rnb.execute.sex.removal}}, \code{\link{rnb.section.sex.removal}}, \code{\link{rnb.run.filtering}}
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.sex.removal <- function(rnb.set, report) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	result <- rnb.step.sex.removal.internal(class(rnb.set), integer(), report,
		annotation(rnb.set, add.names = inherits(rnb.set, "RnBeadSet")))
	return(list(dataset = remove.sites(rnb.set, result$filtered), report = result$report))
}

rnb.step.sex.removal.internal <- function(dataset.class, sites2ignore, report, anno.table) {
	logger.start("Removal of Sites on Sex Chromosomes")
	filtered <- rnb.execute.sex.removal.internal(sites2ignore, anno.table)
	logger.status(c("Removed", length(filtered), "site(s) on sex chromosomes"))
	report <- rnb.section.sex.removal.internal(report, dataset.class, filtered, anno.table)
	logger.status("Added a corresponding section to the report")
	logger.completed()
	list(report = report, filtered = filtered)
}

########################################################################################################################
########################################################################################################################

#' rnb.execute.low.coverage.masking
#'
#' Replaces all low coverage sites by \code{NA}.
#'
#' @param rnb.set        Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param covg.threshold Threshold for minimal acceptable coverage, given as a non-negative \code{integer} value. All
#'                       methylation measurements with lower coverage than this threshold are set to \code{NA}. If this
#'                       parameter is \code{0}, calling this method has no effect.
#' @return List of three elements:
#'         \describe{
#'           \item{\code{"dataset.before"}}{Copy of \code{rnb.set}.}
#'           \item{\code{"dataset"}}{The (possibly) modified dataset after retaining sites on autosomes only.}
#'           \item{\code{"mask"}}{A logical matrix of dimension \code{meth(rnb.set,type="sites")} indicating which
#' 				   methylation values have been masked}
#'         }
#'
#' @author Fabian Mueller
#' @export
rnb.execute.low.coverage.masking <- function(rnb.set, covg.threshold = rnb.getOption("filtering.coverage.threshold")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (is.double(covg.threshold) && all(covg.threshold == as.integer(covg.threshold))) {
		covg.threshold <- as.integer(covg.threshold)
	}
	if (!(is.integer(covg.threshold) && length(covg.threshold) == 1 && isTRUE(0 <= covg.threshold))) {
		stop("invalid value for covg.threshold")
	}
	mask <- rnb.execute.low.coverage.masking.internal(rnb.set, integer(), covg.threshold)
	dataset <- rnb.set
	if (any(mask)) {
		dataset@meth.sites[,][mask] <- NA
		dataset <- updateRegionSummaries(dataset)
	}
	list(dataset.before = rnb.set, dataset = dataset, mask = mask)
}

rnb.execute.low.coverage.masking.internal <- function(rnb.set, sites2ignore, covg.threshold) {
	coverage.matrix <- covg(rnb.set)
	if (!(is.matrix(coverage.matrix) && all(dim(coverage.matrix) != 0))) {
		return(NULL)
	}
	mask <- (coverage.matrix < covg.threshold) & (!is.na(coverage.matrix)) & (!is.na(meth(rnb.set)))
	if (length(sites2ignore) != 0) {
		mask[sites2ignore, ] <- FALSE
	}
	return(mask)
}

########################################################################################################################

#' rnb.section.low.coverage.masking
#'
#' Adds a section on masking low coverage sites with NAs to the specified report.
#'
#' @param report Report to contain the new section. This must be an object of type \code{\linkS4class{Report}}.
#' @param stats a list as outputted from \code{rnb.execute.low.coverage.masking}
#' @return The modified report.
#'
#' @seealso \code{\link{rnb.execute.low.coverage.masking}}, \code{\link{rnb.run.filtering}}
#'
#' @author Fabian Mueller
#' @noRd
rnb.section.low.coverage.masking <- function(report, stats,
	covg.threshold = rnb.getOption("filtering.coverage.threshold")) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	## TODO: Validate stats
	if (is.double(covg.threshold) && all(covg.threshold == as.integer(covg.threshold))) {
		covg.threshold <- as.integer(covg.threshold)
	}
	if (!(is.integer(covg.threshold) && length(covg.threshold) == 1 && (!is.na(covg.threshold)) &&
		  0 <= covg.threshold)) {
		stop("invalid value for covg.threshold")
	}
	rnb.section.low.coverage.masking.internal(report, class(stats$dataset), stats$mask,
		annotation(stats$dataset.before, add.names = inherits(stats$dataset, "RnBeadSet")), covg.threshold)
}

rnb.section.low.coverage.masking.internal <- function(report, dataset.class, mask, anno.table, covg.threshold) {
	txt.site <- rnb.get.row.token(dataset.class, )
	txt.sites <- rnb.get.row.token(dataset.class, plural = TRUE)
	num.masked <- colSums(mask)
	num.masked.total <- sum(num.masked)
	if (num.masked.total < 1) {
		txt <- c("No ", txt.sites, " were masked")
	} else {
		masked.table <- data.frame(Sample = names(num.masked), masked = num.masked, check.names = FALSE)
		colnames(masked.table)[2] <- paste("Masked", txt.sites)
		fname <- "masked_sites_coverage.csv"
		fname.full <- file.path(rnb.get.directory(report, "data", absolute = TRUE), fname)
		utils::write.csv(masked.table, file = fname.full, row.names = FALSE)
		rnb.status(c("Saved numbers of masked sites per sample to", fname.full))
		fname.relative <- paste(rnb.get.directory(report, "data"), fname, sep = "/")
		txt <- c("A total of ",num.masked.total," ",txt.sites, " with coverage less than ",covg.threshold," were masked by NA in the methylation table")
		txt <- c(txt, paste(" The numbers of masked",txt.sites,"per sample"))
		txt <- c(txt, " are available in a <a href=\"", fname.relative, "\">dedicated table</a> accompanying this report.")
	}
	report <- rnb.add.section(report, paste("Masking of", capitalize(txt.sites), "with Low Coverage"), txt)
	return(report)
}

########################################################################################################################

#' rnb.step.low.coverage.masking
#'
#' Sets methylation values of sites with low coverage to NA
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return List of up to three elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The (possibly modified) dataset after replacing low coverage sites}
#'           \item{\code{"report"}}{The modified report.}
#'         }
#'
#' @seealso \code{\link{rnb.run.filtering}},  \code{\link{rnb.execute.low.coverage.masking}},  \code{\link{rnb.section.low.coverage.masking}}
#'
#' @author Fabian Mueller
#' @noRd
rnb.step.low.coverage.masking <- function(rnb.set, report, covg.threshold = rnb.getOption("filtering.coverage.threshold")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (is.double(covg.threshold) && all(covg.threshold == as.integer(covg.threshold))) {
		covg.threshold <- as.integer(covg.threshold)
	}
	if (!(is.integer(covg.threshold) && length(covg.threshold) == 1 && (!is.na(covg.threshold)) &&
		  	0 <= covg.threshold)) {
		stop("invalid value for covg.threshold")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	rnb.step.low.coverage.masking.internal(rnb.set, integer(), report,
		annotation(rnb.set, add.names = inherits(rnb.set, "RnBeadSet")), covg.threshold)
}

rnb.step.low.coverage.masking.internal <- function(rnb.set, sites2ignore, report, anno.table, covg.threshold) {
	logger.start("Replacing Low Coverage Sites by NA")
	mask <- rnb.execute.low.coverage.masking.internal(rnb.set, sites2ignore, covg.threshold)
	logger.status(c("Masked ", sum(colSums(mask)), " site(s) based on coverage threshold ", covg.threshold))
	report <- rnb.section.low.coverage.masking.internal(report, class(rnb.set), mask, anno.table, covg.threshold)
	logger.status("Added a corresponding section to the report")
	logger.completed()
	list(report = report, mask = mask)
}

########################################################################################################################

#' rnb.execute.high.coverage.removal
#'
#' Removes methylation sites with a coverage larger than 100 times the 95-percentile of coverage in each sample.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBiseqSet}}.
#' @return \code{list} of two elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The (possibly) modified dataset after retaining sites on autosomes only.}
#'           \item{\code{"filtered"}}{\code{integer} vector storing the indices of all removed sites.}
#'         }
#'
#' @author Fabian Mueller
#' @export
rnb.execute.high.coverage.removal <- function(rnb.set) {
	if (!inherits(rnb.set, "RnBiseqSet")) {
		stop("invalid value for rnb.set")
	}
	filtered <- rnb.execute.high.coverage.removal.internal(rnb.set, integer())
	if (is.null(filtered)) {
		stop("invalid value for rnb.set; no coverage information present")
	}
	if (length(filtered) != 0) {
		dataset <- remove.sites(rnb.set, filtered)
	} else {
		dataset <- rnb.set
	}
	return(list(dataset.before = rnb.set, dataset = dataset, filtered = filtered))
}

rnb.execute.high.coverage.removal.internal <- function(rnb.set, sites2ignore) {
	cover <- covg(rnb.set)
	if (is.null(cover)) {
		return(NULL)
	}
	cover[cover<1] <- NA
	filtered <- matrix(FALSE, nrow = nrow(cover), ncol = ncol(cover)) # allocate indication matrix
	for (i in 1:ncol(cover)) {
		qqs <- ceiling(quantile(cover[, i], probs = HIGH.COVER.OUTLIER.QUANTILE, na.rm = TRUE)) * HIGH.COVER.OUTLIER.FACTOR
		filtered[, i] <- (cover[, i] > qqs)
	}
	filtered <- which(apply(filtered, 1, any, na.rm = TRUE))
	setdiff(filtered, sites2ignore)
}

########################################################################################################################

#' rnb.section.high.coverage.removal
#'
#' Adds a section on masking low coverage sites with NAs to the specified report.
#'
#' @param report Report to contain the new section. This must be an object of type \code{\linkS4class{Report}}.
#' @param masking.result a list as outputted from \code{rnb.execute.high.coverage.removal}
#' @return The modified report.
#'
#' @seealso \code{\link{rnb.execute.high.coverage.removal}}, \code{\link{rnb.run.filtering}}
#'
#' @author Fabian Mueller
#' @noRd
rnb.section.high.coverage.removal <- function(report, stats) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	validate.stats(stats)
	rnb.section.high.coverage.removal.internal(report, class(stats$dataset), stats$filtered, annotation(stats$dataset.before))
}

rnb.section.high.coverage.removal.internal <- function(report, dataset.class, filtered, anno.table) {
	txt.site <- rnb.get.row.token(dataset.class)
	txt.sites <- rnb.get.row.token(dataset.class, plural = TRUE)
	if (is.null(filtered)) {
		txt <- c("No high coverage outlier ", txt.sites, " were found.")
	} else {
		fname <- "removed_sites_high_coverage.csv"
		fname <- rnb.save.removed.sites(anno.table[filtered, ], report, fname)

		txt <- length(filtered)
		txt <- c(ifelse(txt == 1, paste("One", txt.site), paste(txt, txt.sites)), " was detected as a high coverage ",
			"outlier in at least one sample and removed at this step. An outlier site is defined as one whose ",
			"coverage exceeds ", HIGH.COVER.OUTLIER.FACTOR, " times the ", HIGH.COVER.OUTLIER.QUANTILE, "-quantile ",
			"of coverage values in its sample. The list of removed ", txt.sites, " is available in a <a href=\"", fname,
			"\">dedicated table</a> accompanying this report.")
	}

	report <- rnb.add.section(report, "Removal of High Coverage Outlier Sites", txt)
	return(report)
}

########################################################################################################################

#' rnb.step.high.coverage.removal
#'
#' Removes methylation sites with a coverage larger than 100 times the 95-percentile of coverage in each sample.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBiseqSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return List of up to three elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The (possibly modified) dataset after removing high coverage outlier sites}
#'           \item{\code{"report"}}{The modified report.}
#'         }
#'
#' @seealso \code{\link{rnb.execute.high.coverage.removal}}, \code{\link{rnb.section.high.coverage.removal}}
#'
#' @author Fabian Mueller
#' @noRd
rnb.step.high.coverage.removal <- function(rnb.set, report) {
	if (!inherits(rnb.set, "RnBiseqSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	result <- rnb.step.high.coverage.removal.internal(rnb.set, integer(), report, annotation(rnb.set))
	return(list(dataset = remove.sites(rnb.set, result$filtered), report = result$report))
}

rnb.step.high.coverage.removal.internal <- function(rnb.set, sites2ignore, report, anno.table) {
	logger.start("Removal of High Coverage (Outlier) Sites")
	filtered <- rnb.execute.high.coverage.removal.internal(rnb.set, sites2ignore)
	if (is.null(filtered)) {
		filtered <- integer()
		report <- rnb.section.high.coverage.removal.internal(report, class(rnb.set), NULL, anno.table)
		logger.warning("No coverage information present")
	} else {
		logger.status(c("Removed", length(filtered), "high coverage outlier sites"))
		report <- rnb.section.high.coverage.removal.internal(report, class(rnb.set), filtered, anno.table)
		logger.status("Added a corresponding section to the report")
	}
	logger.completed()
	list(report = report, filtered = filtered)
}

########################################################################################################################

#' rnb.execute.na.removal
#'
#' Removes all probes with missing value (if such exists) from the given dataset.
#'
#' @param rnb.set   Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param threshold Maximum quantile of \code{NA}s allowed per site. This must be a value between 0 and 1.
#' @return List of four or five elements:
#'         \describe{
#'           \item{\code{"dataset.before"}}{Copy of \code{rnb.set}.}
#'           \item{\code{"dataset"}}{The (possibly modified) dataset after performing the missing value removal.}
#'           \item{\code{"filtered"}}{\code{integer} vector storing the indices (in beta matrix of the unfiltered
#'                dataset) of all removed sites.}
#' 			 \item{\code{"threshold"}}{Copy of \code{threshold}.}
#'         }
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.set.filtered <- rnb.execute.na.removal(rnb.set.example, 0)$dataset
#' identical(meth(rnb.set.example), meth(rnb.set.filtered)) # TRUE
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.na.removal <- function(rnb.set, threshold = rnb.getOption("filtering.missing.value.quantile")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!(is.double(threshold) && length(threshold) == 1 && (!is.na(threshold)))) {
		stop("invalid value for threshold")
	}
	if (!(0 <= threshold && threshold <= 1)) {
		stop("invalid value for threshold; expected a value between 0 and 1")
	}
	filtered <- which(rowMeans(is.na(meth(rnb.set))) > threshold)
	list(dataset.before = rnb.set, dataset = remove.sites(rnb.set, filtered), filtered = filtered,
		threshold = threshold)
}

rnb.execute.na.removal.internal <- function(mm, sites2ignore, threshold) {
	setdiff(which(rowMeans(is.na(mm)) > threshold), sites2ignore)
}

########################################################################################################################

#' rnb.section.na.removal
#'
#' Adds a section on removing sites or probes with missing values to the specified report.
#'
#' @param report Report to summarize the outcome of the procedure. This must be an object of type
#'               \code{\linkS4class{Report}}.
#' @param stats  Statistics on filtering based on missing values, as returned by \code{\link{rnb.execute.na.removal}}.
#'               See the documentation of the function for more details.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.na.removal <- function(report, stats) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	validate.stats(stats)
	if (!(is.double(stats$threshold) && length(stats$threshold) == 1 && (!is.na(stats$threshold)))) {
		stop("invalid value for stats$threshold")
	}
	if (!(0 <= stats$threshold && stats$threshold <= 1)) {
		stop("invalid value for stats$threshold; expected a value between 0 and 1")
	}
	rnb.section.na.removal.internal(report, class(stats$dataset), meth(stats$dataset.before),
		stats$filtered, stats$threshold, annotation(stats$dataset.before, add.names = TRUE))
}

rnb.section.na.removal.internal <- function(report, dataset.class, mm, filtered, threshold, anno.table) {
	txt.site <- rnb.get.row.token(dataset.class)
	txt.sites <- rnb.get.row.token(dataset.class, plural = TRUE)
	threshold.abs <- as.integer(floor(threshold * ncol(mm)))
	txt.title <- paste("Removal of", capitalize(txt.sites), "with (Many) Missing Values")
	na.counts <- as.integer(rowSums(is.na(mm)))
	removed.count <- length(filtered)
	if (removed.count == 0) {
		txt <- "No sites with too many missing values were found in the methylation table."
		report <- rnb.add.section(report, txt.title, txt)
		return(report)
	}

	## Create a table of sites that contain missing values
	na.c <- na.counts[filtered]
	fname <- "removed_sites_na.csv"
	fname <- rnb.save.removed.sites(anno.table[filtered, ], report, fname, `Number of Missing Values` = na.c)
	txt <- c(removed.count, " ", ifelse(removed.count == 1, paste(txt.site, "was"), paste(txt.sites, "were")),
		" removed because ", ifelse(removed.count == 1, "it contains", "they contain"), " more than ", threshold.abs,
		" missing values in the methylation table. This threshold corresponds to ", round(threshold * 100, 1),
		"% of all samples. The total number of missing values in the methylation table before this filtering step was ",
		sum(na.counts), ". A <a href=\"", fname, "\">dedicated table of all removed ", txt.sites,
		"</a> is attached to this report.")
	report <- rnb.add.section(report, txt.title, txt)
	rm(na.c, fname, txt)

	if (removed.count != 1) {
		## Create a histogram of number of NAs per site
		values <- list("all" = list(values = na.counts))
		if (sum(na.counts == 0) > 1/3 * length(na.counts)) {
			values[["positive"]] <- list(values = na.counts[filtered])
		}
		for (i in names(values)) {
			values[[i]][["fname"]] <- paste("histogram_na_counts", i, sep = "_")
		}
		report.plots <- lapply(values, function(x) {
				dframe <- data.frame(x = x$values)
				binwidth <- range(x$values)
				binwidth <- max((binwidth[2] - binwidth[1]) / 40, 1)
				rplot <- createReportPlot(x$fname, report, width = 5, height = 5)
				pp <- ggplot(dframe, aes_string(x = "x")) + labs(x = "Number of missing values", y = "Frequency") +
					geom_histogram(aes_string(y = "..count.."), binwidth = binwidth) 
				if (0 < threshold && threshold < 1) {
					pp <- pp + geom_vline(xintercept = threshold * ncol(mm), linetype = "dotted")
				}
				print(pp)
				return(off(rplot))
			})
		setting.names <- list()
		if (length(report.plots) > 1) {
			ename <- paste(capitalize(txt.sites), "to include")
			setting.names[[ename]] <- c("all" = "all", "positive" = paste(txt.sites, "with missing values"))
		}
		txt <- c("The figure below shows the distribution of missing values per ", txt.site, ".")
		rnb.add.paragraph(report, txt)
		txt <- c("Histogram of number of ", txt.sites, " that contain missing values.")
		if (0 < threshold && threshold < 1) {
			txt <- c(txt, " The vertical line, if visible, denotes the applied threshold.") 
		}
		report <- rnb.add.figure(report, txt, report.plots, setting.names)
	}

	return(report)
}

########################################################################################################################

#' rnb.step.na.removal
#'
#' Performs the procedure for missing value removal (if applicable) to filter the given methylation dataset and adds a
#' corresponding section to the provided report.
#'
#' @param rnb.set   Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report    Report to summarize the outcome of this procedure. This must be an object of type
#'                  \code{\linkS4class{Report}}.
#' @param threshold Maximum quantile of \code{NA}s allowed per site. This must be a value between 0 and 1.
#' @return List up to three elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The (possibly modified) dataset after performing the missing value removal.}
#'           \item{\code{"report"}}{The modified report.}
#'           \item{\code{"betas"}}{\code{matrix} of the beta values for sites that contain missing values. These values
#'                were essentially removed from the dataset. This element is added to the list only when \code{rnb.set}
#'                contains such sites.}
#'         }
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.na.removal <- function(rnb.set, report, threshold = 0) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	if (!(is.double(threshold) && length(threshold) == 1 && (!is.na(threshold)))) {
		stop("invalid value for threshold")
	}
	if (!(0 <= threshold && threshold <= 1)) {
		stop("invalid value for threshold; expected a value between 0 and 1")
	}
	result <- rnb.step.na.removal.internal(class(rnb.set), meth(rnb.set), report, annotation(rnb.set, add.names = TRUE),
		threshold)
	return(list(dataset = remove.sites(rnb.set, result$filtered), report = result$report))
}

rnb.step.na.removal.internal <- function(dataset.class, mm, sites2ignore, report, anno.table, threshold) {
	logger.start("Missing Value Removal")
	logger.status(c("Using a sample quantile threshold of", threshold))
	filtered <- rnb.execute.na.removal.internal(mm, sites2ignore, threshold)
	logger.status(c("Removed", length(filtered), "site(s) with too many missing values"))
	report <- rnb.section.na.removal.internal(report, dataset.class, mm, filtered, threshold, anno.table)
	logger.status("Added a corresponding section to the report")
	logger.completed()
	list(report = report, filtered = filtered)
}

########################################################################################################################

#' rnb.execute.variability.removal
#'
#' Removes all sites or probes with low variability from the given dataset.
#'
#' @param rnb.set       Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param min.deviation Threshold for standard deviation per site. This must be a scalar between 0 and 1. All sites, for
#'                      which the standard deviation of methylation values (for all samples in \code{rnb.set}) is lower
#'                      than this threshold, will be filtered out.
#' @return List of four elements:
#'         \describe{
#'           \item{\code{"dataset.before"}}{Copy of \code{rnb.set}.}
#'           \item{\code{"dataset"}}{The (possibly modified) dataset after removing sites with low variability.}
#'           \item{\code{"filtered"}}{\code{integer} vector storing the indices (in beta matrix of the unfiltered
#'                dataset) of all removed sites.}
#' 			 \item{\code{"threshold"}}{The value of the given parameter \code{min.deviation}.}
#'         }
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.set.filtered <- rnb.execute.variability.removal(rnb.set.example, 0.01)
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.variability.removal <- function(rnb.set, min.deviation = rnb.getOption("filtering.deviation.threshold")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (is.integer(min.deviation)) {
		ovalue <- as.double(min.deviation)
	}
	if (!(is.double(min.deviation) && length(min.deviation) == 1 && (!is.na(min.deviation)))) {
		stop("invalid value for min.deviation")
	}
	min.deviation <- min.deviation[1]
	if (min.deviation < 0 || min.deviation > 1) {
		stop("invalid value for min.deviation; expected a number between 0 and 1")
	}
	filtered <- which(apply(meth(rnb.set), 1, sd, na.rm = TRUE) < min.deviation)
	list(dataset.before = rnb.set, dataset = remove.sites(rnb.set, filtered), filtered = filtered,
		threshold = min.deviation)
}

rnb.execute.variability.removal.internal <- function(mm, sites2ignore, threshold) {
	site.deviations <- apply(mm, 1, sd, na.rm = TRUE)
	filtered <- setdiff(which(site.deviations < threshold), sites2ignore)
	list(deviations = site.deviations, filtered = filtered)
}

########################################################################################################################

#' rnb.section.variability.removal
#'
#' Adds a section on removing low-variable probes to the specified report.
#'
#' @param report        Report to summarize the outcome of the procedure. This must be an object of type
#'                      \code{\linkS4class{Report}}.
#' @param dataset.class Type of the dataset object, e.g. \code{"RnBiseqSet"}.
#' @param deviations    Vector of type \code{double} storing the standard deviations of all sites.
#' @param filtered      Vector of type \code{integer} storing the indices of all removed sites.
#' @param threshold     Threshold for minimal deviation; sites with lower deviations should have been filtered out.
#' @param anno.table    Annotation information for all sites (before filtering) in the form of a \code{data.frame}.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.variability.removal <- function(report, dataset.class, deviations, filtered, threshold, anno.table) {

	txt.site <- rnb.get.row.token(dataset.class)
	txt.sites <- rnb.get.row.token(dataset.class, plural = TRUE)
	txt.title <- paste("Removal of Consistent", capitalize(txt.sites))
	removed.count <- length(filtered)
	if (removed.count == 0) {
		txt <- c("No ", txt.sites, " in the methylation table were found that exhibit standard deviation lower than ",
			threshold, ".")
		report <- rnb.add.section(report, txt.title, txt)
		return(report)
	}

	## Create a table of probes that have low methylation variability
	fname <- "removed_sites_variability.csv"
	fname <- rnb.save.removed.sites(anno.table[filtered, ], report, fname,
		`Standard Deviation` = deviations[filtered])
	plural <- (removed.count != 1)
	txt <- c(removed.count, ifelse(plural, paste(txt.sites, " were"), paste(txt.site, " was")), " removed because ",
		ifelse(plural, "their", "its"), " beta values exhibit standard deviation lower than ", threshold,
		". A <a href=\"", fname, "\">dedicated table of all removed ", txt.sites, "</a> is attached to this report.")
	report <- rnb.add.section(report, txt.title, txt)
	rm(fname, plural, txt)

	## Create a density estimation plot
	dframe <- data.frame(x = deviations)
	xlim <- range(c(range(deviations, na.rm = TRUE), threshold))
	rplot <- createReportPlot("site_deviations", report, width = 6, height = 5)
	pp <- ggplot(dframe, aes_string(x = "x")) + coord_cartesian(xlim = xlim) +
		labs(x = "Standard deviation", y = "Density") + geom_vline(xintercept = threshold, color = "#FF0000") +
		geom_density(kernel = "gaussian", color = "#000080") + theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in"))
	print(pp)
	rplot <- off(rplot)
	txt <- c("Density plot of observed standard deviations of beta values per ", txt.site, ". The red vertical line ",
		"shows the applied threshold for ", txt.site, " variability.")
	report <- rnb.add.figure(report, txt, rplot)

	return(report)
}

########################################################################################################################

#' rnb.step.variability.removal
#'
#' Performs the procedure for removal of probes with low variability (if applicable) to filter the given methylation
#' dataset and adds a corresponding section to the specified report.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return List of up to three elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The (possibly modified) dataset after removing probes with low variability in
#'                methylation.}
#'           \item{\code{"report"}}{The modified report.}
#'           \item{\code{"betas"}}{\code{matrix} of the beta values for sites that exhibit low variability. These values
#'                were essentially removed from the dataset. This element is added to the list only when \code{rnb.set}
#'                contains such sites.}
#'         }
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.variability.removal <- function(rnb.set, report) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	threshold <- rnb.getOption("filtering.deviation.threshold")
	rnb.step.variability.removal.internal(rnb.set, report,
		annotation(rnb.set, add.names = inherits(rnb.set, "RnBeadSet")))
}

rnb.step.variability.removal.internal <- function(dataset.class, mm, sites2ignore, report, anno.table, threshold) {
	logger.start("Site Removal Based on Standard Deviation")
	logger.info(c("Using standard deviation threshold of", threshold))
	result <- rnb.execute.variability.removal.internal(mm, sites2ignore, threshold)
	logger.status(c("Removed", length(result$filtered), "site(s) with variance lower than the threshold"))
	report <- rnb.section.variability.removal(report, dataset.class, result$deviations, result$filtered,
		threshold, anno.table)
	logger.status("Added a corresponding section to the report")
	list(report = report, filtered = result$filtered)
}

## E N D ###############################################################################################################
