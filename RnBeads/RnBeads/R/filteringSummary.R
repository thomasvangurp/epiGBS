########################################################################################################################
## filteringSummary.R
## created: 2013-12-12
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the summary section after all filtering steps have been performed.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

## rnb.get.filtered.sites.samples
##
## Validates the given datasets before and after filtering and extracts the lists of filtered sites and samples.
##
## @param old.set Methylation dataset before filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
## @param new.set Methylation dataset after filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
## @return \code{list} of three elements: \code{"mm"}, \code{"samples"} and \code{"sites"}.
##
## @author Yassen Assenov
rnb.get.filtered.sites.samples <- function(old.set, new.set) {

	## Validate that new.set is a subset of old.set
	mm.old <- meth(old.set, row.names = TRUE)
	mm.new <- meth(new.set, row.names = TRUE)
	sites.old <- rownames(mm.old)
	sites.new <- rownames(mm.new)
	samples.old <- colnames(mm.old)
	samples.new <- colnames(mm.new)
	rm(mm.new)
	if (length(setdiff(sites.new, sites.old)) != 0) {
		stop("inconsistent sites in old.set and new.set")
	}
	if (length(setdiff(samples.new, samples.old)) != 0) {
		stop("inconsistent samples in old.set and new.set")
	}

	removed.samples <- which(!(samples.old %in% samples.new))
	removed.sites <- which(!(sites.old %in% sites.new))
	return(list(mm = mm.old, samples = removed.samples, sites = removed.sites))
}

########################################################################################################################

#' rnb.execute.filter.summary
#'
#' Calculates a table summarizing the effect of the applied filtering procedures.
#'
#' @param old.set Methylation dataset before filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param new.set Methylation dataset after filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @return \code{matrix} summarizing the number of removed and retained sites, samples, and (optionally) reliable and
#'         unreliable measurements.
#'
#' @details
#' This function expects that the sites and samples in \code{new.set} are subsets of the sites and samples in
#' \code{old.set}, respectively. If this is not the case, it exists with an error.
#'
#' @seealso \code{\link{rnb.run.preprocessing}} for running the whole preprocessing module
#'
#' @author Yassen Assenov
#' @export
rnb.execute.filter.summary <- function(old.set, new.set) {
	if (!inherits(old.set, "RnBSet")) {
		stop("invalid value for old.set")
	}
	if (!inherits(new.set, "RnBSet")) {
		stop("invalid value for new.set")
	}

	removed <- rnb.get.filtered.sites.samples(old.set, new.set)
	relm <- rnb.get.reliability.matrix(old.set)
	rnb.execute.filter.summary.internal(class(old.set), removed$mm, relm, removed$samples, removed$sites)
}

rnb.execute.filter.summary.internal <- function(dataset.class, mm, relm, removed.samples, removed.sites) {

	cont.matrix <- rbind(
		"Sites" = as.double(c(nrow(mm) - length(removed.sites), length(removed.sites))),
		"Samples" = as.double(c(ncol(mm) - length(removed.samples), length(removed.samples))))
	rownames(cont.matrix)[1] <- capitalize(rnb.get.row.token(dataset.class, plural = TRUE))
	colnames(cont.matrix) <- c("Retained", "Removed")

	if (!is.null(relm)) {
		count.measurements <- function(x) {
			result <- table(x)[c("TRUE", "FALSE")]
			result[is.na(result)] <- 0L
			return(as.double(result))
		}
		cmatrix <- matrix(0, nrow = 2, ncol = 2)
		rownames(cmatrix) <- paste(c("Reliable", "Unreliable"), "measurements")
		for (i in 1:ncol(relm)) {
			if (i %in% removed.samples) {
				cmatrix[, 2] <- cmatrix[, 2] + count.measurements(relm[, i])
			} else if (length(removed.sites) != 0) {
				cmatrix[, 2] <- cmatrix[, 2] + count.measurements(relm[removed.sites, i])
				cmatrix[, 1] <- cmatrix[, 1] + count.measurements(relm[-removed.sites, i])
			} else {
				cmatrix[, 1] <- cmatrix[, 1] + count.measurements(relm[, i])
			}
		}
		cont.matrix <- rbind(cont.matrix, cmatrix)
	}
	cbind(cont.matrix, "Total" = rowSums(cont.matrix))
}

########################################################################################################################

#' rnb.step.filter.summary
#'
#' Calculates a table summarizing the effect of the applied filtering procedures and adds a corresponding section to the
#' given report.
#'
#' @param old.set Methylation dataset before filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param new.set Methylation dataset after filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return The report, possibly modified.
#'
#' @details
#' This function expects that the sites and samples in \code{new.set} are subsets of the sites and samples in
#' \code{old.set}, respectively. If this is not the case, it exists with an error. If \code{old.set} and \code{new.set}
#' are identical, no information is added to the given report.
#'
#' @seealso \code{\link{rnb.execute.filter.summary}}, \code{\link{rnb.run.filtering}}
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.filter.summary <- function(old.set, new.set, report) {
	if (!inherits(old.set, "RnBSet")) {
		stop("invalid value for old.set")
	}
	if (!inherits(new.set, "RnBSet")) {
		stop("invalid value for new.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}

	removed <- rnb.get.filtered.sites.samples(old.set, new.set)
	relm <- rnb.get.reliability.matrix(old.set)
	rnb.step.filter.summary.internal(class(old.set), removed$mm, relm, removed$samples, removed$sites, report,
		TRUE)
}

rnb.step.filter.summary.internal <- function(dataset.class, mm, relm, removed.samples, removed.sites, report,
	log.section = FALSE, section.name="Filtering Summary", section.order=1) {

	if (log.section) {
		logger.start("Summary of Filtering Procedures")
	}

	## Create a summary table of removed sites, samples and unreliable measurements
	table.summary <- rnb.execute.filter.summary.internal(dataset.class, mm, relm, removed.samples, removed.sites)
	if (all(table.summary[, "Removed"] == 0)) {
		if (log.section) logger.completed()
		return(report)
	}

	## Save table and create figure
	fname <- sprintf("summary%d.csv", section.order)
	utils::write.csv(table.summary, file = file.path(rnb.get.directory(report, "data", TRUE), fname))
	dframe <- table.summary[, -3]
	colnames(dframe) <- tolower(colnames(dframe))
	dframe <- data.frame(
		x = as.vector(dframe),
		item = factor(rep(rownames(dframe), ncol(dframe)), levels = rev(rownames(dframe))),
		group = factor(rep(colnames(dframe), each = nrow(dframe)), levels = colnames(dframe)))
	pp <- ggplot(dframe, aes_string("item", weight = "x", fill = "factor(group, levels = rev(levels(group)))")) +
		scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = percent_format()) +
		labs(x = NULL, y = "Fraction", fill = "Group") + theme(axis.ticks.y = element_blank()) +
		ggplot2::geom_bar(position = "fill", width = 0.7) + coord_flip() +
		theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in"))
	pname <- sprintf("summary%d_barchart", section.order)
	rplot <- createReportPlot(pname, report, width = 7, height = 0.55 + 0.4 * nrow(table.summary))
	suppressWarnings(print(pp))
	rplot <- off(rplot)
	fname <- paste(rnb.get.directory(report, "data"), fname, sep = "/")
	txt.site <- rnb.get.row.token(dataset.class)
	txt.sites <- rnb.get.row.token(dataset.class, plural = TRUE)
	rems <- table.summary[c(capitalize(txt.sites), "Samples"), "Removed"]
	txt <- c("As a final outcome of the filtering procedures, ", rems[1], " ",
		ifelse(rems[1] == 1, txt.site, txt.sites), " and ", rems[2], " sample", ifelse(rems[2] != 1, "s", ""), " were ",
		"removed. These statistics are presented in <a href=\"", fname, "\">a dedicated table</a> that accompanies ",
		"this report and visualized in the figure below.")
	report <- rnb.add.section(report, section.name, txt)
	txt <- "Fractions of removed values in the dataset after applying filtering procedures."
	report <- rnb.add.figure(report, txt, rplot)
	rm(table.summary, fname, dframe, pp, pname, rplot, rems, txt)
	logger.status("Added summary table of removed and retained items")

	## Construct vectors of removed and retained betas
	if (length(removed.samples) != 0) {
		betas.removed <- as.vector(mm[, removed.samples])
		mm <- mm[, -removed.samples]
	} else {
		betas.removed <- double()
	}
	if (length(removed.sites) != 0) {
		betas.removed <- c(betas.removed, mm[removed.sites, ])
		mm <- mm[-removed.sites, ]
	}
	beta.values <- list(Removed = betas.removed, Retained = as.vector(mm))
	for (i in 1:length(beta.values)) {
		beta.values[[i]] <- beta.values[[i]][!is.na(beta.values[[i]])]
	}
	logger.status("Constructed sequences of removed and retained methylation values")
	rm(betas.removed, mm, i)

	## Compare removed vs. retained methylation beta values
	if ((!is.null(beta.values)) && min(sapply(beta.values, length)) >= 501) {
		report.plots <- rnb.plot.beta.comparison(beta.values, sprintf("summary%d_betas", section.order), report)
		setting.names <- list("Plot type" =
			c("density" = "density estimation", "histogram" = "histograms", "qq" = "quantile-quantile plot"))
		txt <- c("The figure below compares the distributions of the removed methylation &beta; values and of the ",
			"retained ones.")
		rnb.add.paragraph(report, txt)
		txt <- "Comparison of removed and retained &beta; values."
		txt <- c(txt, add.text.subsampling(attr(report.plots, "subsampled"), paste(names(beta.values), "betas")))
		report <- rnb.add.figure(report, txt, report.plots, setting.names)
		logger.status("Added comparison between removed and retained beta values")
	}

	if (log.section) {
		logger.completed()
	}
	return(report)
}

## E N D ###############################################################################################################