########################################################################################################################
## greedycut.R
## created: 2012-03-27
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the iterative procedure for filtering probe and samples based on the given scores.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.greedycut.table
#'
#' Extracts the indication matrix to be used by Greedycut. This function also prints the appropriate messages to the
#' log.
#'
#' @param rnb.set Methylation dataset of interest.
#' @param report  Report to contain the Greedycut section. This parameter is only used in case rnb.set does not contain
#'                the required data for estimating beta value quality. If this is set to \code{NULL} (default), adding a
#'                section is disabled.
#' @return \code{list} of two items, named \code{"matrix"} and \code{"report"}. The first item is a \code{logical}
#'         matrix of the same dimensions as the beta value matrix of \code{rnb.set}. Every row signifies a site, and
#'         every column - a sample in the dataset. Unreliable measurements are encoded by \code{TRUE}. The second item
#'         in the returned list is the possibly modified report.
#' @author Yassen Assenov
#' @noRd
rnb.greedycut.table <- function(rnb.set, report = NULL) {
	if (inherits(rnb.set, "RnBeadSet")) {
		result <- dpval(rnb.set)
		if (is.null(result)) {
			rnb.info("Omitting Greedycut because detection p-values are missing")
			if (!is.null(report)) {
				txt <- "Greedycut was not executed because detection p-values are missing."
				report <- rnb.add.section(report, "Greedycut", txt)
			}
		} else {
			pval.threshold <- rnb.getOption("filtering.greedycut.pvalue.threshold")
			logger.info(c("Working with a p-value threshold of", pval.threshold))
			result[is.na(result)] <- 1
			result <- (result >= pval.threshold)
		}
	} else { # inherits(rnb.set, "RnBiSeq")
		result <- covg(rnb.set)
		if (is.null(result))	{
			rnb.info("Omitting Greedycut because read coverage data are missing")
			if (!is.null(report)) {
				txt <- "Greedycut was not executed because read coverage data are missing."
				report <- rnb.add.section(report, "Greedycut", txt)
			}
		} else {
			min.coverage <- rnb.getOption("filtering.coverage.threshold")
			rnb.info(c("Working with a minimal acceptable coverage of", min.coverage))
			result[is.na(result)] <- 0L
			result <- (result < min.coverage)
		}
	}
	return(list(matrix = result, report = report))
}

########################################################################################################################

#' rnb.execute.greedycut
#'
#' Executes the Greedycut procedure for probe and sample filtering based on the detection p-values, and calculates
#' statistics on its iterations.
#'
#' @param rnb.set HumanMethylation450K dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param rc.ties Flag indicating what the behaviour of the algorithm should be in case of ties between values of rows
#'                (probes) and columns (samples). See the corresponding parameter in
#'                \code{\link{greedycut.filter.matrix}} for more details.
#' @return \code{NULL} if \code{rnb.set} does not contain a matrix of detection p-values, or if all p-values denote
#'         reliable measurements. Otherwise, a list of the following elements:
#'         \describe{
#'           \item{"infos"}{Table summarizing the iterations of the algorithm, as returned by
#'                \code{\link{greedycut.filter.matrix}}.}
#'           \item{"statistics"}{Additional statistics on all iterations, as returned by
#'                \code{\link{greedycut.get.statistics}}.}
#'           \item{"iteration"}{Number of Greedycut iterations + \code{1} applied to the dataset, that is,
#'                a value of 1 indicates that the dataset was not modified.}
#'           \item{"sites"}{Indices of all sites to be removed.}
#'           \item{"samples"}{Indices of all samples to be removed.}
#'         }
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' greedy.result <- rnb.execute.greedycut(rnb.set.example)
#' # Number of applied iterations
#' greedy.result$iteration
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.greedycut <- function(rnb.set, rc.ties = rnb.getOption("filtering.greedycut.rc.ties")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!(is.character(rc.ties) && length(rc.ties) == 1 && (rc.ties %in% c("row", "column", "any")))) {
		stop("invalid value for rc.ties; expected one of: row, column, any")
	}
	rc.ties <- rc.ties[1]
	if (is.null(dpval(rnb.set))) {
		return(NULL)
	}
	beta.state <- rnb.greedycut.table(rnb.set)$matrix
	if (is.null(beta.state)) {
		return(NULL)
	}

	rnb.execute.greedycut.internal(beta.state, integer(), rc.ties)
}

rnb.execute.greedycut.internal <- function(beta.state, sites2ignore, rc.ties) {
	
	## Run the iterative algorithm for filtering out probes and samples
	iter.infos <- greedycut.filter.matrix(beta.state, sites2ignore, rc.ties)
	if (nrow(iter.infos) == 1) {
		rnb.info("No unreliable measurements found")
		return(NULL) # all(beta.state == 0)
	}
	rnb.status(c("Calculated a total of", nrow(iter.infos), "iterations"))

	## Calculate FPR and sensitivity
	iter.statistics <- greedycut.get.statistics(iter.infos)

	## Determine number of iterations to execute
	ind.iteration <- which.max(iter.statistics[["D"]])
	rnb.info(c("Optimal number of iterations is", ind.iteration))

	if (ind.iteration != 1) {
		## Update the matrices of beta and detection p-values
		rcindices <- iter.infos[1:ind.iteration, "Index"]
		rctypes <- iter.infos[1:ind.iteration, "Type"]
		filtered.sites <- rcindices[rctypes == "r"]
		filtered.samples <- rcindices[rctypes == "c"]
	} else {
		filtered.sites <- integer()
		filtered.samples <- integer()
	}

	return(list(infos = iter.infos, statistics = iter.statistics, iteration = ind.iteration,
		sites = filtered.sites, samples = filtered.samples))
}

########################################################################################################################

#' rnb.step.greedycut
#'
#' Performs the Greedycut procedure (if applicable) to filter the given methylation dataset and adds a
#' corresponding section to the specified report.
#'
#' @param rnb.set Dataset to be filtered. This must be an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @param rc.ties Flag indicating what the behaviour of the algorithm should be in case of ties between values of rows
#'                (sites or probes) and columns (samples). See the corresponding parameter in
#'                \code{\link{greedycut.filter.matrix}} for more details.
#' @return        List of three elements:
#'                \describe{
#'                  \item{\code{"dataset"}}{The (possibly modified) dataset after performing the Greeducyt procedure.}
#'                  \item{\code{"report"}}{The modified report.}
#'                  \item{\code{"filtered"}}{Indices of the sites in \code{rnb.set} that were filtered out.}
#'                }
#'
#' @seealso \code{\link{greedycut.filter.matrix}} for running all Greedycut iterations
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.greedycut <- function(rnb.set, report, rc.ties = rnb.getOption("filtering.greedycut.rc.ties")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (!(is.character(rc.ties) && length(rc.ties) == 1 && (rc.ties %in% c("row", "column", "any")))) {
		stop("invalid value for rc.ties; expected one of: row, column, any")
	}
	rc.ties <- rc.ties[1]
	rnb.step.greedycut.internal(rnb.set, integer(), annotation(rnb.set, add.names = inherits(rnb.set, "RnBeadSet")),
		rc.ties)
}

rnb.step.greedycut.internal <- function(rnb.set, sites2ignore, report, anno.table,
	rc.ties = rnb.getOption("filtering.greedycut.rc.ties")) {

	result <- rnb.greedycut.table(rnb.set, report)
	beta.state <- result$matrix
	report <- result$report
	rm(result)
	if (is.null(beta.state)) {
		return(list(sites = integer(), samples = integer(), report = report))
	}
	logger.start("Greedycut")

	result <- rnb.execute.greedycut.internal(beta.state, sites2ignore, rc.ties)
	if (is.null(result)) {
		logger.info("All measurements are reliable based on the specified threshold")
		txt <- "Greedycut was not executed because all measurements are reliable based on the specified threshold"
		report <- rnb.add.section(report, "Greedycut", txt)
		logger.completed()
		return(list(sites = integer(), samples = integer(), report = report))
	}
	unreliable.counts <- list("persite" = rowSums(beta.state), "persample" = colSums(beta.state))
	iter.infos <- result$infos
	iter.statistics <- result$statistics
	ind.iteration <- result$iteration
	rm(beta.state)

	## ----------------------------------------
	## Generate text and figures for the report

	txt.site <- rnb.get.row.token(rnb.set)
	txt.sites <- rnb.get.row.token(rnb.set, plural = TRUE)
	txt <- ifelse(grepl("Bead", class(rnb.set)), "detection p-value", "read coverage")
	txt <- c("The Greedycut algorithm iteratively removes from the dataset ", txt.sites, " and samples of highest ",
		"impurity. These correspond to the rows and columns in the ", txt, " table that contain the largest fraction ",
		"of unreliable measurements. This section summarizes the results of applying Greedycut on the analyzed ",
		"dataset.")
	report <- rnb.add.section(report, "Greedycut", txt)
	if (inherits(rnb.set, "RnBeadSet")) {
		txt <- c("We considered every &beta; value to be unreliable when its corresponding detection p-value is ",
			"not below the threshold <i>T</i>:")
		txt.formula <- c("<i>p</i> &ge; <i>T</i> = ", rnb.getOption("filtering.greedycut.pvalue.threshold"))
	} else {
		txt <- c("We considered every methylation value to be unreliable when its corresponding read coverage is ",
			"below the threshold <i>T</i>:")
		txt.formula <- c("<i>c</i> &lt; <i>T</i> = ", rnb.getOption("filtering.coverage.threshold"))
	}
	report <- rnb.add.section(report, "Unreliable Measurements", txt, level = 2)
	rnb.add.paragraph(report, txt.formula, "centered")
	rm(txt.formula)

	## Generate text and figures for the report using unreliable.counts, iter.infos, iter.statistics, ind.iteration
	report.plots <- list()
	for (iname in names(unreliable.counts)) {
		fname <- paste0("greedycut_cdf_", iname)
		report.plots <- c(report.plots, plotcdf(fname, report, unreliable.counts[[iname]]))
	}
	txt <- "The figure below summarizes the observed number of unreliable measurements per probe and per sample."
	rnb.add.paragraph(report, txt)
	txt <- c("Cumulative distribution function of number of unreliable values per ", txt.site, "/sample.")
	setting.names <- list("Number of values per" = c(txt.site, "sample"))
	names(setting.names[[1]]) <- names(unreliable.counts)
	report <- rnb.add.figure(report, txt, report.plots, setting.names)

	report.plots <- list()
	mdims <- iter.infos[, c("Columns", "Rows")]
	mdims[, 2] <- mdims[, 2] / 1000
	colnames(mdims) <- c("Samples", paste(capitalize(txt.sites), "(in thousand)"))
	iter.statistics[["Iteration"]] <- 1:nrow(iter.statistics)

	## Plot ROC curve
	rpoints <- unique(iter.statistics[, c("False positive rate", "Sensitivity")])
	colnames(rpoints) <- c("x", "y")
	plotroc <- function(fname, report, rocpoints, ind.iteration = NA) {
		rplot <- createReportPlot(fname, report, width = 6.2, height = 6.2)
		pp <- ggplot(data.frame(x = c(0, 1), y = c(0, 1)), aes_string(x = "x", y = "y")) +
			labs(x = "False positive rate", y = "Sensitivity") + coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
			geom_line(col = "#A0A0A0") + geom_line(data = rocpoints, col = "#000080", lwd = 2)
		if (!is.na(ind.iteration)) {
			pp <- pp + geom_point(data = rocpoints[ind.iteration, ], col = "#FF0000", pch = 20)
		}
		pp <- pp + theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in"))
		print(pp)
		return(off(rplot))
	}
	report.plots <- c(report.plots, plotroc("greedycut_roc_all", report, rpoints, ind.iteration))
	i2plot <- list("all" = 1:nrow(iter.infos))
	if (ind.iteration != 1 && ind.iteration != nrow(iter.infos)) {
		i2plot[["stop"]] <- 1:ind.iteration
		report.plots <- c(report.plots, plotroc("greedycut_roc_stop", report, rpoints[1:ind.iteration, ]))
	}
	logger.status("Created ROC plot")

	## Create line plots of matrix dimensions, FPR, sensitivity and more
	plotline <- function(fname, dframe, ind.stop = NA, col = "#000080", col.stop = "#800000", ...) {
		rplot <- createReportPlot(fname, report, width = 6.2, height = 6.2)
		vnames <- paste0("`", colnames(dframe), "`")
		pp <- ggplot(dframe, aes_string(x = vnames[1], y = vnames[2])) + geom_path(col = col)
		if (!is.na(ind.stop)) {
			pp <- pp + geom_point(data = dframe[ind.stop, ], col = col.stop, pch = 20)
		}
		pp <- pp + theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in"))
		print(pp)
		return(off(rplot))
	}
	for (i in 1:length(i2plot)) {
		ind.stop <- ifelse(length(i2plot) == 2 && i == 1, ind.iteration, NA)
		dframe <- iter.statistics[i2plot[[i]], ]
		icolumn <- which(colnames(dframe) == "Iteration")
		if (nrow(dframe) > 10000) {
			dframe[, icolumn] <- dframe[, icolumn] / 1000
			colnames(dframe)[icolumn] <- "Iteration (in thousand)"
		}
		icolumn <- colnames(dframe)[icolumn]
		## Plot matrix dimensions
		fname <- paste("greedycut_stairs_", names(i2plot)[i], sep = "")
		rplot <- plotline(fname, mdims[i2plot[[i]], ], ind.stop)
		report.plots <- c(report.plots, rplot)
		## Plot False positive rate
		fname <- paste("greedycut_fpr_", names(i2plot)[i], sep = "")
		rplot <- plotline(fname, dframe[, c(icolumn, "False positive rate")], ind.stop)
		report.plots <- c(report.plots, rplot)
		## Plot sensitivity
		fname <- paste("greedycut_sensitivity_", names(i2plot)[i], sep = "")
		rplot <- plotline(fname, dframe[, c(icolumn, "Sensitivity")], ind.stop)
		report.plots <- c(report.plots, rplot)
		## Plot distance from the diagonal
		fname <- paste("greedycut_ddiagonal_", names(i2plot)[i], sep = "")
		rplot <- plotline(fname, dframe[, c(icolumn, "D")], ind.stop)
		report.plots <- c(report.plots, rplot)
	}
	rm(mdims, rpoints, plotroc, plotline, i2plot, i, ind.stop, dframe, icolumn, fname, rplot)
	txt <- c("RnBeads executed Greedycut using the threshold given above and applied all its steps. ",
		"Briefly, Greedycut is an iterative algorithm that filters out the ", txt.site, " or sample with the highest ",
		"fraction of unreliable measurements one at a time. Note that every iteration of the algorithm produces a ",
		"matrix of retained measurements and a set of removed ones.")
	report <- rnb.add.section(report, paste("Filtered", capitalize(txt.sites), "and Samples"), txt, level = 2)
	txt <- c("We calculated false positive rate (<i>&alpha;</i>) and sensitivity (<i>s</i>) when the retained ",
		"measurements are considered as prediction for the reliable ones. Among all matrices produced by Greedycut, ",
		"we selected the one that maximizes the value of the expression <i>s</i> + 1 - <i>&alpha;</i>, thereby giving ",
		"equal weights to the sensitivity and specificity. Presented geometrically on a ROC curve, this is the point ",
		"that is furthest from the diagonal. The results of the Greedycut procedure and the selected iteration are ",
		"presented in the figure below.")
	rnb.add.paragraph(report, txt)
	txt <- c("Change of table dimensions / metric related to accuracy as Greedycut progressively removes ", txt.sites,
		" and samples. Accuracy is calculated by treating the retained entries as predictive of reliable ",
		"measurements. The red circle, if present, marks the last iteration that was executed.")
	setting.names <- list(
		"Metric" = c(
			"stairs" = "table dimensions",
			"fpr" = "false positive rate",
			"sensitivity" = "sensitivity",
			"roc" = "ROC curve",
			"ddiagonal" = "distance to the diagonal in a ROC plot"),
		"Iterations to show" = c("all" = "all", "stop" = "executed only"))
	report <- rnb.add.figure(report, txt, report.plots, setting.names)
	logger.status("Created line plots for matrix dimensions and other statistics")

	## Generate tables with the filtered sites and samples
	r.removed <- length(result$sites)
	c.removed <- length(result$samples)
	if (length(r.removed) != 0 || length(c.removed) != 0) {
		tbl.removed <- data.frame(
			Type = c(capitalize(txt.sites), "Samples"),
			Removed = c(r.removed, c.removed),
			Table = as.character(c(NA, NA)), stringsAsFactors = FALSE)

		if (r.removed != 0) {
			fname <- "removed_sites_greedycut.csv"
			fname.full <- rnb.save.removed.sites(anno.table, report, fname)
			tbl.removed[1, "Table"] <- paste0('<a href="', fname.full, '">', fname, '</a>')
		}
		if (c.removed != 0) {
			fname <- "removed_samples_greedycut.csv"
			fname.full <- file.path(rnb.get.directory(report, "data", TRUE), fname)
			write.csv(pheno(rnb.set)[result$samples, , drop = FALSE], fname.full, row.names = FALSE)
			fname <- paste0('<a href="', rnb.get.directory(report, "data"), '/', fname, '">', fname, '</a>')
			tbl.removed[2, "Table"] <- fname
		}

		txt <- c("Based on the criteria described above, ", r.removed, " ", ifelse(r.removed == 1, txt.site, txt.sites),
			" and ", c.removed, " sample", ifelse(c.removed == 1, "", "s"), " were filtered out. Links to the lists of ",
			"removed items are given below.")
		rnb.add.paragraph(report, txt)
		rnb.add.table(report, tbl.removed, row.names = FALSE, na = "")
	}

	logger.completed()
	return(list(sites = result$sites, samples = result$samples, report = report))
}

########################################################################################################################

#' greedycut.filter.matrix
#'
#' Performs all iterations of the Greedycut algorithm for removing rows and columns from the given matrix.
#'
#' @param mm          Numeric matrix to filter.
#' @param rows2ignore \code{integer} vector containing indices of rows in \code{mm} to be ignored by this  function.
#' @param rc.ties     Flag indicating what the behaviour of the algorithm should be in case of ties between values of
#'                    rows and columns. The value of this parameter must be one of \code{"row"}, \code{"column"} or
#'                    \code{"any"} (the last one indicating random choice).
#' @return Table summarizing the iterations of the algorithm in the form of a \code{data.frame} with the following
#'         columns : Index, Type, Score, Normalized score, Rows, Columns.
#'
#' @seealso \code{\link{greedycut.get.submatrix}} for extracting the resulting matrix after filtering
#' @author Yassen Assenov
#' @export
greedycut.filter.matrix <- function(mm, rows2ignore = integer(), rc.ties = "row") {
	if (any(mm < 0)) {
		stop("Negative values found")
	}
	if (!(rc.ties %in% c("row", "column", "any"))) {
		stop("Invalid value for rc.ties")
	}

	mm[rows2ignore, ] <- 0
	N <- nrow(mm); M <- ncol(mm)
	indices <- rep.int(as.integer(0), times = N + M)
	names(indices) <- rep.int("x", times = N + M)
	i <- 1
	fullscores <- c(rowSums(mm), colSums(mm))
	normscores <- c(rowMeans(mm), colMeans(mm))
	dimensions <- matrix(0, nrow = N + M, ncol = 2,
		dimnames = list(NULL, c("Rows", "Columns")))
	dimensions[1, ] <- c(N, M)

	maxfullscores <- rep(0, times = N + M)
	maxnormscores <- rep(0, times = N + M)

	repeat {
		## Find the maximum row or column score
		maxscore <- max(normscores)
		if (maxscore <= 0) {
			break
		}
		i <- i + 1
		max.ind <- which(normscores == maxscore)
		if (length(max.ind) != 1) {
			## Make a (random) selection in case of ties
			if (rc.ties == "row") {
				i1 <- max.ind[max.ind <= N]
				i2 <- max.ind[max.ind > N]
			} else if (rc.ties == "column") {
				i1 <- max.ind[max.ind > N]
				i2 <- max.ind[max.ind <= N]
			} else { # rc.ties = "any"
				i1 <- max.ind
				i2 <- c()
			}
			l <- length(i1)
			max.ind <- base::ifelse(l > 0, base::ifelse(l > 1, sample(i1, size = 1), i1), sample(i2, size = 1))
			rm(i1, i2, l)
		}

		## Save the maximum score and row or column number
		maxnormscores[i] <- maxscore
		maxfullscores[i] <- fullscores[max.ind]
		fullscores[max.ind] <- 0
		normscores[max.ind] <- 0
		if (max.ind <= N) {
			indices[i] <- max.ind
			names(indices)[i] <- "r"
			dimensions[i, ] <- dimensions[i - 1, ] - c(1, 0)
			fscores <- abs(fullscores[(N + 1):length(fullscores)] - mm[max.ind, ])
			fullscores[(N + 1):length(fullscores)] <- fscores
			normscores[(N + 1):length(normscores)] <- fscores / dimensions[i, 1]
			mm[max.ind, ] <- 0
		} else { # nrow(mm) < max.ind <= nrow(mm) + ncol(mm)
			max.ind <- max.ind - N
			indices[i] <- max.ind
			dimensions[i, ] <- dimensions[i - 1, ] - c(0, 1)
			names(indices)[i] <- "c"
			fscores <- abs(fullscores[1:N] - mm[, max.ind])
			fullscores[1:N] <- fscores
			normscores[1:N] <- fscores / dimensions[i, 2]
			mm[, max.ind] <- 0
		}
	}

	return(data.frame(
			"Index" = indices[1:i],
			"Type" = names(indices)[1:i],
			"Score" = maxfullscores[1:i],
			"Normalized score" = maxnormscores[1:i],
			"Rows" = dimensions[1:i, "Rows"] - length(rows2ignore),
			"Columns" = dimensions[1:i, "Columns"], check.names = FALSE, stringsAsFactors = TRUE))
}

########################################################################################################################

#' greedycut.get.statistics
#'
#' Calculates various statistics on the iterations of Greedycut.
#'
#' @param filterinfo Information on the filtering iterations as a \code{data.frame} returned by
#'                   \code{\link{greedycut.filter.matrix}}.
#' @return Additional statistics on the iterations in the form of a \code{data.frame} with the following columns:
#'         \code{"Elements retained"}, \code{"Elements removed"}, \code{"Mismatches retained"},
#'         \code{"Mismatches removed"}, \code{"False Positive Rate"}, \code{"Sensitivity"}, \code{"D"}. The last column
#'         signifies distance from the diagonal in a ROC curve. 
#'
#' @author Yassen Assenov
#' @export
greedycut.get.statistics <- function(filterinfo) {
	elements.retained <- filterinfo[, "Rows"] * filterinfo[, "Columns"]
	elements <- elements.retained[1]
	elements.removed <- elements - elements.retained
	mismatches.removed <- cumsum(filterinfo[["Score"]])
	mismatches <- mismatches.removed[length(mismatches.removed)]
	mismatches.retained <- abs(mismatches - mismatches.removed)
	fprs <- mismatches.retained / mismatches
	sensitivities <- (elements.retained - mismatches.retained) / (elements - mismatches)
	return(data.frame(
			"Elements retained" = elements.retained,
			"Elements removed" = elements.removed,
			"Mismatches retained" = mismatches.retained,
			"Mismatches removed" = mismatches.removed,
			"False positive rate" = fprs,
			"Sensitivity" = sensitivities,
			"D" = sqrt((fprs - sensitivities)^2 / 2),
			check.names = FALSE))
}

########################################################################################################################

#' greedycut.get.submatrix
#'
#' Filters a data matrix executing the given number of iterations of Greedycut.
#'
#' @param mm          Data \code{matrix} to be filtered.
#' @param filter.info Information on the filtering iterations as a \code{data.frame} returned by
#'                    \code{\link{greedycut.filter.matrix}}.
#' @param it.num      Number of iterations to execute. Defaults to all iterations.
#' @return Data matrix containing subsets of the rows and columns of \code{mm}.
#'
#' @author Yassen Assenov
#' @export
greedycut.get.submatrix <- function(mm, filter.info, it.num = nrow(filter.info) - as.integer(1)) {
	if (!is.matrix(mm)) {
		stop("invalid value for mm; expected matrix")
	}
	if (!(is.data.frame(filter.info) &&
			identical(colnames(filter.info), c("Index", "Type", "Score", "Normalized score", "Rows", "Columns")))) {
		stop("invalid value for filterinfo")
	}
	if (is.double(it.num) && all(it.num == as.integer(it.num))) {
		it.num <- as.integer(it.num)
	}
	if (!(is.integer(it.num) && length(it.num) == 1 && (!is.na(it.num)))) {
		stop("invalid value for it.num; expected single integer")
	}
	it.num <- it.num[1]
	if (!(1 <= it.num && it.num < nrow(filter.info))) {
		stop("invalid value for it.num")
	}
	i2remove <- filter.info[2:(it.num + 1), "Index"]
	what2remove <- filter.info[2:(it.num + 1), "Type"]
	ri2remove <- i2remove[what2remove == "r"]
	ci2remove <- i2remove[what2remove == "c"]
	if (any(i2remove < 1) || any(ri2remove > nrow(mm)) || any(ci2remove > ncol(mm))) {
		stop("unexpected dimensions of the data matrix")
	}
	ris <- setdiff(1:nrow(mm), ri2remove)
	cis <- setdiff(1:ncol(mm), ci2remove)
	return(mm[ris, cis])
}

## E N D ###############################################################################################################
