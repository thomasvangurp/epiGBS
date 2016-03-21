########################################################################################################################
## batch.R
## created: 2012-06-21
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of batch effects analysis and correction.
########################################################################################################################

## G L O B A L S #######################################################################################################

## Minimum percentage of variance explained; the first n principal components that explain the variance are considered
VAR.EXPLAINED <- 95

## Color codes used to represent tests for association between traits
COLORS.TEST <- c("Correlation" = "#F8766D", "Fisher" = "#7CAE00", "Wilcoxon" = "#00BFC4", "Kruskal-Wallis" = "#C77CFF")

## Color codes used to represent non-significant and significant p-values
COLORS.SIGNIFICANT <- c("FALSE" = "#9C9CCC", "TRUE" = "#F1B1B1")

## Point shapes used to represent reasons for failure of an association test
SHAPES.FAILURES <- c("not enough common values" = 22, "not enough categories" = 24, "test failed" = 21)

## F U N C T I O N S ###################################################################################################

#' mds
#'
#' Performs Kruskal's non-metric multidimensional scaling using the supplied distance metric.
#'
#' @param m           Data matrix to be scaled. Rows should denote observations, and columns - dimensions.
#' @param dist.metric Distance metric. Must be one of \code{"euclidean"} or \code{"manhattan"}.
#' @return            Two-column matrix of transformed points.
#' @author Yassen Assenov
#' @noRd
mds <- function(m, dist.metric) {
	dist.matrix <- dist(m, method = dist.metric)
	if (any(dist.matrix == 0)) {
		## Get rid of zeros in the distance matrix by adding an epsilon
		dist.matrix <- dist.matrix + min(dist.matrix[dist.matrix != 0]) / 1000
	}
	isoMDS(dist.matrix, k = 2, maxit = 100, trace = FALSE)$points
}

########################################################################################################################

#' validate.pcoordinates
#'
#' Validates that the given object is returned by \code{\link{rnb.execute.dreduction}}.
#'
#' @param pcoordinates Object to validate.
#' @return If \code{pcoordinates} is valid - the number of samples in the dataset that were processed; otherwise - an
#'         error message in the form of a one-element \code{character}.
#'
#' @author Yassen Assenov
#' @noRd
validate.pcoordinates <- function(pcoordinates) {
	## Validate the data structure
	if (!is.list(pcoordinates)) {
		return("invalid value for pcoordinates")
	}
	if (!setequal(names(pcoordinates), c("pca", "mds"))) {
		return("invalid value for pcoordinates")
	}
	if (!(is.null(pcoordinates[["pca"]]) || class(pcoordinates[["pca"]]) == "prcomp")) {
		return("invalid value for pcoordinates$pca")
	}
	pc.mds <- pcoordinates[["mds"]]
	if (!is.list(pc.mds)) {
		return("invalid value for pcoordinates$mds")
	}
	if (length(pc.mds) != 0) {
		if (!(setequal(names(pc.mds), c("euclidean", "manhattan")) &&
			  	all(sapply(pc.mds, is.matrix)) && all(sapply(pc.mds, is.double)))) {
			return("invalid value for pcoordinates$mds")
		}
		nsamples <- sapply(pc.mds, dim)
		if (any(nsamples[2,] != 2)) {
			return("invalid value for pcoordinates$mds; expected 2 dimensions")
		}
		if (any(nsamples[1, ] != nsamples[1, 1])) {
			return("invalid value for pcoordinates$mds; inconsistent number of samples")
		}
		nsamples <- as.integer(nsamples[1, 1])
	} else {
		nsamples <- 0L
	}

	## Validate the number of samples is consistent
	if (!is.null(pcoordinates[["pca"]])) {
		if (nsamples == 0) {
			nsamples <- nrow(pcoordinates[["pca"]]$x)
		} else if (!isTRUE(nsamples == nrow(pcoordinates[["pca"]]$x))) {
			return("invalid value for pcoordinates; inconsistent number of samples")
		}
	}
	return(nsamples)
}

########################################################################################################################

## validate.pcoordinates.all
##
## Validates that the given object is either returned by \code{\link{rnb.execute.dreduction}}, or is list of objects
## returned by that functions, applied on the same samples.
##
## @param pcoordinates Object to validate.
## @return If \code{pcoordinates} is valid - the number of samples in the dataset that were processed; otherwise - an
##         error message in the form of a one-element \code{character}.
##
## @author Yassen Assenov
validate.pcoordinates.all <- function(pcoordinates) {
	if (!(is.list(pcoordinates) && length(pcoordinates) != 0)) {
		return("invalid value for pcoordinates")
	}
	if (setequal(names(pcoordinates), c("mds", "pca"))) {
		nsamples <- validate.pcoordinates(pcoordinates)
	} else {
		nsamples <- lapply(pcoordinates, validate.pcoordinates)
		i.error <- which(sapply(nsamples, is.character))
		if (length(i.error) != 0) {
			nsamples <- nsamples[[i.error[1]]]
		} else {
			nsamples <- sapply(nsamples, identity)
			nsamples.max <- max(nsamples)
			if (!all(nsamples[nsamples != 0] == nsamples.max)) {
				return("invalid value for pcoordinates; inconsistent number of samples")
			}
			nsamples <- nsamples.max
		}
	}
	return(nsamples)
}

########################################################################################################################

## Checks if the given argument is a valid sample index permutation matrix.
##
## @param perm.matrix Matrix to be validated.
## @return \code{TRUE} if \code{perm.matrix} is a valid index permutation matrix, that is, an integer matrix in which
##         every column denotes a permutation; \code{FALSE} otherwise.
## @author Yassen Assenov
is.valid.permutations <- function(perm.matrix) {
	if (!(is.integer(perm.matrix) && is.matrix(perm.matrix) && nrow(perm.matrix) > 0 && ncol(perm.matrix) > 0)) {
		return(FALSE)
	}
	all.indices <- 1:nrow(perm.matrix)
	if (!identical(perm.matrix[, 1], all.indices)) {
		return(FALSE)
	}
	if (ncol(perm.matrix) > 1) {
		if (!all(apply(perm.matrix[, -1], 2, function(x) { identical(sort(x), all.indices) }))) {
			return(FALSE)
		}
	}
	return(TRUE)
}

########################################################################################################################

#' Tests for association between two traits.
#'
#' @param x           Sample values for the first trait. This must be a vector of type \code{factor}, \code{integer} or
#'                    \code{numeric}.
#' @param y           Sample values for the second trait. This must be a vector of type \code{factor}, \code{integer} or
#'                    \code{numeric}.
#' @param perm.matrix Matrix of sample permutations to be used in case none of the traits is a \code{factor}, and thus
#'                    permutation-based p-value from correlations is computed. If this parameter is \code{NULL} and
#'                    both \code{x} and \code{y} are sequences of numbers, no p-value is calculated.
#' @return            List of four elements:
#'                    \describe{
#'                      \item{error}{Error, if any, that prevented this function from computing a p-value for trait
#'                           association.}
#'                      \item{test}{Type of test performed. This is one of \code{"Fisher"}, \code{"Wilcoxon"},
#'                           \code{"Kruskal-Wallis"}, \code{"Correlation"} or \code{NA}. The last value indicates that
#'                           the traits cannot be tested for association.}
#'                      \item{correlation}{Value of the pearson correlation coefficient between \code{x} and \code{y},
#'                           or \code{NA} if any of them is \code{factor}.}
#'                      \item{pvalue}{Calculated p-value, or \code{NA} if the traits cannot be tested for association.}
#'                    }
#' @author Yassen Assenov
#' @noRd
test.traits <- function(x, y, perm.matrix = NULL) {
	result <- list(
			"error" = as.character(NA),
			"test" = as.character(NA),
			"correlation" = as.double(NA),
			"pvalue" = as.double(NA))
	
	## Focus on common values
	if (class(x) == "Date") { x <- as.integer(x) }
	if (class(y) == "Date") { y <- as.integer(y) }
	inds <- which(!(is.na(x) | is.na(y)))
	if (length(inds) < 2) {
		result[["error"]] <- "not enough common values"
		return(result)
	}
	x <- x[inds]
	if (is.factor(x)) {
		x <- as.factor(as.character(x))
		if (nlevels(x) < 2) {
			## Not enough categories in y
			result[["error"]] <- "not enough categories"
			return(result)
		}
	}
	y <- y[inds]
	if (is.factor(y)) {
		y <- as.factor(as.character(y))
		if (nlevels(y) < 2) {
			## Not enough categories in y
			result[["error"]] <- "not enough categories"
			return(result)
		}
	}
	
	## Perform a test or compute correlation
	get.p <- function(expr) { tryCatch(suppressWarnings(expr$p.value), error = function(er) { as.double(NA) }) }
	if (is.factor(x)) {
		if (is.factor(y)) {
			simulate <- (nlevels(x) > 2 || nlevels(y) > 2)
			result[["test"]] <- "Fisher"
			result[["pvalue"]] <- get.p(fisher.test(x, y, conf.int = FALSE, simulate.p.value = simulate, B = 50000))
		} else if (nlevels(x) == 2) {
			result[["test"]] <- "Wilcoxon"
			values <- tapply(y, x, identity)
			result[["pvalue"]] <- get.p(wilcox.test(values[[1]], values[[2]], alternative = "two.sided"))
		} else {
			result[["test"]] <- "Kruskal-Wallis"
			result[["pvalue"]] <- get.p(kruskal.test(y, x))
		}
	} else if (is.factor(y)) {
		if (nlevels(y) == 2) {
			result[["test"]] <- "Wilcoxon"
			values <- tapply(x, y, identity)
			result[["pvalue"]] <- get.p(wilcox.test(values[[1]], values[[2]], alternative = "two.sided"))
		} else {
			result[["test"]] <- "Kruskal-Wallis"
			result[["pvalue"]] <- get.p(kruskal.test(x, y))
		}
	} else {
		result[["test"]] <- "Correlation"
		if (is.null(perm.matrix)) {
			result[["correlation"]] <- cor(x, y)
		} else {
			N <- length(inds)
			values <- apply(perm.matrix, 2, function(i) { cor(x[i[i <= N]], y) })
			result[["correlation"]] <- values[1]
			values <- abs(values)
			result[["pvalue"]] <- mean(values[1] <= values)
		}
	}
	
	if (is.na(result[["pvalue"]])) {
		result[["error"]] <- "test failed"
	}
	return(result)
}

########################################################################################################################

## plot.heatmap.pc.correlations
##
## Creates a heatmap displaying a table of correlation values.
##
## @param report Report to contain the generated heatmap.
## @param tbl    Table of correlation values to be plotted.
## @param fname  File name for the plot.
## @param width  Width of the plot in inches. If this is set to \code{NULL} (default), the width is calculated based on
##               the number of columns in \code{tbl}.
## @param height Height of the plot in inches. If this is set to \code{NULL} (default), the height is calculated based
##               on the number of rows in \code{tbl}.
## @return \code{list} with two elements: \code{"plot"} and \code{"description"} containing the \code{\link{ReportPlot}}
##         object and plot's description, respectively.
## @author Yassen Assenov
plot.heatmap.pc.correlations <- function(report, tbl, fname, width = NULL, height = NULL) {
	tbl.melt <- melt(tbl, varnames = c("x", "y"))
	colnames(tbl.melt)[3] <- "correlation"
	tbl.melt[[1]] <- factor(as.character(tbl.melt[[1]]), levels = rev(rownames(tbl)))
	tbl.melt[[2]] <- factor(as.character(tbl.melt[[2]]), levels = colnames(tbl))
	xlab <- names(dimnames(tbl))[2]
	if ((!is.null(xlab)) && (is.na(xlab) || identical(xlab, ""))) { xlab <- NULL }
	ylab <- names(dimnames(tbl))[1]
	if ((!is.null(ylab)) && (is.na(ylab) || identical(ylab, ""))) { ylab <- NULL }
	xsize <- ifelse(is.null(xlab), 2, 0.25)
	ysize <- ifelse(is.null(ylab), 2, 0.25)
	if (is.null(width)) {
		width <- ysize + ncol(tbl) * 0.28 + 1.6
	}
	if (is.null(height)) {
		height <- xsize + nrow(tbl) * 0.28 + 0.2
	}
	colors.g <- rnb.getOption("colors.3.gradient")
	pp <- ggplot(tbl.melt) + aes_string("y", "x") + labs(x = xlab, y = ylab) + coord_fixed() +
		geom_tile(aes_string(fill = "correlation"), color = "white") +
		scale_fill_gradient2(limits = c(-1, 1), low = colors.g[1], mid = colors.g[2], high = colors.g[3]) +
		scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
		theme(axis.ticks = element_blank(), legend.justification = c(0, 0.5), legend.position = c(1, 0.5)) +
		theme(panel.border = element_blank(), plot.margin = unit(c(0.1, 1.5, 0.1, 0.1), "in"))

	rplot <- createReportPlot(fname, report, width = width, height = height)
	pp <- suppressWarnings(ggplot_gtable(ggplot_build(pp)))
	pp$widths[[3]] <- unit(ysize, "in")
	pp$heights[[length(pp$heights) - 2L]] <- unit(xsize, "in")
	grid.draw(pp)
	rplot <- off(rplot)
	txt <- "Heatmap presenting a table of correlations. Grey cells, if present, denote missing values."
	return(list("plot" = rplot, "description" = txt))
}

########################################################################################################################

plot.heatmap.pc.pvalues <- function(report, tbl, fname, width = NULL, height = NULL) {
	tbl.melt <- melt(tbl, varnames = c("x", "y"))
	colnames(tbl.melt)[3] <- "pvalue"
	tbl.melt[[1]] <- factor(as.character(tbl.melt[[1]]), levels = rev(rownames(tbl)))
	tbl.melt[[2]] <- factor(as.character(tbl.melt[[2]]), levels = colnames(tbl))
	tbl.melt[["significant"]] <-
		as.character(tbl.melt[, "pvalue"] < rnb.getOption("exploratory.correlation.pvalue.threshold"))
	tbl.melt[is.na(tbl.melt[["significant"]]), "significant"] <- "NA"
	tbl.melt[["pvaltext"]] <- toupper(format(tbl.melt[, "pvalue"], digits = 2, scientific = TRUE))
	tbl.melt[["pvaltext"]] <- ifelse(tbl.melt[["significant"]] == "TRUE", tbl.melt[["pvaltext"]], "")
	xlab <- names(dimnames(tbl))[2]
	if ((!is.null(xlab)) && (is.na(xlab) || identical(xlab, ""))) { xlab <- NULL }
	ylab <- names(dimnames(tbl))[1]
	if ((!is.null(ylab)) && (is.na(ylab) || identical(ylab, ""))) { ylab <- NULL }
	if (is.null(width)) {
		width <- ifelse(is.null(ylab), 2, 1) + ncol(tbl) * 0.62
	}
	if (is.null(height)) {
		height <- ifelse(is.null(xlab), 2, 1) + nrow(tbl) * 0.25
	}
	rplot <- createReportPlot(fname, report, width = width, height = height)
	pp <- ggplot(tbl.melt, aes_string("y", "x", label = "pvaltext")) + labs(x = xlab, y = ylab) +
		coord_fixed(ratio = 0.25 / 0.62) + geom_tile(aes_string(fill = "significant"), color = "white") +
		scale_fill_manual(values = c("NA" = "#D0D0D0", "FALSE" = "#9C9CCC", "TRUE" = "#F1B1B1")) +
		geom_text(size = 3) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
		theme(panel.border = element_blank(), panel.grid = element_blank()) +
		theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in")) +
		theme(axis.ticks = element_blank(), legend.position = "none") 
	## Fix the areas for x and y axis labels
	pp <- suppressWarnings(ggplot_gtable(ggplot_build(pp)))
	if (is.null(ylab)) {
		pp$widths[[3]] <- unit(2, "in")
	}
	if (is.null(xlab)) {
		pp$heights[[length(pp$heights) - 2L]] <- unit(2, "in")
	}
	grid.draw(pp)
	txt <- paste0("Heatmap presenting a table of p-values. Significant p-values (less than ",
		rnb.getOption("exploratory.correlation.pvalue.threshold"), ") are printed in pink boxes. Non-significant ",
		"values are represented by blue boxes. Bright grey cells, if present, denote missing values.")
	return(list("plot" = off(rplot), "description" = txt))
}

########################################################################################################################

## plot.heatmap.symm
##
## Creates a lower-diagonal heatmap representing the selected tests applied between pairs of traits, or the resulting
## p-values.
##
## @param report       Report to contain the generated plot(s).
## @param tbl.symm     Symmetric \code{matrix} of type \code{character} storing test names, or of type \code{double}
##                     storing resulting p-values. This function expects the test names to be among the following:
##                     \code{"Correlation"}, \code{"Fisher"}, \code{"Wilcoxon"}, \code{"Kruskal-Wallis"}.
## @param tbl.failures Optional, a symmetric \code{matrix} of type \code{character} storing explanations for test
##                     failures. If provided, this matrix must have the same dimensions and dimension names as
##                     \code{tbl.symm} and its values must be among the following: \code{"not enough common values"},
##                     \code{"not enough categories"} and \code{"test failed"}.
## @param fname        File name for the plot.
## @return Generated plot(s) as an object of type \code{\linkS4class{ReportPlot}} or a \code{list} of such objects.
##
## @details
## The type of heatmap to be draw is determined based on the type of the provided symmetric matrix:
## \code{typeof(tbl.symm)}. Test names are assumed when this type is \code{character}, and p-values otherwise.
##
## @author Yassen Assenov
plot.heatmap.symm <- function(report, tbl.symm, tbl.failures = NULL, fname) {
	do.tests <- (typeof(tbl.symm) == "character")
	# Note that tbl.symm must be symmetric of size at least 2x2; this is not checked

	## Create a data frame with tests/p-values
	tbl <- tbl.symm
	tbl[upper.tri(tbl)] <- NA
	tbl <- tbl[-1, -ncol(tbl)]
	if (nrow(tbl.symm) == 2) {
		tbl <- matrix(tbl, nrow = 1, ncol = 1, dimnames = list(rownames(tbl.symm)[2], colnames(tbl.symm)[1]))
	}
	tbl.melt <- melt(tbl, varnames = c("x", "y"))
	colnames(tbl.melt)[3] <- "test"
	if (do.tests) {
		col.mapping <- COLORS.TEST
		tbl.melt[["test"]] <- factor(as.character(tbl.melt[["test"]]), levels = names(col.mapping))
		tbl.melt[["pvaltext"]] <- ""
	} else {
		col.mapping <- COLORS.SIGNIFICANT
		tbl.melt[["significant"]] <-
			as.character(tbl.melt[, "test"] < rnb.getOption("exploratory.correlation.pvalue.threshold"))
		tbl.melt[["pvaltext"]] <- toupper(format(tbl.melt[, "test"], digits = 2, scientific = TRUE))
		tbl.melt[["pvaltext"]] <- ifelse(tbl.melt[["significant"]] == "TRUE", tbl.melt[["pvaltext"]], "")
		tbl.melt[["test"]] <- tbl.melt[["significant"]]
	}
	tbl.melt[[1]] <- factor(as.character(tbl.melt[[1]]), levels = rev(rownames(tbl)))
	tbl.melt[[2]] <- factor(as.character(tbl.melt[[2]]), levels = colnames(tbl))

	if (is.matrix(tbl.failures) && (!all(is.na(tbl.failures)))) {
		## Create a data frame with test failures
		tbl <- tbl.failures
		tbl[upper.tri(tbl)] <- NA
		tbl <- tbl[-1, -ncol(tbl)]
		tbl.f.melt <- melt(tbl, varnames = c("x", "y"))
		colnames(tbl.f.melt)[3] <- "failure"
		tbl.f.melt[[1]] <- factor(as.character(tbl.f.melt[[1]]), levels = levels(tbl.melt[[1]]))
		tbl.f.melt[[2]] <- factor(as.character(tbl.f.melt[[2]]), levels = levels(tbl.melt[[2]]))
		tbl.f.melt[[3]] <- factor(as.character(tbl.f.melt[[3]]), levels = names(SHAPES.FAILURES))
		tbl.f.melt <- tbl.f.melt[!is.na(tbl.f.melt[[3]]), ]
	}

	## Create the heatmap
	height <- 2.2 + nrow(tbl) * 0.25
	width <- 4 + nrow(tbl) * 0.62
	pp <- ggplot(tbl.melt, aes_string("y", "x", label = "pvaltext")) + labs(x = NULL, y = NULL) +
		coord_fixed(ratio = 0.25 / 0.62) + geom_tile(aes_string(fill = "test"), color = "white") +
		scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
		scale_fill_manual(na.value = "white", values = col.mapping)
	if (exists("tbl.f.melt")) {
		pp <- pp + geom_point(aes_string(shape = "failure", label = NULL), data = tbl.f.melt, color = "#808080",
				fill = "#808080") + scale_shape_manual(values = SHAPES.FAILURES)
	}
	if (do.tests) {
		pp <- pp + theme(legend.justification = c(0, 1), legend.position = c(1, 1))
	} else {
		pp <- pp + geom_text(size = 3) + theme(legend.position = "none")
	}
	pp <- pp + theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
		theme(panel.grid.major = element_blank(), panel.background = element_blank()) +
		theme(panel.border = element_blank(), plot.margin = unit(c(0.1, 1.9, 0.1, 0.1), "in"))
	rplot <- createReportPlot(fname, report, width = width, height = height)
	## Fix the areas for x and y axis labels
	pp <- suppressWarnings(ggplot_gtable(ggplot_build(pp)))
	pp$widths[[3]] <- unit(2, "in")
	pp$heights[[length(pp$heights) - 2L]] <- unit(2, "in")
	grid.draw(pp)

	return(off(rplot))
}

########################################################################################################################
########################################################################################################################

#' rnb.execute.dreduction
#'
#' Performs principal component analysis (PCA) and multi-dimensional scaling (MDS) of the samples in the given
#' methylation dataset.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}. This dataset
#'                must contain at least four samples.
#' @param target  \code{character} singleton specifying the level of DNA methylation infromation. If this is
#' 				  \code{"sites"}, the DNA methylation information for the individual sites or probes is analyzed.
#' 				  Otherwise, this should be one of the supported region types, as returned by
#'                \code{\link{rnb.region.types}}. 
#' @return Results of the dimension reduction in the form of a list with the following elements:
#'         \describe{
#'           \item{\code{pca}}{Results of the PCA as returned by the function \code{\link{prcomp}}.}
#'           \item{\code{mds}}{List of two elements - \code{"manhattan"} and \code{"euclidean"}, each of which is a
#'                two-column \code{matrix} storing the coordinates of the samples in a two-dimensional space. The
#'                matrices are computed using the function \code{\link{isoMDS}}.}
#'         }
#'
#' @details
#' Row names in the returned matrices are sample identifiers, determined based on the package option
#' \code{"identifiers.column"}. See \emph{\link[=rnb.options]{RnBeads Options}} for more information on this option.
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' regs <- c("sites", summarized.regions(rnb.set.example))
#' dreduction <- function(x) rnb.execute.dreduction(rnb.set.example, x)
#' pcoordinates <- lapply(regs, dreduction)
#' names(pcoordinates) <- regs
#' str(pcoordinates)
#' }
#' @seealso \code{\link{rnb.run.exploratory}} for running the whole exploratory analysis module
#' @author Yassen Assenov
#' @export
rnb.execute.dreduction <- function(rnb.set, target = "sites") {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!(is.character(target) && length(target) == 1 && (!is.na(target)))) {
		stop("invalid value for target")
	}
	X <- t(meth(rnb.set, target))
	if (nrow(X) < 4) {
		stop("too few samples")
	}
	rownames(X) <- samples(rnb.set)

	pcoordinates <- list(mds = list(), pca = NULL)
	top.sites <- rnb.getOption("exploratory.top.dimensions")
	site.nas <- colMeans(is.na(X))
	i <- which(site.nas == 1)
	if (length(i) != 0) {
		if (ncol(X) - length(i) < 4) {
			stop("too many missing values")
		}
		## Remove sites that contain NAs only
		X <- X[, -i]
		site.nas <- site.nas[-i]
	}
	dframe.selected <- data.frame( # table with statistics on used sites
		target = rep(target, 2),
		technique = c("MDS", "PCA"),
		all = rep(ncol(X), 2),
		na = rep(length(i), 2),
		selected = rep(ncol(X) - length(i), 2),
		explained = c(1, 1), row.names = c("mds", "pca"), stringsAsFactors = FALSE)
	site.nas <- which(site.nas != 0)
	if (top.sites != 0) {
		select.top <- function(X, dr.method) {
			site.vars <- colVars(X, na.rm = TRUE)
			i.sites <- order(site.vars, decreasing = TRUE)[1:top.sites]
			dframe.selected[dr.method, "selected"] <<- top.sites
			dframe.selected[dr.method, "explained"] <<- sum(site.vars[i.sites]) / sum(site.vars, na.rm=TRUE)
			rnb.status(c("Selected the", top.sites, "dimensions with highest variance"))
			X[, i.sites, drop = FALSE]
		}
	}

	## Run MDS
	rnb.logger.start("MDS")
	if (0 != top.sites && top.sites < ncol(X)) {
		X.mds <- select.top(X, "mds")
		if (length(site.nas) == 0) {
			X <- X.mds
			dframe.selected["pca", "selected"] <- dframe.selected["mds", "selected"]
			i <- gc()
		}
	} else {
		X.mds <- X
	}
	for (dist.metric in c("manhattan", "euclidean")) {
		pcoordinates$mds[[dist.metric]] <- tryCatch(mds(X.mds, dist.metric), error = function(er) { NULL })
		rnb.status(c("Calculated MDS coordinates using", dist.metric, "distance"))
	}
	rnb.logger.completed()
	rm(X.mds)

	## Run PCA
	rnb.logger.start("PCA")
	if (length(site.nas) != 0) {
		## Remove sites that contain NA in any sample
		dframe.selected["pca", "na"] <- dframe.selected["pca", "na"] + length(site.nas)
		dframe.selected["pca", "selected"] <- dframe.selected["pca", "selected"] - length(site.nas)
		if (length(site.nas) >= ncol(X) - nrow(X)) {
			attr(pcoordinates, "selected") <- dframe.selected
			rnb.info("Skipped due to too many missing values")
			rnb.logger.completed()
			return(pcoordinates)
		}
		X <- X[, -site.nas]
		rnb.info(c("Removed", length(site.nas), "loci (", target, ") because they contain missing values"))
		if (0 != top.sites && top.sites < ncol(X)) {
			X <- select.top(X, "pca")
		}
	}
	suppressWarnings(rm(site.nas, top.sites, i, select.top))
	pcoordinates$pca <- prcomp(X, center = TRUE, scale. = FALSE)
	pcoordinates$pca$rotation <- NULL # remove rotation matrix to reduce memory load
	attr(pcoordinates, "selected") <- dframe.selected
	rnb.logger.completed()
	return(pcoordinates)
}

########################################################################################################################

#' rnb.section.dreduction
#'
#' Creates a report section dedicated to dimension reduction.
#'
#' @param report            Report to contain the section. This must be an object of type \code{\linkS4class{Report}}.
#' @param pcoordinates      Coordinates of data points, as an object returned by \code{\link{rnb.execute.dreduction}},
#'                          or a list of such objects.
#' @param sample.phenotypes Table of sample phenotype information in the form of a \code{data.frame}. The number of rows
#'                          in this table must be consistent with the samples in \code{pcoordinates}. This parameter is
#'                          used in determining point colors and/or shapes when plotting the samples in the reduced
#'                          dimensional representations. Set this parameter to \code{NULL} if no phenotype data is to be
#'                          visualized.
#' @return The modified report.
#'
#' @seealso \code{\link{rnb.step.dreduction}}
#' @author Yassen Assenov
#' @noRd
rnb.section.dreduction <- function(report, pcoordinates, sample.phenotypes = NULL) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	nsamples <- validate.pcoordinates.all(pcoordinates)
	
	if (is.character(nsamples)) {
		stop(nsamples)
	}
	if (!is.null(sample.phenotypes)) {
		if (!is.data.frame(sample.phenotypes)) {
			stop("invalid value for sample.phenotypes")
		}
		if (nrow(sample.phenotypes) != nsamples) {
			stop("invalid value for sample.phenotypes; unexpected number of samples")
		}
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	if (setequal(names(pcoordinates), c("mds", "pca"))) {
		pcoordinates <- list("unknown" = pcoordinates)
	}
	rnb.section.dreduction.internal(report, pcoordinates, sample.phenotypes, "Low-dimensional Representation")
}

rnb.section.dreduction.internal <- function(report, pcoordinates, sample.phenotypes, section.title) {

	targets <- names(pcoordinates)
	names(targets) <- 1:length(targets)

	## Define point colors and types based on the phenotypic information
	col.visual <- names(rnb.sample.groups(sample.phenotypes, rnb.getOption("exploratory.columns")))
	use.colors <- (length(col.visual) != 0)
	logger.info(c("Mapping", length(col.visual), "traits to point colors and types"))
	
	setting.names <- list("Location type" = targets,
		"Principal components" = c("1" = "first and second", "2" = "second and third"),
		"Distance" = c("manhattan" = "Manhattan", "euclidean" = "Euclidean"))
	snames <- c("identifiers", "circles")
	if (use.colors) {
		snames <- c(snames, paste("based on", col.visual))
	}
	names(snames) <- 1:length(snames)
	setting.names[["Sample representation"]] <- snames
	if (use.colors) {
		snames <- c("black", paste("based on", col.visual))
		names(snames) <- 0:length(col.visual)
		setting.names[["Sample color"]] <- snames
	}
	rm(snames)
	
	create.scatters <- function(dpoints, fprefix) {
		report.plots <- list()
		dframe <- data.frame(x = dpoints[, 1], y = dpoints[, 2], id = rownames(dpoints))
		for (ptype in 1:(length(col.visual) + 2)) {
			plot.labs <- list("x" = colnames(dpoints)[1], "y" = colnames(dpoints)[2])
			if (ptype > 2) {
				sample.categories <- length(unique(sample.phenotypes[, col.visual[ptype - 2]]))
				dframe[["pvalues"]] <- as.factor(sample.phenotypes[, col.visual[ptype - 2]])
				plot.labs[["shape"]] <- col.visual[ptype - 2]
			}
			for (i in 0:length(col.visual)) {
				fname <- paste(fprefix, ptype, sep = "_")
				if (use.colors) {
					fname <- paste(fname, i, sep = "_")
				}
				rplot <- createReportPlot(fname, report, width = 5.4 + 2.6 * use.colors, height = 5.4)
				if (i != 0) {
					sample.categories <- length(unique(sample.phenotypes[, col.visual[i]]))
					dframe[["color"]] <- as.factor(sample.phenotypes[, col.visual[i]])
					cvalues <- rep(rnb.getOption("colors.category"), length.out = sample.categories)
					plot.labs[["color"]] <- col.visual[i]
				}
				if (i != 0) {
					pp <- ggplot(dframe, aes_string(x = "x", y = "y", color = "color", label = "id")) +
							scale_colour_manual(na.value = "#C0C0C0", values = cvalues)
				} else {
					pp <- ggplot(dframe, aes_string(x = "x", y = "y", label = "id"))
				}
				pp <- pp + do.call(labs, plot.labs)
				if (ptype == 1) { # plot sample IDs
					pp <- pp + geom_text(size = 2)
				} else if (ptype == 2) { # all samples are circles
					pp <- pp + geom_point()
				} else {
					ptvalues <- rep(rnb.getOption("points.category"), length.out = nlevels(dframe[, "pvalues"]))
					if (ptype - 2 == i) {
						pp <- pp + geom_point(aes_string(shape = "color"))
					} else {
						pp <- pp + geom_point(aes_string(shape = "pvalues"))
					}
					pp <- pp + scale_shape_manual(na.value = 1L, values = ptvalues)
				}
				if (use.colors) {
					pp <- pp + theme(plot.margin = unit(0.1 + c(0, 2.6, 0, 0), "in"),
							legend.justification = c(0, 0.5), legend.position = c(1, 0.5))
				}
				suppressWarnings(print(pp)) # "Removed ... rows containing missing values (geom_point)"
				report.plots <- c(report.plots, off(rplot))
			}
		}
		return(report.plots)
	}

	## Start the section on dimension reduction
	stext <- c("Dimension reduction is used to visually inspect the dataset for a strong signal in the methylation ",
		"values that is related to samples\' clinical or batch processing annotation. RnBeads implements two methods ",
		"for dimension reduction - principal component analysis (PCA) and multidimensional scaling (MDS).")
	if (is.null(section.title)) {
		rnb.add.paragraph(report, stext)
	} else {
		report <- rnb.add.section(report, "Low-dimensional Representation", stext)
	}

	## Mention that only some sites and/or regions are selected
	if (!is.null(attr(pcoordinates, "selected"))) {
		selected <- attr(pcoordinates, "selected")
#		save(selected, file = "selected.RData", compression_level = 9L)
		top.selected <- isTRUE(any(selected$explained != 1))
		if (top.selected || isTRUE(any(selected$na != 0))) {
			if (top.selected) {
				stext <- c("The analyses in the following sections are based on selected sites and/or regions with ",
					"highest variability in methylation across all samples. The following table shows the maximum ",
					"dimensionality and the selected dimensions in each setting (column names <i>Dimensions</i> ",
					"and <i>Selected</i>, respectively).")
			} else {
				stext <- c(ifelse(length(pcoordinates) == 1, "The methylation matrix",
						"One or more of the methylation matrices"), " was augmented before applying the dimension ",
					"reduction techniques because it contains missing values.")
			}
			selected <- data.frame(
				"Sites/regions" = selected$target,
				"Technique" = selected$technique,
				"Dimensions" = selected$all,
				"Missing" = selected$na,
				"Selected" = selected$selected,
				"Variance explained" = sprintf("%.1f", selected$explained * 100),
				check.names = FALSE, stringsAsFactors = FALSE)
			if (!top.selected) {
				selected <- selected[, setdiff(colnames(selected), "Variance explained")]
			}
			if (isTRUE(any(selected[, "Missing"] != 0))) {
				stext <- c(stext, " The column <i>Missing</i> lists the number of dimensions ignored due to ",
					"missing values. In the case of MDS, dimensions are ignored only if they contain missing ",
					"values for all samples. In contrast, sites or regions with missing values in any sample ",
					"are ignored prior to PCA.")
			} else {
				selected <- selected[, setdiff(colnames(selected), "Missing")]
			}
			rnb.add.paragraph(report, stext)
			rnb.add.table(report, selected, row.names = FALSE)
		}
		rm(top.selected)
	}

	## Create scatter plots of all samples in the reduced representation by MDS
	stext <- c("The scatter plot below visualizes the samples transformed into a two-dimensional space using MDS.")
	report <- rnb.add.section(report, "Multidimensional Scaling", stext, level = 2)
	report.plots <- list()

	for (target.id in 1:length(targets)) {
		for (dist.metric in names(pcoordinates[[target.id]][["mds"]])) {
			dpoints <- pcoordinates[[target.id]][["mds"]][[dist.metric]]
			colnames(dpoints) <- paste("MDS dimension", 1:ncol(dpoints))
			fprefix <- paste("scatter", "mds", target.id, dist.metric, sep = "_")
			report.plots <- c(report.plots, create.scatters(dpoints, fprefix))
		}
	}
	stext <- "Scatter plot showing samples after performing Kruskal\'s non-metric mutidimensional scaling."
	sn.indices <- match("Principal components", names(setting.names))
	tryCatch(
		report <- rnb.add.figure(report, stext, report.plots, setting.names[-sn.indices],
			selected.image = 2 + length(col.visual)),
		error = function(err) {
			print(setting.names[-sn.indices])
			print(sapply(report.plots, slot, 'fname'))
		}
	)

	## Create scatter plots of all samples in the first few principle components
	max.component <- 3L
	report <- rnb.add.section(report, "Principal Component Analysis", NULL, level = 2)
	txt <- "Similarly, the figure below shows the values of selected principal components in a scatter plot."
	report.plots <- list()
	for (target.id in 1:length(targets)) {
		dpoints <- pcoordinates[[target.id]][["pca"]]
		if (length(dpoints) != 0) {
			dpoints <- dpoints$x[, 1:max.component]
			colnames(dpoints) <- paste("Principal component", 1:ncol(dpoints))
			for (k in 1:(max.component - 1)) {
				fprefix <- paste("scatter", "pca", target.id, k, sep = "_")
				rplot <- create.scatters(dpoints[, c(k, k + 1)], fprefix)
				report.plots <- c(report.plots, rplot)
			}
		} else if (length(txt) == 1) {
			txt <- c(txt, " Note that PCA plots are not available for all site and region types; in some cases PCA ",
					"was not performed due to too many missing values.")
		}
	}
	if (length(report.plots) == 0) {
		txt <- "PCA was skipped in all studied site and/or region types due to too many missing values."
		rnb.add.paragraph(report, txt)
		return(report)
	}
	rnb.add.paragraph(report, txt)

	txt <- "Scatter plot showing the samples\' coordinates on principal components."
	sn.indices <- match("Distance", names(setting.names))
	report <- rnb.add.figure(report, txt, report.plots, setting.names[-sn.indices])
	suppressWarnings(rm(max.component, txt, report.plots, target.id, dpoints, k, fprefix, rplot, sn.indices))

	## Calculate variance explained by the top principle components
	var.tbl <- NULL
	rplots <- list()
	for (target.id in 1:length(targets)) {
		var.fraction <- pcoordinates[[target.id]][["pca"]]$sdev^2
		if (length(var.fraction) == 0) {
			next
		}
		var.fraction <- 100 * var.fraction / sum(var.fraction)
		var.explained <- cumsum(var.fraction)
		max.comp <- min(which(var.explained >= VAR.EXPLAINED))
		txt <- c("Principal components that explain at least", VAR.EXPLAINED, "% of the total variance:", max.comp)
		logger.info(txt)
		var.percentage <- data.frame(
			"Principal Component" = 1:length(var.fraction),
			"Percentage of Total Variance" = var.fraction,
			"Cumulative Variance Explained" = var.explained,
			check.names = FALSE, stringsAsFactors = FALSE)

		fname <- sprintf("pca_variance_explained_%s.csv", target.id)
		fname.full <- file.path(rnb.get.directory(report, "data", absolute = TRUE), fname)
		utils::write.csv(var.percentage, file = fname.full, row.names = FALSE)
		logger.info(c("Saved percentage of total variance to", fname))

		fname.full <- paste0('<a href="', rnb.get.directory(report, "data"), '/', fname, '">csv</a>')
		var.tbl <- rbind(var.tbl, c(targets[target.id], as.character(max.comp), fname.full))

		if (max.comp == 1) {
			rplot <- list(NULL)
		} else { # max.comp > 1

			## Plot a CDF of variance explained
			dframe <- data.frame(x = 1:max.comp, y = var.explained[1:max.comp])
			ylab <- paste0(ifelse(is.null(selected), "Total v", "V"), "ariance explained (%)")
			fname <- paste0("cdf_pca_variability_regions_", target.id)
			rplot <- createReportPlot(fname, report, width = 6.2, height = 6.2)
			pp <- ggplot(dframe, aes_string(x = "x", y = "y")) + labs(x = "Number of components", y = ylab) +
				geom_line(color = "blue") + geom_point(color = "blue") +
				theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in"))
			print(pp)
			rplot <- off(rplot)
			rm(dframe, ylab, pp)
		}
		rm(var.explained, max.comp, txt, var.percentage, fname, fname.full)
		rplots <- c(rplots, rplot)
		rm(rplot)
	}
	rm(target.id, var.fraction)

	## Create a figure of variance explained
	i.targets <- which(!sapply(rplots, is.null))
	if (length(i.targets) != 0) {
		stext <- c("The figure below shows the cumulative distribution functions of variance explained by the ",
			"principal components.")
		rnb.add.paragraph(report, stext)
		stext <- ifelse(is.null(selected), "total ", "")
		stext <- c("Cumulative distribution function of percentange of ", stext, "variance explained.")
		tryCatch(
			report <- rnb.add.figure(report, stext, rplots[i.targets], setting.names[1L]),
			error = function(err) {
				print(setting.names[1])
				print(sapply(rplots, slot, 'fname'))
			}
		)
	}

	## Add links to the tables of principal component variances
	if (!is.null(var.tbl)) {
		if (nrow(var.tbl) == 1) {
			txt <- c("table", "is", "a ", "file")
		} else { # nrow(var.tbl) > 1
			txt <- c("tables", "are", "", "files")
		}
		txt <- c("The table below gives for each location type a number of principal components that explain at least ",
			VAR.EXPLAINED, " percent of the total variance. The full ", txt[1], " of variances explained by all ",
			"components ", txt[2], " available in ", txt[3], "comma-separated values ", txt[4], " accompanying this ",
			"report.")
		rnb.add.paragraph(report, txt)
		
		colnames(var.tbl) <- c("Location Type", "Number of Components", "Full Table File")
		rnb.add.table(report, var.tbl, row.names = FALSE)
	}

	return(report)
}

########################################################################################################################

#' rnb.step.dreduction
#'
#' Executes the dimension reduction step and adds a dedicated section to the given report.
#'
#' @param rnb.set            Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report             Report to contain the dimension reduction section. This must be an object of type
#'                           \code{\linkS4class{Report}}.
#' @param return.coordinates Flag indicating if the computed sample coordinates in the prinicipal component space should
#'                           also be returned.
#' @return If \code{return.coordinates} is \code{FALSE}: the modified report. Otherwise, a list of two elements:
#'         \describe{
#'           \item{report}{The modified report.}
#'           \item{pcoordinates}{Coordinates of the samples in the dataset, as returned by
#'                \code{\link{rnb.execute.dreduction}}.}
#'         }
#'
#' @seealso \code{\link{rnb.execute.dreduction}}, \code{\link{rnb.section.dreduction}}
#' @author Yassen Assenov
#' @noRd
rnb.step.dreduction <- function(rnb.set, report, return.coordinates = FALSE) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnbSet")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (!parameter.is.flag(return.coordinates)) {
		stop("invalid value for return.coordinates; expected TRUE or FALSE")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}

	## Extract list of site and region types to analyze
	section.title <- "Low-dimensional Representation"
	targets <- rnb.step.analyze.targets(rnb.set, report, section.title)
	if (length(targets) == 0) {
		if (return.coordinates) {
			return(list(report = report, pcoordinates = list()))
		}
		return(report)
	}
	if (!is.null(attr(targets, "warning"))) {
		txt <- attr(targets, "warning")
		report <- rnb.add.section(report, section.title, txt)
		section.title <- NULL
	}

	## Map the samples to a low-dimensional space
	logger.start("Dimension Reduction Techniques")
	names(targets) <- targets
	if (parallel.isEnabled()) {
		pcoordinates <- foreach(target = targets) %dopar%
			tryCatch(RnBeads::rnb.execute.dreduction(rnb.set, target = target), error = function(e) { e$message } )
		names(pcoordinates) <- names(targets)
	} else {
		pcoordinates <- lapply(targets, function(target) {
			tryCatch(rnb.execute.dreduction(rnb.set, target = target), error = function(e) { e$message } )
		})
	}

	selected <- lapply(pcoordinates, attr, which = "selected")
	selected <- selected[!sapply(selected, is.null)]
	if (length(selected) != 0) {
		attr(pcoordinates, "selected") <- do.call(rbind, selected)
		for (i in 1:length(pcoordinates)) {
			attr(pcoordinates[[i]], "selected") <- NULL
		}
		rm(i)
	}
	rm(selected)

	i.failed <- which(sapply(pcoordinates, is.character))
	if (length(i.failed) != 0) {
		if (length(i.failed) == length(pcoordinates)) {
			if (identical(pcoordinates[[1]], "too few samples")) {
				txt <- "Dimension reduction techniques were not applied because the dataset contains too few samples."
				logger.warning("Skipped due to too few samples")
			} else { # too many sites with missing values
				txt <- c("Dimension reduction techniques were not applied.", "")
				logger.warning("Skipped due to too many missing values")
			}
		} else {
			txt <- paste0("<b>", names(pcoordinates)[i.failed], "</b>", collapse = ", ")
			txt <- c("Dimension reduction techniques were not applied for the following types: ", txt, ".")
		}
		if (length(txt) != 1) {
			txt <- c(txt, " The most likely reason for this error is inability to calculate a distance between a pair ",
				"of samples due to (too many) missing values in the respective methylation ",
				ifelse(length(i.failed) == 1, "matrix", "matrices"), ".")
		}
		if (is.null(section.title)) {
			rnb.add.paragraph(report, txt)
		} else {
			report <- rnb.add.section(report, section.title, txt)
			section.title <- NULL
		}
		pcoordinates <- pcoordinates[setdiff(1:length(pcoordinates), i.failed)]
	}

	if (length(pcoordinates) != 0) {
		## Create plots of PCA and MDS
		report <- rnb.section.dreduction.internal(report, pcoordinates, pheno(rnb.set), section.title)
		logger.status("Created scatter plots and CDFs summarizing the reduced dimensional representations")
	}
	logger.completed()
	if (return.coordinates) {
		return(list(report = report, pcoordinates = pcoordinates))
	}
	return(report)
	
}

########################################################################################################################

#' rnb.execute.batcheffects
#'
#' Performs tests for association between traits and principal components.
#'
#' @param rnb.set      Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param pcoordinates Coordinates of the samples of \code{rnb.set} in the principal components space, as returned by
#'                     \code{\link{rnb.execute.dreduction}}.
#' @return Results of attempted tests for associations in the form of a list with up to three elements:
#'         \describe{
#'           \item{\code{"permutations"}}{\code{integer} matrix of index permutations. The number of rows in the matrix
#'                is \emph{N} - the number of samples in \code{rnb.set}. Every column in this matrix denotes a sample
#'                permutation; the first column is the sequence 1 to \emph{N}. This element is included only when
#'                \code{rnb.getOption("exploratory.correlation.permutations")} is non-zero and there are numeric traits
#'                to be tested.}
#'           \item{\code{"pc"}}{List of four matrices named \code{"failures"}, \code{"tests"}, \code{"correlations"}
#'                and \code{"pvalues"}. The rows in each of these matrices correspond to the first several principal
#'                components, and the columns - to selected traits. This element is not included in the returned list
#'                when \code{pcoordinates} is \code{NULL}.}
#'           \item{\code{"traits"}}{List of four square symmetric matrices named \code{"failures"}, \code{"tests"},
#'                \code{"correlations"} and \code{"pvalues"}, containing information about the performed tests for
#'                pairwise trait association. This element is included only if two or more traits were tested.}
#'         }
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' regs <- c("sites", summarized.regions(rnb.set.example))
#' dreduction <- function(x) rnb.execute.dreduction(rnb.set.example, x)
#' pcoordinates <- lapply(regs, dreduction)
#' names(pcoordinates) <- regs
#' result <- rnb.execute.batcheffects(rnb.set.example, pcoordinates)
#' }
#' @seealso \code{\link{rnb.run.exploratory}} for running the whole exploratory analysis module
#' @author Yassen Assenov
#' @export
rnb.execute.batcheffects <- function(rnb.set, pcoordinates = NULL) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	pheno.table <- pheno(rnb.set)
	if (is.null(pheno.table)) {
		stop("invalid value for rnb.set; missing phenotype data")
	}
	if (!is.null(pcoordinates)) {
		nsamples <- validate.pcoordinates.all(pcoordinates)
		if (is.character(nsamples)) {
			stop(nsamples)
		}
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	logger.start("Tests for Associations")
	
	## Get a list of comparable traits
	predefined.columns <- rnb.getOption("exploratory.columns")
	if (!is.null(predefined.columns)) {
		if (is.character(predefined.columns)) {
			predefined.columns <- intersect(predefined.columns, colnames(pheno.table))
		} else { # is.integer(columns)
			predefined.columns <- intersect(predefined.columns, 1:ncol(pheno.table))
		}
		pheno.table <- pheno.table[, predefined.columns]
	}
	traits <- lapply(pheno.table, function(x) {
				if (length(unique(na.omit(x))) < 2) {
					return(NULL)
				}
				if (is.character(x)) {
					x <- as.factor(x)
				}
				if (is.logical(x)) {
					x <- as.factor(x)
				}
				if (is.factor(x)) {
					y <- na.omit(x)
					if (anyDuplicated(y) == 0) {
						return(NULL)
					}
				}
				return(x)
			}
	)
	traits <- traits[sapply(traits, is.null) == FALSE]
	NT <- length(traits)
	if (NT == 0) {
		logger.info("No suitable traits found")
		logger.completed()
		return(NULL)
	}
	logger.info(c("Testing the following traits for associations:", paste(names(traits), collapse = "; ")))
	
	result <- list()
	
	## Create sample permutations if necessary
	perm.matrix <- NULL
	perm.count <- rnb.getOption("exploratory.correlation.permutations")
	if ((!is.null(pcoordinates)) && perm.count != 0 && sum(!sapply(traits, is.factor)) >= 2) {
		perm.matrix <- mapply(sample, rep(nrow(pheno.table), times = perm.count))
		perm.matrix[, 1] <- 1:nrow(perm.matrix)
		result[["permutations"]] <- perm.matrix
		logger.status(c("Created", ncol(perm.matrix), "sample permutations"))
	}

	pc.association.count <- rnb.getOption("exploratory.principal.components")
	if ((!is.null(pcoordinates)) && pc.association.count != 0) {
		## Compute correlations between principle components and traits
		result[["pc"]] <- list()
		for (target in names(pcoordinates)) {
			dpoints <- pcoordinates[[target]]$pca$x
			if(!is.null(dpoints)){
				if (ncol(dpoints) > pc.association.count) {
					dpoints <- dpoints[, 1:pc.association.count]
				}
	
				init.matrix <- function(x) {
					matrix(x, nrow = NT, ncol = ncol(dpoints),
							dimnames = list(names(traits), "Principal component" = 1:ncol(dpoints)))
				}
				table.test.failures <- init.matrix(as.character(NA))
				table.test.names <- init.matrix(as.character(NA))
				table.correlations <- init.matrix(as.double(NA))
				table.pvalues <- init.matrix(as.double(NA))
	
				for (i in 1:NT) {
					for (j in 1:ncol(dpoints)) {
						t.result <- test.traits(traits[[i]], dpoints[, j], perm.matrix)
						table.test.failures[i, j] <- t.result[["error"]]
						table.test.names[i, j] <- t.result[["test"]]
						table.correlations[i, j] <- t.result[["correlation"]]
						table.pvalues[i, j] <- t.result[["pvalue"]]
					}
				}
				result[["pc"]][[target]] <- list(
						"failures" = table.test.failures,
						"tests" = table.test.names,
						"correlations" = table.correlations,
						"pvalues" = table.pvalues)
			}
		}
		rm(dpoints, init.matrix, table.test.failures, table.test.names, table.correlations, table.pvalues)
		rm(i, j, t.result)
		logger.status("Computed correlations between principal components and traits.")
	}
	
	## Compute correlations and perform tests
	if (NT > 1) {
		init.matrix <- function(x) {
			matrix(x, nrow = NT, ncol = NT, dimnames = list(names(traits), names(traits)))
		}
		table.test.failures <- init.matrix(as.character(NA))
		table.test.names <- init.matrix(as.character(NA))
		table.correlations <- init.matrix(as.double(NA))
		table.pvalues <- init.matrix(as.double(NA))
		for (i in 1:(NT - 1)) {
			for (j in (i + 1):NT) {
				t.result <- test.traits(traits[[i]], traits[[j]], perm.matrix)
				table.test.failures[i, j] <- table.test.failures[j, i] <- t.result[["error"]]
				table.test.names[i, j] <- table.test.names[j, i] <- t.result[["test"]]
				table.correlations[i, j] <- table.correlations[j, i] <- t.result[["correlation"]]
				table.pvalues[i, j] <- table.pvalues[j, i] <- t.result[["pvalue"]]
			}
		}
		rm(perm.matrix, NT, i, j, t.result)
		result[["traits"]] <- list(
				"failures" = table.test.failures,
				"tests" = table.test.names,
				"correlations" = table.correlations,
				"pvalues" = table.pvalues)
		logger.status("Computed pairwise correlations between traits.")
	}
	logger.completed()
	return(result)
}

########################################################################################################################

#' rnb.section.batcheffects
#'
#' Creates a report section summarizing associations between traits and principal components, as well as trait pairwise
#' associations.
#'
#' @param report       Report to contain the section. This must be an object of type \code{\linkS4class{Report}}.
#' @param batcheffects Results of attempted tests for associations, as returned by
#'                     \code{\link{rnb.execute.batcheffects}}.
#' @return The modified report.
#'
#' @seealso \code{\link{rnb.step.batcheffects}}, \code{\link{rnb.run.exploratory}}
#' @author Yassen Assenov
#' @noRd
rnb.section.batcheffects <- function(report, batcheffects) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (!is.null(batcheffects)) {
		if (!(is.list(batcheffects) && all(names(batcheffects) %in% c("permutations", "pc", "traits")))) {
			stop("invalid value for batcheffects")
		}
		if ("permutations" %in% names(batcheffects)) {
			be <- batcheffects[names(batcheffects) != "permutations"]
			if (!is.valid.permutations(batcheffects$permutations)) {
				stop("invalid value for batcheffects$permutations")
			}
		} else {
			be <- batcheffects
		}
		check.batcheffect.matrices <- function(x) {
			is.list(x) && identical(names(x), c("failures", "tests", "correlations", "pvalues")) &&
				all(sapply(x, is.matrix)) && all(sapply(x[1:2], is.character)) && all(sapply(x[3:4], is.numeric))
		}
		if (!all(sapply(be[names(be) != "pc"], check.batcheffect.matrices))) {
			stop("invalid value for batcheffects")
		}
		if (!all(sapply(be[["pc"]], check.batcheffect.matrices))) {
			stop("invalid value for batcheffects")
		}
		rm(be, check.batcheffect.matrices)
	} else {
		txt <- c("Batch effects were not studied because none of the traits can be tested for association with sample ",
			"coordinates in the principal components space.")
		report <- rnb.add.section(report, "Batch Effects", txt)
		return(report)
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	
	txt <- c("In this section, different properties of the dataset are tested for significant associations. The ",
			"properties can include sample coordinates in the principal component space, phenotype traits and ",
			"intensities of control probes. The tests used to calculate a p-value given two properties depend on the ",
			"essence of the data:")
	report <- rnb.add.section(report, "Batch Effects", txt)
	nperm <- rnb.getOption("exploratory.correlation.permutations")
	txt <- list(
		c("If both properties contain categorical data (e.g. tissue type and sample processing date), the test of ",
				"choice is a two-sided Fisher's exact test."),
		c("If both properties contain numerical data (e.g. coordinates in the first principal component and age of ",
				"individual), the correlation coefficient between the traits is computed. ",
				ifelse(nperm > 0, paste("A p-value is estimated using permutation tests with", nperm, "permutations."),
						"No p-value is estimated in this case because permutation tests are disabled.")),
		c("If property <i>A</i> is categorical and property <i>B</i> contains numeric data, p-value for association ",
				"is calculated by comparing the values of <i>B</i> for the different categories in <i>A</i>. The test ",
				"of choice is a two-sided Wilcoxon rank sum test (when <i>A</i> defines two categories) or ",
				"a Kruskal-Wallis one-way analysis of variance (when <i>A</i> separates the samples into three or more ",
				"categories)."))
	rnb.add.list(report, txt)
	txt <- "Note that the p-values presented in this report are <em>not corrected</em> for multiple testing."
	rnb.add.paragraph(report, txt)
	data.dir.abs <- rnb.get.directory(report, "data", TRUE)
	data.dir.rel <- rnb.get.directory(report, "data")

	## -----------------------------------------------------------------------------------------------------------------
	## Summarize the results of association between principal components and traits

	if ("pc" %in% names(batcheffects)) {
		
		targets <- names(batcheffects[["pc"]])
		names(targets) <- 1:length(targets)
		setting.names <- list("Region type" = targets)

		txt <- ifelse(is.null(rnb.getOption("exploratory.columns")), "available", "specified")
		txt <- c("The computed sample coordinates in the principal component space were tested for association with ",
			"the ", txt, " traits. Below is a list of the traits and the tests performed.")
		report <- rnb.add.section(report, "Associations between Principal Components and Traits", txt, level = 2)

		## Display summary table of tests performed
		pc.tables <- batcheffects$pc[[1]]
		tbl <- matrix(c(rownames(pc.tables$tests), pc.tables$tests[, 1]), ncol = 2)
		colnames(tbl) <- c("Trait", "Test")
		rnb.add.table(report, tbl)

		append.table <- function(tbl, fnames, heatmap.fun) {
			mytbl <- cbind("Trait \\ Principal Component" = rownames(tbl), as.data.frame(tbl, check.names = FALSE))
			fnames.full <- file.path(data.dir.abs, fnames)
			fnames[1] <- paste(data.dir.rel, fnames[1], sep = "/")
			utils::write.csv(mytbl, file = fnames.full[1], na = "", row.names = FALSE)
			sprintf("<a href=\"%s\">csv</a>", fnames[1])
			hmap <- heatmap.fun(report, tbl, fnames[2])
			list(link = sprintf("<a href=\"%s\">csv</a>", fnames[1]), plot = hmap$plot, description = hmap$description)
		}

		## Create heatmaps of correlations and save the values to CSV files
		file.tbl <- matrix(character(), nrow = 0, ncol = 2, dimnames = list(NULL, c("Location Type", "File Name")))
		rplots <- list()
		for (target.id in 1:length(targets)) {
			tbl <- batcheffects$pc[[target.id]]$correlations
			if (!all(is.na(tbl))) {
				fnames <- paste0(c("correlations_pc_trait_", "heatmap_pc_correlations_"), target.id, c(".csv", ""))
				result <- append.table(tbl, fnames, plot.heatmap.pc.correlations)
				file.tbl <- rbind(file.tbl, c(targets[target.id], result$link))
				rplots[[targets[target.id]]] <- result$plot
				description <- result$description
			}
		}
		tbl.correlations.present <- (length(rplots) != 0)
		if (tbl.correlations.present) {
			txt <- c("The next figure shows the computed correlations between the first ",
				rnb.getOption("exploratory.principal.components"), " principal components and the sample traits.")
			rnb.add.paragraph(report, txt)
			report <- rnb.add.figure(report, description, rplots, setting.names)
			colnames(file.tbl) <- c("Location type", "Table file")
			rnb.add.exported.tables(report, file.tbl)
		}

		## Create heatmaps of p-values and save the values to CSV files
		file.tbl <- matrix(character(), nrow = 0, ncol = 2, dimnames = list(NULL, c("Location Type", "File Name")))
		rplots <- list()
		for (target.id in 1:length(targets)) {
			pc.tables <- batcheffects$pc[[target.id]]
			tbl <- pc.tables$pvalues
			if (!all(is.na(tbl))) {
				fnames <- paste0(c("pvalues_pc_trait_", "heatmap_pc_pvalues_"), target.id, c(".csv", ""))
				result <- append.table(tbl, fnames, plot.heatmap.pc.pvalues)
				file.tbl <- rbind(file.tbl, c(targets[target.id], result$link))
				rplots[[targets[target.id]]] <- result$plot
				description <- result$description
			}
		}
		tbl.pvalues.present <- (length(rplots) != 0)
		if (nperm == 0) {
			if (tbl.correlations.present) {
				txt <- c("Since permutation tests were disabled, no p-values were estimated from the correlations ",
					"shown above.")
				rnb.add.paragraph(report, txt)
			}
		} else if (tbl.pvalues.present) {
			txt <- c("The heatmap below summarizes the results of permutation tests performed for associations. ",
				"Significant p-values (values less than ", rnb.getOption("exploratory.correlation.pvalue.threshold"),
				") are displayed in pink background.")
			rnb.add.paragraph(report, txt)
			report <- rnb.add.figure(report, description, rplots, setting.names)
			txt <- c("The full tables of p-values for each location type are available in CSV (comma-separated ",
				"value) files below.")
			rnb.add.paragraph(report, txt)
			rnb.add.table(report, file.tbl)
		} else {
			## TODO: All p-value tables are full of NAs
		}

		rm(pc.tables)
	}

	## -----------------------------------------------------------------------------------------------------------------
	## Summarize the results of trait association

	if ("traits" %in% names(batcheffects)) {
		trait.tables <- batcheffects$traits
		txt <- c("This section summarizes the associations between pairs of traits.")
		report <- rnb.add.section(report, "Associations between Traits", txt, level = 2)

		## Create a triangular heatmap of the performed tests for associations
		plots.associations <- list(NULL, NULL)
		txt <- c("The figure below visualizes the tests that were performed on trait pairs based on the description ",
			"provided above.")
		tbl <- trait.tables$tests
		tbl.failures <- trait.tables$failures
		if (!all(is.na(tbl.failures))) {
			txt <- c(txt, " In some cases, pairs of traits could not be tested for associations. These scenarios are ",
				"marked by grey shapes, and the underlying reason is given in the figure legend.")
		} else {
			tbl.failures <- NULL
		}
		plots.associations[[1]] <- plot.heatmap.symm(report, tbl, tbl.failures, "heatmap_traits_tests")

		## Create a triangular heatmap of the calculated p-values for associations
		tbl <- trait.tables$pvalues
		plots.associations[[2]] <- plot.heatmap.symm(report, tbl, NULL, "heatmap_traits_pvalues")
		
		## Attach the table of the calculated p-values for associations
		tbl <- cbind("Traits" = rownames(tbl), as.data.frame(tbl, check.names = FALSE))
		fname <- "pvalues_traits.csv"
		utils::write.csv(tbl, file = file.path(data.dir.abs, fname), na = "", row.names = FALSE)

		## Add a figure of the triangular heatmaps
		txt <- c(txt, " In addition, the calculated p-values for associations between traits are shown. ",
			"Significant p-values (values less than ", rnb.getOption("exploratory.correlation.pvalue.threshold"),
			") are displayed in pink background. The full table of p-values is available in a <a href=\"",
			data.dir.rel, "/", fname, "\">dedicated file</a> that accompanies this report.")
		rnb.add.paragraph(report, txt)
		txt <- c("(1) Table of performed tests on pairs of traits. Test names (Correlation + permutation test, ",
			"Fisher's exact test, Wilcoxon rank sum test and/or Kruskal-Wallis one-way analysis of variance) are ",
			"color-coded according to the legend given above.<br />")
		txt <- c(txt, "(2) Table of resulting p-values from the performed tests on pairs of traits. Significant ",
			"p-values (less than ", rnb.getOption("exploratory.correlation.pvalue.threshold"), ") are printed in pink boxes ",
			"Non-significant values are represented by blue boxes. White cells, if present, denote missing values.")
		setting.names <- list("Heatmap of" = c("tests" = "tests performed", "pvalues" = "p-values")) 
		report <- rnb.add.figure(report, txt, plots.associations, setting.names,  selected.image = 2)
		rm(plots.associations, txt, tbl, tbl.failures, fname, setting.names)
		
		if (!all(is.na(trait.tables$correlations))) {
			## Attach the table of correlations to the report
			txt <- "In some cases, a correlation was computed between a pair of traits."
			if (nperm > 0) {
				txt <- c(txt, " As described earlier in this report, these correlation values are used as the basis ",
					"for a permutation-based test.")
			}
			tbl <- trait.tables$correlations
			tbl <- cbind("Traits" = rownames(tbl), as.data.frame(tbl, check.names = FALSE))
			fname <- "correlations_traits.csv"
			utils::write.csv(tbl, file = file.path(data.dir.abs, fname), na = "", row.names = FALSE)
			txt <- c(txt, " The table of computed correlations is available as a <a href=\"",
				data.dir.rel, "/", fname, "\">comma-separated file</a> accompanying this report.")
			rnb.add.paragraph(report, txt)
		}
	}
	
	return(report)
}

########################################################################################################################

#' rnb.step.batcheffects
#'
#' Executes the batch effects step and adds a dedicated section to the given report.
#'
#' @param rnb.set             Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report              Report to contain the dimension reduction section. This must be an object of type
#'                            \code{\linkS4class{Report}}.
#' @param pcoordinates        Coordinates of the samples of \code{rnb.set} in the principal components space, as returned
#'                            by \code{\link{rnb.execute.dreduction}}.
#' @param return.permutations Flag indicating if the generated sample permutations should also be returned.
#' @return If \code{return.permutations} is \code{FALSE}: the modified report. Otherwise, a list of two elements:
#'         \describe{
#'           \item{report}{The modified report.}
#'           \item{permutations}{\code{integer} matrix of sample index permutations. See the return value of
#'                \code{\link{rnb.execute.batcheffects}} for more information. This element is \code{NULL} if no
#'                permutations were generated.}
#'         }
#'
#' @seealso \code{\link{rnb.execute.batcheffects}}, \code{\link{rnb.section.batcheffects}},
#'   \code{\link{rnb.run.exploratory}}
#' @author Yassen Assenov
#' @noRd
rnb.step.batcheffects <- function(rnb.set, report, pcoordinates, return.permutations = FALSE) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	nsamples <- validate.pcoordinates.all(pcoordinates)
	if (is.character(nsamples)) {
		stop(nsamples)
	}
	rm(nsamples)
	if (!parameter.is.flag(return.permutations)) {
		stop("invalid value for return.permutations; expected TRUE or FALSE")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	
	if (is.null(pheno(rnb.set))) {
		txt <- "Batch effects were not studied because the dataset contains no phenotype data."
		report <- rnb.add.section(report, "Batch Effects", txt)
		return(report)
	}

	batcheffects <- rnb.execute.batcheffects(rnb.set, pcoordinates)
	report <- rnb.section.batcheffects(report, batcheffects)

	if (return.permutations) {
		return(list(report = report, permutations = batcheffects$permutations))
	}
	return(report)
}
