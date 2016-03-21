########################################################################################################################
## gender.R
## created: 2014-02-28
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions to infer gender from Infinium 450k data, and to visualize the results.
########################################################################################################################

## G L O B A L S #######################################################################################################

## Columns added to the sample annotation table when gender is predicted
RNB.COLUMNS.PREDICTED.GENDER <- c("Predicted Male Probability", "Predicted Gender")

## F U N C T I O N S ###################################################################################################

#' rnb.contains.sex
#'
#' Checks if the given dataset contains sites on sex chromosomes.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{RnBSet}.
#' @return \code{TRUE} if the dataset contains sites on the X or Y chromosome; \code{FALSE} otherwise.
#' @author Yassen Assenov
#' @noRd
rnb.contains.sex <- function(rnb.set) {
	ci.sex <- which(rnb.get.chromosomes(rnb.set@assembly) %in% c("X", "Y"))
	any(rnb.set@sites[, 2] %in% ci.sex)
}

########################################################################################################################

#' rnb.get.XY.shifts
#'
#' Calculates the increase of signal on the sex chromosomes in the given dataset.
#'
#' @param rnb.set     Dataset of interest. Currently only \code{\linkS4class{RnBeadRawSet}} is supported.
#' @param signal.type Matrix or matrices to use in calculating the shifts. Currently ignored.
#' @return Calculated shifts in the form of a matrix in which every row corresponds to a sample in \code{rnb.set} and
#'         the columns denote the shifts in the X and Y chromosomes.
#'
#' @author Yassen Assenov
#' @noRd
rnb.get.XY.shifts <- function(rnb.set, signal.type = "raw") {

	## Identify the indices of sites on the sex chromosomes and autosomes
	ci.X <- which(rnb.get.chromosomes(rnb.set@assembly) == "X")
	ci.Y <- which(rnb.get.chromosomes(rnb.set@assembly) == "Y")
	site.chroms <- rnb.set@sites[, 2]
	chrom2probes <- list(
		X = as.integer(which(site.chroms == ci.X)),
		Y = as.integer(which(site.chroms == ci.Y)))
	chrom2probes$autosome <- setdiff(1:length(site.chroms), unlist(chrom2probes, use.names = FALSE))

	## Validate that not too many sites are missing
	probes.max <- rnb.annotation.size(rnb.set@target, assembly = rnb.set@assembly)
	probes.max <- c(X = sum(probes.max[ci.X]), Y = sum(probes.max[ci.Y]), autosome = sum(probes.max[-c(ci.X, ci.Y)]))
	fr.available <- sapply(chrom2probes, length) / probes.max
	if (any(fr.available < 0.2)) {
		return(NULL)
	}
	rm(ci.X, ci.Y, site.chroms)
	rm(probes.max, fr.available)

	## Calculate signal shifts
	shifts <- t(apply(rnb.set@M[, , drop = FALSE] + rnb.set@U[, , drop = FALSE], 2, function(x) {
			t.signals <- sapply(chrom2probes, function(i) { mean(x[i], na.rm = TRUE) })
			c(t.signals[1], t.signals[2]) - t.signals[3]
		})
	)
	return(shifts / 1000)
}

########################################################################################################################

#' rnb.set.update.predicted.gender
#'
#' Adds two columns (RNB.COLUMNS.PREDICTED.GENDER) to the sample annotation of the given methylation dataset, based
#' on the calculated methylation shifts. Also sets the gender inferred covariate to \code{TRUE}.
#'
#' @param rnb.set         Dataset of interest. Currently only RnBeadRawSet is supported.
#' @param shifts          Matrix of calculated mean signal increases, as returned by \code{\link{rnb.get.XY.shifts}}.
#' @param pr.coefficients Coefficients of the logistic regression model used for the prediction of gender.
#' @return The possibly modified dataset with two columns added to its sample annotation table. If \code{shifts} is
#'         \code{NULL}, the returned dataset is \code{rnb.set}.
#'
#' @author Yassen Assenov
#' @noRd
rnb.set.update.predicted.gender <- function(rnb.set, shifts, pr.coefficients = c(2, -1)) {
	if (!is.null(shifts)) {
		male.probabilities <- as.vector(1 / (1 + exp(shifts %*% pr.coefficients)))
		rnb.set@pheno[, RNB.COLUMNS.PREDICTED.GENDER[1]] <- male.probabilities
		p.genders <- factor(ifelse(male.probabilities > 0.5, "male", "female"), levels = c("female", "male", "unknown"))
		p.genders[is.na(p.genders)] <- "unknown"
		rnb.set@pheno[, RNB.COLUMNS.PREDICTED.GENDER[2]] <- p.genders
		rnb.set@inferred.covariates$gender <- TRUE
	}
	rnb.set
}

########################################################################################################################

#' rnb.execute.gender.prediction
#'
#' Infers the gender of every sample in the given Infinium 450k dataset, based on average signal intensity values on
#' the autosomes and the sex chromosomes.
#'
#' @param rnb.set Methylation dataset as an object of type \code{\linkS4class{RnBeadRawSet}}.
#' @return The possibly modified dataset. If gender could be predicted, the sample annotation table is enriched with
#' two more columns - \code{"Predicted Male Probability"} and \code{"Predicted Gender"}.
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.set.example <- rnb.execute.gender.prediction(rnb.set.example)
#' table(rnb.set.example[, "Predicted Gender"])
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.gender.prediction <- function(rnb.set) {
	if (!inherits(rnb.set, "RnBeadRawSet")) {
		stop("invalid value for rnb.set")
	}
	if (rnb.set@target != "probes450") {
		stop("unsupported platform")
	}

	shifts <- rnb.get.XY.shifts(rnb.set)
	if (is.null(shifts)) {
		rnb.warning("Prediction skipped due to incomplete data")
	}
	rnb.set.update.predicted.gender(rnb.set, shifts)
}

########################################################################################################################

#' rnb.section.gender.prediction
#'
#' Adds a dedicated section on gender prediction results to the given report.
#'
#' @param rnb.set Methylation dataset after running the gender prediction step, as an object of type 
#'                \code{\linkS4class{RnBSet}}.
#' @param shifts  Matrix of calculated mean signal increases, as returned by \code{\link{rnb.get.XY.shifts}}.
#' @param report  Report on annotation inferrence to contain the gender prediction section. This must be an object of
#'                type \code{\linkS4class{Report}}.
#' @return The modified report.
#'
#' @seealso \code{\link{rnb.execute.gender.prediction}} for performing gender prediction
#' @author Yassen Assenov
#' @noRd
rnb.section.gender.prediction <- function(rnb.set, shifts, report) {
	if (inherits(rnb.set, "RnBeadRawSet")) {
		if (is.null(shifts)) {
			if (rnb.contains.sex(rnb.set)) {
				txt <- c("Gender prediction was not performed because the dataset contains too few sites. Please note ",
					"that gender prediction is started when the dataset contains reliable measurements for at least ",
					"20% of the sites on autosomes and on each sex chromosome.")
			} else {
				txt <- c("Gender prediction was not performed because the dataset contains no sites on sex ",
					"chromosomes. The prediction is based on signal intensities on the sex chromosomes relative to ",
					"the ones on autosomes.")
				if (rnb.getOption("filtering.sex.chromosomes.removal")) {
					txt <- c(txt, " Please disable sex chromosomes removal (analysis option ",
						"<code>filtering.sex.chromosomes.removal</code>) in order to enable gender prediction.")
				}
			}
		} else {
			pred.genders <- rnb.set@pheno[, RNB.COLUMNS.PREDICTED.GENDER[2]]
			colors.gender <- c(muted("pink"), muted("blue"), "#808080")
			names(colors.gender) <- levels(pred.genders)
			pred.genders <- table(pred.genders)
			txt <- c("RnBeads predicted the gender of the samples in the dataset using a logistic regression model. ",
				"The results are summarized in the table below.")
			pred.genders <- data.frame("Gender" = names(pred.genders), "Samples" = as.integer(pred.genders),
				check.names = FALSE, stringsAsFactors = FALSE)
		}
	} else {
		txt <- c("Gender prediction is skipped because this operation is not supported for the dataset. Currently, ",
			"gender prediction can be performed only on raw Infinium 450k datasets. These are, for example, datasets ",
			"loaded from IDAT files.")
		return(report)
	}
	report <- rnb.add.section(report, "Gender Prediction", txt)

	if (exists("colors.gender", inherits = FALSE)) {
		rnb.add.table(report, pred.genders)

		## Display the shifts and the separation line
		txt <- c("Gender was predicted based on the increase (or decrease) of mean signal intensities in the sex ",
			"chromosomes w.r.t. the corresponding value in autosomes. The figure below displays these characteristics ",
			"of the samples.")
		rnb.add.paragraph(report, txt)
		s.colorings <- c("prob" = "predicted male probability", "gend" = "predicted gender")
		dframe <- data.frame(
			prob = rnb.set@pheno[, RNB.COLUMNS.PREDICTED.GENDER[1]],
			gend = rnb.set@pheno[, RNB.COLUMNS.PREDICTED.GENDER[2]], X = shifts[, 1], Y = shifts[, 2])
		rplots <- list()
		for (s.coloring in names(s.colorings)) {
			a.labels <- paste("Mean signal increase in the", c("X", "Y"), "chromosome")
			c.label <- gsub("predicted ", "", s.colorings[s.coloring], fixed = TRUE)
			pp <- ggplot2::ggplot(dframe, aes_string(x = 'X', y = 'Y', color = s.coloring)) + ggplot2::geom_point() +
				ggplot2::coord_fixed() + ggplot2::labs(x = a.labels[1], y = a.labels[2], color = c.label)
			if (s.coloring == "prob") {
				pp <- pp + ggplot2::scale_color_gradient2(limits = c(0, 1), low = colors.gender[1], mid = "white",
						high = colors.gender[2], midpoint = 0.5, na.value = colors.gender[3])
			} else { # s.coloring == "gend"
				pp <- pp + ggplot2::scale_color_manual(na.value = colors.gender[3], values = colors.gender)
			}
			pp <- pp + ggplot2::geom_abline(intercept = 0, slope = -(-2) / 1) +
				ggplot2::theme(plot.margin = unit(0.1 + c(0, 1, 0, 0), "in")) +
				ggplot2::theme(legend.position = c(1, 0.5), legend.justification = c(0, 0.5))
			fname <- paste0("gender_prediction_signals_", s.coloring)
			rplot <- createReportPlot(fname, report, width = 8.2, height = 7.2)
			print(pp)
			rplots <- c(rplots, off(rplot))
		}
		txt <- c("Gender prediction based on mean signal increase. The decision boundary between the two genders is ",
			"visualized by a black line. Sample colors denote predicted male probability / gender.")
		report <- rnb.add.figure(report, txt, rplots, list("Colors denote" = s.colorings))
		rm(pred.genders, txt, s.colorings, dframe, rplots, s.coloring, pp, fname, rplot)
	}
	report
}

########################################################################################################################

#' rnb.step.gender.prediction
#'
#' Executes the gender prediction step and adds a dedicated section to the report.
#' 
#' @param rnb.set Methylation dataset as an object of type \code{\linkS4class{RnBSet}}.
#' @param report  Report on annotation inferrence to contain the gender prediction section. This must be an object of
#'                type \code{\linkS4class{Report}}.
#' @return List of two elements:
#'         \describe{
#'           \item{\code{"dataset"}}{The dataset, possibly enriched with predicted gender annotation.}
#'           \item{\code{"report"}}{The modified report.}
#'         }
#'
#' @note This function is not used because gender prediction is performed in the data import module, whereas results are
#'       moved to a different module - quality control.
#' @author Yassen Assenov
#' @noRd
rnb.step.gender.prediction <- function(rnb.set, report) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}

	if (inherits(rnb.set, "RnBeadRawSet")) {
		shifts <- rnb.get.XY.shifts(rnb.set)
	} else {
		shifts <- NULL
	}
	rnb.set <- rnb.set.update.predicted.gender(rnb.set, shifts)
	report <- rnb.section.gender.prediction(rnb.set, shifts, report)
	return(list(dataset = rnb.set, report = report))
}
