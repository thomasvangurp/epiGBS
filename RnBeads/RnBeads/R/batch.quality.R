########################################################################################################################
## batch.quality.R
## created: 2012-06-28
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Batch effects associated with quality probe information.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.execute.batch.qc
#'
#' Computation of correlations and permutation-based p-values for detecting quality-associated batch effects.
#'
#' @param rnb.set      HumanMethylation450K dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param pcoordinates Coordinates of the samples of \code{rnb.set} in the principal components space, as returned by
#'                     \code{\link{rnb.execute.dreduction}}.
#' @param permutations Matrix of sample index permutations, as returned by \code{\link{rnb.execute.batcheffects}}. If
#'                     this parameter is \code{NULL}, permutation-based p-values are not calculated.
#' @return \code{NULL} if no principal components for batch analysis are specified (
#'         \code{rnb.getOption("exploratory.principal.components") == 0}); otherwise, a hierarchical structure of
#'         matrices in the form of a nested list. The root branches are represented by the elements
#'         \code{"correlations"} and \code{"pvalues"}. Every element is a list of control probe types; each type is in
#'         turn a list of up to two matrices of correlations between probe values and principal components - one for the
#'         probes on the green channel and one for the red channel. Note that the \code{"pvalues"} branch is not
#'         returned when \code{permutations} is \code{NULL}.
#'
#' @author Pavlo Lutsik
#' @export
rnb.execute.batch.qc <- function(rnb.set, pcoordinates, permutations = NULL) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (is.null(qc(rnb.set))) {
		stop("invalid value for rnb.set; quality control data is missing")
	}
	nsamples <- validate.pcoordinates.all(pcoordinates)
	if (is.character(nsamples)) {
		stop(nsamples)
	}
	rm(nsamples)
	if (!(is.null(permutations) || is.valid.permutations(permutations))) {
		stop("invalid value for permutations")
	}

	pc.association.count <- rnb.getOption("exploratory.principal.components")
	if (pc.association.count == 0) {
		return(NULL)
	}

	## Initialize data points, control probe types and control probe information
	if(rnb.set@target=="probes450"){
		CONTROL.TYPES <- rnb.infinium.control.targets("probes450")[c(1:4, 6, 11:14)]
		id.col<-"ID"
		type.col<-"Target"
		ctrls<-"controls450"
	}else if(rnb.set@target=="probes27"){
		CONTROL.TYPES <- rnb.infinium.control.targets("probes27")[c(1:4, 6, 9:11)]
		id.col<-"Address"
		type.col<-"Type"
		ctrls<-"controls27"
	}
	control.probe.infos <- rnb.get.annotation(ctrls)[, c(id.col, type.col)]
	channels <- list("green" = qc(rnb.set)$Cy3, "red" = qc(rnb.set)$Cy5)

	targets <- names(pcoordinates)
	names(targets) <- 1:length(targets)
	results<-lapply(1:length(targets), function(target.id){
		
		dpoints <- pcoordinates[[target.id]]$pca$x
		if (ncol(dpoints) > pc.association.count) {
			dpoints <- dpoints[, 1:pc.association.count]
		}
		
		result <- list("correlations" = list())
		if (!is.null(permutations)) {
			result[["pvalues"]] <- list()
		}
	
		get.probes.channel <- function(channel, ids) {
			i <- intersect(ids, rownames(channel))
			if (length(i) == 0) {
				return(NULL)
			}
			result <- channel[i, ]
			if (length(i) == 1) {
				result <- t(result)
				rownames(result) <- i
			}
			return(result)
		}
		init.matrix <- function(N) {
			matrix(as.double(NA), nrow = N, ncol = ncol(dpoints),
				dimnames = list("Probe" = 1:N, "Principal component" = 1:ncol(dpoints)))
		}
	
		## Compute correlations between QC probe values and principal components
		target2ids <- tapply(control.probe.infos[[id.col]], control.probe.infos[[type.col]], as.character)
		for (ci in 1:length(CONTROL.TYPES)) {
			ctype <- names(CONTROL.TYPES)[ci]
			ids <- list(target2ids[[CONTROL.TYPES[ci]]])
			p.channels <- mapply(get.probes.channel, channels, ids, SIMPLIFY = FALSE)
			if (all(sapply(p.channels, is.null))) {
				next
			}
			result[["correlations"]][[ctype]] <- list()
			if (!is.null(permutations)) {
				result[["pvalues"]][[ctype]] <- list()
			}
			for (c.channel in names(p.channels)) {
				matrix.channel <- p.channels[[c.channel]]
				if (is.null(matrix.channel)) {
					next
				}
				table.correlations <- init.matrix(nrow(matrix.channel))
				table.pvalues <- init.matrix(nrow(matrix.channel))
	
				for (i in 1:nrow(matrix.channel)) {
					for (j in 1:ncol(dpoints)) {
						t.result <- test.traits(matrix.channel[i, ], dpoints[, j], permutations)
						table.correlations[i, j] <- t.result[["correlation"]]
						table.pvalues[i, j] <- t.result[["pvalue"]]
					}
				}
				result[["correlations"]][[ctype]][[c.channel]] <- table.correlations
				if (!is.null(permutations)) {
					result[["pvalues"]][[ctype]][[c.channel]] <- table.pvalues
				}
			}
			rnb.status(c("Computed associations between", ncol(dpoints), "principal component(s) and", ctype, "for",
					targets[target.id]))
		}
		attr(result, "Target") <- targets[target.id]
		result
	})

	return(results)
}

########################################################################################################################

#' rnb.section.batch.qc
#'
#' Adds a section on quality-associated batch effects to the given report.
#'
#' @param report         Report to contain the quality-associated section. This must be an object of type
#'                       \code{\linkS4class{Report}}.
#' @param qccorrelations Nested lists of correlation and p-value matrices, as returned by
#'                       \code{\link{rnb.execute.batch.qc}}.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.batch.qc <- function(report, qccorrelations) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	## TODO: Validate qccorrelations
	
	targets<-sapply(qccorrelations, attr, "Target")
	report.descriptions <- c("correlations" = as.character(NA), "pvalues" = as.character(NA))
	
	channels<-NULL
	ctypes<-NULL
	pvals.present<-NULL
	table.pvalue.links<-list()
	
	all.plots<-lapply(1:length(targets), function(target.id){
	
		## Create heatmaps of correlations and p-values
		pvals.present <<- ("pvalues" %in% names(qccorrelations[[target.id]]))
		ctypes <<- names(qccorrelations[[target.id]][["correlations"]])
		channels <<- sort(unique(unlist(lapply(qccorrelations[[target.id]][["correlations"]], names))))
		probecount <- max(sapply(unlist(qccorrelations[[target.id]][["correlations"]], recursive = FALSE), nrow))
		imgheight <- list("correlations" = 1 + probecount * 0.4, "pvalues" = 1 + probecount * 0.25)
		report.plots <- list()

		if (pvals.present) {
			## Create a table of links to the p-value matrices
			table.pvalue.links[[target.id]] <<- matrix(as.character(NA), nrow = length(ctypes), ncol = length(channels),
				dimnames = list(ctypes, channels))
		}
	
		heatmap.functions <- list("correlations" = plot.heatmap.pc.correlations, "pvalues" = plot.heatmap.pc.pvalues)
		for (tblname in names(qccorrelations[[target.id]])) {
			heatmap.function <- heatmap.functions[[tblname]]
			report.plots[[tblname]] <- list()
			for (i in 1:length(ctypes)) {
	#i <- 1
				tables <- qccorrelations[[target.id]][[tblname]][[i]]
				for (cchannel in channels) {
	#cchannel <- "green"
					fname <- paste("heatmap", tblname, "pc", target.id, cchannel, i, sep = "_")
					if (cchannel %in% names(tables)) {
						## Create a heatmap of the table of correlations or pvalues
						rplot <- heatmap.function(report, tables[[cchannel]], fname, height = imgheight[[tblname]])
						if (is.na(report.descriptions[tblname])) {
							report.descriptions[tblname] <<- rplot$description
						}
						rplot <- rplot$plot
	
						## Save the table of p-values to a file
						if (tblname == "pvalues") {
							fname <- paste("pvalues_pc_", target.id, "_", cchannel, "_", i, ".csv", sep = "")
							fname.full <- file.path(rnb.get.directory(report, "data", absolute = TRUE), fname)
							tbl <- cbind("Probe \\ Principal component" = rownames(tables[[cchannel]]),
								as.data.frame(tables[[cchannel]], check.names = FALSE, stringsAsFactors = FALSE))
							utils::write.csv(tbl, file = fname.full, row.names = FALSE)
							table.pvalue.links[[target.id]][i,grep(cchannel, colnames(table.pvalue.links[[target.id]]))] <<- paste("<a href=\"", rnb.get.directory(report, "data"), "/",
								fname, "\">csv</a>", sep = "")
							
							rm(fname.full, tbl)
						}
					} else {
						## TODO: Create an image "n/a"
					}
					report.plots[[tblname]] <- c(report.plots[[tblname]], rplot)
				}
			}
		}
#		rm(imgheight, tblname, i, tables, cchannel, fname, rplot)
		report.plots
	})
	setting.names <- list("Location type" = targets, "Channel" = channels, "Probe group" = ctypes)
	names(setting.names[["Location type"]]) <- 1:length(targets)
	names(setting.names[["Channel"]]) <- channels
	names(setting.names[["Probe group"]]) <- 1:length(ctypes)
	
	txt <- c("The heatmaps below visualize the Pearson correlation coefficients between the principal ",
		"components and the signal levels of selected quality control probes.")
	rnb.add.paragraph(report, txt)
	
	report.plots.corrs<-unlist(lapply(all.plots, el, "correlations"), recursive = F)
	report.plots.pvals<-unlist(lapply(all.plots, el, "pvalues"), recursive = F)
	
	report <- rnb.add.figure(report, report.descriptions["correlations"], report.plots.corrs, setting.names)

	if (!is.null(report.plots.pvals) && all(!sapply(report.plots.pvals, is.null))) {
		txt <- c("In a complete analogy to the heatmaps above, the figure below visualizes the p-values calculated ",
			"using permutation tests.")
		rnb.add.paragraph(report, txt)
		report <- rnb.add.figure(report, report.descriptions["pvalues"], report.plots.pvals, setting.names)

		txt <- c("All computed p-values for associations are available as comma-separated files that accompany this ",
			"report. The links to the dedicated files are provided in the table below.")
		rnb.add.paragraph(report, txt)

		table.pvalue.links<-lapply(table.pvalue.links, function(tpl){
					rn<-rownames(table.pvalue.links[[1]])
					tpl<-cbind("Control probe type \\ Target - Channel" = rn, tpl)
					colnames(tpl) <- capitalize(colnames(tpl))
					tpl
				})
		names(table.pvalue.links) <- paste("batchQcPvalueTable", 1:length(targets), sep = "_")

		table.header <- c("<colgroup>", paste0("\t<col width=\"", c(72, 14, 14), "%\" />"), "</colgroup>")
		report <- rnb.add.tables(report, table.pvalue.links, setting.names[1], row.names = FALSE,
			first.col.header = TRUE, thead = table.header)
	}

	return(report)
}

########################################################################################################################

#' rnb.step.batch.qc
#'
#' Detection of quality-associated batch effects.
#'
#' @param rnb.set      HumanMethylation450K dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param report       Report to contain the quality-associated section. This must be an object of type
#'                     \code{\linkS4class{Report}}.
#' @param pcoordinates Coordinates of the samples of \code{rnbSet} in the principal components space, as returned by
#'                     \code{\link{rnb.execute.dreduction}}.
#' @param permutations Matrix of sample index permutations, as returned by \code{\link{rnb.execute.batcheffects}}.
#' @return The modified report.
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.step.batch.qc <- function(rnb.set, report, pcoordinates, permutations = NULL) {
	if (!inherits(rnb.set, "RnBeadSet")){
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	logger.start("Quality-associated Batch Effects")
	if (!is.null(pcoordinates)) {
		nsamples <- validate.pcoordinates.all(pcoordinates)
		if (is.character(nsamples)) {
			stop(nsamples)
		}
		rm(nsamples)
	}
	if (!(is.null(permutations) || is.valid.permutations(permutations))) {
		stop("invalid value for permutations")
	}

	if (is.null(qc(rnb.set))) {
		## No quality information is available; skip this step
		logger.warning("No quality information available")
		txt <- c("The studied dataset does not include quality probe signal values, therefore, it cannot be examined ",
			"for quality-associated batch effects.")
		report <- rnb.add.section(report, "Quality-associated Batch Effects", txt)
		logger.completed()
		return(report)
	}
	if (is.null(pcoordinates)) {
		txt <- c("No low-dimensional representation of the dataset is available, therefore, it cannot be examined for ",
			"quality-associated batch effects.")
		report <- rnb.add.section(report, "Quality-associated Batch Effects", txt)
		logger.completed()
		return(report)
	}
	txt <- "This section examines the methylation values of the dataset for quality-associated batch effects."
	report <- rnb.add.section(report, "Quality-associated Batch Effects", txt)

	## Compute correlations between QC probes and principal components
	qccorrelations <- rnb.execute.batch.qc(rnb.set, pcoordinates, permutations)

	## Add heatmaps of correlation between QC probes and principal components
	report <- rnb.section.batch.qc(report, qccorrelations)

	logger.status("Added heatmaps of correlation between principal components and signal levels of Q probes")

	logger.completed()
	return(report)
}
