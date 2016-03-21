########################################################################################################################
## clustering.R
## created: 2012-06-18
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Definition and implementation of clustering algorithms and visualization.
########################################################################################################################

## G L O B A L S #######################################################################################################

## Dissimilarity metrics used in clustering
DIST.METRICS <- c(
	"correlation" = "correlation-based",
	"manhattan" = "Manhattan distance",
	"euclidean" = "Euclidean distance")

## Clustering algorithms applied
ALGORITHMS <- c("hierarchical" = "hierarchical")

## Linkage methods used in hierarchical clustering
AGGLOMERATIONS <- c("average" = "average", "complete" = "complete", "median" = "median-based")

## C L A S S ###########################################################################################################

#' RnBeadClustering Class
#'
#' Storage class for the results of a clustering algorithm applied on an \code{\linkS4class{RnBSet}} dataset.
#'
#' @section Slots:
#' \describe{
#'   \item{dissimilarity}{Dissimilarity metric used in the form of a one-element \code{character} vector.}
#'   \item{dimensionality}{Dimensionality of the clustered points in the form of a one-element \code{integer} vector.}
#'   \item{algorithm}{Clustering algorithm (and optionally, type) as a \code{character} vector of length 1 or 2.}
#'   \item{result}{Resulting object after applying the clustering algorithm on a dataset.}
#'   \item{assignments}{Cluster assignments for the samples in the dataset as a matrix. Row names in this matrix are
#'        sample identifiers, and each column is dedicated to partitioning into \emph{k} clusters for a fixed \emph{k}.}
#'   \item{silhouettes}{\code{numeric} vector of mean silhouette values for each tested value of \emph{k}.}
#' }
#'
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{samples}}{Gets the identifiers of all samples used in the clustering.}
#' }
#'
#' @name RnBeadClustering-class
#' @rdname RnBeadClustering-class
#' @aliases initialize,RnBeadClustering-method
#' @aliases samples,RnBeadClustering-method
#' @author Yassen Assenov
#' @exportClass RnBeadClustering
setClass("RnBeadClustering",
	representation(
		dissimilarity = "character",
		dimensionality = "integer",
		algorithm = "character",
		result = "list",
		assignments = "matrix",
		silhouettes = "numeric"),
	package = "RnBeads")

## M E T H O D S #######################################################################################################

setValidity("RnBeadClustering",
	function(object) {
		x <- object@dissimilarity
		if (!(is.character(x) && length(x) == 1 && (!is.na(x)))) {
			return("invalid value for dissimilarity; expected single character")
		}
		if (!(x %in% names(DIST.METRICS))) {
			return("unsupported value for dissimilarity")
		}
		x <- object@dimensionality
		if (!(is.integer(x) && length(x) == 1 && isTRUE(x > 0))) {
			return("invalid value for dimensionality; expected single positive integer")
		}
		x <- object@algorithm
		if (!(is.character(x) && (length(x) %in% c(1, 2)) && (!any(is.na(x))))) {
			return("invalid value for algorithm; expected character")
		}
		if (!(x[1] %in% names(ALGORITHMS))) {
			return("unsupported value for algorithm")
		}
		if (length(x) == 2) {
			if (!(x[2] %in% names(AGGLOMERATIONS))) {
				return("unsupported value for algorithm; unsupported linkage method")
			}
		}
		if (!(length(object@result) != 0)) {
			return("invalid value for result")
		}
		x <- object@assignments
		if (!(is.integer(x) && is.matrix(x) && (!any(is.na(x))))) {
			return("invalid value for assignments")
		}
		## TODO: Validate that colnames(x) are "integer" values in the range 2:nrow(x)
		y <- object@silhouettes
		if (!(is.numeric(y) && length(y) == ncol(x))) {
			return("invalid value for silhouettes")
		}
		if (!identical(colnames(x), names(y))) {
			return("mismatch between assignments and silhouettes")
		}
		TRUE
	}
)

########################################################################################################################

setMethod("initialize", "RnBeadClustering",
	function(.Object, dissimilarity, dimensionality, algorithm, result, assignments, silhouettes) {
		.Object@dissimilarity <- dissimilarity
		.Object@dimensionality <- dimensionality
		.Object@algorithm <- algorithm
		.Object@result <- result
		.Object@assignments <- assignments
		.Object@silhouettes <- silhouettes
		validObject(.Object)
		.Object
	}
)

########################################################################################################################

if (!isGeneric("samples")) {
	setGeneric("samples", function(object) standardGeneric("samples"))
}

#' @rdname samples-methods
#' @export
setMethod("samples", signature(object = "RnBeadClustering"),
	function(object) {
		rownames(object@assignments)
	}
)

## F U N C T I O N S ###################################################################################################

## dist.correlation
##
## Computes correlation-based dissimilarities between the row vectors of the given matrix.
##
## @param x Matrix to use in calculating dissimilarities.
## @return Object of class \code{"dist"}.
## @author Yassen Assenov
dist.correlation <- function(x) {
	res <- as.dist(1 - cor(t(x), use = "pairwise.complete.obs"))
	#take care of the following type of warnings that occur for invariable rows:
	#	Warning message:
	#	In cor(t(x), use = "pairwise.complete.obs") :
	#	the standard deviation is zero
	res[is.na(res)] <- 1
	return(res)
}

########################################################################################################################

## get.silhouette
##
## Computes the average silhouette value of the specified clustering using the given distance metric.
##
## @param x    \code{integer} vector with k different integer cluster codes.
## @param dist Distance metric to use.
## @return Average among the silhouette values of all points in \code{x}.
## @author Yassen Assenov
get.silhouette <- function(x, dist) {
	mean(silhouette(x, dist = dist)[, "sil_width"])
}

########################################################################################################################

## plot.heatmap.rand
##
## Creates a heatmap visualizing the given Rand indices.
##
## @param tbl    Matrix of computed adjusted Rand indices. All values are expected to be between -1 and 1.
## @param fname  File name to contain the generated plot. If this is \code{NA} (default), the heatmap is plotted to the
##               currently active graphics device.
## @param report Report to host the generated image. This parameter is used only when \code{fname} is not \code{NA}. If
##               provided, this parameter must be an instance of class \code{Report}.
##
## @author Yassen Assenov
plot.heatmap.rand <- function(tbl, fname = NA, report = NULL) {
	tbl.melt <- melt(tbl, varnames = c("x", "y"))
	tbl.melt[[1]] <- factor(as.character(tbl.melt[[1]]), levels = rev(rownames(tbl)))
	tbl.melt[[2]] <- factor(as.character(tbl.melt[[2]]), levels = colnames(tbl))
	colnames(tbl.melt)[3] <- "rindex"
	if (!is.na(fname)) {
		width <- 2 + ncol(tbl) * 0.4 + 1.5
		height <- 3 + nrow(tbl) * 0.4
		rplot <- createReportPlot(fname, report, width = width, height = height)
	}
	colors.g <- rnb.getOption("colors.3.gradient")
	pp <- ggplot(tbl.melt) + aes_string("y", "x") + coord_fixed() + labs(x = NULL, y = NULL, fill = "Rand index") +
		geom_tile(aes_string(fill = "rindex"), color = "white") +
		scale_fill_gradient2(limits = c(-1, 1), low = colors.g[1], mid = colors.g[2], high = colors.g[3]) +
		scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
		theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
			plot.margin = unit(c(0.5, 1.5, 0.5, 0.5), c("lines", "in", "lines", "lines")),
			legend.justification = c(0, 0.5), legend.position = c(1, 0.5))
	print(pp)
	if (!is.na(fname)) {
		return(invisible(off(rplot)))
	}
}

########################################################################################################################

#' rnb.execute.clustering
#'
#' Performs hierarchical clustering on the samples of the given dataset using multiple distance metrics and
#' agglomeration methods for a single given region type.
#'
#' @param rnb.set     Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param region.type the clustering is performed on methylation levels from regions of that type. 
#' 					  see \code{\link{rnb.region.types}} for possible values.
#' @return List of clustering results, whereby each element is an object of type \code{\linkS4class{RnBeadClustering}}.
#'         In case clustering cannot be performed, the return value is \code{NULL}. Reasons for a failure include, among
#'         others, the case when \code{rnb.set} contains less than 3 samples, or undefined distances between a pair of
#'         samples due to (too many) missing values in the respective methylation matrix.
#'         
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' results <- rnb.execute.clustering(rnb.set.example, "promoters")
#' # List applied dissimilarity metrics
#' sapply(results, slot, "dissimilarity")
#' # List applied clustering algorithms
#' str(lapply(results, slot, "algorithm"))
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.clustering <- function(rnb.set, region.type = "sites") {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}

	X <- t(meth(rnb.set, type = region.type))
	if (nrow(X) < 3) {
		## Too few samples; cannot compute silhouette values
		rnb.status(c("Skipped clustering on", region.type, ": too few samples"))
		return(NULL)
	}
	if (rnb.getOption("exploratory.clustering") == "top") {
		N <- max(rnb.getOption("exploratory.clustering.top.sites"))
		if (N < ncol(X)) {
			X <- X[, order(colVars(X, na.rm = TRUE), decreasing = TRUE)[1:N], drop = FALSE]
		}
		rm(N)
	}
	results <- list()
	for (dist.metric in names(DIST.METRICS)) {
		if (dist.metric == "correlation") {
			dmatrix <- dist.correlation(X)
		} else {
			dmatrix <- dist(X, method = dist.metric)
		}
		for (agglomeration in names(AGGLOMERATIONS)) {
			result <- tryCatch(hclust(dmatrix, method = agglomeration), error = function(e) { NULL })
			if (is.null(result)) {
				logger.warning(c("Failed clustering on", region.type, "( using", dist.metric, ")"))
				return(NULL)
			}
			cluster.assignments <- as.matrix(cutree(result, k = c(2:(nrow(X) - 1))))
			silhouette.values <- apply(cluster.assignments, 2, get.silhouette, dist = dmatrix)
			attr(result, "class") <- NULL # remove S3 class information in order to store it in RnBeadClustering
			result <- new("RnBeadClustering", dist.metric, ncol(X), c("hierarchical", agglomeration), result,
				cluster.assignments, silhouette.values)
			results <- c(results, result)
		}
		rnb.status(c("Performed clustering on", region.type, "using", dist.metric, "as a distance metric"))
	}
	return(results)
}

########################################################################################################################

#' rnb.execute.clustering.all
#'
#' Performs hierarchical clustering on the samples of the given dataset using multiple distance metrics and
#' agglomeration methods for all suggested site and region types.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @return List of list of clustering results; each element corresponds to one region type and is a list of objects
#'         of type \code{\linkS4class{RnBeadClustering}}.
#'
#' @seealso \code{\link{rnb.execute.clustering}} for performing clustering using a single site or region type.
#'
#' @author Fabian Mueller
#' @export
rnb.execute.clustering.all <- function(rnb.set) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}

	result <- list()
	if (rnb.getOption("analyze.sites")) {
		result[["sites"]] <- rnb.execute.clustering(rnb.set, "sites")
	}
	for (rtype in rnb.region.types.for.analysis(rnb.set)) {
		clust.res <- rnb.execute.clustering(rnb.set, rtype)
		if (is.null(clust.res)){
			result[rtype] <- list(NULL)
		} else {
			result[[rtype]] <- clust.res
		}
	}

	return(result)
}

########################################################################################################################

### Just a helper function to not copy and paste code
rnb.section.clustering.add.heatmap <- function(report, X, fname, cluster.rows, clust.result, sample.ids, locus.colors,
	sample.colors = NULL) {

	dist.metric <- clust.result@dissimilarity
	agglomeration <- clust.result@algorithm[2]
	rplot <- createReportPlot(fname, report, width = 7.2, height = 7.2,
		create.pdf = rnb.getOption("exploratory.clustering.heatmaps.pdf"), high.png = 200)
	cresult <- clust.result@result
	attr(cresult, "class") <- "hclust"
	cresult <- as.dendrogram(cresult)

	heatmap.parameters <- list(x = X, Rowv = FALSE, Colv = cresult, dendrogram = "column", scale = "none",
		col = get.methylation.color.panel(), trace = "none")
	if (cluster.rows) {
		## Plotting sites/regions, attempt to cluster rows
		heatmap.parameters$labRow <- NA
		heatmap.parameters$margins <- c(4, 1)
		if (dist.metric == "correlation") {
			distfun <- dist.correlation
		} else {
			distfun <- function(x) { dist(x, method = dist.metric) }
		}
		cc <- tryCatch(hclust(distfun(X), method = agglomeration), error = function(err) {
				## Row dendrogram cannot be drawn
				logger.warning(c("Could not cluster rows of the heatmap with method", dist.metric, "and linkage",
						agglomeration, ". Discarding row dendrogram"))
				return(NULL)
			}
		)
		if (!is.null(cc)) {
			heatmap.parameters$Rowv <- as.dendrogram(cc)
			heatmap.parameters$dendrogram <- "both"
		}
	} else {
		## Plotting quantiles
		heatmap.parameters$labRow <- rownames(X)
		heatmap.parameters$margins <- c(4, 3)
	}
	heatmap.parameters$labCol <- sample.ids
	heatmap.parameters$RowSideColors <- locus.colors
	heatmap.parameters$ColSideColors <- sample.colors
	tryCatch(do.call(heatmap.2, heatmap.parameters), error = function(err) {
			print(rnb.message.plot("Heatmap could not be created"))
			list()
		}
	)
	return(off(rplot))
}

########################################################################################################################

#' rnb.section.clustering
#'
#' Creates a report section dedicated to sample clustering.
#'
#' @param report          Report on methylation profiles to contain the clustering section. This must be an object of
#'                        type \code{\linkS4class{Report}}.
#' @param rnb.set         Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param clust.results   Results of applying hierarchical clustering as an object of type
#'                        \code{\linkS4class{RnBeadClustering}}, or a non-empty list of such objects.
#' @param site.reg.info   Information about traits to be visualized and tested for associations with clusters, as
#'                        returned by \code{\link{get.site.and.region.types}}.
#' @param section.started Flag indicating if the clustering section in the given report is already created.
#' @return The modified report.
#'
#' @author Yassen Assenov, Fabian Mueller
#' @noRd
rnb.section.clustering <- function(report, rnb.set, clust.results, site.reg.info, section.started = FALSE) {

	rnb.require("mclust")
	beta.values <- meth.matrices(rnb.set)

	if (inherits(clust.results, "RnBeadClustering")) {
		clust.results <- list(clust.results)
	}
	if (!(is.list(clust.results) && length(clust.results) != 0 &&
			all(sapply(unlist(clust.results,recursive=FALSE), inherits, "RnBeadClustering")))) {
		stop("invalid value for clust.results; expected list of list of RnBeadClustering objects")
	}
	if (!all(sapply(names(clust.results), function(rr) { all(sapply(clust.results[[rr]],
						FUN=function(x){identical(samples(x), colnames(beta.values[[rr]]))})) }))) {
		stop("invalid value for clust.results; does not match rnb.set")
	}
	dist.metrics <- unique(sapply(unlist(clust.results,recursive=FALSE), slot, "dissimilarity"))
	algorithms <- lapply(clust.results, FUN=function(x){
		lapply(x, slot, "algorithm")
	})
	if (!all(sapply(unlist(algorithms,recursive=FALSE), length) == 2)) {
		stop("unsupported clustering algorithms")
	}
	agglomerations <- unique(sapply(unlist(algorithms,recursive=FALSE), "[", 2))
	if (any(sapply(clust.results,FUN=function(x){length(dist.metrics) * length(agglomerations) != length(x)}))) {
		stop("invalid value for clust.results; list is incomplete")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	logger.start("Clustering Section")

	## Define sample IDs and coloring
	sample.ids <- samples(rnb.set)
	pheno.columns <- names(rnb.sample.groups(rnb.set,rnb.getOption("exploratory.columns")))
	if (length(pheno.columns) != 0) {
		pheno.colors <- rnb.pheno2colors(pheno(rnb.set)[, pheno.columns])
		if (length(pheno.columns) == 1) {
			pheno.colors$colors <- as.matrix(pheno.colors$colors)
			colnames(pheno.colors$colors) <- pheno.columns
			pheno.colors$mapping <- list(pheno.colors$mapping)
		}
		scolors <- 1:length(pheno.columns)
	} else {
		pheno.colors <- NULL
		scolors <- 0L
	}

	locus.columns <- unique(unlist(lapply(site.reg.info,FUN=function(x){colnames(x)})))
	lcolors <- 0L
	if (length(locus.columns) != 0) {
		lcolors <- 1:length(locus.columns)
	}
	
	reg.types <- intersect(names(clust.results),names(beta.values))

	## Create heatmaps
	description <- "The figure below shows clustering of samples using several algorithms and distance metrics."
	if (section.started) {
		rnb.add.paragraph(report, description)
	} else {
		report <- rnb.add.section(report, "Clustering", description)
	}
	probe.order <- lapply(beta.values,FUN=function(x){ 
		order(rowVars(x,na.rm = TRUE), decreasing = TRUE)
	})
	logger.start("Generating Heatmaps")
	report.plots <- list(quant=list(),loc=list())
	topprobes <- c(0L, rnb.getOption("exploratory.clustering.top.sites"))
	locus.colors <- list()
	for (ri in 1:length(reg.types)) {
		rt <- reg.types[ri]
		if (rt %in% names(site.reg.info)) {
			cur.annot <- colnames(site.reg.info[[rt]])
			locus.colors[[rt]] <- lapply(locus.columns,FUN=function(aa){
				if (aa %in% cur.annot){
					get.site.and.region.types.colors(site.reg.info,rt,aa)
				} else {
					NULL
				}
			})
			names(locus.colors[[rt]]) <- locus.columns
		} else {
			locus.colors[[rt]] <- NULL
		}
		
		logger.start(c("Region type:",rt))
		clust.res <- clust.results[[rt]]
		for (tprobes in topprobes){
			tprobes.trunc <- tprobes #truncate if tprobes exceeds the number of sites/regions
			if (tprobes == 0) {
				X <- apply(beta.values[[rt]], 2, quantile, probs = seq(0, 1, length = 101), na.rm = TRUE)
				rownames(X)[1:nrow(X) %% 10 != 1] <- ""
			} else {
				tprobes.trunc <- min(c(tprobes,length(probe.order[[rt]])))
				X <- beta.values[[rt]][probe.order[[rt]][1:tprobes.trunc], ]
			}
			for (i in 1:length(clust.res)){
				clust.result <- clust.res[[i]]
				dist.metric <- clust.result@dissimilarity
				agglomeration <- clust.result@algorithm[2]
				fname.base <- paste(ri, dist.metric, agglomeration, sep = "_")
				for (si in scolors) {
					
					sample.colors <- NULL
					locus.colors.cur <- NULL
					fname <- fname.base
					if (si != 0) {
						fname <- paste(fname.base, si, sep = "_")
						sample.colors <- pheno.colors$colors[, si]
					}
					if (tprobes != 0){
						fname.base.2 <- paste("heatmapLoc",fname,sep="_")
						fname <- fname.base.2
						for (li in lcolors) {
							if (li != 0) {
								fname <- paste(fname.base.2, li, sep = "_")
								locus.colors.cur <- locus.colors[[rt]][[locus.columns[li]]][probe.order[[rt]][1:tprobes.trunc]]
							}
							fname <- paste(fname, tprobes, sep = "_")
							rplot <- rnb.section.clustering.add.heatmap(report,X,fname,TRUE,clust.result,sample.ids,locus.colors.cur,sample.colors)

							report.plots[["loc"]] <- c(report.plots[["loc"]],rplot)
						}
					} else {
						fname <- paste("heatmapQuant",fname,sep="_")
						rplot <- rnb.section.clustering.add.heatmap(report,X,fname,FALSE,clust.result,sample.ids,NULL,sample.colors)
						report.plots[["quant"]] <- c(report.plots[["quant"]],rplot)
					}
				}
			}
		}
		logger.completed()
	}
	setting.names.common <- list(
		"Site/region level" = reg.types,
		"Dissimilarity metric" = DIST.METRICS[dist.metrics],
		"Agglomeration strategy (linkage)" = AGGLOMERATIONS[agglomerations])
	names(setting.names.common[["Site/region level"]]) <- 1:length(reg.types)

	setting.names <- setting.names.common
	if (max(scolors) != 0) {
		sam.colors <- colnames(pheno.colors$colors)
		names(sam.colors) <- 1:ncol(pheno.colors$colors)
		setting.names[["Sample color based on"]] <- sam.colors
		rm(sam.colors)
	}
	txt.clustering <- paste0("Hierarchical clustering of samples based on ",
		ifelse(rnb.getOption("exploratory.clustering") == "all", "all methylation values",
		paste(max(rnb.getOption("exploratory.clustering.top.sites")), "most variable loci")), ". The heatmap displays ")
	description <- c(txt.clustering, "methylation percentiles per sample. ",
		"The legend for sample coloring can be found in the figure below.")
	report <- rnb.add.figure(report, description, report.plots[["quant"]], setting.names)

	if (max(lcolors) != 0) {
		lcolors <- locus.columns
		names(lcolors) <- 1:length(locus.columns)
		setting.names[["Site/region color based on"]] <- lcolors
	}
	setting.names[["Visualize"]] <- paste(topprobes[-1], "most variable loci")
	names(setting.names[["Visualize"]]) <- topprobes[-1]
	description <- c(txt.clustering, "only selected sites/regions with the highest variance across all samples. ",
		"The legend for locus and sample coloring can be found in the figure below.")
	report <- rnb.add.figure(report, description, report.plots[["loc"]], setting.names)

	logger.status(c("Created", length(report.plots[["quant"]]) + length(report.plots[["loc"]]),
		"heatmaps based on the clustering results"))
	logger.completed()

	## Add color legends
	logger.start("Adding Color Legends")
	if (!is.null(pheno.colors) || length(locus.columns) != 0) {
		add.legend <- function(s.index, legend.colors) {
			screen(s.index)
			par(mar = c(0, 0, 0, 0) + 0.1)
			plot.new()
			legend("center", legend = names(legend.colors), fill = legend.colors, bty = "n", ncol = 1)
		}
		ns <- length(scolors)
		nl <- length(lcolors)
		## TODO: Simplify this weird parameter computations
		param.mat <- c()
		if (max(scolors) == 0) {
			param.mat <- cbind(0,1:nl)
		} else {
			for (i in 1:ns){
				if (nl == 0) {
					param.mat <- rbind(param.mat,cbind(i,0))
				} else {
					param.mat <- rbind(param.mat,cbind(i,1:nl))
				}
			}
		}

		report.plots <- list()
		for (ri in 1:length(reg.types)) {
			rt <- reg.types[ri]
			for (i in 1:(dim(param.mat)[1])) {
				fname <- paste("heatmapLegend",ri, param.mat[i,1],param.mat[i,2], sep = "_")
				rplot <- createReportPlot(fname, report, width = 7.2, height = 1.8, high.png = 200)
				screens <- rbind(c(0, 0.5, 0, 1), c(0.5, 1, 0, 1))
				s.indices <- split.screen(screens)
				if (param.mat[i,2] != 0){
					leg.loc <-  attr(site.reg.info[[rt]],"legend")[[lcolors[param.mat[i,2]]]]
					if (!is.null(leg.loc)){
						#convert color gradient
						if (!is.character(leg.loc)){
							if (leg.loc == 2){
								colors.grad <- rnb.getOption("colors.gradient")
								vals <- round(site.reg.info[[rt]][,lcolors[param.mat[i,2]]],2)
								leg.loc <- colors.grad[1:2]
								names(leg.loc) <- c(min(vals),max(vals))
							} else {
								stop("Invalid legend")
							}
						}
						add.legend(s.indices[1],leg.loc)
					}

				}
				if (param.mat[i,1] != 0) {
					add.legend(s.indices[2], pheno.colors$mapping[[param.mat[i,1]]])
				}
				close.screen(all.screens = TRUE)
				off(rplot)
				report.plots <- c(report.plots,rplot)
			}
		}
		if (max(lcolors) == 0) {
			lcolors <- "n.a."
			names(lcolors) <- 0
		}
		if (max(scolors) == 0) {
			scolors <- "n.a."
			names(scolors) <- 0
		} else {
			scolors <- colnames(pheno.colors$colors)
			names(scolors) <- 1:ncol(pheno.colors$colors)
		}
		regs <- reg.types
		names(regs) <- 1:length(reg.types)
		setting.names.leg <- list(
				"Site/region level" = regs,
				"Sample color based on" = scolors,
				"Site/region color based on" = lcolors)
		
		description <- c("Probe and sample colors used in the heatmaps in the previous figures.")
		report <- rnb.add.figure(report, description, report.plots, setting.names.leg)
	}
	rm(scolors, description, probe.order, topprobes, tprobes, X)
	logger.completed()

	## Create plots of silhouette values
	logger.start("Estimating Optimal Numbers of Clusters")
	refText <- c("Rousseeuw, P. J. (1987) Silhouettes: A graphical aid to the interpretation and validation of ",
		"cluster analysis. <i>Journal of Computational and Applied Mathematics</i>, <b>20</b>, 53-65")
	report <- rnb.add.reference(report, refText)
	description <- c("Using the average silhouette value as a measure of cluster assignment ",
		rnb.get.reference(report, refText), ", it is possible to infer the number of clusters produced by each of the ",
		"studied methods. The figure below shows the corresponding mean silhouette value for every observed ",
		"separation into clusters.")
	report <- rnb.add.section(report, "Identified Clusters", description, level = 3)
	## TODO: Support also different lengths of silhouettes
	colnames.sil <- lapply(reg.types,FUN=function(rt){
		as.character(AGGLOMERATIONS[sapply(algorithms[[rt]], "[", 2)])
	})
	names(colnames.sil) <- reg.types
	silhouettes <- lapply(reg.types,FUN=function(rt){
			sils <- as.matrix(sapply(clust.results[[rt]], slot, "silhouettes"))
			colnames(sils) <- colnames.sil[[rt]]
			return(sils)
	})
	names(silhouettes) <- reg.types
	algorithm.names <- lapply(algorithms,FUN=function(algs){
		sapply(algs, function(a.name) {
			base::ifelse(length(a.name) == 2, paste0(a.name[1], " (", a.name[2], " linkage)"), a.name)
		})
	})
	clusters.identified.inner <- data.frame(
		"Metric" = character(0),
		"Algorithm" = character(0),
		"Clusters" = integer(0),
		stringsAsFactors = FALSE)
	clusters.identified <- list()
	report.plots.all <- list()
	for (rt in 1:length(reg.types)){
		rt <- reg.types[ri]
		clusters.identified.cur <- clusters.identified.inner
		report.plots <- lapply(dist.metrics, function(dist.metric) {
			cinds <- (sapply(clust.results[[rt]], slot, "dissimilarity") == dist.metric)
			for (i in which(cinds)) {
				cn <- as.integer(rownames(silhouettes[[rt]])[which.max(silhouettes[[rt]][, i])])
				clusters.identified.cur[nrow(clusters.identified.cur) + 1, "Metric"] <<- DIST.METRICS[dist.metric]
				clusters.identified.cur[nrow(clusters.identified.cur), "Algorithm"] <<- algorithm.names[[rt]][i]
				clusters.identified.cur[nrow(clusters.identified.cur), "Clusters"] <<- cn
			}
			dframe <- data.frame(
				"clusters" = rep(as.integer(rownames(silhouettes[[rt]])), sum(cinds)),
				"silhouette" = as.vector(silhouettes[[rt]][, cinds]),
				"linkage" = rep(colnames(silhouettes[[rt]])[cinds], each = nrow(silhouettes[[rt]])), check.names = FALSE)
			fname <- paste("silhouette", ri, dist.metric, sep = "_")
			rplot <- createReportPlot(fname, report, width = 7, height = 5)
			pp <- ggplot(dframe, aes_string(x = "clusters", y = "silhouette", colour = "linkage", group = "linkage")) +
				scale_colour_manual(values = rnb.getOption("colors.category")) +
				labs(x = paste("Number of clusters"), y = paste("Silhouette value")) + geom_line()
			print(pp)
			off(rplot)
			rplot
		})
		report.plots.all <- c(report.plots.all,report.plots)
		clusters.identified <- c(clusters.identified,list(clusters.identified.cur))
	}
	names(clusters.identified) <- 1:length(reg.types)
	description <- c("Line plot visualizing mean silhouette values of the clustering algorithm outcomes for each ",
		"applicable value of <i>K</i> (number of clusters).")
	report <- rnb.add.figure(report, description, report.plots.all, setting.names.common[1:2])
	description <- "The table below summarizes the number of clusters identified by the algorithms."
	rnb.add.paragraph(report, description)
	report <- rnb.add.tables(report, clusters.identified, setting.names.common[1], row.names = FALSE)
	logger.status("Estimated number of clusters based on mean silhouette value")
	logger.completed()

	logger.start("Overlapping Clusters with Sample Traits")
	if (!is.null(pheno.colors)) {
		report <- rnb.add.section(report, "Clusters and Traits", NULL, level = 3)

		## Compute Rand indices to test for associations between clusterings and traits	
		randIndices.all <- list()
		csv.files <- sapply(reg.types, function(x) { as.character(NA) })
		for (ri in 1:length(reg.types)){
			rt <- reg.types[ri]
			randIndices <- matrix(as.double(NA), nrow = nrow(clusters.identified[[ri]]), ncol = ncol(pheno.colors$colors))
			dimnames(randIndices) <- list(clusters.identified[[ri]][, "Algorithm"], colnames(pheno.colors$colors))
			dist.metrics <- DIST.METRICS[sapply(clust.results[[rt]], slot, "dissimilarity")]
			for (i in 1:nrow(clusters.identified[[ri]])) {
				j <- which(dist.metrics == clusters.identified[[ri]][i, "Metric"] &
						algorithm.names[[rt]] == clusters.identified[[ri]][i, "Algorithm"])
				cl.assignments <- clust.results[[rt]][[j]]@assignments[, as.character(clusters.identified[[ri]][i, "Clusters"])]
				for (j in colnames(randIndices)) {
					tr.assignments <- as.integer(as.factor(pheno(rnb.set)[, j]))
					valids <- (!is.na(tr.assignments))
					randIndices[i, j] <- mclust::adjustedRandIndex(cl.assignments[valids], tr.assignments[valids])
				}
			}
			randIndices.all <- c(randIndices.all,list(randIndices))
			rm(dist.metrics, i, j, cl.assignments, tr.assignments, valids)
	
			## Save the table to an external file
			rIndices <- cbind(clusters.identified[[ri]][, c("Metric", "Algorithm")], randIndices)
			csv.files[rt] <- paste("adjusted_rand_indices_",ri,".csv",sep="")
			fname.full <- file.path(rnb.get.directory(report, "data", absolute = TRUE), csv.files[rt])
			utils::write.csv(rIndices, file = fname.full, row.names = FALSE)
			csv.files[rt] <- paste(rnb.get.directory(report, "data"), csv.files[rt], sep = "/")
			logger.status(c("Computed adjusted rand indices and saved to", csv.files[rt]))
		}
		names(randIndices.all) <- reg.types
		refText <- c("Hubert, L. and Arabie, P. (1985) Comparing partitions. ",
			"<i>Journal of Classification</i>, <b>2</b>(1), 193-218")
		report <- rnb.add.reference(report, refText)
		txt <- c("The figure below shows associations between clusterings and the examined traits. Associations are ",
			"quantified using the adjusted Rand index ", rnb.get.reference(report, refText), ". Rand indices near 1 ",
			"indicate high agreement while values close to -1 indicate seperation. The full table of all computed ",
			"indices is stored in the following comma separated files:")
		rnb.add.paragraph(report, txt)
		## Add a list of links to the created CSV files
		txt <- as.list(paste("<a href=\"", csv.files, "\">", names(csv.files), " (csv)</a>", sep = ""))
		rnb.add.list(report, txt)
		rm(csv.files, ri, rt, randIndices, rIndices, fname.full, refText, txt)

		## Create heatmaps of adjusted Rand indices
		report.plots.all <- list()
		for (ri in 1:length(reg.types)){
			rt <- reg.types[ri]
			report.plots <- lapply(names(DIST.METRICS), function(dist.metric) {
				fname <- paste("heatmap", "rindex", ri, dist.metric, sep = "_")
				cinds <- which(DIST.METRICS[dist.metric] == clusters.identified[[ri]][, "Metric"])
				return(plot.heatmap.rand(t(randIndices.all[[rt]][cinds, ]), fname, report))
			})
			report.plots.all <- c(report.plots.all,report.plots)
		}
		description <- c("Heatmap visualizing Rand indices computed between sample traits (rows) and clustering ",
			"algorithm outcomes (columns).")
		report <- rnb.add.figure(report, description, report.plots.all, setting.names.common[1:2])
	}
	logger.completed()
	logger.completed()
	return(report)
}

########################################################################################################################

#' rnb.step.clustering
#'
#' Performs hierarchical clustering on the samples of the given dataset and generates a corresponding section in the
#' report.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report on methylation profiles to contain the clustering section. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return The modified report.
#'
#' @seealso \code{\link{rnb.execute.clustering}}, \code{\link{rnb.section.clustering}}
#' @author Yassen Assenov
#' @noRd
rnb.step.clustering <- function(rnb.set, report) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	rnb.step.clustering.internal(rnb.set, report, get.site.and.region.types(rnb.set))
}

rnb.step.clustering.internal <- function(rnb.set, report, rinfo) {
	logger.start("Sample Clustering")

	## Perform clustering algorithms
	logger.start("Agglomerative Hierarchical Clustering")
	clust.results <- rnb.execute.clustering.all(rnb.set)
	logger.completed()

	## Handle the cases when clusterings on some or all of the methylation matrices failed
	clust.failed <- sapply(clust.results, is.null)
	clust.edited <- sum(clust.failed) != 0
	if (clust.edited) {
		failed.single <- sum(clust.failed) == 1
		failed.all <- sum(clust.failed) == length(clust.results)
		txt <- "The clustering algorithm was not completed"
		if (failed.all) {
			txt <- c(txt, ifelse(failed.single, ".", " for all requested levels of methylation."))
		} else {
			txt <- c(txt, " for the following requested level", ifelse(failed.single, "", "s"), " of methylation - ")
			txt <- c(txt, paste("<b>", names(clust.results)[clust.failed], "</b>", sep = "", collapse = ", "), ".")
		}
		if (length(samples(rnb.set)) < 3) {
			txt <- c(txt, " Note that clustering is performed only when the dataset contains at least 3 samples.")
		} else {
			txt <- c(txt, " The most likely reason for this error is inability to calculate a distance between a pair ",
				"of samples due to (too many) missing values in the respective methylation ",
				ifelse(failed.single, "matrix", "matrices"), ".")
		}
		report <- rnb.add.section(report, "Clustering", txt)
		if (failed.all) {
			logger.completed()
			return(report)
		}
		clust.results <- clust.results[!clust.failed]
		rm(failed.single, failed.all, txt)
	}
	rm(clust.failed)

	## Generate clustering section in the report
	report <- rnb.section.clustering(report, rnb.set, clust.results, rinfo, clust.edited)
	logger.completed()
	return(report)
}
