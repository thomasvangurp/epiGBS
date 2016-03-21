########################################################################################################################
## plottingUtils.R
## created: 2012-05-12
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Constants and functions related to plotting.
########################################################################################################################

## G L O B A L S #######################################################################################################

## Full legends for site and region annotations that are supported by RnBeads
LEGENDS.ALL <- list(
	"CGI Relation" = c("Open Sea" = "#0000FF", "Shelf" = "#00FFFF", "Shore" = "#B200FF", "Island" = "#FF0000"),
	"Design" = c("I" = "#FF8080", "II" = "#8080FF"),
	"CpG Density" = 2L)

## color palettes to be used when there are multiple density estimates in one plot
## (e.g. for scatterplots or differential methylation)
## colors based on ColorBrewer color palette "Paired", pairs 1 and 3
DENS.COLORS.LOW  <- c("#1F78B4","#E31A1C") 
DENS.COLORS.HIGH <- c("#A6CEE3","#FB9A99")

## F U N C T I O N S ###################################################################################################

## Applies the given coloring scheme to the vector of values.
##
## @param val    a vector of values
## @param rng	a 2-tuple descibing the possible range of values
## @param colscheme.col.val the colorscheme that is applied
## @return a vector of color values
## @author Fabian Mueller
colorize.value <- function(val, rng=c(min(val),max(val)),
	colscheme.col.val=c(gplots::colorpanel(100,"blue","white"),gplots::colorpanel(100,"white","red")))
{
	#require(gplots)
	return(colscheme.col.val[round((val - rng[1]) / (rng[2] -rng[1])  * (length(colscheme.col.val)-1),0)+1])
}

## Returns a color panel of color values for methylation levels from the "colors.meth option"
##
## @author Fabian Mueller
get.methylation.color.panel <- function()
{
	meth.color.base <- rnb.getOption("colors.meth")
	if (length(meth.color.base) == 3) {
		meth.color.panel <- gplots::colorpanel(100,meth.color.base[1],meth.color.base[2],meth.color.base[3])
	} else { # length(meth.color.base) == 2
		meth.color.panel <- gplots::colorpanel(100,meth.color.base[1],meth.color.base[2])
	}
	return(meth.color.panel)
}

########################################################################################################################

## Applies the given coloring scheme to the vector of values.
##
## @param x      Non-empty vector of values.
## @param cols   Coloring scheme to be applied. This must be either a non-empty \code{character} vector storing
##               different colors, or a function that takes one integer argument \code{n} and returns a vector of
##               \code{n} different colors.
## @param na.col Color to be used for \code{NA} elements of \code{x}. If the value vector contains \code{NA}s, this
##               color is appended to the mapping and used to indicate a missing value.
## @return List of two elements:
##         \itemize{
##           \item{"colors" }{Vector of length \code{length(x)} containing the assigned colors for the values in
##                \code{x}.}
##           \item{"mapping" }{\code{character} vector storing the applied mapping from value to color.}
##         }
## @author Yassen Assenov
rnb.get.cols <- function(x, cols, na.col = "#C0C0C0") {
	contains.na <- any(is.na(x))
	xvalues <- unique(sort(x))
	if (typeof(cols) == "closure") {
		color.map <- as.character(cols(length(xvalues)))
	} else if (is.character(cols) && length(cols) != 0) {
		if (length(cols) < length(xvalues)) {
			cols <- rep(cols, ceiling(length(xvalues) / length(cols)))
		}
		color.map <- cols[1:length(xvalues)]
	} else {
		stop("invalid value for cols")
	}
	names(color.map) <- as.character(xvalues)
	result <- color.map[as.character(x)]
	if (any(is.na(x))) {
		color.map <- c(color.map, na.col)
		names(color.map)[length(color.map)] <- NA
		result[is.na(result)] <- na.col
	}
	return(list("colors" = result, "mapping" = color.map))
}

########################################################################################################################

#' get.site.and.region.types
#'
#' Initializes color legends and extracts annotation values of site and/or region types covered in the given methylation
#' dataset.
#'
#' @param rnb.set Methylation dataset of interest.
#' @return List of mappable annotation values, one \code{data.frame} per site or region type. The columns of this
#'         \code{data.frame} are the supported annotations. Every table also contains an attribute named \code{"legend"}
#'         that contains a \code{list} of the legends for all columns in the respective annotation table. Every legend
#'         item is either a \code{character} vector storing the mapping from values to colors, or a single
#'         \code{integer} value, signifying the number of colors in a gradient.
#' @author Yassen Assenov, Fabian Mueller
#' @noRd
get.site.and.region.types <- function(rnb.set) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	result <- list()
	rtypes <- rnb.region.types.for.analysis(rnb.set)
	if (rnb.getOption("analyze.sites")) {
		rtypes <- c("sites", rtypes)
	}
	for (rtype in rtypes) {
		anno.table <- tryCatch(annotation(rnb.set, type = rtype), error = function(e) { NULL })
		if (is.null(anno.table)) {
			next
		}
		## Add a column for CpG Density
		if (!("CpG Density" %in% colnames(anno.table)) && ("CpG" %in% colnames(anno.table))) {
			if (rtype == "sites") {
				LENGTH.NEIGHBORHOOD <- 100L # TODO: Get this value from the data package
				anno.table[, "CpG Density"] <- anno.table[, "CpG"] / LENGTH.NEIGHBORHOOD
			} else {
				anno.table[, "CpG Density"] <- anno.table[, "CpG"] / (anno.table[, "End"] - anno.table[, "Start"] + 1)
			}
		}

		## Extract suitable columns from the annotation table
		c.names <- intersect(names(LEGENDS.ALL), colnames(anno.table))
		if (length(c.names) == 0) {
			next
		}
		if (length(c.names) != 1) {
			anno.table <- anno.table[, c.names]
		} else {
			anno.table <- data.frame("1" = anno.table[, c.names])
			colnames(anno.table) <- c.names
		}
		## Adjust the CGI Relation values if necessary
		if (("CGI Relation" %in% c.names) &&
			(!identical(levels(anno.table[, "CGI Relation"]), names(LEGENDS.ALL[["CGI Relation"]])))) {
			cgi.relations <- as.character(anno.table[, "CGI Relation"])
			cgi.relations[cgi.relations %in% c("North Shelf", "South Shelf")] <- "Shelf"
			cgi.relations[cgi.relations %in% c("North Shore", "South Shore")] <- "Shore"
			anno.table[, "CGI Relation"] <- factor(cgi.relations, levels = names(LEGENDS.ALL[["CGI Relation"]]))
		}

		## Construct color legends
		legend <- LEGENDS.ALL[c.names]
		for (c.name in c.names) {
			lg <- legend[[c.name]]
			if (is.character(lg)) {
				lg <- lg[intersect(names(lg), unique(anno.table[, c.name]))]
				legend[[c.name]] <- lg
			}
		}
		attr(anno.table, "legend") <- legend
		result[[rtype]] <- anno.table
	}
	return(result)
}

########################################################################################################################

## get.site.and.region.types.colors
##
## Takes the returned value from \code{\link{get.site.and.region.types}} and returns a vector containing the color
## values for the given locus type.
##
## @param snrt      List of site and region types, as returned by \code{get.site.and.region.types}.
## @param type      Site or region type.
## @param annot.ind Column index or name for \code{snrt[[type]]}, i.e. the index of the annotation in the object.
## @return a vector containing a color value for each locus
get.site.and.region.types.colors <- function(snrt,type,annot.ind){
	vals <- snrt[[type]][,annot.ind]
	leg <- attr(snrt[[type]],"legend")[[annot.ind]]
	if (is.character(leg)){
		return(leg[vals])
	}
	else if (leg == 2){
		colors.grad <- rnb.getOption("colors.gradient")
		return(colorize.value(vals,colscheme.col.val=c(gplots::colorpanel(100,colors.grad[1],colors.grad[2]))))	
	} else {
		stop("Error in converting values to colors: invalid legend")
	}
}

########################################################################################################################

## rnb.pheno2colors
##
## Maps the given sample phenotypic information to colors.
##
## @param sample.table Phenotypical information of a single trait, stored in a non-empty \code{vector} or \code{factor};
##                     or of multiple traits, stored in a non-empty table (\code{matrix} or \code{data.frame}). In case
##                     this is a table, rows are treated as samples, and columns - as traits. Note that missing values
##                     (\code{NA}), if present in a trait, are assigned a hard-coded color irrespective of the provided
##                     scheme.
## @param cols         Color scheme(s) to be mapped to the values in \code{sample.table}. This may be specified in one
##                     of three forms: \code{character} vector, function, or a list of schemes, each element in the list
##                     being a vector or a function. If the number of the given color schemes is smaller than the number
##                     of traits in \code{sample.table}, the schemes are recycled.
## @return     List of two elements:
##             \itemize{
##                \item{"colors" }{\code{character} vector or matrix of the same dimensions as \code{sample.table},
##                     containing the assigned colors based on the given values.}
##                \item{"mapping" }{\code{character} vector storing the applied mapping from value to color. or a list
##                     of such vector, one for each column of \code{colors}.}
##             }
## @author Yassen Assenov
rnb.pheno2colors <- function(sample.table, cols = rnb.getOption("colors.category")) {
	if ((is.vector(sample.table) || is.factor(sample.table)) && length(sample.table) != 0) {
		if (is.list(cols) && length(cols) != 0) {
			cs <- cols[[1]]
		} else if (is.character(cols) || typeof(cols) == "closure") {
			cs <- cols
		} else {
			stop("invalid value for cols")
		}
		return(rnb.get.cols(sample.table, cs))
	}
	if (!(is.matrix(sample.table) || is.data.frame(sample.table))) {
		stop("invalid value for sample.table")
	}
	if (nrow(sample.table) == 0 || ncol(sample.table) == 0) {
		stop("invalid value for sample.table")
	}
	if (is.list(cols) && length(cols) != 0) {
		ci <- 0
	} else if (is.character(cols) || typeof(cols) == "closure") {
		cs <- cols
	} else {
		stop("invalid value for cols")
	}
	result <- list(
		"colors" = matrix("", nrow = nrow(sample.table), ncol = ncol(sample.table), dimnames = dimnames(sample.table)),
		"mapping" = list())
	for (i in 1:ncol(sample.table)) {
		if (is.list(cols)) {
			ci <- ci %% length(cols) + 1
			cs <- cols[[ci]]
		}
		ri <- rnb.get.cols(sample.table[, i], cs)
		result$colors[, i] <- ri$colors
		result$mapping[[i]] <- ri$mapping
	}
	return(result)
}

########################################################################################################################

#' rnb.plot.beta.comparison
#'
#' Draws plots that compare two distributions of beta values.
#'
#' @param beta.values      Two beta value sequences in the form of a named \code{list} of two non-empty vectors of type
#'                         \code{double}. If any of the vectors contains \code{NA}s, this method may exit with an error.
#' @param fprefix          File name prefix for the plots. This function appends the suffixes \code{"_density"},
#'                         \code{"_histogram"} and \code{"_qq"} to this prefix.
#' @param report           Report to which the plots are to be added.
#' @param qq.length        Positive \code{integer} value showing the number of quantiles to be calculated and presented
#'                         in the generated Q-Q plot.
#' @param points.per.group Maximum number of values to use in plotting a group's distribution. Groups that contain more
#'                         observations than this threshold are subsampled. Setting this parameter to a value less than
#'                         \code{2} disables subsampling.
#' @return List of all generated plots, each being an object ot type \code{\linkS4class{ReportPlot}}.
#'
#' @author Yassen Assenov
#' @export
rnb.plot.beta.comparison <- function(beta.values, fprefix, report = NULL, qq.length = 501L,
	points.per.group = rnb.getOption("distribution.subsample")) {

	if (!(is.list(beta.values) && length(beta.values) == 2 && all(sapply(beta.values, is.double)))) {
		stop("invalid value for beta.values")
	}
	if (any(sapply(beta.values, length) == 0)) {
		stop("invalid value for beta.values; both sequences must be non-empty")
	}
	if (!(is.character(fprefix) && length(fprefix) == 1 && (!is.na(fprefix)))) {
		stop("invalid value for fprefix")
	}
	if (!(is.null(report) || inherits(report, "Report"))) {
		stop("invalid value for report")
	}
	if (is.double(qq.length) && isTRUE(all(as.integer(qq.length) == qq.length))) {
		qq.length <- as.integer(qq.length)
	}
	if (!(is.integer(qq.length) && length(qq.length) == 1 && isTRUE(qq.length > 0))) {
		stop("invalid value for qq.length; expected positive integer")
	}
	if (is.double(points.per.group) && isTRUE(all(as.integer(points.per.group) == points.per.group))) {
		points.per.group <- as.integer(points.per.group)
	}
	if (!(is.integer(points.per.group) && length(points.per.group) == 1 && (!is.na(points.per.group)))) {
		stop("invalid value for points.per.group")
	}

	## Perform subsampling, if required
	subsampled <- rep(FALSE, length(beta.values))
	if (points.per.group > 1) {
		for (i in 1:length(beta.values)) {
			if (length(beta.values[[i]]) > points.per.group) {
				beta.values[[i]] <- sample(beta.values[[i]], size = points.per.group)
				subsampled[i] <- TRUE
			}
		}
	}

	## Construct a data.frame to plot
	dframe <- data.frame(
	 	"value" = unlist(beta.values, use.names = FALSE),
	 	"vtype" = factor(rep(1:length(beta.values), sapply(beta.values, length))), check.names = FALSE)
	levels(dframe[[2]]) <- names(beta.values)

	report.plots <- list()
	fnames <- paste(fprefix, c("density", "histogram", "qq"), sep = "_")
	pmargins <- unit(c(0.1, 1.6, 0.1, 0.1), "in")

	## Create a density estimation plot
	rplot <- createReportPlot(fnames[1], report, width = 6.8, height = 5.2)
	pp <- ggplot(dframe, aes_string(x = "value")) +
		labs(x = expression(beta), y = "Density", color = "Values") +
		geom_density(aes_string(color = "vtype"), kernel = "gaussian") +
		theme(plot.margin = pmargins, legend.justification = c(0, 0.5), legend.position = c(1, 0.5))
	print(pp)
	off(rplot)
	report.plots <- c(report.plots, rplot)
	
	## Create histograms
	rplot <- createReportPlot(fnames[2], report, width = 6.8, height = 5.2)
	pp <- ggplot(dframe, aes_string(x = "value")) + coord_cartesian(xlim = c(0, 1)) +
		labs(x = expression(beta), y = "Density") +
		geom_histogram(aes_string(y = "..density.."), binwidth = 0.02) +
		facet_grid(vtype ~ .) + theme(plot.margin = pmargins)
	print(pp)
	off(rplot)
	report.plots <- c(report.plots, rplot)
	
	## Create a Q-Q plot
	v.quantiles <- as.data.frame(lapply(beta.values, quantile, probs = seq(0, 1, length.out = qq.length)))
	v.diagonal <- data.frame("x" = c(0, 1), "y" = c(0, 1))
	colnames(v.quantiles) <- colnames(v.diagonal) <- names(beta.values)
	cnames <- paste0("`", colnames(v.quantiles), "`")
	rplot <- createReportPlot(fnames[3], report, width = 6.8, height = 5.2)
	pp <- ggplot(v.quantiles, aes_string(x = cnames[1], y = cnames[2])) +
		coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + coord_fixed() +
		geom_path(data = v.diagonal, color = "#A0A0A0") + geom_path(lineend = "round", size = 2) +
		theme(plot.margin = pmargins)
	print(pp)
	off(rplot)
	report.plots <- c(report.plots, rplot)

	attr(report.plots, "subsampled") <- subsampled
	return(report.plots)
}

########################################################################################################################

## add.text.subsampling
##
## Generates, if necessary, text about subsampling. This text should be added to a figure description.
##
## @param subsampled       \code{logical} vector of length \code{2} containing information if each of the targeted two
##                         value sequences was subsampled.
## @param group.names      Names of the two groups of values.
## @param points.per.group Size of sample, if sampling was applied.
## @return \code{character} vector containing descriptive text about subsampling; an empty vector if subsampling was not
##         performed.
## @author Yassen Assenov
add.text.subsampling <- function(subsampled, group.names, points.per.group = rnb.getOption("distribution.subsample")) {
	if (any(subsampled)) {
		if (all(subsampled)) {
			txt <- "Both distributions are"
		} else {
			txt <- paste(" The distribution of", tolower(group.names[which(subsampled)]), "is")
		}
		txt <- c(txt, " estimated by randomly sampling ", points.per.group, " values",
			ifelse(all(subsampled), " in each group.", "."))
	} else {
		txt <- character()
	}
	return(txt)
}

########################################################################################################################

#' rnb.message.plot
#'
#' Creates a plot, using \pkg{ggplot2}, with a single text message.
#'
#' @param txt Text to be plotted.
#' @return The newly initialized \code{ggplot} instance.
#'
#' @examples
#' \dontrun{
#' x11(width = 5, height = 5)
#' rnb.message.plot("Missing data")
#' }
#' @author Yassen Assenov
#' @export
rnb.message.plot <- function(txt) {
	if (!(is.character(txt) && length(txt) == 1 && (!is.na(txt)))) {
		stop("invalid value for txt")
	}
	ggplot(data.frame(x = 1, y = 1, labeltext = txt), aes_string("x", "y", label = "labeltext")) +
		geom_text(color = "grey50") +
		theme(axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
			axis.ticks = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
			panel.background = element_blank(), plot.background = element_blank())
}

########################################################################################################################

#' rnb.color.legends
#'
#' Creates a figure in the given report that contains one or more color legends.
#'
#' @param report        Report to contain the legend figure. This must be an object of type \code{\linkS4class{Report}}.
#' @param legends       Color legend in the form of a non-empty \code{character} vector. Element names denote legend
#'                      labels, and the elements themselves specify colors. This parameter can also be a \code{list} of
#'                      color legends. Special restrictions apply to the names of the list elements, see \emph{Details}.
#' @param fprefix       File name or prefix for the plot files.
#' @param description   Text of the figure description. See the correponding parameter in
#'                      \code{\link{rnb.add.figure}} for more details.
#' @param setting.names One-element list containing a plot file descriptor, when \code{legends} is a list. See the
#'                      corresponding parameter in \code{\link{rnb.add.figure}} for more details. If this is set to
#'                      \code{NULL} (default), the list is automatically created using \code{names(legends)} (when
#'                      \code{legends} is a list), or as an empty list (when \code{legends} is a vector).
#' @param size.factor   Relative size, in inches of the plots. Legends are displayed in columns of up to 10 items; each
#'                      column is effectively a square with the specified size.
#' @return The modified report.
#'
#' @details In case \code{legends} specifies multiple legends in the form of a list, \code{names(legends)} are appended
#'          to \code{fprefix} to generate file names. In order to ensure independence of the operating system, there are
#'          strong restrictions on these names. They can consist of the following symbols only: Latin letters, digits,
#'          dot (\code{.}), dash (\code{-}) and underline (\code{_}).
#'
#' @author Yassen Assenov
#' @export
rnb.color.legends <- function(report, legends, fprefix = ifelse(is.character(legends), "legend", "legend_"),
	description = "Color legend.", setting.names = NULL, size.factor = 3) {
	if (!(is.character(legends) && length(legends) != 0)) {
		if (!(is.list(legends) && all(sapply(legends, is.character)) && all(sapply(legends, length) != 0))) {
			stop("invalid value for legends; expected character or list")
		}
		lnames <- names(legends)
		if (!(is.character(lnames) && all(sapply(lnames, is.valid.relative)))) {
			stop("invalid value for legends; missing or invalid legend names")
		}
		if (anyDuplicated(lnames) != 0) {
			stop("invalid value for legends; duplicated legend names")
		}
		if (!all(grepl("^[A-Za-z0-9._-]+$", lnames))) {
			stop("invalid value for legends; unsupported legend names")
		}
		if (is.null(setting.names)) {
			setting.names <- list("Based on" = lnames)
			names(setting.names[[1]]) <- lnames
		}
	} else {
		setting.names <- list()
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report; expected Report")
	}
	if (!(is.vector(fprefix) && length(fprefix) == 1 && (!is.na(fprefix)))) {
		stop("invalid value for fprefix; expected single-element vector")
	}
	fprefix <- as.character(fprefix)
	if (!((is.double(size.factor) || is.integer(size.factor)) && length(size.factor) == 1 && (!is.na(size.factor)))) {
		stop("invalid value for size.factor; expected scalar of type integer or double")
	}

	## Create plots for the legends
	plotlegend <- function(fname, report, legendinfo) {
		legend.columns <- (length(legendinfo) - 1) %/% 10 + 1
		rplot <- createReportPlot(fname, report, width = size.factor * legend.columns, height = size.factor)
		par(mar = c(0, 0, 0, 0) + 0.1)
		plot.new()
		legend("center", legend = names(legendinfo), fill = legendinfo, bty = "n", ncol = legend.columns)
		off(rplot)
		return(rplot)
	}
	if (is.character(legends)) {
		report.plots <- plotlegend(fprefix, report, legends)
	} else { # is.list(legends) && all(sapply(legends, is.character))
		report.plots <- lapply(names(legends), function(lname) {
			plotlegend(paste(fprefix, lname, sep = ""), report, legends[[lname]])
		})
	}

	## Add a figure to the report
	return(rnb.add.figure(report, description, report.plots, setting.names))
}

########################################################################################################################

#' rnb.plot.pheno.categories
#'
#' Generates bar charts summarizing the categorical traits in a sample annotation table.
#'
#' @param annotations  Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}, or its sample
#'                     annotations in the form of a \code{data.frame}. If this parameter is a dataset, the annotation
#'                     information is extracted using the method \code{\link[=pheno,RnBSet-method]{pheno}}.
#' @param columns      Optional; predefined column names (in the form of a \code{character} vector) or indices (an
#'                     \code{integer} vector) to consider. All other columns in the annotation table will be ignored.
#' @param fileprefix   \code{character} vector with one element storing the file name prefix of the output files,
#'                     without the extension. Only a limited set of symbols is allowed to be used in this prefix.
#' @param report       Report to contain the generated plots. If specified, this must be an object of type
#'                     \code{\linkS4class{Report}}.
#' @param color.values Non-empty \code{character} vector containing the color scheme to be mapped to the categories
#'                     defined in the annotation table. Colors are recycled if necessary, that is, if the length of this
#'                     vector is smaller than the number of categories in a trait.
#' @return List of report plots. The names in this list are the column names in the annotation table that were selected
#'         for visualization. In case no suitable categorical traits are found among the provided annotations, this
#'         function returns an empty list.
#' 
#' @details This function identifies the traits that define sample subgroups and then generates one report plot per
#' trait. Every report plot consists of two files. File names are formed by appending an index and file extension to
#' \code{fileprefix}. Thus, the suffixes appended are \code{"_1.pdf"}, \code{"_1.png"}, \code{"_2.pdf"},
#' \code{"_2.png"}, ... Existing files with the generated filenames are overwritten.
#' 
#' @seealso \code{\link{rnb.sample.groups}} for identifying traits in the annotation table that define sample subgroups;
#'          \code{\link{createReportPlot}} for the allowed symbols to be used in \code{fileprefix}
#' 
#' @author Yassen Assenov
#' @export
rnb.plot.pheno.categories <- function(annotations, columns = NULL, fileprefix = "barchart_pheno", report = NULL,
	color.values = rnb.getOption("colors.category")) {

	## Validate the function's parameters
	pheno.groups <- lapply(rnb.sample.groups(annotations, columns), sapply, length)
	if (!(is.character(fileprefix) && length(fileprefix) == 1 && (!is.na(fileprefix)))) {
		stop("invalid value for fileprefix")
	}
	if (!grepl("^[A-Za-z0-9._-]+$", fileprefix)) {
		stop("invalid value for fileprefix")
	}
	if (!(is.null(report) || inherits(report, "Report"))) {
		stop("invalid value for report")
	}
	if (!(is.character(color.values) && length(color.values) != 0 && all(!is.na(color.values)))) {
		stop("invalid value for color.values")
	}

	result <- list()
	if (length(pheno.groups) == 0) {
		names(result) <- character()
		return(result)
	}

	## Recycle colors if necessary
	color.values <- rep_len(color.values, max(sapply(pheno.groups, length)))

	## Create bar charts
	for (i in 1:length(pheno.groups)) {
		x <- pheno.groups[[i]]
		x <- data.frame(x = 1L, t = factor(names(x), levels = names(x)), count = as.integer(x))
		pp <- ggplot(x, aes_string(x = "x", y = "count", fill = "t")) +
			ggplot2::geom_bar(stat = "identity", position = "stack") +
			labs(x = NULL, y = "Samples", fill = names(pheno.groups)[i]) +
			scale_fill_manual(na.value = "C0C0C0", values = color.values, guide = guide_legend(reverse = TRUE)) +
			theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
				legend.justification = c(0, 0.5), legend.position = c(1, 0.5),
				plot.margin = unit(c(0.1, 2.1, 0.1, 0.1), "in"))
		fname <- paste(fileprefix, i, sep = "_")
		rplot <- createReportPlot(fname, report, width = 3.1, height = 5.2)
		print(pp)
		result[[names(pheno.groups)[i]]] <- off(rplot)
	}
	result 
}

########################################################################################################################

#' rnb.plot.dreduction
#' 
#' Creates a dimension reduction plot based on the methylation values of the given dataset.
#' 
#' @param rnb.set         Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}. This dataset
#'                        must contain at least four samples.
#' @param plot.type       Type of plot to be created. This must be one of \code{"pca"} (projection to two principal
#'                        components) or \code{"mds"} (multidimensional scaling to two dimensions). The section
#'                        \emph{Details} provides more details on how the dimension reduction techniques are applied.
#' @param dimensions      Vector of two positive \code{integer} values giving the principle components to be shown in
#'                        the horizontal and vertical axis of the plot. This parameter is considered only when
#'                        \code{plot.type} is \code{"pca"}.
#' @param distance.metric Distance metric to be applied when reducing the dimensionality of the methylation data. This
#'                        must be one of \code{"eucledian"} or \code{"manhattan"}. The second metric is supported only
#'                        in multidimensional scaling.
#' @param target          Site or region type to be used in the dimension reduction technique. This must be either
#'                        \code{"sites"} (individual CpGs) or one of the region types summarized in \code{rnb.set}.
#' @param point.types     Trait, specified as column name or index in the sample annnotation table of \code{rnb.set}, to
#'                        be used to define point types in the plot. Setting this parameter to zero (default) or to a
#'                        trait that does not define categories results in all samples being displayed as filled
#'                        circles. If this parameter specifies a column that can be used as sample identifiers, the plot
#'                        displays the samples as identifiers instead of points.
#' @param point.colors    Trait, specified as column name or index in the sample annnotation table of \code{rnb.set}, to
#'                        be used to define sample colors in the plot. Setting this parameter to zero (default) or to a
#'                        trait that does not define categories results in all samples being displayed in black.
#' @param legend.space    Width, in inches, of the space dedicated for legends that will be assigned on the right side
#'                        of the plot. This parameter is considered only if legends are actually included, that is, if
#'                        sample traits are mapped to point types and/or colors.
#' @return The generated plot as an object of type \code{\link[ggplot2:ggplot]{ggplot}}. The object also contains an
#'         attribute \code{"info"}, which is a list with the following elements:
#'         \describe{
#'           \item{\code{"Target" }}{Targeted sites or regions; the value of the parameter \code{target}.}
#'           \item{\code{"Technique" }}{Dimension reduction technique applied; one of \code{"PCA"} or \code{"MDS"}.}
#'           \item{\code{"All" }}{Total number of sites or regions defining the high dimensional methylation space.}
#'           \item{\code{"Missing" }}{Number of dimensions ignored because they contain (only) missing values.}
#'           \item{\code{"Selected" }}{Number of dimensions used when applying a dimension reduction technique.}
#'           \item{\code{"Explained" }}{Value between \code{0} and \code{1} showing the variance explained by the
#'             selected dimensions, as a fraction of the total variance of all dimensions.}
#'         }
#' 
#' @details
#' The analysis option \code{"exploratory.top.dimensions"} controls whether dimension reduction is applied on all
#' probes, sites or regions available in the given dataset, or only on the most variable ones. In case a trait is mapped
#' to point types, the shapes to use are taken from the option \code{"points.category"}. Similary, the option
#' \code{"colors.category"} determines which colors are used when mapping to color is applied. See
#' \emph{\link[=rnb.options]{RnBeads Options}} for more information on these options.
#' 
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' pdf("PCA.pdf", width = 7.2, height = 5.2)
#' print(rnb.plot.dreduction(rnb.set.example, point.colors="Sample_Group"))
#' dev.off()
#' }
#' @seealso \code{\link[=summarized.regions,RnBSet-method]{summarized.regions}} for listing all region types summarized
#'          in a dataset
#' 
#' @author Yassen Assenov
#' @export
rnb.plot.dreduction <- function(rnb.set, plot.type = "pca", dimensions = 1:2, distance.metric = "euclidean",
	target = "sites", point.types = 0L, point.colors = 0L, legend.space = 2) {

	## Validate parameters
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!(is.character(plot.type) && length(plot.type) == 1 && tolower(plot.type) %in% c("mds", "pca"))) {
		stop("invalid value for plot.type")
	}
	plot.type <- tolower(plot.type)
	if (plot.type == "pca") {
		if (is.double(dimensions) && isTRUE(all(dimensions == as.integer(dimensions)))) {
			dimensions <- as.integer(dimensions)
		}
		if (!(is.integer(dimensions) && length(dimensions) == 2 && isTRUE(all(dimensions > 0)))) {
			stop("invalid value for dimensions")
		}
	}
	if (!(is.character(distance.metric) && length(distance.metric) == 1 &&
		  	distance.metric %in% c("manhattan", "euclidean"))) {
		stop("invalid value for distance.metric")
	}
	if (plot.type != "mds" && distance.metric != "euclidean") {
		stop(paste("invalid distance metric for", plot.type))
	}
	if (!(is.character(target) && length(target) == 1 && (!is.na(target)))) {
		stop("invalid value for target")
	}
	if (!(target %in% c("sites", summarized.regions(rnb.set)))) {
		stop("unsupported region type")
	}
	validate.column <- function(x, column.names) {
		if (!(length(x) == 1 && !is.na(x))) {
			return(-1L)
		}
		if (is.double(x) && isTRUE(x == as.integer(x))) {
			x <- as.integer(x)
		}
		if (is.integer(x)) {
			if (x < 0) {
				x <- -1L
			} else if (x > length(column.names)) {
				x <- integer()
			}
		} else if (is.character(x)) {
			x <- which(column.names == x)
		} else {
			x <- -1L
		}
		x
	}
	s.annotation <- pheno(rnb.set)
	point.types <- validate.column(point.types, colnames(s.annotation))
	if (length(point.types) == 0) {
		stop("unsupported column in point.types")
	}
	if (point.types == -1L) {
		stop("invalid value for point.types")
	}
	point.colors <- validate.column(point.colors, colnames(s.annotation))
	if (length(point.colors) == 0) {
		stop("unsupported column in point.colors")
	}
	if (point.colors == -1L) {
		stop("invalid value for point.colors")
	}
	if (!((is.double(legend.space) || is.integer(legend.space)) && length(legend.space) == 1 && (!is.na(legend.space)) &&
		  	0 <= legend.space)) {
		stop("invalid value for legend.space")
	}

	## Extract high-dimensional methylation values
	X <- t(meth(rnb.set, target))
	if (nrow(X) < 4) {
		stop("too few samples")
	}

	site.nas <- colMeans(is.na(X))
	if (plot.type == "mds") {
		i <- which(site.nas == 1)
	} else { # plot.type == c("pca", "pca")
		i <- which(site.nas != 0)
	}
	info <- list("Target" = target, "Technique" = toupper(substr(plot.type, 1, 3)),
		"All" = ncol(X), "Missing" = length(i), "Selected" = ncol(X) - length(i), "Explained" = 1)
	if (length(i) != 0) {
		if (ncol(X) - length(i) < 4) {
			stop("too many missing values")
		}
		## Remove sites that contain (only) NAs
		X <- X[, -i]
		site.nas <- site.nas[-i]
	}
	top.sites <- rnb.getOption("exploratory.top.dimensions")
	site.nas <- which(site.nas != 0)
	if (0 < top.sites && top.sites < ncol(X)) {
		if (top.sites < 4) {
			stop("too few sites selected")
		}
		site.vars <- colVars(X, na.rm = TRUE)
		i <- order(site.vars, decreasing = TRUE)[1:top.sites]
		info$selected <- top.sites
		info$explained <- sum(site.vars[i]) / sum(site.vars)
		X <- X[, i]
	}

	## Perform dimension reduction
	if (plot.type == "mds") {
		X <- mds(X, distance.metric)
		colnames(X) <- paste("MDS dimension", 1:2)
	} else {
		X <- prcomp(X, center = TRUE, scale. = FALSE)$x
		if (ncol(X) < max(dimensions)) {
			stop(paste("unsupported value for dimensions;", ncol(X), "principal components available"))
		}
		X <- X[, dimensions]
		colnames(X) <- paste("Principal component", dimensions)
	}

	## Define the data and mappings to plot
	dframe <- data.frame(x = X[, 1], y = X[, 2])
	p.aes <- list("x" = "x", "y" = "y")
	plot.labs <- list("x" = colnames(X)[1], "y" = colnames(X)[2])
	if (point.colors != 0) {
		i <- s.annotation[, point.colors]
		if (is.character(i) || is.factor(i) || is.integer(i)) {
			dframe$pcolor <- as.factor(i)
			p.aes[["color"]] <- "pcolor"
			plot.labs[["color"]] <- colnames(s.annotation)[point.types]
		} else { # the trait cannot be mapped to point colors
			rnb.warning("Cannot map point.colors to categories")
			point.colors <- 0L
		}
	}
	if (point.types != 0) {
		i <- s.annotation[, point.types]
		if ((is.character(i) || is.factor(i) || is.integer(i)) && (!any(is.na(i))) && anyDuplicated(i) == 0) {
			dframe$id <- as.character(i)
			p.aes[["label"]] <- "id"
		} else if (point.types == point.colors) {
			p.aes[["shape"]] <- p.aes[["color"]]
			plot.labs[["shape"]] <- plot.labs[["color"]]
		}  else if (is.character(i) || is.factor(i) || is.integer(i)) {
			dframe$ptype <- as.factor(i)
			p.aes[["shape"]] <- "ptype"
			plot.labs[["shape"]] <- colnames(s.annotation)[point.types]
		} else { # the trait cannot be mapped to point types
			rnb.warning("Cannot map point.types to categories")
		}
	}

	## Create the ggplot object
	pp <- ggplot2::ggplot(dframe) + do.call(ggplot2::aes_string, p.aes) + do.call(ggplot2::labs, plot.labs)
	if ("id" %in% colnames(dframe)) {
		pp <- pp + ggplot2::geom_text()
	} else {
		pp <- pp + ggplot2::geom_point()
	}
	if ("pcolor" %in% colnames(dframe)) {
		v2color <- rep_len(rnb.getOption("colors.category"), length.out = nlevels(dframe$pcolor))
		pp <- pp + ggplot2::scale_color_manual(na.value = "#C0C0C0", values = v2color)
	}
	if ("shape" %in% names(p.aes)) {
		v2shape <- rep(rnb.getOption("points.category"), length.out = nlevels(dframe[, p.aes[["shape"]]]))
		pp <- pp + ggplot2::scale_shape_manual(na.value = 1L, values = v2shape)
	}
	if (length(intersect(c("pcolor", "ptype"), colnames(dframe))) == 0) {
		pp <- pp + theme(plot.margin = unit(rep(0.1, 4), "in"))
	} else {
		pp <- pp + theme(plot.margin = unit(0.1 + c(0, legend.space, 0, 0), "in"),
			legend.justification = c(0, 0.5), legend.position = c(1, 0.5))
	}

	attr(pp, "info") <- info
	pp
}

########################################################################################################################

#' plotcdf
#'
#' Plots the cumulative distribution function estimated from the given sequence of values.
#'
#' @param fname  Base file name for the generated plot.
#' @param report Report that will contain the plot.
#' @param values Sequence of values sampled from a distribution.
#' @param width  Width, in inches, of the generated PDF image file.
#' @param height Height, in inches, of the generated PDF image file.
#' @param main   Title of the plot.
#' @param xlab   Label of the x axis.
#' @return Newly generated plots as an object of type \code{\linkS4class{ReportPlot}}.
#'
#' @author Yassen Assenov
#' @noRd
plotcdf <- function(fname, report, values, width = 6.2, height = 6.2, main = NA, xlab = "Count", ...) {
	rplot <- createReportPlot(fname, report, width = width, height = height)
	cd.function <- ecdf(values)
	xs <- knots(cd.function)
	dframe <- data.frame("x" = xs, "y" = cd.function(xs), "xend" = c(xs[-1], xs[length(xs)]))
	pp <- ggplot(dframe, aes_string(x = "x", y = "y")) + labs(x = xlab, y = "Fn(x)") +
		geom_segment(aes_string(xend = "xend", yend = "y"), color = "blue") +
		coord_cartesian(xlim = range(values), ylim = c(0, 1)) + theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in"))
	print(pp)
	return(off(rplot))
}

########################################################################################################################

corrHeatmap<-function(cor.mat, ncols, min.corr,
                      leg=TRUE, sp=c(.87,0.9,0.25,0.75), rv=NA, cv=NA, xa=1, ...){


  rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
  color.function<-function(ncols){

    #cm.colors(ncols)
    #rainbow(ncols)
    #rgb.palette(ncols)
    if(ncols<=11) brewer.pal(ncols,"RdGy")
    else rainbow(ncols, start=1/12, end=7/12)
    #topo.colors(ncols)
  }

  color.space.merge.up<-ceiling(ncols*(range(na.omit(as.numeric(cor.mat)))[2]-min.corr)/(1-min.corr))
  color.space.merge.bot<-floor(ncols*(1-(1-range(na.omit(as.numeric(cor.mat)))[1])/(1-min.corr)))


  heatmap.mod(cor.mat,
              Rowv=rv,Colv=cv,
              col=color.function(ncols)[color.space.merge.bot:color.space.merge.up],
              scale="none",
              xaxis=xa, ylab.side=2,...)

  if(leg) image.plot(legend.only=T,
                     graphics.reset=T, smallplot= sp,
                     horizontal=F,
                     zlim=c(min.corr,1), col=color.function(ncols),
                     legend.shrink = 0.5,
                     axis.args=list(cex.axis=0.67), ...)

}


# Heatmap with optional axis' positions
#
###############################################################################

heatmap.mod<-function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
                       distfun = dist, hclustfun = hclust, reorderfun = function(d,
                                                                                 w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv,
                                                                                                                                            "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE,
                       margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 +
                         1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
                       labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE,
                       verbose = getOption("verbose"),
                       xaxis=1, xaxt="nn", ylab.side=4,...)
{
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("'x' must be a numeric matrix")
  nr <- di[1L]
  nc <- di[2L]
  if (nr <= 1 || nc <= 1)
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2L)
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (!doRdend && identical(Colv, "Rowv"))
    doCdend <- FALSE
  if (is.null(Rowv))
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv))
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram"))
      ddr <- Rowv
    else {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv)
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr)))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1L:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram"))
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      hcc <- hclustfun(distfun(if (symm)
        x
                               else t(x)))
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv)
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc)))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1L:nc
  x <- x[rowInd, colInd]
  labRow <- if (is.null(labRow))
    if (is.null(rownames(x)))
      (1L:nr)[rowInd]
  else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol))
    if (is.null(colnames(x)))
      (1L:nc)[colInd]
  else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = na.rm)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2L, sd, na.rm = na.rm)
    x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0,
            4)
  if (!missing(ColSideColors)) {
    if (!is.character(ColSideColors) || length(ColSideColors) !=
      nc)
      stop("'ColSideColors' must be a character vector of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], 0.2, lhei[2L])
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) !=
      nr)
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1L], 0.2, lwid[2L])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei,
        "; lmat=\n")
    print(lmat)
  }
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1L], 0, 0, 0.5))
    image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  if (xaxis==3){
    par(mar = c(margins[1L],margins[2L], margins[1L], 0))
  }else{
    par(mar = c(margins[1L], 0, 0, margins[2L]))
  }


  if (!symm || scale != "none")
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend)
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
    c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
  if(xaxt!="n") axis(xaxis, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
                     cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1L] - 1.25)
  axis(2, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)

  if (!is.null(ylab))
    mtext(ylab, side = ylab.side, line = margins[2L] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  par(mar = c(margins[1L], 0, 0, 0))
  if (doRdend)
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  if(xaxis==3) par(mar = c(0, margins[2L], if (!is.null(main)) 1 else 0, 0))
  else par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
  if (doCdend)
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main))
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&
    doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}

########################################################################################################################

## Creates a deviation plot using the supplied statistics.
##
## @param stats           ...
## @param ylim            ...
## @param additional.info ...
## @param line.col        ...
## @param fill.col        ...
## @param border.col      ...
## @param subtitle        ...
## @param xlab            ...
## @param mar             ...
## @param ...             ...
## @return
## @author Yassen Assenov
deviation.plot <- function(stats, ylim = range(stats), additional.info = NULL,
	line.col = "blue", fill.col = "#FFFF00", border.col = "#BEBEBE",
	subtitle = NA, xlab = NA, mar = c(3, 4, 1, 1) + 0.1, ...) {
	s.indices <- NA
	N <- ncol(stats)
	if (is.matrix(additional.info)) {
		xlabcols <- (ncol(additional.info) == N)
	} else if (is.vector(additional.info)) {
		xlabcols <- (length(additional.info) == N)
	} else {
		xlabcols <- FALSE
	}
	if (xlabcols) {
		screens <- rbind(c(0, 1, 0.20, 1), c(0, 1, 0, 0.20))
		s.indices <- split.screen(screens)
		screen(s.indices[1])
		mar.bottom <- 0.5
	} else {
		mar.bottom <- mar[1]
	}
	par(mar = c(mar.bottom, mar[2:4]), las = 1)
	med.index <- ifelse(nrow(stats) == 1, 1, 2)
	plot(x = c(1, N), y = range(stats), type = "n",
		bty = "n", ylim = ylim, xaxs = "i", xaxt = "n", lwd = 2,
		col = line.col, xlab = NA, ...)
	if (!is.na(subtitle)) {
		mtext(subtitle)
	}
	if (med.index == 2) {
		xx <- c(1:N, N:1)
		yy <- c(stats[1, 1:N], stats[3, N:1])
		polygon(xx, yy, col = fill.col, border = border.col, lwd = 1)
	}
	lines(x = 1:N, stats[med.index, 1:N], lwd = 2, col = line.col)
	polygon(x = c(1, N, N, 1), y = rep(ylim, each = 2))

	if (xlabcols) {
		## Add color bar below the deviation plot
		screen(s.indices[2])
		par(mar = c(mar[1], mar[2], 0, mar[4]))
		if (is.matrix(additional.info)) {
			plot(c(1, N), c(0, 1), type = "n", main = NA, xlab = NA, ylab = NA, xaxs = "i", axes = FALSE)
			for (i in 1:nrow(additional.info)) {
				if (i == 1) {
					ys.start <- rep(0, N)
				} else {
					ys.start <- additional.info[i - 1, ]
				}
				ys.end <- additional.info[i, ]
				j <- which(ys.start != ys.end)
				if (length(j) != 0) {
					ys <- as.vector(rbind(ys.start[j], ys.end[j], rep(as.double(NA), length(j))))
					lines(rep(j, each = 3), ys, col = rownames(additional.info)[i], lwd = 1)
				}
			}
		} else { # is.vector(additional.info)
			plot(rep(1, N), type = "h", col = additional.info,
				ylim = c(0, 1), main = NA, sub = NA, xlab = NA, ylab = NA,
				xaxs = "i", axes = FALSE)
		}
		if (!is.na(xlab)) {
			mtext(xlab, side = 1, line = 1)
		}
		close.screen(all.screens = TRUE)
	} else if (!is.na(xlab)) {
		mtext(xlab, side = 1, line = 1)
	}

	if (med.index == 2) {
		return(mean(stats[3, ] - stats[1, ]))
	}
	return(0)
}

########################################################################################################################

## deviation.plot.beta.get.cuts
##
## Determines the grouping of sites or regions (if necessary) for a deviation plot.
## 
## @param bstats Total number of sites, probes or regions that will be summarized in a deviation plot.
## @return \code{integer} vector of length \code{N}, storing the group indices for every site or region; \code{NULL} if
##         no grouping is necessery.
## @author Yassen Assenov
deviation.plot.beta.get.cuts <- function(N) {
	if (N <= 4000L) {
		return(NULL)
	}
	if (N <= 40000L) {
		return(cut(1:N, breaks = 2000L, labels = FALSE))
	}
	# N > 40000
	return(cut(1:N, breaks = 4000L, labels = FALSE))
}

########################################################################################################################

#' deviation.plot.beta
#'
#' Creates a deviation plot based on the methylation beta values of a population.
#'
#' @param betas    Non-empty numeric \code{matrix} of methylation beta values. Rows in this matrix must denote sites or
#'                 regions, and columns - samples. If a locus (row in the matrix) contains missing values only, it is
#'                 not included in the plot.
#' @param c.values Vector (usually a \code{factor}) storing category or quantitative values for each site or region. The
#'                 length of this vector must be equal to \code{nrow(betas)}, the \emph{i}-th element storing the
#'                 property values for the \emph{i}-th locus in \code{betas}. Note that this vector's names, if present,
#'                 are ignored.
#' @param c.legend If \code{c.values} stores categories, this parameter specifies the mapping from property values to
#'                 colors. The mapping is in the form of a named \code{character} vector. All values that appear in
#'                 \code{c.values} must be present among the names of this vector. The order of the values in this
#'                 mapping determines in which order the colors are stacked (when the number of loci is large). If
#'                 \code{c.values} denotes a quantitative measure, this parameter is a singleton \code{integer},
#'                 specifying the color scheme for visualizing the values. Currently, the only supported values are
#'                 \code{2} and \code{3}. See \code{\link{rnb.options}} for more details.
#' @return Methylation variability as a number between \code{0} and \code{1}, invisibly. This number denotes the relative
#'         area of variation in the generated plot.
#'
#' @author Yassen Assenov
#' @export
deviation.plot.beta <- function(betas, c.values = NULL, c.legend = NULL) {
	if (!(is.numeric(betas) && is.matrix(betas))) {
		stop("invalid value for betas; numeric matrix expected")
	}
	if (nrow(betas) * ncol(betas) == 0) {
		stop("invalid value for betas; unexpected dimensions")
	}
	if (is.null(c.values)) {
		if (!is.null(c.legend)) {
			stop("invalid value for c.legend; expected NULL")
		}
	} else { # !is.null(c.values)
		if (!is.vector(c.values)) {
			stop("invalid value for c.values")
		}
		if (is.character(c.legend)) {
			val.occurring <- unique(c.values)
			val.legend <- names(c.legend)
			if (is.null(val.legend)) {
				stop("invalid value for c.legend; expected named character")
			}
			if (length(setdiff(val.legend, val.occurring)) != 0) {
				stop("invalid value for c.legend; mapping is incomplete")
			}
			c.legend <- c.legend[intersect(val.legend, val.occurring)]
		} else if (!(is.integer(c.legend) && length(c.legend) == 1 && c.legend %in% c(2L, 3L))) {
			stop("invalid value for c.legend")
		}
	}
	bstats <- apply(betas, 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
	site.order <- order(bstats[2, ], bstats[1, ], bstats[3, ], na.last = NA)
	bstats <- bstats[, site.order]
	c.values <- c.values[site.order]
	invisible(deviation.plot.beta.internal(bstats, c.values, c.legend))
}

deviation.plot.beta.internal <- function(bstats, c.values, c.legend, cuts = deviation.plot.beta.get.cuts(ncol(bstats))) {
	N <- ncol(bstats)
	additional.info <- NULL
	if (is.null(cuts)) {
		if (!is.null(c.values)) {
			additional.info <- c.legend[c.values]
		}
	} else {
		bstats <- sapply(tapply(1:N, cuts, function(x) {
					if (length(x) != 1) { return(rowMeans(bstats[, x])) }
					bstats[, x]
				}), identity)
		if (!is.null(c.values)) {
			if (is.character(c.legend)) {
				value.stats <- rep.int(0L, length(c.legend))
				names(value.stats) <- names(c.legend)
				additional.info <- sapply(tapply(c.values, cuts, function(x) {
							result <- value.stats
							values <- table(x)
							result[names(values)] <- values
							return(result / sum(result))
						}), cumsum)
				if (is.matrix(additional.info)) {
					rownames(additional.info) <- c.legend
				}
			} else { # is.integer(c.legend)
				if (c.legend == 2) {
					colors.grad <- rnb.getOption("colors.gradient")
					color.scheme <- gplots::colorpanel(20, low = colors.grad[1], high = colors.grad[2])
				} else { # c.legend == 3
					colors.grad <- rnb.getOption("colors.3.gradient")
					color.scheme <- gplots::colorpanel(20, low = colors.grad[1], mid = colors.grad[2], high = colors.grad[3])
				}
				v.range <- range(c.values, na.rm = TRUE)
				v.sections <- seq(v.range[1], v.range[2], length.out = 21)[-1]
				get.fractions <- function(x) { ecdf(x)(v.sections) }
				additional.info <- sapply(tapply(c.values, cuts, get.fractions), identity)
				if (is.matrix(additional.info)) {
					rownames(additional.info) <- color.scheme
				}
			}
			
		}
	}

	deviation.plot(bstats, ylim = c(0, 1), additional.info, mar = c(1, 4, 1, 1) + 0.1, ylab = expression(beta))
}


#### modified CGH profile plotting function from GLAD package


plotCGHProfile <- function(profileCGH, variable="LogRatio", Chromosome=NULL,
		Smoothing=NULL, GNL="ZoneGNL", Bkp=FALSE,
		labels=TRUE, plotband=TRUE, unit=0,
		colDAGLAD=c("black","blue","red","green","yellow"),
		pchSymbol=c(20,13),
		colCytoBand=c("white","darkblue"),
		colCentro="red", text=NULL, cytoband = NULL,
		main="", ylim=NULL, ...)
{
	
	if(is.null(cytoband))
	{
		stop("Error: cytoband must be provided")
	}
	
	if (length(intersect(names(profileCGH$profileValues),"PosBase"))<1)
	{
		stop("Error in plotProfile.profileCGH: PosBase is not available")
	}
	
	if (!is.null(Smoothing))
	{
		if (length(intersect(names(profileCGH$profileValues),Smoothing))<1)
		{
			print(paste("Warning in plotProfile.profileCGH:", Smoothing," is not available"))
		}
	}
	
	if (Bkp)
	{
		if (length(intersect(names(profileCGH$profileValues),"Breakpoints"))<1)
		{
			print("Warning in plotProfile.profileCGH: Breakpoints is not available")
		}
	}
	
	profileCGH$profileValues$VarToPlot <- profileCGH$profileValues[,variable]
	profileCGH$profileValues$Chromosome <- ChrNumeric(profileCGH$profileValues$Chromosome)
	
	ChrNum <- TRUE
	
	
	indexna <- attr(na.omit(profileCGH$profileValues[variable]),"na.action")
	if(!is.null(indexna))
	{
		profileCGH$profileValues <- profileCGH$profileValues[-indexna,]
	}
	
#    cytoband <- NULL
#    data("cytoband")
	
	if (!is.null(Chromosome))
	{
		ind <- NULL
		for (Chr in Chromosome)
		{
			indChr <- which(profileCGH$profileValues$Chromosome==Chr)
			ind <- c(ind,indChr)
		}
		profileCGH$profileValues <- profileCGH$profileValues[ind,]
		profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
	}
	
	
	LabelChr <- unique(na.omit(profileCGH$profileValues$Chromosome))
	NbChr <- length(LabelChr)
	LabelChr <- data.frame(Chromosome=LabelChr)
	
	### Information dans les cytobandes
	genomeInfo <- aggregate(cytoband$End, list(Chromosome=cytoband$Chromosome, ChrNumeric=cytoband$ChrNumeric), max, na.rm=TRUE)
	names(genomeInfo) <- c("Chromosome", "ChrNumeric", "Length")
	genomeInfo$Chromosome <- as.character(genomeInfo$Chromosome)
	genomeInfo$ChrNumeric <- as.integer(as.character(genomeInfo$ChrNumeric))
	
	if (ChrNum)
	{        
		LabelChr <- merge(LabelChr, genomeInfo[,c("ChrNumeric","Length")], by.x="Chromosome", by.y="ChrNumeric", all.x=TRUE)
		LabelChr <- LabelChr[order(LabelChr$Chromosome),]
	}
	else
	{
		LabelChr <- merge(LabelChr, genomeInfo, by="Chromosome", all.x=TRUE)
		LabelChr <- LabelChr[order(LabelChr$ChrNumeric),]
	}
	
	
	
	if (NbChr > 1)
	{
		gap <- 100000000/3
		LabelChr$Length <- LabelChr$Length + gap
		
		cumulLength <- cumsum(LabelChr$Length)
		LabelChr$Length <- c(0, cumulLength[1:(NbChr-1)])                                                             
		LabelChr$Length <- LabelChr$Length/(10^unit)
	}
	
	else
	{
		LabelChr$Length <- 0
	}
	
	
	if (ChrNum)
	{
		cytobandNew <- subset(cytoband, select=setdiff(names(cytoband),"Chromosome"))
		cytobandNew <- merge(LabelChr, cytobandNew, by.x="Chromosome", by.y="ChrNumeric")
	}
	else
	{
		cytobandNew <- subset(cytoband, select=setdiff(names(cytoband),"ChrNumeric"))
		cytobandNew <- merge(LabelChr, cytobandNew, by="Chromosome") 
	}
	
	cytobandNew$Start <- cytobandNew$Start/(10^unit)
	cytobandNew$End <- cytobandNew$End/(10^unit)
	
	cytobandNew$Start <- cytobandNew$Start + cytobandNew$Length
	cytobandNew$End <- cytobandNew$End +  cytobandNew$Length
	
	profileCGH$profileValues <- merge(profileCGH$profileValues, LabelChr, by="Chromosome")
	
	profileCGH$profileValues$NewPosBase <- profileCGH$profileValues$PosBase + profileCGH$profileValues$Length
	
	def.par <- par(no.readonly = TRUE)
	
	##    if (plotband)
	##       {
	##         layout(c(1,2), heights=c(4,1))
	##         par(mar=c(0,4,5,2))
	##       }
	
	if (plotband)
	{
		### Cytobandes
		
		layout(c(1,2), heights=c(1,4))
		
		par(mar=c(0,4,4,2))
		
		plot(0, type="n", xlim=c(0, max(cytobandNew$End)),
				ylim=c(-1.5,1.5), xaxt="n", yaxt="n", ylab="", xlab="")
		
		LabelChrCyto <- unique(cytobandNew$Chromosome)
		
		
		for (i in 1:NbChr)
		{
			GLAD::plotCytoBand(cytobandNew, Chromosome=LabelChrCyto[i], labels=labels, y=0, height=2, colCytoBand=colCytoBand, colCentro=colCentro)
		}
		
		par(mar=c(4,4,0,2))
	}
	
	if (!is.null(Smoothing))
	{
		profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$Chromosome,profileCGH$profileValues$PosBase),]
		NbPos <- length(profileCGH$profileValues[,1])
		PosMax <- max(profileCGH$profileValues$NewPosBase) + 1
		Pos <- profileCGH$profileValues$NewPosBase[1:(NbPos-1)]
		PosNext <- profileCGH$profileValues$NewPosBase[2:NbPos]
		InterPos <- Pos + (PosNext-Pos)/2
		InterPos <- c(0, InterPos, PosMax)
		
		SmtStart <- profileCGH$profileValues[,Smoothing][1]
		SmtEnd <- profileCGH$profileValues[,Smoothing][NbPos]
		
		Smt1 <- profileCGH$profileValues[,Smoothing][1:(NbPos-1)]
		Smt1 <- c(SmtStart, Smt1, SmtEnd)
		
		Smt2 <- profileCGH$profileValues[,Smoothing][2:NbPos]
		Smt2 <- c(SmtStart, Smt2, SmtEnd)
		
		
		datasmt <- data.frame(PosBase=c(InterPos,InterPos),Smoothing=c(Smt1,Smt2))
		datasmt <- unique(datasmt)
		datasmt <- datasmt[order(datasmt$PosBase),]
	}
	
	if (length(intersect(names(profileCGH$profileValues),GNL))>=1)
	{
		
		col <- rep(colDAGLAD[5],length(profileCGH$profileValues$PosOrder))
		col[which(profileCGH$profileValues[GNL]==-1)] <- colDAGLAD[4]
		col[which(profileCGH$profileValues[GNL]==1)] <- colDAGLAD[3]
		col[which(profileCGH$profileValues[GNL]==2)] <- colDAGLAD[2]
		col[which(profileCGH$profileValues[GNL]==-10)] <- colDAGLAD[1]
		
		
		outliers <- rep(pchSymbol[1],length(profileCGH$profileValues$PosOrder))
		outliers[which(profileCGH$profileValues$OutliersTot!=0)] <- pchSymbol[2]
		
		if (plotband)
		{
			
			plot(VarToPlot ~ NewPosBase, data=profileCGH$profileValues,
					pch=outliers, col=col, xaxt="n", xlab=main,
					ylab=variable, ylim=ylim, xlim=c(0,max(cytobandNew$End)),...)
		}
		else
		{
			plot(VarToPlot ~ NewPosBase, data=profileCGH$profileValues,
					pch=outliers, col=col, xaxt="n", xlab="", ylab=variable, main=main, ylim=ylim, ...)
		}
		
		
		if (!is.null(Smoothing))
		{
			lines(datasmt$Smoothing ~ datasmt$PosBase, col="black")
		}
		
		if (Bkp)
		{
			if (is.data.frame(profileCGH$BkpInfo))
			{
				profileCGH$BkpInfo <- merge(profileCGH$BkpInfo, LabelChr, by="Chromosome")
				profileCGH$BkpInfo$NewPosBase <- profileCGH$BkpInfo$PosBase + profileCGH$BkpInfo$Length
				abline(v=profileCGH$BkpInfo$NewPosBase+0.5, col="red", lty=2)
			}
		}
	}
	else
	{
		if (plotband)
		{
			plot(VarToPlot ~ NewPosBase, data=profileCGH$profileValues,
					pch=20, xaxt="n", xlab=main, ylab=variable, ylim=ylim, xlim=c(0,max(cytobandNew$End)), ...)
		}
		else
		{
			plot(VarToPlot ~ NewPosBase, data=profileCGH$profileValues,
					pch=20, xaxt="n", xlab="", ylab=variable, main=main, ylim=ylim, ...)
		}
		
		if (Bkp)          
		{
			if (is.data.frame(profileCGH$BkpInfo))
			{
				profileCGH$BkpInfo <- merge(profileCGH$BkpInfo, LabelChr, by="Chromosome")
				profileCGH$BkpInfo$NewPosBase <- profileCGH$BkpInfo$PosBase + profileCGH$BkpInfo$Length
				abline(v=profileCGH$BkpInfo$NewPosBase+0.5, col="red", lty=2)
			}
		}
		
		if (!is.null(Smoothing))
		{
			lines(datasmt$Smoothing ~ datasmt$PosBase, col="red")
		}
	}
	
	if (!is.null(text))
	{
		text(text$x, text$y, labels=text$labels, cex=text$cex)
	}
	
#    par(def.par)
}

########################################################################################################################

#' densRanks
#'
#' Rank the points accordind to density of the region they fall in. Densities are computed
#' as Kernel Density estimates. The method and parameters are implemented in analogy to
#' \code{grDevices::densCols}
#' @param x x-coordinate
#' @param y y-coordinate
#' @param nbin number of bins
#' @param bandwidth bandwidth
#' @author Fabian Mueller
densRanks <- function (x, y = NULL, nbin = 128, bandwidth) 
{
    xy <- xy.coords(x, y)
    select <- is.finite(xy$x) & is.finite(xy$y)
    x <- cbind(xy$x, xy$y)[select, ]
    map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
    mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
    xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
    ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
    dens <- map$fhat[cbind(xbin, ybin)]
    dens[is.na(dens)] <- 0

    res <- rep(NA_integer_, length(select))
    rrs <- rank(dens,ties.method="max")
    res[select] <- rrs
    res
}

########################################################################################################################

#' create.densityScatter
#' 
#' Creates a density scatterplot highlighting points in sparsely populated plot regions
#' as well as points marked as special in a seperate color
#' @param df2p \code{data.frame} to be plotted. Only the fist two columns are taken into account as
#' 			x and y coordinates respectively
#' @param is.special boolean vector of length equal to the number of rows in \code{df2p}. Specifies
#' 			which points should be highlighed seperately in a different color
#' @param dens.subsample if the number of points exceeds this number, subsample the number of points for the
#'			density estimation to that number. Any non-numeric value disables subsampling.
#' @param dens.special Flag indicating whether the points of the special population should be colored
#' 			according to their density
#' @param sparse.points Either percentage (\code{<=1,>=0}) or the absolute number 
#'          of points in the sparsely populated area that should be drawn seperately. A value of 0 means that these points
#'			will not be drawn.
#' @param dens.n passed on to \code{ggplot2::stat_density2d}: argument: \code{n}
#' @param add.text.cor flag indicating whether a text token with the correlation coefficient should be included in the lower
#'          right corner of the plot
#' @return \code{ggplot} object
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' d <- data.frame(x=rnorm(1000),y=rnorm(1000))
#' s <- rep(FALSE,1000)
#' s[sample(1:length(s),100)] <- TRUE
#' create.densityScatter(d,s)
#' }
create.densityScatter <- function(df2p,is.special=NULL,dens.subsample=FALSE,dens.special=TRUE,
		sparse.points=0.01,dens.n=100,add.text.cor=FALSE){
	if (!(is.numeric(sparse.points) && sparse.points>=0)) {
		stop("Invalid parameter value: sparse.points")
	}
	if (!is.null(is.special)) is.special[is.na(is.special)] <- FALSE
	if (sum(is.special)<1){
		is.special <- NULL
	}
	if (!is.null(is.special)){
		df2p$is.special <- is.special
	}
	if (is.null(df2p) || nrow(df2p)<1){
		logger.warning(c("Could not create density scatterplot"))
		pp <- rnb.message.plot("Could not create plot")
		return(pp)
	}
	df2p <- na.omit(df2p)
	if (is.null(df2p) || nrow(df2p)<1){
		logger.warning(c("Could not create density scatterplot (NA omission removed all entries)"))
		pp <- rnb.message.plot("Could not create plot")
		return(pp)
	}
	df2p.sub <- df2p
	dens.ranks <- NULL
	tryCatch(
		dens.ranks <- densRanks(x=df2p[,1],y=df2p[,2]),
		error=function(ee){
			logger.warning(c("Could not assess density ranking:",ee$message))
		}
	)
	if (is.numeric(dens.subsample) && dens.subsample>0){
		ss <- as.integer(dens.subsample)
		if (nrow(df2p) > ss) {
			df2p.sub <- df2p[sample(nrow(df2p),ss),]
		}
	}

	#the standard bandwith function of MASS::kde2d is unstable when looking at
	#distributions with very low variance. Here's a more stable version
	stable.bandwidth.fun <- function(x,eps=1e-4){
	    r <- quantile(x, c(0.25, 0.75))
	    h <- (r[2] - r[1])/1.34
	    if (h==0) h <- eps
	    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
	}
	stable.h <- c(stable.bandwidth.fun(df2p.sub[,1]),stable.bandwidth.fun(df2p.sub[,2]))

	if (is.null(dens.ranks)){
  		pp <- rnb.message.plot("Could not assess density")
  	} else {		
  		pp <- ggplot(df2p.sub) + aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2]) + 
		  stat_density2d(geom="tile", fill=DENS.COLORS.LOW[1], aes(,alpha=..density..^0.25), contour=FALSE, n=dens.n, h=stable.h) +
		  scale_alpha(range = c(0.0, 1))
		if (sparse.points > 0){
			if (sparse.points <= 1){
				thres <- ceiling(nrow(df2p)*sparse.points)
			} else {
				thres <- sparse.points
			}
			df2p.loose <- df2p[dens.ranks<=thres,]#the sub data.frame in of the least dens points
			pp <- pp + geom_point(data=df2p.loose,aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2]),colour=DENS.COLORS.LOW[1],size=0.4)
		}
		if (!is.null(is.special)){
			df2p.special <- df2p[df2p$is.special,]
			colors.dmp <- DENS.COLORS.LOW[2]
			if (dens.special){
				tryCatch(
					colors.dmp   <- densCols(x=df2p.special[,1],y=df2p.special[,2],colramp = colorRampPalette(c(DENS.COLORS.LOW[2],DENS.COLORS.HIGH[2]))),
					error=function(ee){
						logger.warning(c("Could not assess density colors using densCols:",ee$message))
					}
				)
			}
			df2p.special$color <- colors.dmp

			pp <- pp + geom_point(data=df2p.special,aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2],colour="color"),size=1) + scale_color_identity()
		}
		if (add.text.cor) {
			cc <- cor(df2p[,1],df2p[,2],use="pairwise.complete.obs")
			txt.cor <- paste0('rho',paste0("==",round(cc,4)))
			pp <- pp + annotate("text", x=max(df2p[,1],na.rm=TRUE),y=min(df2p[,2],na.rm=TRUE),label=txt.cor,parse=TRUE,hjust=1,vjust=1,size=4)
		}
		pp <- pp + theme(legend.position="none")
  	}
	return(pp)
}

########################################################################################################################
#' create.scatter.dens.points
#' 
#' Creates a scatterplot containing all points in a given data.frame. Points are colored according to point
#' density. Optionally, a selection of points are shown in a different color
#' @param df2p \code{data.frame} to be plotted. Only the fist two columns are taken into account as
#' 			x and y coordinates respectively
#' @param is.special boolean vector of length equal to the number of rows in \code{df2p}. Specifies
#' 			which points should be highlighed seperately in a different color
#' @param dens.special Flag indicating whether the points of the special population should be colored
#' 			according to their density
#' @param mock Should only the axis be plotted? useful when exporting scatterplots with lots of points
#'        as immage and the corresponding axis as vector graphics.
#' @return \code{ggplot} object
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' d <- data.frame(x=rnorm(1000),y=rnorm(1000))
#' s <- rep(FALSE,1000)
#' s[sample(1:length(s),100)] <- TRUE
#' create.scatter.dens.points(d,s)
#' }
create.scatter.dens.points <- function(df2p,is.special=NULL,dens.special=TRUE,mock=FALSE){
	if (!is.null(is.special)) is.special[is.na(is.special)] <- FALSE
	if (sum(is.special)<1){
		is.special <- NULL
	}
	if (!is.null(is.special)){
		df2p$is.special <- is.special
	}
	df2p <- na.omit(df2p)

	#plot order: plot DMRs last
	n.points <- nrow(df2p)
	df2p$plotOrder <- NA
	num.not.special <- sum(!df2p$is.special)
	df2p$plotOrder[!df2p$is.special] <- seq_len(num.not.special)
	df2p$plotOrder[df2p$is.special] <- seq((num.not.special+1),n.points)
	
	df2p$color <- NA
	if (sum(!df2p$is.special)>1){
		colors.nodmr <- DENS.COLORS.LOW[1]
		tryCatch(
			colors.nodmr <- densCols(x=df2p[!df2p$is.special,1],y=df2p[!df2p$is.special,2],colramp = colorRampPalette(c(DENS.COLORS.LOW[1],DENS.COLORS.HIGH[1]))),
			error=function(ee){
				logger.warning(c("Could not assess density colors using densCols:",ee$message))
			}
		)
		df2p[!df2p$is.special,"color"] <- colors.nodmr
	} else if (sum(!df2p$is.special)==1){
		df2p[!df2p$is.special,"color"] <- DENS.COLORS.LOW[1]
	}
	if (sum(df2p$is.special)>1){
		colors.dmr <- DENS.COLORS.LOW[2]
		if (dens.special) {
			tryCatch(
				colors.dmr   <- densCols(x=df2p[ df2p$is.special,1],y=df2p[ df2p$is.special,2],colramp = colorRampPalette(c(DENS.COLORS.LOW[2],DENS.COLORS.HIGH[2]))),
				error=function(ee){
					logger.warning(c("Could not assess density colors using densCols:",ee$message))
				}
			)
		}
		df2p[df2p$is.special,"color"] <- colors.dmr
		
	} else if (sum(df2p$is.special)==1){
		df2p[df2p$is.special,"color"] <- DENS.COLORS.LOW[2]
	}
	pp <- ggplot(df2p) + aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2])
	if (mock){
		pp <- pp + geom_blank()
	} else {
		pp <- pp + geom_point(aes_string(color="color",order="plotOrder")) + scale_color_identity()
	}
	
	return(pp)
}

########################################################################################################################
#' create.hex.summary.plot
#' 
#' Creates a summary plot binning the data given by a certain quantity in heagonal bins
#' @param df2p \code{data.frame} to be plotted.
#' @param x name of the variable in \code{df2p} considered as x-axis
#' @param y name of the variable in \code{df2p} considered as y-axis
#' @param q name of the variable in \code{df2p} considered as quantity to be summarized over bins
#' @param bins,fun,... arguments to be passed on to \code{stat_summary_hex}
#' @return \code{ggplot} object
#' @author Fabian Mueller
create.hex.summary.plot <- function(df2p,x=colnames(df2p)[1],y=colnames(df2p)[2],q=colnames(df2p)[3],
	bins=128,fun=median,...){
	pp <- ggplot(df2p) + aes_string(x=x,y=y,z=q) +
			stat_summary_hex(bins=bins,fun=fun,...) + scale_fill_gradientn(colours=rev(rnb.getOption("colors.gradient")))
	return(pp)
}

