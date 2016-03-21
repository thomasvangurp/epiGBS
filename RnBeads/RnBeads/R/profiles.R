########################################################################################################################
## profiles.R
## created: 2012-05-28
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the methylation profiles step of the analysis pipeline.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.section.sample.groups
#'
#' Creates leading section ("Sample Groups") in the report on methylation profiles.
#'
#' @param rnb.set Methylation dataset of interest.
#' @param report  Report to contain the Sample Groups section. 
#' @return List of two elements:
#'         \describe{
#'           \item{\code{"sample.inds"}}{List of traits defining subgroups in the dataset, as returned by
#'                \code{\link{rnb.sample.groups}}.}
#'           \item{\code{"report"}}{The modified report.}
#'         }
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.sample.groups <- function(rnb.set, report) {
	predefined.columns <- rnb.getOption("exploratory.columns")
	sample.inds <- rnb.sample.groups(rnb.set, columns = predefined.columns)
	pheno.columns <- data.frame(
		"Trait" = names(sample.inds),
		"Number of groups" = sapply(sample.inds, length), check.names = FALSE, stringsAsFactors = FALSE)
	if (is.null(predefined.columns)) {
		if (is.null(pheno(rnb.set))) {
			txt <- "Since no phenotypic data are available, no predefined sample groups were examined."
		} else if (nrow(pheno.columns) == 0) {
			txt <- c("The table of phenotypic information contains no traits that meet the criteria for defining ",
				"sample groups.")
		} else {
			txt <- c("Traits in the table of phenotypic information were automatically selected based on criteria for ",
				"defining sample groups.")
		}
	} else {
		if (is.null(pheno(rnb.set))) {
			txt <- "Since no phenotypic data are available, no predefined sample groups were examined."
		} else if (nrow(pheno.columns) == 0) {
			txt <- c("None of the predefined traits meets the criteria for defining sample groups.")
		} else {
			txt <- c("The specified traits were tested based on criteria for defining sample groups.")
		}
	}
	
	if (nrow(pheno.columns) != 0) {
		## Display a summary table of traits
		txt <- c(txt, " The table below summarizes these traits.")
		report <- rnb.add.section(report, "Sample Groups", txt)
		rnb.add.table(report, pheno.columns, row.names = FALSE)
	} else {
		report <- rnb.add.section(report, "Sample Groups", txt)
	}
	rnb.status("Designed color mappings for probe type and CGI status")
	return(list(sample.inds = sample.inds, report = report))
}

########################################################################################################################

## Creates a density estimation plot of beta values for the specified sample and/or probe group.
##
## @param dframe        \code{data.frame} containing at least 2 variables:
##						\code{beta:} Vector of beta values of the distributions to be displayed.
##                      \code{group:} Vector of group labels for the values.
## @param annotation    Name of the annotation being visualized, in the form of a \code{character} vector of length 1.
## @param cvalues       Vector storing the color mapping to be applied based on the groups. This is passed to the
##                      function \code{\link{scale_color_manual}}.
## @param legend.margin Length, in inches, of the space reserved for the plot legend. Legend is placed to the right of
##                      the plot.
## @return The generated ggplot object.
## @author Fabian Mueller
rnb.plot.beta.density.group <- function(dframe, annotation, cvalues, legend.margin = 2.2) {
	
	if (length(cvalues) > 1) {
		pp <- ggplot(dframe, aes(x = beta, color = group)) +
			labs(x = expression(beta), y = "Density", color = annotation) +
			geom_density() + scale_color_manual(values = cvalues) +
			theme(plot.margin = unit(c(1, legend.margin, 0.5, 0.5), c("lines", "in", "lines", "lines")),
				legend.justification = c(0, 0.5), legend.position = c(1, 0.5))
	} else {
		pp <- ggplot(dframe, aes(x = beta)) + labs(x = expression(beta), y = "Density") + geom_density() +
			theme(plot.margin = unit(c(1, legend.margin, 0.5, 0.5), c("lines", "in", "lines", "lines")))
	}
	return(pp)
}

########################################################################################################################

## transform a list of disjoint integer vectors to one factor vector. The result is a vector
## of length n where the indices specified in the input list
## correspond to the level specified by the list element name
## Useful for transforming result from \code{rnb.sample.groups} to factor vectors
##
## @param ll   index list
## @param n    length of the output vector
## @return factor vector with levels corresponding to names(ll)
## @author Fabian Mueller
index.list.2.factor <- function(ll,n){
	res <- rep(NA,n)
	ggs <- names(ll)
	if(is.null(ggs)) ggs <- 1:length(ll)
	inds.names.na <- which(is.na(ggs)||(ggs==""))
	ggs[inds.names.na] <- inds.names.na
	for(i in 1:length(ll)){
		if (any(!is.na(res[ll[[i]]]))) stop("Non-disjoint index sets")
		res[ll[[i]]] <- ggs[i]
	}
	return(factor(res,levels=ggs))
}

########################################################################################################################

## Given a metrix of beta values and a vector of sample or site group memberships,
## get a data.frame that can be passed on to \code{rnb.plot.beta.density.group}. 
## NA values are ommitted.
## If subsampling, is enabled (i.e. \code{points.per.group}>0), 
## observations per group are subsampled according to the following procedure:
## Given a data.frame with K groups and number of observations
## N_1,...,N_K, target number of points per group T, the total number of points N = sum(N_1,...,N_K) is computed
## Afterwards the proportions p_k = N_k/N is computed and from each group, S_k = p_k*(K*T) observations
## are randomly selected from all observations belonging to group k.
## @param beta.matrix matrix containing the values to be plotted. in this package, mostly 
##                    beta methylation values
## @param groups      group memberships as factor vector for either the columns or rows of \code{beta.matrix}
## @param byrow		  dows the grouping corrospond to the rows of beta.matrix. If so, beta.matrix will be transposed
##					  and the rest of function is handled in complete analogy to the column case
## @param log.name    character to print as name for the logger. If NULL, logging is ommited
## @param points.per.group the targeted number of points (T) per group. Set this to a value < 1 to disable subsampling
## @return \code{data.frame} containing at least 2 variables:
##				\code{beta:} Vector of beta values of the distributions to be displayed.
##              \code{group:} Vector of group labels for the values.
## @author Fabian Mueller
get.density.dframe <- function(beta.matrix,groups,byrow=FALSE,log.name=NULL,points.per.group=0){
	if (byrow){
		beta.matrix <- t(beta.matrix)
	}
	if (length(groups)!=ncol(beta.matrix)){
		stop("Grouping does not match beta.matrix")
	}
	ggs <- levels(groups)
	n.grps <- length(ggs)
	grp.ncol <- table(groups)
	group.mats <- lapply(ggs,FUN=function(gg){
		beta.matrix[,groups==gg]
	})
	names(group.mats) <- ggs
	grp.n.obs <- sapply(group.mats,FUN=function(x){
		sum(!is.na(x))
	})
	N <- sum(grp.n.obs)
	
	if (!is.null(log.name)){
		grp.na.str <- paste(paste0(ggs,":",grp.n.obs,"/",grp.ncol*nrow(beta.matrix)),collapse="; ")
		logger.info(c("Density estimation (",log.name,"): Groupwise retained observations after missing value removal:",grp.na.str))
	}
	
	N.sub <- points.per.group*n.grps
	do.subsample <- points.per.group>0 && N.sub < N
	
	grp.beta.list <- lapply(ggs,FUN=function(gg){
		#select observations based on not missing and subsampling
		obs.sel <- which(!is.na(group.mats[[gg]]))
		if (do.subsample){
			p.sub <- N.sub/N
			target.num <- round(grp.n.obs[gg]*p.sub)
			if (target.num < grp.n.obs[gg]) {
				obs.sel <- sample(obs.sel,target.num)
			}
		}
		as.vector(group.mats[[gg]][obs.sel])
	})
	names(grp.beta.list) <- ggs
	grp.n.obs.subsample <- sapply(grp.beta.list,length)
	if (!is.null(log.name) && do.subsample){
		grp.ss.str <- paste(paste0(ggs,":",grp.n.obs.subsample,"/",grp.n.obs),collapse="; ")
		logger.info(c("Density estimation (",log.name,"): Groupwise retained observations after subsampling:",grp.ss.str))
	}
	return(data.frame(beta=unlist(grp.beta.list,use.names=FALSE),
				group=factor(rep(ggs,grp.n.obs.subsample),levels=ggs)))
}

########################################################################################################################

#' rnb.plot.betadistribution.sampleGroups
#'
#' Plots beta value distrubions given a sample grouping.
#'
#' @param beta.matrix       Beta values in the form of a non-empty \code{matrix} of type \code{double}. Rows in this
#'                          matrix must correspond to Infinium probes, and columns - to samples.
#' @param sample.group.inds Named \code{list} that contains indices for the samples contained in the groups in
#' 						    \code{beta.matrix}. The number of groups is determined by the length of the list, and its
#'                          names are used as group names.
#' @param annotation        Name of the annotation being visualized, in the form of a \code{character} vector of length 1.
#' @param log.str		    string specifying more details for the log file
#' @param points.per.group the targeted number of points per group. Set this to a value < 1 to disable subsampling. More
#' 				  information in the Details section of \code{\link{rnb.step.betadistribution}}
#' @return the plot as a \code{ggplot2} object
#'
#' @seealso rnb.plot.betadistribution.probeCategories
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' meth.mat <- meth(rnb.set.example)
#' sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
#' rnb.plot.betadistribution.sampleGroups(meth.mat,sample.groups)
#' }
rnb.plot.betadistribution.sampleGroups <- function(beta.matrix, sample.group.inds, annotation = "Group",
		log.str=NULL,points.per.group=rnb.getOption("distribution.subsample")) {
	if (!(is.matrix(beta.matrix) && is.double(beta.matrix) && nrow(beta.matrix) * ncol(beta.matrix) != 0)) {
		stop("invalid value for beta.matrix")
	}
	## TODO: Validate sample.group.inds
	if (!(is.character(annotation) && length(annotation) == 1 && (!is.na(annotation)))) {
		stop("invalid value for annotation")
	}
	grps <- index.list.2.factor(sample.group.inds,ncol(beta.matrix))
	dframe <- get.density.dframe(beta.matrix,grps,byrow=FALSE,log.name=paste(c(annotation,log.str),collapse="--"),
				points.per.group=points.per.group)
	cvalues <- rep(rnb.getOption("colors.category"), length.out = length(sample.group.inds)) #adjust length of the color vector to the appropriate size
	rnb.plot.beta.density.group(dframe, annotation, cvalues)
}

########################################################################################################################

#' rnb.plot.betadistribution.probeCategories
#'
#' plot beta value distrubions given probe categories
#'
#' @param beta.matrix  Beta values in the form of a non-empty \code{matrix} of type \code{double}. Rows in this matrix
#'                     must correspond to Infinium probes, and columns - to samples.
#' @param probe.cat    \code{factor} vector of length \code{nrow(beta.matrix)} corresponding to the
#'                     probe categories.
#' @param annotation   Name of the annotation being visualized, in the form of a \code{character} vector of length 1.
#' @param color.legend Color legend to use in the form of a \code{character} vector with element names. The values in
#'                     this vector should encode colors. All values in \code{probe.cat} must be present in the names of
#'                     this color legend. If this parameter is \code{NULL}, a default color legend is be constructed.
#' @param log.str	   string specifying more details for the log file
#' @param points.per.group the targeted number of points per group. Set this to a value < 1 to disable subsampling. More
#' 				  information in the Details section of \code{\link{rnb.step.betadistribution}}
#' @return The plot as a \code{ggplot2} object.
#'
#' @seealso rnb.plot.betadistribution.sampleGroups
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' meth.mat <- meth(rnb.set.example)
#' probe.types <- annotation(rnb.set.example)[, "Design"]
#' rnb.plot.betadistribution.probeCategories(meth.mat,probe.types,annotation="Infinium probe type")
#' }
rnb.plot.betadistribution.probeCategories <- function(beta.matrix, probe.cat, annotation = "Group", color.legend = NULL,
		log.str=NULL,points.per.group=rnb.getOption("distribution.subsample")) {
	if (!(is.matrix(beta.matrix) && is.double(beta.matrix) && nrow(beta.matrix) * ncol(beta.matrix) != 0)) {
		stop("invalid value for beta.matrix")
	}
	if (!(is.factor(probe.cat) && length(probe.cat)==nrow(beta.matrix))) {
		stop("invalid value for probe.cat")
	}
	if (!(is.character(annotation) && length(annotation) == 1 && (!is.na(annotation)))) {
		stop("invalid value for annotation")
	}
	if (is.null(color.legend)) {
		## Use default palette if none is provided
		color.legend <- rep(rnb.getOption("colors.category"), length.out = length(unique(probe.cat)))
	} else if (!(is.character(color.legend) && length(color.legend) != 0 && (!any(is.na(color.legend))))) {
		stop("invalid value for color.legend")
	}
	dframe <- get.density.dframe(beta.matrix,probe.cat,byrow=TRUE,log.name=paste(c(annotation,log.str),collapse="--"),
					points.per.group=points.per.group)
	rnb.plot.beta.density.group(dframe,annotation,color.legend)
}
	
########################################################################################################################

#' rnb.step.betadistribution
#'
#' Computes the distributions of beta values across various sample groups and adds a corresponding section to the
#' report.
#'
#' @param rnb.set HumanMethylation450K dataset as an object of type \code{\linkS4class{RnBSet}}.
#' @param report  Report to contain the methylation deviation section. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @param columns Optional; predefined column names (in the form of a \code{character} vector) or indices (an
#'                \code{integer} vector) in the sample annotation table. Only these columns are considered for grouping
#'                samples and defining profiles. All other columns in the phenotype table are ignored.
#' @param points.per.group the targeted number of points (T) per group. Set this to a value < 1 to disable subsampling. More
#' 				  information in the Details section
#' @return The modified report.
#' @section Details:
#' If subsampling is enabled (i.e. \code{points.per.group}>0), 
#' observations per group are subsampled according to the following procedure:
#' Given K groups and numbers of observed beta values per group
#' N_1,...,N_K, and the target number of points per group T: the total number of points N = sum(N_1,...,N_K) is computed
#' Afterwards the proportions p_k = N_k/N is computed and from each group, S_k = p_k*(K*T) observations
#' are randomly selected from all observations belonging to group k.
#'
#' @author Fabian Mueller
#' @export
rnb.step.betadistribution <- function(rnb.set, report, columns = rnb.getOption("exploratory.columns"),
		points.per.group=rnb.getOption("distribution.subsample")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	## rnb.sample.groups validates the value of columns
	sample.inds <- rnb.sample.groups(rnb.set, columns = columns)
	if (rnb.getOption("analyze.sites") == FALSE) {
		return(report)
	}
	logger.start("Methylation Value Distribution Plots")
	pinfos <- get.site.and.region.types(rnb.set)[[1]]
	return(rnb.step.betadistribution.internal(rnb.set, report, sample.inds, pinfos, points.per.group))
}

rnb.step.betadistribution.internal <- function(rnb.set, report, sample.inds, pinfos,
		points.per.group=rnb.getOption("distribution.subsample")) {
#	init.plot <- function(fname) {
#		createReportPlot(fname, report, width = 7.8, height = 5.6, high.png = 200L)
#	}
	do.ggplot <- function(pp,fname) {
		res <- createReportGgPlot(pp, fname, report, width = 7.8, height = 5.6, high.png = 200L)
		res <- off(res,handle.errors=TRUE)
		return(res)
	}
	## Keep only site annotations that define categories
	is.category <- sapply(attr(pinfos, "legend"), is.character)
	pinfos.legend <- attr(pinfos, "legend")[is.category]
	pinfos <- as.data.frame(pinfos[, is.category])
	colnames(pinfos) <- names(pinfos.legend)
	attr(pinfos, "legend") <- pinfos.legend
	rm(is.category, pinfos.legend)

	X <- meth.matrices(rnb.set)
	sample.inds.ext <- c(list("all samples" = list("all" = 1:ncol(X[[1]]))), sample.inds)

	txt.site <- rnb.get.row.token(rnb.set)
	txt.covered <- paste(txt.site, ifelse(length(X) != 1, " and region", ""), sep = "")
	txt <- c("Methylation value distributions were assessed based on selected sample groups. This was done on ",
			txt.covered, " levels. This section contains the generated density plots.")
	report <- rnb.add.section(report, "Methylation Value Distributions", txt)

	## Beta value distributions for different sample groups
	logger.start("Methylation Value Distributions - Sample Groups")
	param.combinations <- mapply(c, rep(1:length(sample.inds.ext), each = length(X)),
		rep(1:length(X), length(sample.inds.ext)), SIMPLIFY = FALSE)
	create.plot <- function(params) {
		i <- params[1]
		j <- params[2]
		pname <- paste("beta_density_samples", i, sep = "_")
		if (length(X) != 1) {
			pname <- paste(pname, j, sep = "_")
		}
		logger.info(c("processing",pname))
#		rplot <- init.plot(pname)
#		print(rnb.plot.betadistribution.sampleGroups(X[[j]], sample.inds.ext[[i]], names(sample.inds.ext)[i]))
#		off(rplot)
		pp <- rnb.plot.betadistribution.sampleGroups(X[[j]], sample.inds.ext[[i]], names(sample.inds.ext)[i],
				log.str=names(X)[j],points.per.group=points.per.group)
		rplot <- do.ggplot(pp,pname)
		rnb.cleanMem()
		return(rplot)
	}
	if (parallel.isEnabled()) {
		figure.plots <- foreach(params = param.combinations,
			.export = c("logger.info", "rnb.plot.betadistribution.sampleGroups", "rnb.cleanMem")) %dopar% create.plot(params)
	} else {
		figure.plots <- lapply(param.combinations, create.plot)
	}
	txt <- c("The plots below compare the distributions of methylation values in different sample groups, as defined by the ",
		"traits listed above.")
	report <- rnb.add.section(report, "Methylation Value Densities of Sample Groups", txt, level = 2)
	groupings <- names(sample.inds.ext)
	names(groupings) <- 1:length(groupings)
	setting.names <- list("Sample trait" = groupings)
	if (length(X) != 1) {
		i <- length(setting.names) + 1
		setting.names[["Methylation of"]] <- c(rnb.get.row.token(rnb.set, plural = TRUE), names(X)[-1])
		names(setting.names[[i]]) <- 1:length(setting.names[[i]])
	}
	description <- "Beta value density estimation according to sample grouping."
	report <- rnb.add.figure(report, description, figure.plots, setting.names)
	logger.completed()

	## Beta value distributions for different site types
	logger.start(c("Methylation Value Distributions -", capitalize(txt.site), "Categories"))
	site.categories <- colnames(pinfos)
	if(ncol(pinfos)>0){
		names(site.categories) <- 1:length(site.categories)
	}

	#skip this section if all the site annotations are either not categorical or only have 1 category
	skip <- all(sapply(pinfos,FUN=function(x){
		!is.factor(x) || length(unique(x)) < 2
	}))
	if (skip){
		logger.info(c("Site categories are non-categorical. --> skipped"))
	} else {
		if (length(sample.inds) != 0) {
			sample.inds.ext <- unlist(sample.inds, recursive = FALSE, use.names = FALSE)
			names(sample.inds.ext) <- paste(unlist(lapply(sample.inds, names), use.names = FALSE),
					rep(paste("(based on ", names(sample.inds), ")", sep = ""), sapply(sample.inds, length)))
		} else {
			sample.inds.ext <- list()
		}
		sample.inds.ext <- c(list("all samples" = 1:ncol(X[[1]])), sample.inds.ext)

		figure.plots <- list()
		for (i in 1:length(sample.inds.ext)){
			bbs <- as.matrix(X[[1]][, sample.inds.ext[[i]]])
			if (parallel.isEnabled()) {
				figure.plots <- c(figure.plots, foreach(j = 1:length(site.categories)) %dopar% {
					fname <- paste("beta_density_sites", i, j, sep = "_")
#					rplot <- init.plot(fname)
					pp <- rnb.plot.betadistribution.probeCategories(bbs, pinfos[,j],
						annotation = colnames(pinfos)[j], attr(pinfos, "legend")[[j]],
						log.str=names(sample.inds.ext)[i], points.per.group=points.per.group)
#					print(pp)
#					off(rplot)
					rplot <- do.ggplot(pp,fname)
					rplot
				})
			} else {
				for (j in 1:length(site.categories)){
					fname <- paste("beta_density_sites", i, j, sep = "_")
#					rplot <- init.plot(fname)
					pp <- rnb.plot.betadistribution.probeCategories(bbs, pinfos[,j],
						annotation = colnames(pinfos)[j], attr(pinfos, "legend")[[j]],
						log.str=names(sample.inds.ext)[i], points.per.group=points.per.group)
#					print(pp)
#					off(rplot)
					rplot <- do.ggplot(pp,fname)
					figure.plots <- c(figure.plots,list(rplot))
				}
			}
		}

		txt <- c("In a similar fashion, the plot below compares the distributions of beta values in different ",
			txt.site, " types.")
		report <- rnb.add.section(report, paste("Methylation Value Densities of", capitalize(txt.site), "Categories"),
			txt, level = 2)

		groups <- names(sample.inds.ext)
		names(groups) <- 1:length(groups)
		setting.names <- list("Sample group" = groups, "category" = site.categories)
		names(setting.names)[2] <- paste(capitalize(txt.site), names(setting.names)[2])
		description <- c("Methylation value density estimation according to sample grouping and ", txt.site, " category.")
		report <- rnb.add.figure(report, description, figure.plots, setting.names)
	}

	logger.completed()
	return(report)
}

########################################################################################################################

#' rnb.step.intersample
#'
#' Generates one or more deviation plots for the variation step of the methylation profiles module and adds a dedicated
#' section to the given report.
#'
#' @param rnb.set Methylation dataset as an object of type \code{\linkS4class{RnBSet}}.
#' @param report  Report on methylation profiles to contain the methylation deviation section. This must be an object of
#'                type \code{\linkS4class{Report}}.
#' @param columns Optional; predefined column names (in the form of a \code{character} vector) or indices (an
#'                \code{integer} vector) in the sample annotation table. Only these columns are considered for grouping
#'                samples and defining profiles. All other columns in the phenotype table are ignored.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.intersample <- function(rnb.set, report, columns = rnb.getOption("exploratory.columns")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	## rnb.sample.groups validates the value of columns
	sample.inds <- rnb.sample.groups(rnb.set, columns = columns)
	rinfos <- get.site.and.region.types(rnb.set)
	rnb.step.intersample.internal(rnb.set, report, sample.inds, rinfos)
}

rnb.step.intersample.internal <- function(rnb.set, report, sample.inds, rinfos) {
	X <- meth.matrices(rnb.set)
	if (length(X) == 0) {
		return(report)
	}

	## Extract sample group definitions
	logger.start("Scatter Plots of Mean Beta vs Variance")
	if (length(sample.inds) != 0) {
		logger.info(c("Sample subgroups are defined by:", paste(names(sample.inds), collapse = "; ")))
		sample.group.traits <- rep(names(sample.inds), sapply(sample.inds, length))
		sample.group.names <- unlist(lapply(sample.inds, names), use.names = FALSE)
		sample.inds <- unlist(sample.inds, recursive = FALSE, use.names = FALSE)
		names(sample.inds) <- paste(sample.group.names, " (based on ", sample.group.traits, ")", sep = "")
	} else {
		logger.info("None of the traits passes the criteria for sample subgroup definition")
		sample.group.traits <- character()
		sample.group.names <- character()
	}
	sample.inds <- c(list("all samples" = 1:ncol(X[[1]])), sample.inds)
	sample.group.traits <- c("", sample.group.traits)
	sample.group.names <- c(names(sample.inds)[1], sample.group.names)

	## Measure site and/or region variance
	sample.inds.2 <- sample.inds[sapply(sample.inds, length) >= 2]
	sites.supported <- rnb.getOption("analyze.sites")
	regions.supported <- character(0)
	if (sites.supported) {
		if (length(X) == 1) {
			txt.position <- "locus"
			txt.positions <- "locations"
		} else { # length(X) > 1
			txt.position <- "locus/region"
			txt.positions <- "locations/regions"
			regions.supported <- names(X)[-1]
		}
	} else {
		txt.position <- "region"
		txt.positions <- "regions"
		regions.supported <- names(X)
	}
	txt <- c("The variability of the methylation values is measured in two aspects: (1) intra-sample variance, ",
		"that is, differences of methylation between genomic ", txt.positions, " within the same sample, and (2) ",
		"inter-sample variance, i.e. variability in the methylation degree at a specific ", txt.position,
		" across a group of samples.")
	report <- rnb.add.section(report, "Inter-sample Variability", txt, level = 2)

	txt.site <- rnb.get.row.token(rnb.set)
	txt.sites <- rnb.get.row.token(rnb.set, plural = TRUE)
	txt <- c("The following figure shows the relationship between average methylation and methylation variability of ",
		"a ", ifelse(sites.supported, txt.site, "region"), ".")
	if (!identical(sample.inds.2, sample.inds)) {
		txt <- c(txt, " Only groups of 2 or more samples are considered here.")
	}
	rnb.add.paragraph(report, txt)

	scatter.var <- function(fprefix, beta.values, c.values, c.legends) {
		report.plots <- list()
		for (i in 1:length(sample.inds.2)) {
			x.vars <- data.frame(
				"mean" = rowMeans(beta.values[, sample.inds.2[[i]]], na.rm = TRUE),
				"var" = rowVars(beta.values[, sample.inds.2[[i]]], na.rm = TRUE))
			if (!is.null(c.values)) {
				x.vars <- cbind(x.vars, c.values)
			}
			na.rows <- is.na(x.vars[, 1]) | is.na(x.vars[, 2])
			if (any(na.rows)) {
				x.vars <- x.vars[!na.rows, ]
				logger.warning(c("Removed", sum(na.rows), "loci with missing values in group", names(sample.inds.2)[i]))
			}
			for (k in 0:length(c.legends)) {
				fname <- paste(fprefix, i, sep = "_")
				if (length(c.legends) != 0) {
					fname <- paste(fname, k, sep = "_")
				}
				rplot <- createReportPlot(fname, report, create.pdf = FALSE, high.png = 200)
				if (nrow(x.vars) == 0) {
					print(rnb.message.plot("Data not available"))
					pp <- NULL
				} else if (k == 0) {
					pp <- ggplot(x.vars, aes_string(x = "mean", y = "var"))
				} else if (c.legends[k] %in% colnames(c.values)) {
					cl <- paste("`", c.legends[k], "`", sep = "")
					pp <- ggplot(x.vars, aes_string(x = "mean", y = "var", colour = cl))
					cl <- attr(c.values, "legend")[[c.legends[k]]]
					if (is.character(cl)) {
						pp <- pp + scale_colour_manual(values = cl)
					} else { # is.integer(cl) && length(cl) == 1
						if (cl == 2L) {
							g.colors <- rnb.getOption("colors.gradient")
						} else { # cl == 3L
							g.colors <- rnb.getOption("colors.3.gradient")
						}
						pp <- pp + scale_color_gradientn(colours = g.colors, na.value = "grey50")
					}
				} else {
					print(rnb.message.plot("Annotation not available"))
					pp <- NULL
				}
				if (!is.null(pp)) {
					pp <- pp + geom_point() + labs(x = "Mean methylation", y = "Variance")
					pp <- pp + theme(legend.justification = c(0, 0.5), legend.position = c(1, 0.5),
						plot.margin = unit(c(0.5, 1.5, 0.5, 0.5), c("lines", "in", "lines", "lines")))
					print(pp)
				}
				report.plots <- c(report.plots, off(rplot))
			}
		}
		return(report.plots)
	}

	setting.names <- list()
	setting.names.common <- list("Sample group" = names(sample.inds.2))
	names(setting.names.common[[1]]) <- 1:length(sample.inds.2)
	description <- paste("Scatter plot showing the correlation betweeen %s mean methylation and the variance across",
		"a group of samples. Every point corresponds to one %s.")
	if (sites.supported) {
		lg <- names(attr(rinfos[[1]], "legend"))
		fplots <- scatter.var("scatter_meanvariance_sites", X[[1]], rinfos[[1]], lg)
		names(lg) <- 1:length(lg)
		lg <- c("0" = "nothing (all black)", lg)
		setting.names$sites <- setting.names.common
		setting.names$sites[["Point color based on"]] <- lg
		f.description <- sprintf(description, txt.site, txt.site)
		report <- rnb.add.figure(report, f.description, fplots, setting.names$sites)
		logger.status(c("Created", length(fplots), "scatter plots of methylation variance on the site level"))
		rm(lg, fplots, f.description)
	}
	if (length(regions.supported) != 0) {
		ri <- rinfos[regions.supported]
		lg <- unique(unlist(lapply(ri, function(x) { names(attr(x, "legend")) }), use.names = FALSE))
		if (is.null(lg)) { lg <- character(0) }
		fplots <- list()
		for (i in 1:length(regions.supported)) {
			fprefix <- paste("scatter_meanvariance_regions", i, sep = "_")
			fplots <- c(fplots, scatter.var(fprefix, X[[regions.supported[i]]], ri[[i]], lg))
		}
		setting.names$regions <- c(list("Regions" = regions.supported), setting.names.common)
		names(setting.names$regions[[1]]) <- 1:length(setting.names$regions[[1]])
		if (length(lg) != 0) {
			names(lg) <- 1:length(lg)
			lg <- c("0" = "nothing (all black)", lg)
			setting.names$regions[["Point color based on"]] <- lg
		}
		f.description <- sprintf(description, "region", "region")
		if (sites.supported) {
			txt <- c("In a complete analogy to the plots above, the figure below shows the relationship between ",
				"average methylation and methylation variability of a genomic region.")
			rnb.add.paragraph(report, txt)
		}
		report <- rnb.add.figure(report, f.description, fplots, setting.names$regions)
		logger.status(c("Created", length(fplots), "scatter plots of methylation variance on the region level"))
		rm(ri, lg, fplots, i, f.description)
	}
	rm(setting.names.common, description)
	logger.completed()

	## -----------------------------------------------------------------------------------------------------------------
	## Create deviation plots

	o.deviation <- rnb.getOption("exploratory.deviation.plots")
	if (!(isTRUE(o.deviation) || (is.null(o.deviation) && inherits(rnb.set, "RnBeadSet")))) {
		return(report)
	}

	logger.start("Deviation Plots")
	description <- paste0("Deviation plot of a sample group. %s are sorted in increasing order of their median ",
		"methylation%s. The horizontal axis in the plot iterates over %s, and the vertical axis measures methylation ",
		"degree. Median &beta; values are depicted by a blue curve. Grey borders mark the 5th and 95th percentiles of ",
		"&beta; values in a %s%s, ensuring that 90 percent of the observed values lie in the yellow area.<br />",
		"Relative frequencies of %s categories%s are color-coded and plotted below the horizontal ",
		"axis. Every segment in the color legend shows the distribution of %s categories ",
		"that underlie the corresponding segment in the deviation plot above it.")
	txt <- c("The figure below shows a methylation deviation plot for all samples in the dataset",
		ifelse(length(sample.inds) == 1, ".",
			", as well as other sample groups inferred from the table of phenotypic information."))
	rnb.add.paragraph(report, txt)

	setting.names <- lapply(setting.names, function(sn) {
			sn[["Sample group"]] <- names(sample.inds)
			names(sn[["Sample group"]]) <- 1:length(sample.inds)
			i <- which(names(sn) == "Point color based on")
			if (length(i) != 0) {
				names(sn)[i] <- "Color legend based on"
				sn[[i]] <- sn[[i]][-1]
			}
			return(sn)
		})
	sample.group.stats <- data.frame(
		"Loci/regions" = "",
		"Sample Group" = sample.group.names,
		"Based on Trait" = sample.group.traits,
		"Size" = sapply(sample.inds, length),
		"Variability" = as.double(NA), check.names = FALSE, stringsAsFactors = FALSE)
	sample.group.stats <- rep(list(sample.group.stats), length(X))
	report.plots <- list()
	bin.sizes <- rep.int(1L, length(X))
	for (j in 1:length(X)) {
		sample.group.stats[[j]][, "Loci/regions"] <- names(X)[j]
		for (i in 1:length(sample.inds)) {
			## Reorder sites or regions according to their median methylation
			if (length(sample.inds[[i]]) > 1) {
				bstats <- t(rowQuantiles(X[[j]][, sample.inds[[i]]], probs = c(0.05, 0.5, 0.95), na.rm = TRUE))
			} else {
				bstats <- matrix(rep(X[[j]][, sample.inds[[i]]], each = 3), nrow = 3)
			}
			site.order <- order(bstats[2, ], bstats[1, ], bstats[3, ], na.last = NA)
			bstats <- bstats[, site.order]
			sample.group.stats[[j]][i, "Variability"] <- mean(bstats[3, ] - bstats[1, ])
			cuts <- deviation.plot.beta.get.cuts(ncol(bstats))
			if (!is.null(cuts)) {
				bin.sizes[j] <- max(bin.sizes[j], max(table(cuts)))
			}
			if (j == 1 && sites.supported) {
				for (k in 1:ncol(rinfos[[1]])) {
					c.values <- rinfos[[1]][site.order, k]
					c.legend <- attr(rinfos[[1]], "legend")[[k]]
					fname <- paste("deviation_sites", i, k, sep = "_")
					rplot <- createReportPlot(fname, report, width = 8, height = 5)
					deviation.plot.beta.internal(bstats, c.values, c.legend, cuts)
					report.plots <- c(report.plots, off(rplot))
				}
				next
			}
			# length(regions.supported) != 0
			if ("Color legend based on" %in% names(setting.names$regions)) {
				knames <- setting.names$regions[["Color legend based on"]]
				kvalues <- 1:length(knames)
				ri <- rinfos[[names(X)[j]]]
			} else {
				knames <- character()
				kvalues <- 0L
				ri <- NULL
			}
			for (k in kvalues) {
				if ((!is.null(ri)) && (knames[k] %in% colnames(ri))) {
					c.values <- ri[site.order, knames[k]]
					c.legend <- attr(ri, "legend")[[knames[k]]]
				} else {
					c.values <- rep.int("1", ncol(bstats))
					c.legend <- c("1" = "white")
				}
				fname <- paste("deviation_regions", j - sites.supported, i, sep = "_")
				if (!identical(kvalues, 0L)) {
					fname <- paste(fname, k, sep = "_")
				}
				rplot <- createReportPlot(fname, report, width = 8, height = 5)
				deviation.plot.beta.internal(bstats, c.values, c.legend, cuts)
				report.plots <- c(report.plots, off(rplot))
			}
			rm(knames, kvalues, ri)
		}
		if (j == 1 && sites.supported) {
			is.binned <- (bin.sizes[1] != 1)
			f.description <- sprintf(description, capitalize(txt.sites),
				ifelse(is.binned, paste(" and are binned in groups of up to", bin.sizes[1]), ""),
				ifelse(is.binned, paste(txt.site, "groups"), txt.sites), txt.site,
				ifelse(is.binned, " (averaged over the group)", ""), txt.site,
				ifelse(is.binned, " in every group", ""), txt.site)
			report <- rnb.add.figure(report, f.description, report.plots, setting.names$sites)
			logger.status(c("Created", length(report.plots), "deviation plots on the site level"))
			report.plots <- list()
		}
	}
	if (length(report.plots) != 0) {
		if (sites.supported) {
			txt <- "In a similar fashion, the figure below shows deviation plots on the region level."
			rnb.add.paragraph(report, txt)
			bin.sizes <- bin.sizes[-1]
		}
		is.binned <- any(bin.sizes != 1)
		txt <- paste(ifelse(all(bin.sizes != 1), " and", ", and in some cases"), " are binned in groups of up to ",
			max(bin.sizes), sep = "")
		f.description <- sprintf(description, "Regions",
			ifelse(is.binned, txt, ""),
			ifelse(is.binned, "region groups", "regions"), "region",
			ifelse(is.binned, " (averaged over the group)", ""), "region",
			ifelse(is.binned, " in every group", ""), "region")
		report <- rnb.add.figure(report, f.description, report.plots, setting.names$regions)
		logger.status(c("Created", length(report.plots), "deviation plots on the region level"))
	}
	rm(report.plots, j, k, c.values, c.legend, fname, rplot, is.binned, f.description)

	sample.group.stats <- do.call(rbind, sample.group.stats)
	variabilities <- sprintf("%1.4f", sample.group.stats[["Variability"]])
	variabilities[is.na(sample.group.stats[["Variability"]])] <- as.character(NA)
	sample.group.stats[["Variability"]] <- variabilities

	## Add group statistics to the report
	txt <- c("The <i>variability</i> of a sample group is the span between 5th and 95th percentile of &beta; values , ",
		"averaged over all valid ", txt.positions, ". This amounts to a number between 0 and 1 and corresponds to the ",
		"relative area of deviation in the plots presented above. The table below lists the variabilities of the ",
		"studied sample groups.")
	rnb.add.paragraph(report, txt)
	rnb.add.table(report, sample.group.stats, row.names = FALSE)
	logger.status("Added summary table of variabilities")

	logger.completed()
	return(report)
}
