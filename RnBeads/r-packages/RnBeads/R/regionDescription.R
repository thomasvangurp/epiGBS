## rnb.section.region.description
##
## Adds a section summarizing the requested regions.
##
## @param report  Report on loading to contain the newly constructed section.
## @param rnb.set Methylation dataset to be analyzed.
## @param regions Non-empty \code{character} vector of region names.
## @return The modified report.
##
## @author Yassen Assenov, Fabian Mueller
rnb.section.region.description <- function(report, rnb.set, regions) {

	#Descriptive table
	msg <- function(txt, e.class = "disabled") {
		paste("<span class=\"", e.class, "\">", txt, "</span>", sep = "")
	}
	reg.counts <- sapply(regions, function(region) {
			tryCatch(as.character(nrow(meth(rnb.set, region))), error = function(e) { msg("not supported") })
		}
	)
	get.r.description <- function(reg) {
		result <- attr(rnb.get.annotation(reg, assembly = assembly(rnb.set)), "description")
		if (is.null(result) || identical("", result)) {
			result <- msg("n.a.")
		} else {
			result <- paste("<p style=\"font-weight:normal;text-align:left\">", result, "</p>", sep = "")
		}
		result
	}
	reg.descriptions <- sapply(regions, function(region) {
			tryCatch(get.r.description(region), error = function(e) { msg("not imported", "outdated") })
		}
	)
	table.statistics <- data.frame(
		"Annotation" = regions, "Description" = reg.descriptions, "Regions in the Dataset" = reg.counts,
		check.names = FALSE, stringsAsFactors = FALSE)
	txt <- c(ifelse(rnb.getOption("analyze.sites"), "In addition to CpG sites, there", "There"),
		ifelse(length(regions) == 1, " is one set", paste(" are", length(regions), "sets")),
		" of genomic regions to be covered in the analysis. The table below gives a summary of these annotations.")
	report <- rnb.add.section(report, "Region Annotations", txt)
	table.header <- c("<colgroup>", paste("\t<col width=\"", c(210, 420, 150), "px\" />", sep = ""), "</colgroup>")
	rnb.add.table(report, table.statistics, row.names = FALSE, thead = table.header)

	names(regions) <- ifelse(is.valid.fname(regions),regions,paste("reg",1:length(regions),sep=""))
	plots.sizes <- list()
	plots.n.sites <- list()
	plots.site.distr <- list()
	for (i in 1:length(regions)){
		reg <- regions[i]
		reg.file.name <- names(regions)[i]
		aa <- annotation(rnb.set,type=reg)
		regions2sites <- regionMapping(rnb.set,reg)
		reg.sizes <- as.integer(aa$End - aa$Start + 1L) #+1 as we use 1-based, right-inclusive regions
		reg.n.sites <- vapply(regions2sites,length,0L)
		df2p <- data.frame(size=reg.sizes,n.sites=reg.n.sites)

		if (max(reg.sizes)==min(reg.sizes)){
			pp <- rnb.message.plot(paste0("Constant region size: ",min(reg.sizes)))
		} else {
			pp <- ggplot(df2p) + aes(x=size) + geom_density(color=rnb.getOption('colors.category')[1],
				fill=rnb.getOption('colors.category')[1], alpha=0.5) + 
				scale_x_continuous(limits = quantile(df2p$size, c(0.05, 0.95)))
		}
		figName <- paste("regionSizes",reg.file.name,sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=TRUE)
		report.plot <- off(report.plot,handle.errors=TRUE)
		plots.sizes <- c(plots.sizes,list(report.plot))

		pp <- ggplot(df2p) + aes(x=n.sites) + geom_density(color=rnb.getOption('colors.category')[1],
			fill=rnb.getOption('colors.category')[1], alpha=0.5)
		figName <- paste("regionNsites",reg.file.name,sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=TRUE)
		report.plot <- off(report.plot,handle.errors=TRUE)
		plots.n.sites <- c(plots.n.sites,list(report.plot))

		pp <- rnb.plot.region.site.density(rnb.set,reg)
		figName <- paste("regionSiteDistr",reg.file.name,sep="_")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=TRUE)
		report.plot <- off(report.plot,handle.errors=TRUE)
		plots.site.distr <- c(plots.site.distr,list(report.plot))
	}

	setting.names <- list('Region type' = regions)

	#length distribution and number of sites per region
	txt <- c("The plots below show region size distributions for the region types above.")
	report <- rnb.add.section(report, "Region length distributions", txt, level = 2)
	
	description <- "Distribution of region lengths"
	report <- rnb.add.figure(report, description, plots.sizes, setting.names)

	txt <- c("The plots below show the distributions of the number of sites per region type.")
	report <- rnb.add.section(report, "Number of sites per region", txt, level = 2)

	description <- "Distribution of the number of sites per region"
	report <- rnb.add.figure(report, description, plots.n.sites, setting.names)

	txt <- c("The plots below show distributions of sites across the different region types.")
	report <- rnb.add.section(report, "Region site distributions", txt, level = 2)
	
	description <- "Distribution of sites across regions. relative coordinates of 0 and 1 corresponds to the start
					and end coordinates of that region respectively. Coordinates smaller than 0 and greater than 1 denote flanking regions
					normalized by region length."
	report <- rnb.add.figure(report, description, plots.site.distr, setting.names)
	
	return(report)
}
