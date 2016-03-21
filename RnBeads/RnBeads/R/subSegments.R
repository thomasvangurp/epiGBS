#methods for subdividing regions into smaller regions of clustered CpG sites

########################################################################################################################
### get.subsegments
###
### For a given region. Subdivide this region into subsegments by
### hierarchical clustering on the site distances in a particular region and then splitting the into subregions consisting
### of these site clusters. The number of clusters is determined in such way that the mean number of sites per cluster
### is given by the \code{ns} parameter.
### @param reg.annot a \code{data.frame} of exactly 1 row as output by \code{annotation(rnb.set,type=...)} specifying the region to be subdivided
### @param site.coords the genomic start coordinates of the sites to be considered for subdivision
### @param ns the mean number of sites per cluster.
### @param site.length length of the site motif. Typically 2 for CpGs.
### @return a data.frame of the same format as \code{reg.annot} with each row corresponding
###				to a subsegment of the original region spanning first site to last site. All columns except
###				for genomic coordinates will be identical to \code{reg.annot}
### @author Fabian Mueller
get.subsegments <- function(reg.annot,site.coords,ns=10,site.length=2){
	res <- reg.annot
	n.sites <- length(site.coords)
	if (n.sites > ns){
		#do the clustering
		site.dists <- dist(site.coords,method="manhattan") #method doesn't matter
		site.clust <- hclust(site.dists,method="single")
		n.clust <- ceiling(n.sites/ns)
		clusts <- cutree(site.clust,k=n.clust)
		#TODO: see if the dynamicTreeCut package can improve the method in resonable time

		#repeat the region annotation for each cluster and just adjust the coordinates
		res <- reg.annot[rep(1,n.clust),]
		clust.coords.min <- tapply(site.coords,clusts,min)
		clust.coords.max <- tapply(site.coords,clusts,max)+(site.length-1)
		res$Start <- clust.coords.min
		res$End <- clust.coords.max
		rownames(res) <- paste0(rownames(reg.annot),".sub",1:n.clust)
	}
	return(res)
}

########################################################################################################################
### add.subsegment.annotation
###
### For a given region annotation data.frame. Subdivide each region into subsegments by
### hierarchical clustering on the site distances in a particular region and then splitting the into subregions consisting
### of these site clusters. The number of clusters is determined in such way that the mean number of sites per cluster
### is given by the \code{ns} parameter.
### @param reg.name the name/region.type of the newly generated subregions
### @param regions a \code{data.frame} as output by \code{annotation(rnb.set,type=...)} specifying the regions to be subdivided
### @param sites a \code{data.frame} as output by \code{annotation(rnb.set,type="sites")} specifying all sites to be considered
###				for subdivision.
### @param add.region.types.to.options Flag indicating whether to add the newly created subregions to the package's
### 				  \code{region.types} option
### @param regions2sites a mapping of regions to sites. A list as output by \code{regionMapping()}
### @param ns the mean number of sites per cluster.
### @param save.dir a directory to save the annotation to for later reloading. (binary \code{RData} format.)
### @param assembly the genome assembly to be used for the annotation
### @return nothing of particular interest
### @author Fabian Mueller
add.subsegment.annotation <- function(reg.name,regions,sites,regions2sites,ns=10,save.dir=NULL,assembly="hg19"){
	if (nrow(regions) != length(regions2sites)){
		stop("incompatible regions and regions2sites")
	}
	n.regions <- nrow(regions)

	logger.start("Subsegmentation")
	#parallel computing here taks up too much memory and serial computing is not great, but ok timewise
	# if (parallel.isEnabled()) {
	# 	subregions <- foreach(i=1:n.regions, .combine='rbind',.multicombine=TRUE,.maxcombine=100) %dopar% {
	# 		get.subsegments(regions[i,],sites[regions2sites[[i]],]$Start,ns=ns)
	# 	}
	# } else {
		subregions <- do.call(rbind,lapply(1:n.regions,FUN=function(i){
			get.subsegments(regions[i,],sites[regions2sites[[i]],]$Start,ns=ns)
		}))
	# }
	logger.completed()
	#renaming and reformating such that the annotation can be set
	names(subregions$Start) <- NULL
	names(subregions$End) <- NULL
	subregions$CpG <- NULL
	subregions$GC <- NULL
	subregions$Start <- as.integer(subregions$Start)
	subregions$End <- as.integer(subregions$End)
	# #TEMPRORARY for testing
	# fn <- file.path(save.dir,paste0("regions_",reg.name,"_tmp.RData"))
	# save(list=ls(),file=fn)
	rnb.cleanMem()
	logger.start("Adding CpG statistics")
	rnb.set.annotation.and.cpg.stats(reg.name,subregions,description=NULL,assembly=assembly)
	logger.completed()
	if (!is.null(save.dir)) {
		logger.start("Saving annotation")
		fn <- file.path(save.dir,paste0("regions_",reg.name,".RData"))
		rnb.save.annotation(fn,reg.name,assembly=assembly)
		logger.completed()
	}
}

########################################################################################################################
#' addRegionSubsegments
#'
#' For the region annotation of a given \code{RnBSet} object. Subdivide each region into subsegments by
#' hierarchical clustering on the site distances in a particular region and then splitting the region into subregions consisting
#' of these site clusters. The number of clusters is determined in such way that the mean number of sites per cluster
#' is given by the \code{ns} parameter.
#' @param rnb.set an \code{RnBSet object}
#' @param annotation.dir a directory to save the annotation to for later reloading. (binary \code{RData} format.)
#' @param region.types the region types to which subsegmentation should be applied. Must be a non-empty
#' 					   subset of \code{summarized.regions(rnb.set)}. Defaults (\code{NULL}) to all region types in
#' 					   \code{rnb.set}
#' @param add.region.types.to.options Flag indicating whether to add the newly created subregions to the package's
#' 				  \code{region.types} option
#' @param ns the mean number of sites per cluster.
#' @return the modified \code{RnBSet} object
#' @author Fabian Mueller
#' @export addRegionSubsegments
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.set.mod <- addRegionSubsegments(rnb.set.example,tempdir(),region.types=c("tiling","genes"))
#' summary(meth(rnb.set.mod,type="tiling.subsegments"))
#' }
addRegionSubsegments <- function(rnb.set,annotation.dir,region.types=NULL,add.region.types.to.options=FALSE,ns=10){
	logger.start("Adding region sub-segmentation to RnBSet")
	if (is.null(region.types)){
		region.types <- summarized.regions(rnb.set)
	} else {
		region.types <- intersect(region.types,summarized.regions(rnb.set))
	}
	if (length(region.types)<1){
		stop("Invalid region types specified for subsegmentation")
	}
	assembly <- assembly(rnb.set)
	#check if all regions are annotated
	if (!all(region.types %in% rnb.region.types(assembly=assembly))){
		stop(paste("Region annotation not found."))
	}

	aa.s <- annotation(rnb.set,"sites")
	for (rt in region.types){
		logger.status(c("Subsegmenting region:",rt))
		aa.r <- annotation(rnb.set,rt)
		regions2sites <- regionMapping(rnb.set,rt)
		subreg.name <- paste(rt,"subsegments",sep=".")
		logger.start("Annotating")
		add.subsegment.annotation(subreg.name,aa.r,aa.s,regions2sites,ns=ns,save.dir=annotation.dir,assembly=assembly)
		logger.completed()
		logger.start("Summarizing regions in RnBSet")
		rnb.set <- summarize.regions(rnb.set,region.type=subreg.name)
		logger.completed()
		if (add.region.types.to.options){
			rnb.options(region.types=c(rnb.getOption("region.types"),subreg.name))
		}
	}
	logger.completed()
	return(rnb.set)
}

########################################################################################################################
#' load.region.subsegment.annotation
#'
#' For the region annotation of a given \code{RnBSet} object. Subdivide each region into subsegments by
#' hierarchical clustering on the site distances in a particular region and then splitting the region into subregions consisting
#' of these site clusters. The number of clusters is determined in such way that the mean number of sites per cluster
#' is given by the \code{ns} parameter.
#' @param rnb.set The \code{RnBSet object} with subsegments specified in the regions
#' @param annotation.dir a directory to load the annotation from. (binary \code{RData} format.)
#' @return invisible \code{TRUE}
#' @author Fabian Mueller
#' @export
load.region.subsegment.annotation <- function(rnb.set,annotation.dir){
	logger.start("Loading region sub-segmentation annotation")
	region.types <- grep("\\.subsegments$",summarized.regions(rnb.set),value=TRUE)

	for (rt in region.types){
		logger.status(c("Loading annotation",rt))
		fn <- file.path(annotation.dir,paste0("regions_",rt,".RData"))
		rnb.load.annotation(fn, rt)
	}
	logger.completed()
	invisible(TRUE)
}

########################################################################################################################
#' rnb.section.region.subsegmentation
#'
#' Given a \code{RnBSet} with subsegmentation information, add a correspong section to the provided RnBeads report.
#' @param rnb.set The \code{RnBSet object} with subsegments specified in the regions
#' @param report Report to summarize the outcome of the region subsegmentation. This must be an object of type
#'               \code{\linkS4class{Report}}.
#' @return The modified \code{report}, invisibly.
#' @author Fabian Mueller
#' @noRd
rnb.section.region.subsegmentation <- function(report,rnb.set){
	logger.start("Adding region sub-segmentation section")
	txt <- c("Region subsegmentation was performed in order to obtained a more fine-grained view of local methylation ",
			 "levels. Each region was subdivided into subregions containing ",rnb.getOption("region.subsegments")," sites ",
			 "on average. These subregions were obtained by hierarchical clustering using ",
			 "single linkage agglomeration based on genomic coordinates. Each cluster then ",
			 "corresponds to a subsegment. This guarantees that neighboring ",
			 "sites end up in the same subregion. The following table provides some key statistics on the ",
			 "subsegmentation.")
	rnb.add.section(report, "Region sub-segmentation", txt)

	reg.types.subseg <- grep("\\.subsegments$",summarized.regions(rnb.set),value=TRUE)
	reg.types.base    <- sub("\\.subsegments$","",reg.types.subseg)
	reg.types.all <- c(reg.types.base,reg.types.subseg)

	reg.types.base.files <- paste0("reg",1:length(reg.types.base))
	names(reg.types.base) <- reg.types.base.files

	get.summary.df.from.list <- function(ll,value.name="value"){
		res <- melt(ll)
		colnames(res) <- c(value.name,"region")
		res$is.subsegmentation <- grepl("\\.subsegments$",res$region)
		res$region <- sub("\\.subsegments$","",res$region)
		res$region <- factor(res$region,levels=reg.types.base)
		return(res)
	}

	list.allreg.size <- lapply(reg.types.all,FUN=function(rn){
		aa <- annotation(rnb.set,type=rn)
		return(aa[,"End"]-aa[,"Start"])
	})
	names(list.allreg.size) <- reg.types.all
	df.allreg.size <- get.summary.df.from.list(list.allreg.size,"region.size")

	list.allreg.num.sites <- lapply(reg.types.all,FUN=function(rn){
		reg.map <- regionMapping(rnb.set,rn)
		return(vapply(reg.map,length,integer(1)))
	})
	names(list.allreg.num.sites) <- reg.types.all
	df.allreg.num.sites <- get.summary.df.from.list(list.allreg.num.sites,"num.sites")
	#add a table containing key statistics on subsegmentation
	num.regs  <- vapply(list.allreg.size,length,integer(1))
	mean.size <- vapply(list.allreg.size,mean,numeric(1),na.rm=TRUE)
	median.size <- vapply(list.allreg.size,median,numeric(1),na.rm=TRUE)
	mean.num.sites <- vapply(list.allreg.num.sites,mean,numeric(1),na.rm=TRUE)
	median.num.sites <- vapply(list.allreg.num.sites,median,numeric(1),na.rm=TRUE)

	tt <- data.frame(num.regs[reg.types.base],
					 num.regs[reg.types.subseg],
					 mean.num.sites[reg.types.base],
					 mean.num.sites[reg.types.subseg],
					 median.num.sites[reg.types.base],
					 median.num.sites[reg.types.subseg],
					 mean.size[reg.types.base],
					 mean.size[reg.types.subseg],
					 median.size[reg.types.base],
					 median.size[reg.types.subseg])
	colnames(tt) <-c("Number of Regions (base)",
					 "Number of Regions (subseg)",
					 "Mean number of sites (base)",
					 "Mean number of sites (subseg)",
					 "Median number of sites (base)",
					 "Median number of sites (subseg)",
					 "Mean region size (base)",
					 "Mean region size (subseg)",
					 "Median region size (base)",
					 "Median region size (subseg)")
	tt <- round(tt,4)
	
	rnb.add.table(report,tt)

	#create some nice plots
	cvalues <- rep(rnb.getOption("colors.category"), length.out = 2)
	figPlots.size <- list()
	figPlots.num.sites <- list()
	for (i in 1:length(reg.types.base)){
		rn  <- reg.types.base[i]
		rnf <- reg.types.base.files[i]
		#plot distribution of region sizes
		figName <- paste("subseg_distribution","size",rnf,sep="_")
		df2p <- df.allreg.size[df.allreg.size$region==rn,]
		pp <- ggplot(df2p, aes(x = log2(region.size), color = is.subsegmentation)) +
			labs(y = "Density", color = "Region subsegmentation") +
			geom_density() + scale_color_manual(values = cvalues) +
			theme(legend.position="bottom")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=TRUE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots.size <- c(figPlots.size,list(report.plot))

		#plot distribution of number of sites per region
		figName <- paste("subseg_distribution","numSites",rnf,sep="_")
		df2p <- df.allreg.num.sites[df.allreg.num.sites$region==rn,]
		pp <- ggplot(df2p, aes(x = log2(num.sites), color = is.subsegmentation)) +
			labs(y = "Density", color = "Region subsegmentation") +
			geom_density() + scale_color_manual(values = cvalues) +
			theme(legend.position="bottom")
		report.plot <- createReportGgPlot(pp,figName, report,create.pdf=TRUE,high.png=200)
		report.plot <- off(report.plot,handle.errors=TRUE)
		figPlots.num.sites <- c(figPlots.num.sites,list(report.plot))
	}
	setting.names <- list(
			'region type' = reg.types.base)
	description <- c('Distribution of the log2 of the region sizes of the base regions ',
					 'and of the regions resulting from sub-segmentation.')
	report <- rnb.add.figure(report, description, figPlots.size, setting.names)
	
	description <- c('Distribution of the log2 of the number of sites in the base regions ',
					 'and in the regions resulting from sub-segmentation.')
	report <- rnb.add.figure(report, description, figPlots.num.sites, setting.names)

	logger.completed()
	invisible(report)
}

########################################################################################################################

#' rnb.step.region.subsegmentation
#'
#' Perform region subsegmentation (see \code{\link{addRegionSubsegments}}) for details and add a corresponding
#' section to the provided report
#' @param rnb.set an \code{RnBSet object}
#' @param report Report to summarize the outcome of this procedure. This must be an object of type
#'               \code{\linkS4class{Report}}. The subsegments will be stored in the 'data/subsegmentation' directory
#' 				 of the report
#' @param region.types the region types to which subsegmentation should be applied. Must be a non-empty
#' 					   subset of \code{summarized.regions(rnb.set)}. Defaults (\code{NULL}) to all region types in
#' 					   \code{rnb.set}
#' @param add.region.types.to.options,ns arguments passed on to \code{\link{addRegionSubsegments}}
#' @return a list with two elements: \code{rnb.set} the modified \code{RnBSet} object with subregions included
#'         and \code{report} the modified report with a section concerning region subsegmentation added
#' @author Fabian Mueller
#' @noRd
rnb.step.region.subsegmentation <- function(rnb.set,report,region.types=NULL,
		add.region.types.to.options=TRUE,ns=rnb.getOption("region.subsegments")){
	logger.start("Region sub-segmentation")
	annot.dir <- file.path(rnb.get.directory(report, dir="data",absolute=TRUE),"subsegmentation")
	create.path(annot.dir, accept.existing=FALSE)
	rnb.set <- addRegionSubsegments(rnb.set,annot.dir,region.types=region.types,add.region.types.to.options=add.region.types.to.options,ns=ns)
	report <- rnb.section.region.subsegmentation(report,rnb.set)
	res <- list(rnb.set=rnb.set,report=report)
	logger.completed()
	invisible(res)
}
