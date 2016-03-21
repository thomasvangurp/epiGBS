library(argparse)
suppressPackageStartupMessages(library(RnBeads))

ap <- ArgumentParser()
ap$add_argument("-x", "--xml", action="store", help="Configuration xml file")
ap$add_argument("-s", "--rnbSet", action="store", dest="rnbSetLoc", help="Location of stored RnBSet")
ap$add_argument("-d", "--diffMethChunkId", action="append", dest="diffMethChunkIds", help="Ids for chunks for which differential methylation has already been computed. Will be used to construct the input file locations.")
ap$add_argument("-o", "--output", action="store", help="Output directory")
ap$add_argument("-c", "--cores", action="store", type="integer", default=1, help="Number of cores used for the analysis")
cmdArgs <- ap$parse_args()
module.name <- "differential_wrapup"

logger.start(fname=NA)
logger.status(c("...Started module:",module.name))

logger.start("Configuring Analysis")
	rnb.settings <- rnb.xml2options(cmdArgs$xml,return.full.structure=TRUE)

	report.dir <- rnb.settings$analysis.params[["dir.reports"]]
	analysis.options <- rnb.settings$options
	dm.chunk.ids <- cmdArgs$diffMethChunkIds

	if ("preanalysis.script" %in% names(rnb.settings)){
		source(rnb.settings$preanalysis.script)
	} 
	## Set options
	if (length(analysis.options) != 0) {
		do.call(rnb.options, analysis.options)
	}

	logger.machine.name()

	if (cmdArgs$cores > 1) {
		parallel.setup(cmdArgs$cores)
	}

	aname <- rnb.getOption("analysis.name")
	if (!(is.null(aname) || is.na(aname) || nchar(aname) == 0)) {
		logger.info(c("Analysis Title:", aname))
	}
	ncores <- parallel.getNumWorkers()
	if (ncores == -1) {
		ncores <- 1L
	}
	logger.info(c("Number of cores:", ncores))
	rm(aname, ncores)
logger.completed()

logger.start("Loading RnBSet")
	rnb.set <- load.rnb.set(cmdArgs$rnbSetLoc)
logger.completed()

logger.start(fname=c(file.path(report.dir,paste0("analysis_",module.name,".log")),NA))

################################################################################
# main script
################################################################################
if (length(dm.chunk.ids)<1){
	stop("Specify at least one id for a differential methylation chunk")
}

disk.dump <- rnb.getOption("disk.dump.big.matrices")

logger.start("Loading and Combining DiffMeth Files")
	chunk.id <- dm.chunk.ids[1]
	diffmeth.path <- file.path(cmdArgs$output,paste0("differential_chunk_",chunk.id,"_rnbDiffMeth"))
	enrichment.file <- file.path(diffmeth.path, "enrichment.RData")
	logger.info(c("Loading:",diffmeth.path))
	diffmeth <- load.rnb.diffmeth(diffmeth.path)
	diffmeth.enrichment <- NULL
	if (file.exists(enrichment.file)){
		logger.info(c("Loading:",enrichment.file))
		load(enrichment.file) #loads 'diffmeth.enrichment' object
	}
	diffmeth.res <- list(report=NULL,diffmeth=diffmeth,dm.enrich=diffmeth.enrichment)

	rm(diffmeth,diffmeth.enrichment) # to make sure a result.diffmeth object is not stored twice if a following file does not contain one
	RnBeads:::rnb.cleanMem()

	if (length(dm.chunk.ids)>1){
		for (chunk.id in dm.chunk.ids[2:length(dm.chunk.ids)]){
			diffmeth.path <- file.path(cmdArgs$output,paste0("differential_chunk_",chunk.id,"_rnbDiffMeth"))
			enrichment.file <- file.path(diffmeth.path, "enrichment.RData")
			logger.info(c("Loading:",diffmeth.path))
			diffmeth <- load.rnb.diffmeth(diffmeth.path)
			diffmeth.enrichment <- NULL
			if (file.exists(enrichment.file)){
				logger.info(c("Loading:",enrichment.file))
				load(enrichment.file) #loads 'diffmeth.enrichment' object
			}
			diffmeth.res.cur <- list(report=NULL,diffmeth=diffmeth,dm.enrich=diffmeth.enrichment)

			ll <- list(diffmeth.res,diffmeth.res.cur)
			diffmeth.res <- RnBeads:::combine.diffMeth.objs(ll)

			destroy(diffmeth)
			rm(diffmeth.res.cur,diffmeth,diffmeth.enrichment) # to make sure a result.diffmeth object is not stored twice if a following file does not contain one
			RnBeads:::rnb.cleanMem()
		}
	}
	if (!is.valid(diffmeth.res$diffmeth,verbose=TRUE)){
		stop("RnBDiffMeth object is invalid")
	}
logger.completed()

logger.start("Saving combined data")
	diffmeth.path <- file.path(cmdArgs$output,paste0("differential_rnbDiffMeth"))
	save.rnb.diffmeth(diffmeth.res$diffmeth, diffmeth.path)
	diffmeth.enrichment <- diffmeth.res$dm.enrich
	if (!is.null(diffmeth.enrichment)){
		save(diffmeth.enrichment, file=file.path(diffmeth.path, "enrichment.RData"))
	}
logger.completed()

logger.start("Differential Methylation")
	logger.start("Report Generation")
		aname <- rnb.getOption("analysis.name")
		if (!(is.null(aname) || is.na(aname) || nchar(aname) == 0)) {
			logger.info(c("Analysis Title:", aname))
		}
		ncores <- parallel.getNumWorkers()
		if (ncores == -1) {
			ncores <- 1L
		}
		logger.info(c("Number of cores:", ncores))
		rm(aname, ncores)

		diffmeth <- diffmeth.res$diffmeth
		dm.enrich <- diffmeth.res$dm.enrich

		init.configuration <- !file.exists(file.path(report.dir, "configuration"))
		report <- RnBeads:::init.pipeline.report("differential_methylation", report.dir, init.configuration)
		optionlist <- rnb.options("analyze.sites","region.types", "differential.permutations", "differential.comparison.columns",
			"differential.comparison.columns.all.pairwise","columns.pairing","differential.site.test.method","covariate.adjustment.columns",
			"differential.adjustment.sva","differential.adjustment.celltype","differential.enrichment")
		report <- RnBeads:::rnb.add.optionlist(report, optionlist)
		
		if (is.null(diffmeth)){
			txt <- "Differential methylation analyis was skipped because no valid grouping information could be found."
			report <- rnb.add.section(report, "Differential Methylation Analysis", txt)
		} else {
			gz <- rnb.getOption("gz.large.files")
			includeSites <- rnb.getOption("analyze.sites")
			report <- RnBeads:::rnb.section.diffMeth.introduction(diffmeth,report)
			if (includeSites){
				report <- RnBeads:::rnb.section.diffMeth.site(rnb.set,diffmeth,report,gzTable=gz)
			}
			if (length(get.region.types(diffmeth))>0){
				report <- RnBeads:::rnb.section.diffMeth.region(rnb.set,diffmeth,report,dm.enrich=dm.enrich,gzTable=gz)
			}
		}

		off(report)

	logger.completed()
logger.completed()

logger.status(c("...Completed module:",module.name))
quit(save='no')
