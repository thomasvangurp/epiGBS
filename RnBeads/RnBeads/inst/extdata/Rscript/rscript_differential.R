library(argparse)
suppressPackageStartupMessages(library(RnBeads))

ap <- ArgumentParser()
ap$add_argument("-x", "--xml", action="store", help="Configuration xml file")
ap$add_argument("-s", "--rnbSet", action="store", dest="rnbSetLoc", help="Location of stored RnBSet")
ap$add_argument("-o", "--output", action="store", help="Output directory")
ap$add_argument("-c", "--cores", action="store", type="integer", default=1, help="Number of cores used for the analysis")
cmdArgs <- ap$parse_args()
module.name <- "differential"

logger.start(fname=NA)
logger.status(c("...Started module:",module.name))

logger.start("Configuring Analysis")
	rnb.settings <- rnb.xml2options(cmdArgs$xml,return.full.structure=TRUE)


	data.source <- rnb.settings$analysis.params[["data.source"]]
	data.type <- rnb.settings$analysis.params[["data.type"]]
	report.dir <- rnb.settings$analysis.params[["dir.reports"]]
	analysis.options <- rnb.settings$options

	if ("preanalysis.script" %in% names(rnb.settings)){
		source(rnb.settings$preanalysis.script)
	} 
	## Set options
	if (length(analysis.options) != 0) {
		do.call(rnb.options, analysis.options)
	}

	if(rnb.getOption("region.subsegments")>1L){
		logger.start("Loading region subsegmentation data")
			subseg.annot.dir <- file.path(report.dir,"preprocessing_data","subsegments")
			load.region.subsegment.annotation(rnb.set,subseg.annot.dir)
		logger.completed()
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

result.diffmeth <- rnb.run.differential(rnb.set, report.dir, close.report = TRUE)

logger.start("Saving")
	diffmeth.path <- file.path(cmdArgs$output,paste0(module.name,"_rnbDiffMeth"))
	save.rnb.diffmeth(result.diffmeth$diffmeth, diffmeth.path)
	diffmeth.enrichment <- result.diffmeth$dm.enrich
	if (!is.null(diffmeth.enrichment)){
		save(diffmeth.enrichment, file=file.path(diffmeth.path, "enrichment.RData"))
	}
logger.completed()

logger.status(c("...Completed module:",module.name))
quit(save='no')
