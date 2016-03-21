library(argparse)
suppressPackageStartupMessages(library(RnBeads))

ap <- ArgumentParser()
ap$add_argument("-x", "--xml", action="store", help="Configuration xml file")
ap$add_argument("-s", "--rnbSet", action="store", dest="rnbSetLoc", help="Location of stored RnBSet")
ap$add_argument("-c", "--cores", action="store", type="integer", default=1, help="Number of cores used for the analysis")
cmdArgs <- ap$parse_args()
module.name <- "qc"

logger.start(fname=NA)
logger.status(c("...Started module:",module.name))

logger.start("Configuring Analysis")
	rnb.settings <- rnb.xml2options(cmdArgs$xml,return.full.structure=TRUE)

	report.dir <- rnb.settings$analysis.params[["dir.reports"]]
	analysis.options <- rnb.settings$options

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

rnb.run.qc(rnb.set, report.dir, close.report = TRUE)

logger.status(c("...Completed module:",module.name))
quit(save='no')
