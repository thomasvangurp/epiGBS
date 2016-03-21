########################################################################################################################
## test_logger.R
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Simple tests on the logging functionality.
########################################################################################################################

test_logger <- function() {
	rnb.options(logging.memory = TRUE)
	logger.close()
	checkTrue(logger.isinitialized() == FALSE)

	## Print a dummy log in a file
	logfiles <- c("testme.log", NA)
	logger.start("Analysis", fname = logfiles)
	checkTrue(logger.isinitialized())
	logger.info(c("Loaded information from", "data.RData"))
	logger.start("Processing Detection P-values")
	logger.info(c("Removed", 3979, "probes that overlap with SNPs"))
	logger.info(c("Completed Greedycut on", 510, "samples"))
	logger.completed()
	logger.warning(c("File not found:", "data2.RData"))
	logger.completed()
	logger.close()
	checkTrue(logger.isinitialized() == FALSE)

	## Delete the generated file
	logfiles <- logfiles[!is.na(logfiles)]
	if (length(logfiles) != 0) {
		file.remove(logfiles)
	}
}
