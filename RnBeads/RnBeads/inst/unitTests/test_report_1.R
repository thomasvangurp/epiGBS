########################################################################################################################
## test_report_1.R
## created: 2012-04-16
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Tests on the Report and ReportPlot class definitions, in the form of an example report generation.
########################################################################################################################

test_report <- function() {
	
	## Initialize the reports locations and the log
	report.directory <- "RnBeads_report_testing"
	if (!rnb.initialize.reports(report.directory)) {
		stop(paste("Could not initialize reports in", report.directory))
	}
	logfile <- file.path(report.directory, "quality.log")
	logger.start("Testing Report Generation", fname = logfile)

	## Initialize the report
	report <- createReport(file.path(report.directory, "example.html"), "Example")
	logger.status("Initialized the example report")
	
	## Define references
	rf1 <- "Marina Bibikova et al. Genome-wide DNA methylation profiling using Infinium assay. Epigenomics 1 (2009) 177-200"
	rf2 <- c("The ENCODE Project Consortium. ", "An integrated encyclopedia of DNA elements in the human genome. ",
		"Nature 489 (2012) 57-74")
	report <- rnb.add.reference(report, rf1)
	report <- rnb.add.reference(report, rf2)

	## Add report sections
	txt <- c("This is example report ", rnb.get.reference(report, rf1), ". It is used in testing.")
	report <- rnb.add.section(report, "Introduction", txt)
	txt <- c("Some background knowledge is required ", rnb.get.reference(report, rf2), ".")
	report <- rnb.add.section(report, "Background", txt, level = 2)

	## Create plots
	doplot <- function(type = "p", ...) {
		plot(x = c(0.4, 0.6, 0.8, 1), y = c(2, 8, 3, 9), type = type, ..., main = NA, xlab = expression(beta),
			ylab = "Measure")
	}

	report.plot1 <- createReportPlot("example_scatterplot_data_zoom", report, high.png = 200)
	doplot(pch = 16, col = c("#000080", "#00FF00"))
	off(report.plot1)

	report.plot2 <- createReportPlot("example_lines_data_zoom", report, high.png = 200)
	doplot(type = "l", lwd = 2, col = "#00FF00")
	off(report.plot2)

	report.plot3 <- createReportPlot("example_both_data_zoom", report, high.png = 200)
	doplot(type = "b", lwd = 2, col = "#00FF00")
	off(report.plot3)

	report.plot4 <- createReportPlot("example_scatterplot_data_full", report, high.png = 200)
	doplot(pch = 16, col = c("#000080", "#00FF00"), xlim = c(0, 1))
	off(report.plot4)

	report.plot5 <- createReportPlot("example_lines_data_full", report, high.png = 200)
	doplot(type = "l", lwd = 2, col = "#00FF00", xlim = c(0, 1))
	off(report.plot5)

	report.plot6 <- createReportPlot("example_both_data_full", report, high.png = 200)
	doplot(type = "b", lwd = 2, col = "#00FF00", xlim = c(0, 1))
	off(report.plot6)

	logger.info("Generated plots")

	## Generate report figure
	report.plots <- list(report.plot1, report.plot2, report.plot3, report.plot4, report.plot5, report.plot6)
	setting.names <- list(
		"plot type" = c("scatterplot" = "scatter plot", "lines" = "line plot", "both" = "line-and-point plot"),
		"values to visualize" = c("data" = "random data", "letters" = "random letters"),
		"methylation value range" = c("full" = "full", "zoom" = "zoomed in"))
	description <- c("Example figure of four data points displayed in different plots.",
		"The horizontal axis depicts methylation &beta; value, and the vertical axis represents a measurement.",
		"All values and colors were selected randomly.")
	description <- paste(description, collapse = " ")
	report <- rnb.add.section(report, "The Figure", "Here comes <a href=\"#fig1image\">Figure 1</a>:")
	report <- rnb.add.figure(report, description, report.plots, setting.names)
	
	## Close the report
	report <- rnb.add.section(report, "Summary", "The generation of this report was successful.")
	off(report)
	logger.info("Closed the example report")

	logger.completed()

	## Remove the generated files
	unlink(report.directory, recursive = TRUE)
}
