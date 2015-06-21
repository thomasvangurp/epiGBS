########################################################################################################################
## test_report_3.R
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Simple tests on creating and adding listing of tables to a report.
########################################################################################################################

test_report_3 <- function() {

	NAME.POOL <- c("Anna Prentice", "Annie Versary", "Ben Zeen", "Bill Ding", "Cole Desice", "Dee Zaster", "Eddie Bull",
		"Emma Fraid", "Flo Tinaway", "Hugh Manatty", "Isabel Ringing", "Joe King", "Justin Case", "Kerry Dowt",
		"Lois Steem", "Markus Absent", "Matt Schtick", "Maya Magination", "Miss Teri Novel", "Nick Rofillia",
		"Otto Matic", "Perry Scope", "Ray Neday", "Robbin Banks", "Shirley Knott", "Terry Torrie", "Trudy Lite",
		"Vic Tree", "Walter Mellon", "Warren Peece")

	## Generate a few random tables

	setting.names <- list(
		"Eye color" = c("blue" = "blue", "brown" = "brown", "black" = "black", "green" = "green"),
		"Hair color" = c("b" = "blond", "o" = "brown", "r" = "red", "k" = "black"))

	celebrities <- list()
	for (n1 in names(setting.names[[1]])) {
		for (n2 in names(setting.names[[2]])) {
			celebrities[[paste(n1, n2, sep = "_")]] <- matrix(sample(NAME.POOL, 12), nrow = 4, ncol = 3,
				dimnames = list(c("Europe", "Asia", "Africa", "Australia"), c("A", "B", "C")))
		}
	}
	rm(n1, n2)

	## Initialize a report
	report.directory <- "RnBeads_report_testing"
	report.file <- file.path(report.directory, "test_tables.html")
	dirs <- c("configuration" = "data", "data" = "data", "pdfs" = "data", "pngs" = "data", "high" = "data")
	report <- createReport(report.file, "Testing tables", dirs = dirs, init.configuration = TRUE)
	report <- rnb.add.section(report, "Celebrity Names", "There they are:")
	report <- rnb.add.tables(report, celebrities, setting.names, row.names = FALSE)
	off(report)

	unlink(report.directory, recursive = TRUE)
}
