########################################################################################################################
## Report-class.R
## created: 2012-05-08
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Report class definition.
########################################################################################################################

## C L A S S ###########################################################################################################

#' Report Class
#'
#' Handler of a generated HTML report. Reports are initialized using the function \code{\link{createReport}}.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{fname}}{Name of the file that contains the HTML report.}
#'   \item{\code{dir.conf}}{Directory that contains configuration files; usually shared between reports.}
#'   \item{\code{dir.data}}{Directory that contains the generated external lists and tables.}
#'   \item{\code{dir.pngs}}{Directory that contains the generated figure image files.}
#'   \item{\code{dir.pdfs}}{Directory that contains the generated figure PDF files.}
#'   \item{\code{dir.high}}{Directory that contains the generated high-resolution image file.}
#'   \item{\code{sections}}{Number of sections and subsections currently added to the report.}
#'   \item{\code{opensections}}{Indices of currently active section and subsections.}
#'   \item{\code{figures}}{Number of figures currently added to the report.}
#'   \item{\code{tables}}{Number of selectable tables added to the report.}
#'   \item{\code{references}}{List of references to be added at the end of the report.}
#' }
#'
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{\link{rnb.get.directory}}}{Gets the location of a given report-specific directory.}
#'   \item{\code{\link{rnb.add.section}}}{Generates HTML code for a new section in the report.}
#'   \item{\code{\link{rnb.add.paragraph}}}{Generates HTML code for a new paragraph in the report.}
#'   \item{\code{\link{rnb.add.list}}}{Generates HTML code for a list in the report.}
#'   \item{\code{\link{rnb.add.table}}}{Generates HTML code for a table in the report.}
#'   \item{\code{\link{rnb.add.tables}}}{Generates HTML code for a listing of tables in the report.}
#'   \item{\code{\link{rnb.add.figure}}}{Generates HTML code for a figure in the report.}
#'   \item{\code{\link{rnb.add.reference}}}{Adds a reference item to the report.}
#'   \item{\code{\link[=off,Report-method]{off}}}{Completes the HTML report by adding a reference section (if needed),
#'        a footer notice and closing the \code{<body>} and \code{<html>} tags.}
#' }
#'
#' @name Report-class
#' @rdname Report-class
#' @aliases initialize,Report-method
#' @author Yassen Assenov
#' @exportClass Report
setClass("Report",
	representation(fname = "character",
		dir.conf = "character",
		dir.data = "character",
		dir.pngs = "character",
		dir.pdfs = "character",
		dir.high = "character",
		sections = "integer",
		opensections = "integer",
		figures = "integer",
		tables = "integer",
		references = "character"),
	prototype = prototype(fname = "",
		dir.conf = "configutation",
		dir.data = "data",
		dir.pngs = "images",
		dir.pdfs = "images/pdf",
		dir.high = "images/high_resolution",
		sections = 0L,
		opensections = rep(0L, 3L),
		figures = 0L,
		tables = 0L,
		references = character()),
	package = "RnBeads")
