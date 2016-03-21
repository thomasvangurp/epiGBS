########################################################################################################################
## loading.R
## created: 30-05-2012
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## High-level functions handling the loading step of the data import module.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.execute.import
#'
#' Loads the data from the specified type and encapsulates it in either an \code{\linkS4class{RnBSet}}-inheriting object
#'
#' @param data.source 		non-empty \code{character} vector or \code{list} specifying the location of the data items. The
#'                    		expected format depends on the \code{data.type} that is given. See the \emph{Details} section.
#' @param data.type   		type of the input data; must be one of \code{"idat.dir"}, \code{"data.dir"}, \code{"data.files"},
#'                    		\code{"GS.report"}, \code{"GEO"} or \code{"rnb.set"}.
#' @param dry.run			if \code{TRUE} and \code{data.type} is \code{"bs.bed.dir"}, only a test data import is performed and first 10,000 lines are read from each BED file   
#' @param verbose			flag specifying whether diagnostic output should be written to the console or to the RnBeads logger 
#' 							in case the latter is initialized
#'  
#' @return Loaded data as an object of type \code{\linkS4class{RnBSet}} (when the input data type is
#'         \code{"data.dir"}, \code{"data.files"} or \code{"GEO"}) or of type \code{\linkS4class{MethyLumiSet}} (when
#'         the data type is \code{"idat.dir"} or \code{"GS.report"}).
#'
#' @details The interpretation of \code{data.source} depends on the value of \code{data.type} and is summarized in the
#' following table:
#' \tabular{lccl}{
#'   \bold{\code{data.type}} \tab \bold{Type of \code{data.source}} \tab \bold{Maximal length of \code{data.source}} \tab
#'     \bold{Interpretation}\cr
#'   \code{"infinium.idat.dir"}  \tab \code{list} or \code{character} \tab \code{2}       \tab
#'     (1) Directory containing IDAT files; (2) a sample annotation table as a \code{data.frame} or the name of the corresponding file\cr
#'   \code{"infinium.data.dir"}  \tab \code{character} \tab \code{1} \tab
#'    Directory containing data tables in plain text format. The directory should contain one file with \code{Sample|sample} token in the filename 
#' for the table of sample annotations, and one file with a token \code{beta} in the filename, with beta-values. It may also contain tables with p-values 
#' (token \code{pval}) and bead counts (\code{bead}). In the latter case beta-value, p-value and bead count tables should have matching columns and rows. 
#' The beta-value, p-value and bead-count tables should contain row names, i.e. the first column should contain the Infinium CG identifiers and not have
#' a column header (for that the first row should have one entry less than all other rows). Sample annotation table should contain as many rows as there are 
#' columns in other tables. The character used as value separator in the text tables can be set using the \code{import.table.separator} option 
#' (see \code{\link{rnb.options}} for details). \cr
#' 	 \code{"infinium.data.files"} \tab \code{character} \tab \code{2..4} \tab
#' 	The character vector should contain at least full paths to the sample annotation file and beta-value table. Detection p-values and bead counts table
#'  can be added as the third and the fourth elements. The table format requirements are the same as for \code{"data.dir"} above.\cr 
#'   \code{"infinium.GS.report"} \tab \code{character}                \tab \code{1}       \tab
#'     Genome Studio report file\cr
#'   \code{"infinium.GEO"}       \tab \code{character}                \tab \code{1}       \tab
#'     \href{http://www.ncbi.nlm.nih.gov/geo/}{GEO} identifier or downloaded series matrix file\cr
#'   \code{"bs.bed.dir"}       \tab \code{list} or \code{character}                \tab \code{1..3}       \tab
#'     (1) Directory with BED files each giving a DNA methylation profile of a sample; (2) a sample annotation table as a \code{data.frame} or the name of the corresponding file; 
#' 		(3) number of the sample annotation sheet column containing the file names. One of the first two elements have to be present. In case only the directory is specified, 
#' 		it should contain a sample annotation file with a token "sample" in the file name. In case only the sample sheet is specified, one column should be giving full absolute paths
#' 		of the BED-like files with sequencing information. If both elements (1) and (2) are specified, the files should reside in the directory, specified as element (1).
#' 		If the third element is absent, an attempt will be made to find the file name containing column automatically. For this reason the file names in the sample annotation sheet 
#' 		should be given with extensions (".bed", ".cov" etc).\cr
#'  \code{"rnb.set"}       \tab \code{\linkS4class{RnBSet}}                \tab \code{1}       \tab
#'     object of class inheriting from \code{\linkS4class{RnBSet}}\cr
#' }
#' 
#' @examples
#' \dontrun{
#' # Directory where your data is located
#' data.dir <- "~/RnBeads/data/Ziller2011_PLoSGen_450K"
#' idat.dir <- file.path(data.dir, "idat")
#' sample.annotation <- file.path(data.dir, "sample_annotation.csv")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
#' }
#' @seealso \code{\link{read.data.dir}}, \code{\link{read.idat.files}}, \code{\link{read.GS.report}},
#'          \code{\link{read.geo}}, \code{\link{read.bed.files}}
#' #' 
#' @author Pavlo Lutsik
#' @export
rnb.execute.import<-function(data.source, data.type=rnb.getOption("import.default.data.type"),
		dry.run=FALSE, verbose=TRUE) {

	if (!((is.character(data.source) || is.list(data.source) || inherits(data.source, "RnBSet")) && length(data.source) != 0)) {
		stop("invalid value for data.source; expected list or character")
	}
	if (!(is.character(data.type) && length(data.type) == 1 && (!is.na(data.type)))) {
		stop("invalid value for data.type; expected one-element character")
	}
	if (!(data.type %in% c("data.dir", "data.files", "idat.dir", "GS.report", "GEO", "bed.dir", "rnb.set",
						"infinium.data.dir", "infinium.data.files", "infinium.idat.dir", "infinium.GS.report", "infinium.GEO", 
						"bs.bed.dir","bs.data.dir","bs.data.files","rnb.set"))) {
		stop("invalid value for data.type; expected one of data.dir, idat.dir, GS.report, GEO, bed.dir or rnb.set")
	}
	
	if(data.type %in% c("idat.dir", "GS.report", "infinium.idat.dir", "infinium.GS.report")){
		if(!suppressPackageStartupMessages(require("IlluminaHumanMethylation450k.db")))
			stop("Missing required package IlluminaHumanMethylation450k.db")
	}

	if(data.type %in% c("data.dir", "infinium.data.dir", "bs.data.dir")) {

		if (!(is.character(data.source) && length(data.source) == 1 && (!any(is.na(data.source))))) {
			stop("invalid value for data.source; expected character of length 1")
		}
		
		if(!file.exists(data.source) || !file.info(data.source)[1,"isdir"]){
			rnb.error("invalid data.source parameter, data directory not found or is a file")
		}
		
		result <- read.data.dir(dir=data.source, sep=rnb.getOption("import.table.separator"))
	    
    }else if(data.type %in% c("data.files", "infinium.data.files", "bs.data.files")){
		
		if (!(is.character(data.source) && length(data.source) %in% 2L:4L && (!any(is.na(data.source))))) {
			stop("invalid value for data.source; expected character of length from 2 to 4")
		}
		
		sapply(data.source, function(dsf){
			if(!file.exists(dsf) || file.info(dsf)[1,"isdir"]){
				msg<-sprintf("invalid data.source parameter, file %s not found, or is a directory", dsf)
				rnb.error(msg)
			}
		})
		names(data.source) <- c("pheno", "betas", "p.values", "bead.counts")[1:length(data.source)]
		result <- do.call(read.data.dir, c(as.list(data.source), list(verbose=verbose)))

	} else if(data.type %in% c("idat.dir", "infinium.idat.dir")){
		
		## TODO: Fix this once read.idat.files can also work with one parameter only
	    if(!(is.character(data.source) || is.list(data.source)) && !length(data.source) %in% 1L:2L) {
			stop("invalid value for data.source; expected list or character of length 1 to 2")
		}
		data.source <- as.list(data.source)

		if(!file.exists(data.source[[1]]) || !file.info(data.source[[1]])[1,"isdir"]){
			rnb.error("invalid data.source parameter, idat.dir is not found, or is not directory")
		}
		
		if(length(data.source)==1L){
			result <- read.idat.files(base.dir=data.source[[1]], verbose=verbose
					,load.chunk=rnb.getOption("import.idat.chunk.size"))
		}else{
			if (!is.data.frame(data.source[[2]])) {
				if(!file.exists(data.source[[2]]) || file.info(data.source[[2]])[1,"isdir"]){
					rnb.error("invalid data.source parameter, sample annotation file not found, or is a directory")
				}			
			}
			result <- read.idat.files(base.dir=data.source[[1]], sample.sheet=data.source[[2]], 
					verbose=verbose, load.chunk=rnb.getOption("import.idat.chunk.size"))
		}
		gc()
		rnb.cleanMem()

	} else if(data.type %in% c("GS.report", "infinium.GS.report")) {

		if (!((is.character(data.source) || is.list(data.source)) && length(data.source) %in% c(1L,2L) && (!is.na(data.source[1])))) {
			stop("invalid value for data.source; expected one-element character")
		}
		if(!file.exists(data.source[[1]]) || file.info(data.source)[1,"isdir"]){
			rnb.error("invalid data.source parameter, GS report file not found, or is a directory")
		}		
		if(is.character(data.source) && length(data.source)==2){
			pd<-data.source[2]
		}else if(is.list(data.source) && length(data.source)==2){
			pd<-data.source[[2]]
		}else{
			pd<-NULL
		}
		
		if(rnb.getOption("import.table.separator")!="\t"){
			rnb.warning("The table separator is not tab: this is unusual for Genome Studio reports ")			
		}
		
		result <- read.GS.report(data.source[1], pd=pd, verbose=verbose)

	} else if(data.type %in% c("GEO", "infinium.GEO")) {

		if(!(is.character(data.source) && length(data.source) == 1 && (!is.na(data.source[1])))) {
			stop("invalid value for data.source; expected one-element character")
		}
		if(file.exists(data.source) && ! file.info(data.source)[1,"isdir"]) {
			result <- read.geo(filename = data.source)
		}else{
			result <- read.geo(accession = data.source)
		}
		
	} else if(data.type %in% c("bed.dir", "bs.bed.dir")){
		
		## TODO: Fix this once read.idat.files can also work with one parameter only
		if(!(is.character(data.source) || is.list(data.source)) && !length(data.source) %in% 1L:3L) {
			stop("invalid value for data.source; expected list or character of length 1 to 3")
		}
		data.source <- as.list(data.source)
		
		if(length(data.source)==1L && is.null(data.source[[1]])){
			msg<-"invalid data.source parameter, data directory should be provided"
			rnb.error(msg)
		}
		
		if(!is.null(data.source[[1L]]) && (!file.exists(data.source[[1]]) || !file.info(data.source[[1]])[1,"isdir"])){
			rnb.error("invalid data.source parameter, bed.dir is not found, or is not directory")
		}
		
		if(length(data.source)==1L){
			data.source<-c(data.source, list(NULL))
		}else{
			if(is.character(data.source[[2]])){
				if(!file.exists(data.source[[2]]) || file.info(data.source[[2]])[1,"isdir"]){
					rnb.error("invalid data.source parameter, sample annotation file is not found, or is a directory")
				}
			}
		}
		
		msg<-"No column with file names specified: will try to find one"
		if(length(data.source)<3L){
			filename.column=NA
			if(verbose){
				rnb.info(msg)
			}
		}else{
			if(is.integer(data.source[[3L]])){
				filename.column=data.source[[3L]]
			}else{
				filename.column=NA
				if(verbose){
					rnb.info(msg)
				}
			}
		}
		
		if(dry.run){
			nrows=10000L
		}else{
			nrows=-1L
		}
		
		if (rnb.getOption("import.bed.style")=="EPP"){
			result <- read.bed.files(base.dir=data.source[[1]], sample.sheet=data.source[[2]], file.names.col=filename.column,
					verbose=verbose,
					skip.lines=0,
					coord.shift.plus=0,
					is.epp.style=TRUE,
					nrows=nrows)
		} else if (rnb.getOption("import.bed.style")=="Encode"){
			result <- read.bed.files(base.dir=data.source[[1]], sample.sheet=data.source[[2]], file.names.col=filename.column,
					verbose=verbose,
					chr.col=rnb.getOption("import.bed.columns")["chr"],
					start.col=rnb.getOption("import.bed.columns")["start"],
					end.col=rnb.getOption("import.bed.columns")["end"],
					strand.col=rnb.getOption("import.bed.columns")["strand"],
					mean.meth.col=11L,
					coverage.col=10L,
					c.col=rnb.getOption("import.bed.columns")["c"],
					t.col=rnb.getOption("import.bed.columns")["t"],
					nrows=nrows)
		} else if (rnb.getOption("import.bed.style")=="bismarkCytosine"){
			result <- read.bed.files(base.dir=data.source[[1]], sample.sheet=data.source[[2]], file.names.col=filename.column,
					verbose=verbose,
					coord.shift.plus=1,
					skip.lines=0,
					chr.col=1L,
					start.col=2L,
					end.col=NA,
					c.col=4L,
					t.col=5L,
					strand.col=3L,
					mean.meth.col=NA,
					coverage.col=NA,
					coord.shift = 0L,
					nrows=nrows)
		} else if (rnb.getOption("import.bed.style")=="bismarkCov"){
			result <- read.bed.files(base.dir=data.source[[1]], sample.sheet=data.source[[2]], file.names.col=filename.column,
					verbose=verbose,
					skip.lines=0,
					chr.col=1L,
					start.col=2L,
					end.col=NA,
					c.col=5L,
					t.col=6L,
					strand.col=NA,
					mean.meth.col=NA,
					coverage.col=NA,
					coord.shift = 0L,
					nrows=nrows)
		} else {
			skip.lines <- 1
			coord.shift.plus <- 1
			result <- read.bed.files(base.dir=data.source[[1]], sample.sheet=data.source[[2]], file.names.col=filename.column,
					verbose=verbose,
					coord.shift.plus=coord.shift.plus,
					skip.lines=skip.lines,
					chr.col=rnb.getOption("import.bed.columns")["chr"],
					start.col=rnb.getOption("import.bed.columns")["start"],
					end.col=rnb.getOption("import.bed.columns")["end"],
					strand.col=rnb.getOption("import.bed.columns")["strand"],
					mean.meth.col=rnb.getOption("import.bed.columns")["meth"],
					coverage.col=rnb.getOption("import.bed.columns")["coverage"],
					c.col=rnb.getOption("import.bed.columns")["c"],
					t.col=rnb.getOption("import.bed.columns")["t"],
					nrows=nrows)
		}
		
	} else if (data.type == "rnb.set") {
		
		if (!inherits(data.source, "RnBSet")) {
			stop("invalid value for data.source: expected object of type RnBSet")
		}
		result <- data.source
			
	} else {
		stop("invalid value for data.type; expected one of data.dir, idat.dir, GS.report, GEO or rnb.set")
	}

	if(inherits(result, "MethyLumiSet")){
		result<-as(result, "RnBeadRawSet")
	}

	## TODO: Also extract the useful information from data.type and data.source, and package it in a suitable way to be
	##       passed to rnb.section.import
	return(result)
}

#######################################################################################################################

#' rnb.section.import
#'
#' A data loading section of the analysis report
#'
#' @param report    Report to contain the section. This must be an object of type \code{\linkS4class{Report}}.
#' @param object ...
#' @param data.type a characted vector of length one specifying the type of the input data: one of the \code{"data.dir"}, \code{"idat.dir"},\code{"GS.report"}
#' @param data.source a character vector specifying the location of the data items on disk. The expected length of the vector differs for different values of \code{data.type}:
#'         \describe{
#'           \item{\code{"data.dir"}}{\code{1<=length(data.source)<=4}}
#'           \item{\code{"infinium.idat.dir"}}{\code{length(data.source)=2}}
#'           \item{\code{"infinium.GS.report"}}{\code{length(data.source)=1}}
#'           \item{\code{"GEO"}}{\code{length(data.source)=1}}
#'         }
#' @return The modified report.
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.section.import<-function(report, object, data.source, data.type=rnb.getOption("import.default.data.type")){
	## TODO: Update the method and/or its documentation
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (!(inherits(object, "MethyLumiSet") || inherits(object, "RnBSet"))) {
		stop("invalid value for object; expected MethyLumiSet or RnBSet")
	}
	if (!(is.character(data.source) || is.list(data.source) || inherits(data.source, "RnBSet"))) {
		stop("invalid value for data.source; expected list or character")
	}
	if (!(is.character(data.type) && length(data.type) == 1 && (!is.na(data.type)))) {
		stop("invalid value for data.type; expected one-element character")
	}
	if (!(data.type %in% c("data.dir", "data.files", "idat.dir", "GS.report", "GEO", "bed.dir", "rnb.set",
						"infinium.data.dir", "infinium.data.files", "infinium.idat.dir", "infinium.GS.report", "infinium.GEO", 
						"bs.bed.dir","bs.data.dir","bs.data.files","rnb.set"))) {
		stop("invalid value for data.type; expected one of data.dir, idat.dir, GS.report, GEO, bed.dir or rnb.set")
	}

	pvals.present <- FALSE
	bead.counts.present <- FALSE
	pheno.table <- NULL

	if (inherits(object,"RnBSet")){
		num.samples<-ncol(meth(object))
		num.probes<-nrow(meth(object))
		if(inherits(object,"RnBeadSet"))
			pvals.present<-!is.null(dpval(object)) else pvals.present<-FALSE 
		bead.counts.present<-!is.null(covg(object))
		pheno.table <- pheno(object)
	} else { # class(object) == "MethyLumiSet"
		num.samples<-dim(phenoData(object))[1]
		num.probes<-dim(betas(object))[1]
		pvals.present<-!is.null(pvals(object))
#		bead.counts.present<-!is.null(assayData(object))
		pheno.table <- as(phenoData(object), "data.frame")
	}

	report <- rnb.add.section(report, "Loading", "This section presents a summary of the data import step.")

	if(data.type %in% c("data.dir","infinium.data.dir","bs.data.dir")) dtype="data directory"
	if(data.type %in% c("idat.dir", "infinium.idat.dir")) dtype=".idat files"
	if(data.type %in% c("GS.report", "infinium.GS.report")) dtype="GenomeStudio report"
	if(data.type %in% c("GEO", "infinium.GEO")) dtype="GEO data set"
	if(data.type %in% c("bed.dir", "bs.bed.dir")) dtype="Bisulfite sequencing data set"
	if(data.type=="rnb.set") dtype="RnBead Set object"
	descr<-c(
		   "Data type",
	       "Data source",
		   if(class(object)=="RnBiseqSet") "Genome assembly used" else NULL,
           "Number of samples",
           sprintf("Number of %s",rnb.get.row.token(object, plural=TRUE)),
           "Detection p-values",
		   sprintf("%s information is", rnb.get.covg.token(object)),
           if(inherits(object, "RnBeadSet")) "Quality control information" else NULL)

	vals<-c(dtype,
		if(is.list(data.source)){data.source[[1]] }else if(is.character(data.source)){ data.source[1]} else deparse(substitute(data.source)),
		if(class(object)=="RnBiseqSet") assembly(object) else NULL,
		num.samples,
		num.probes,
		if(pvals.present) "present" else "absent",
		if(bead.counts.present) "present" else "absent",
		if(inherits(object, "RnBeadSet")){if(class(object)=="MethyLumiSet") "present" else "absent"} else NULL)
	table.df<-data.frame(Descr=as.character(descr), Vals=as.character(vals))
	colnames(table.df)<-NULL
	rnb.add.table(report, table.df, row.names=FALSE, first.col.header=TRUE)

	## Study the traits in the sample annotation table
	get.trait.type <- function(x, any.missing = any(is.na(x))) {
		if (length(unique(x)) == 1) {
			if (all(is.na(x))) {
				return("all missing")
			}
			return("all identical")
		}
		if (is.logical(x)) {
			return("yes/no")
		}
		ids <- (any.missing == FALSE) && (anyDuplicated(x) == 0)
		if (is.numeric(x) || is.integer(x)) {
			return(paste0("numbers", ifelse(ids, " or identifiers", "")))
		}
		# is.factor(x) || is.character(x)
		return(paste(ifelse(ids, "identifiers", "categories"), "or text"))
	}
	report <- rnb.add.section(report, "Sample Annotation", NULL, level = 2)
	if (is.null(pheno.table) || (nrow(pheno.table) * ncol(pheno.table) == 0)) {
		rnb.add.paragraph(report, "The loaded dataset does not include sample annotation information.")
	} else {
		txt <- "The loaded dataset includes a table of sample annotations describing the following traits:"
		rnb.add.paragraph(report, txt)
		traits.summary <- data.frame(
			"No." = 1:ncol(pheno.table),
			"Trait" = colnames(pheno.table),
			"Missing Values" = as.character(sapply(pheno.table, function(x) { sum(is.na(x)) })),
			"Type" = as.character(NA), check.names = FALSE, stringsAsFactors = FALSE)
		traits.summary[traits.summary[["Missing Values"]] == "0", "Missing Values"] <- ""
		traits.summary[, "Type"] <- sapply(pheno.table, get.trait.type)
		colnames(traits.summary)[3] <- paste("<span title=\"",
			"number of samples that are not assigned a value for the respective trait",
			"\">Missing Values</span>", sep = "")
		rnb.add.table(report, traits.summary, row.names = FALSE)

		## Export the annotation table
		fname <- "annotation.csv"
		fname <- rnb.write.table(pheno.table, fname, fpath = rnb.get.directory(report, "data", absolute = TRUE),
			format = "csv", gz = FALSE, row.names = FALSE)
		txt <- c("The full sample annotation table is available as a <a href=\"", rnb.get.directory(report, "data"), "/",
			fname, "\">", "comma-separated value file</a> accompanying this report.")
		rnb.add.paragraph(report, txt)
	}

	return(report)
}
########################################################################################################################
#' rnb.step.import
#'
#' Loads the data of the specified type and creates the corresponding report section.
#'
#' @param data.source Non-empty \code{character} vector or \code{list} specifying the location of the data items. The
#'                    expected format depends on the \code{data.type} that is given. See the documentation of
#'                    \code{\link{rnb.execute.import}} for more details.
#' @param data.type   Type of the input data; must be one of \code{"data.dir"}, \code{"data.files"}, \code{"infinium.idat.dir"}, \code{"infinium.GS.report"},
#'                    \code{"infinium.GEO"} or \code{"rnb.set"}.
#' @param report      Report on to contain the Loading section. This must be an object of type
#'                    \code{\linkS4class{Report}}.
#'
#' @return a list with two elements: \describe{
#'
#'                        \item{\code{object}}{Loaded data as \code{\linkS4class{RnBSet}} (in case \code{data.type="data.dir"}) or \code{\linkS4class{MethyLumiSet}} (in case \code{data.type=c("idat.dir", "GS.report")}}
#'                        \item{\code{report}}{the modified report}
#'
#' }
#' @author Pavlo Lutsik
#' @noRd
rnb.step.import <- function(data.source, data.type = rnb.getOption("import.default.data.type"), report) {
	if (!(is.character(data.source) || is.list(data.source) || inherits(data.source, "RnBSet"))) {
		stop("invalid value for data.source; expected list or character")
	}
	if (!(is.character(data.type) && length(data.type) == 1 && (!is.na(data.type)))) {
		stop("invalid value for data.type; expected one-element character")
	}
	if (!(data.type %in% c("data.dir", "idat.dir", "GS.report", "GEO", "bed.dir","rnb.set",
						"infinium.data.dir", "infinium.data.files","infinium.idat.dir", "infinium.GS.report", "infinium.GEO", 
						"bs.bed.dir","bs.data.dir","bs.data.files","rnb.set"))) {
		
		stop("invalid value for data.type; expected one of data.dir, idat.dir, GS.report, GEO, bed.dir or rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}
	
	logger.info(sprintf("Loading data of type \"%s\"",data.type))
	
	## Load the data into an object
	if(data.type %in% c("bed.dir", "bs.bed.dir") && rnb.getOption("import.bed.test")){
		logger.start("Performing loading test")
		logger.info("The first 10000 rows will be read from each data file")
		object <- rnb.execute.import(data.source, data.type, dry.run=TRUE)
		logger.start("Checking the loaded object")
		valid<-check.rnb.biseq.set(object, verbose=TRUE)
		if(valid){
			logger.info("The object loaded during the loading test is valid")
		}else{
			logger.warning(c("The object loaded during the loading test contains invalid information (see details above).",
							"Please check the whether the data source arguments as well as the data import options, like table separator, BED style or BED column assignment, are set correctly"))
		}
		logger.completed()
		rnb.cleanMem()
		logger.completed()
	}
	
	if(!data.type %in% c("bed.dir", "bs.bed.dir") || !rnb.getOption("import.bed.test.only")){
		object <- rnb.execute.import(data.source, data.type)
		rnb.cleanMem()
		if(data.type=="bed.dir"){
			logger.start("Checking the loaded object")
			valid<-check.rnb.biseq.set(object, verbose=TRUE)
			if(valid){
				logger.info("The loaded object is valid")
			}else{
				logger.warning(c("The loaded object contains invalid information (see details above).",
								"Please check the whether the data source arguments as well as the data import options, like table separator, BED style or BED column assignment, are set correctly"))
			}
			logger.completed()
		}
	}
	
	if (inherits(data.source, "RnBSet")) {
		d.source <- paste("object of type", class(data.source))
	} else { # inherits(object, "")
		d.source <- as.character(data.source[1])
	}
	logger.status(c("Loaded data from", d.source))
	
	## Create a section in the report
	report <- rnb.section.import(report, object, data.source, data.type)
	logger.status("Added data loading section to the report")
	if (inherits(object, "RnBSet")) {
		d.source <- paste("object of type", class(data.source))
		nsamples <- length(samples(object))
		nsites <- nrow(meth(object)) # FIXME: See if there is a faster way to get the number of sites
	} else { # inherits(object, "MethyLumiSet")
		nsamples <- dim(phenoData(object))[1]
		nprobes <- dim(betas(object))[1]
	}
	logger.status(c("Loaded", nsamples, "samples and", nsites, "sites"))
	logger.info(c("Output object is of type", class(object)))
	
	return(list(rnb.set=object, report=report))
}
