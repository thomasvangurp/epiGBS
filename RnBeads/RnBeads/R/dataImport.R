########################################################################################################################
## dataImport.R
## created: 2012-04-04
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Data reading routines.
########################################################################################################################

## G L O B A L S #######################################################################################################

## Strings in a sample annotation table that are read as missing values (NAs)
NA.STRINGS <- c("NA","N/A","N.A.","n/a","n.a.","")

## Regular expression for a section header in a CSV file exported from Genome Studio
REGEX.CSV.HEADER <- "^\\[(.+)\\],*$"

## F U N C T I O N S ###################################################################################################

#' check.barcode
#'
#' Validates that the table of sample phenotype information contains a column named \code{"barcode"}, and tries to
#' construct such a column ifit doesn't exist.
#'
#' @param dframe Table of sample phenotype data.
#' @return The possibly modified data frame; \code{NULL} ifthe provided data frame does not a column \code{"barcode"}
#'         and such a column cannot be inferred.
#' @author Yassen Assenov
#' @noRd
check.barcode <- function(dframe){
	if("barcode" %in% colnames(dframe)){
		return(dframe)
	}
	if(all(c("Sentrix_ID", "Sentrix_Position") %in% colnames(dframe))){
		dframe[, "barcode"] <- paste(dframe[, "Sentrix_ID"], dframe[, "Sentrix_Position"], sep = "_")
		return(dframe)
	}
	if(all(c("Sentrix ID", "Sentrix Position") %in% colnames(dframe))){
		dframe[, "barcode"] <- paste(dframe[, "Sentrix ID"], dframe[, "Sentrix Position"], sep = "_")
		return(dframe)
	}
	id.column <- grep("Sentrix_ID", colnames(dframe), fixed = TRUE)[1]
	position.column <- grep("Sentrix_Position", colnames(dframe), fixed = TRUE)[1]
	if(!(is.na(id.column) || is.na(position.column))){
		dframe[, "barcode"] <- paste(dframe[, id.column], dframe[, position.column], sep = "_")
		return(dframe)
	}
	return(NULL)
}

########################################################################################################################

#' read.sample.annotation
#'
#' Reads Illumina Infinium sample annotation.
#'
#' @param fname Name of text file that contains a sample annotation table with a header. This method handles a variety
#'              of file formats, including comma-separated values file exported from Genome Studio.
#' @param sep   One-element \code{character} used as field separator in the tables file.
#' @return Sample annotation table in the form of a \code{data.frame}, in which every row corresponds to a sample, and
#'         every column - to a trait.
#'
#' @author Pavlo Lutsik
#' @examples
#' \dontrun{
#'   TODO add system annotation example
#'   annotation.file<-system.file("")
#'   sa<-read.sample.annotation(annotation.file)
#'   sa
#' }
#' @export 
read.sample.annotation <- function(fname, sep = rnb.getOption("import.table.separator")) {
	if(!(is.character(fname) && length(fname) == 1 && (!is.na(fname)))) {
		stop("invalid value for fname")
	}
	if(!(is.character(sep) && length(sep) == 1 && (!is.na(sep)) && (nchar(sep) != 0))) {
		stop("invalid value for sep")
	}

	file.contents <- scan(fname, what = "", sep = "\n", quiet = TRUE)
	i.headers <- grep(REGEX.CSV.HEADER, file.contents)
	i.first = 1L
	i.last = length(file.contents)
	if(length(i.headers) != 0) {
		headers <- toupper(sub(REGEX.CSV.HEADER, "\\1", file.contents[i.headers]))
		i.data <- which(headers == "DATA")
		if(length(i.data) == 0) {
			rnb.warning(c("File", fname, "may contain sections but [Data] was not found"))
		} else {
			if(length(i.data) > 1) {
				rnb.warning(c("File", fname, "contains multiple sections [Data]; only the first one is read"))
				i.data <- i.data[1]
			}
			i.first <- i.headers[i.data] + 1L
			if(i.data < length(i.headers)) {
				i.last <- i.headers[i.data + 1] - 1L
			}
			if(i.first > i.last) {
				logger.error(c("File", fname, "contains empty section [Data]"))
			}
		}
		rm(i.data)
	}
	rm(i.headers)

	comment.char <- ifelse(sep == "\t", "#", "")
	result <- read.table(header = TRUE, sep = sep, quote = "\"", comment.char = comment.char,
		na.strings = NA.STRINGS, check.names = FALSE, text = file.contents[i.first:i.last])
	if(ncol(result) < 2) {
		rnb.warning("The sample sheet table has only one column. Please check import.table.separator option")
	}
	
	if("" %in% colnames(result)){
		cn<-colnames(result)
		cn[cn==""]<-paste("UnknownCol_",1:length(which(cn=="")))
		colnames(result)<-cn
		rnb.warning("The sample sheet contains columns with empty names, renamed")
	}
	
	result
}

#######################################################################################################################

#' read.data.dir
#'
#' Reads in a directory with Illumina Infinium HumanMethylation450 data. The files shoudl be stored as data
#'
#' @param dir directory containing the table files
#' @param pheno a file containing data sample annotations and phenotypic information
#' @param betas a file containing the beta values. If not supplied, the routine will look in dir for a file containing "beta" token in the filename
#' @param p.values a file containing the detection p values. If not supplied, the routine will look in dir for a file containing "pval" token in the filename
#' @param bead.counts a file containing the bead counts (optional). If not supplied, the routine will look in dir for a file containing "bead" token in the filename
#' @param sep character used as field separator in the tables files. Default value is taken by the call to \code{rnb.getOption("import.table.separator")}
#' @param verbose Flag indicating ifthe messages to the logger should be sent. Note that the logger must be initialized prior to calling this function. 
#' 				  Logging is useful for keeping a record of the downloaded and processed samples. Also, informative messages are stored in case of an error.
#' @details Colnames in all files should match. They will be returned as the samples element of the list.
#' @return Object of type \code{\linkS4class{RnBeadSet}}.
#'
#' @author Pavlo Lutsik
#' @export
read.data.dir<-function(dir,
				pheno,
				betas,
				p.values,
				bead.counts, 
				sep=rnb.getOption("import.table.separator"),
				verbose=TRUE){
			
			if(verbose){
				rnb.logger.start("Loading Data from Table Files")
			}
			
			if(!missingArg(dir)){
				files<-list.files(dir, full.names=T)
			}else{
				files<-character()
			}
			
			if(length(files)>10){
				stop(paste("Too many files in directory:", dir), length(files))
			}
			
			if(missingArg(pheno)){
				pheno<-grep("Sample|sample", files, value=TRUE)
				if(length(pheno)==0){
					stop("Missing file with sample information:", files)
				}else if(length(pheno)!=1){
					stop("More than one file with sample information:", files)
				}
			}
			
			if(missingArg(betas)){
				if(length(grep("beta", files))==1){
					betas<-files[grep("beta", files)] 
				}else{ 
					stop("Missing or more than one file with betas:", files)
				}
			}
			
			if(missingArg(p.values)){
				if(length(grep("pval", files))==1){
					p.values<-files[grep("pval", files)]
				}else{
					p.values<-NULL
				}
				#else stop("Missing or more than one file with p values:", files)
			}
			
			if(missingArg(bead.counts)){
				if(length(grep("bead", files)==1)){
					bead.counts<-files[grep("bead", files)]
				}else{
					bead.counts<-NULL
				}
				#else stop("Missing or more than one file with bead counts:", files)
			}
			
			pheno.table<-read.sample.annotation(pheno, sep=sep)
			if(verbose){
				rnb.info("Read table with sample information")
			}
			
			beta.table<-read.table(betas, sep=sep, header=T, check.names=F)
			if(verbose){
				rnb.info("Read table with betas")
			}
			if(ncol(beta.table)<2){
				rnb.warning("Beta-value table has less than two columns: check the value of sep")
			}
			if(is.null(rownames(beta.table))){
				rnb.error("Beta-value table has no row names: cannot match probes to the annotation")
			}
			
			if(!is.null(p.values)){
				p.value.table<-read.table(p.values, sep=sep, header=T, check.names=F)
				if(verbose){
					rnb.info("Read table with p-values")
				}
				if(ncol(p.value.table)<2){
					rnb.warning("P-value table has less than two columns: check the value of sep")
				}
			}else p.value.table<-NULL  
			
			if(!is.null(bead.counts)){
				bead.count.table<-read.table(bead.counts, sep=sep, header=T, check.names=F)
				if(verbose){
					rnb.info("Read table with bead counts")
				}
				if(ncol(bead.count.table)<2){
					rnb.warning("Bead count table has less than two columns: check the value of sep")
				}
			}else {
				bead.count.table<-NULL
			}
			
			
			if(!is.null(p.value.table)){
				if(sum(colnames(beta.table) != colnames(p.value.table))>0){
					stop("The columns of p-value table do not match the columns of the beta value table ", "")
				}
			}
			if(!is.null(bead.count.table)){
				if(sum(colnames(beta.table) != colnames(bead.count.table))>0){
					stop("The columns of bead counts table do not match the columns of the beta value table ", "")
				}
			}

			data.set<-new("RnBeadSet",
					pheno=pheno.table,
					betas=as.matrix(beta.table),
					p.values=if(is.null(p.value.table)) p.value.table else as.matrix(p.value.table),
					bead.counts=if(is.null(bead.count.table)) bead.count.table else as.matrix(bead.count.table))
			
			if(verbose) {
				rnb.logger.completed()
			}
			return(data.set)
			
}
########################################################################################################################

#' read.idat.files
#'
#' Reads a directory of \code{.idat} files and initializes an object of type \code{\linkS4class{MethyLumiSet}}.
#'
#' @param base.dir       Directory that contains the \code{.idat} files to be read; or a character vector of such
#'                       directories.
#' @param barcodes       Optional non-empty \code{character} vector listing the barcodes of the samples that should be
#'                       loaded. If supplied, this vector must not contain \code{NA} among its elements.
#' @param sample.sheet   Optional file name containing a table of sample annotation data, or the table itself in the
#'                       form of a \code{\link{data.frame}} or \code{matrix}. Only (and all) samples defined in this
#'                       table will be loaded. The table is expected to contain a column named \code{"barcode"} that
#'                       lists the samples' Sentrix  barcodes. If such a column is not present, this function searches
#'                       for columns \code{"Sentrix_ID"} and \code{"Sentrix_Position"} (or similar) that build a
#'                       barcode.
#' @param sep.samples	 character used as field separator in the sample sheet file. 
#' 						 Default value is taken by the call to \code{rnb.getOption("import.table.separator")}
#' @param load.chunk	 \code{integer} of size one, giving the number of IDAT files which should be loaded in 
#' 						 one loading cycle or \code{NULL}, in which case an attempt will be made to load all files 
#' 						 in one go. Should be assigned in case the number of IDATs is more than one thousand.
#' @param keep.methylumi a flag indicating whether the a \code{MethyLumiSet}
#'						 object should be returned instead of a \code{RnBeadRawSet}.
#' @param verbose        Flag indicating ifthe messages to the logger should be sent. Note that the logger
#'                       must be initialized prior to calling this function. Logging is useful for keeping a
#'                       record of the downloaded and processed samples. Also, informative messages are
#'                       stored in case of an error.
#' @return Loaded dataset of HumanMethylation450K samples, encapsulated in an object of type \code{MethyLumiSet}.
#'
#' @details
#' If neither \code{barcodes}, nor \code{sample.sheet} are specified, the function attempts to locate a file in
#' \code{base.dir} containing sample annotation information. It fails ifsuch a file cannot be (unambiguously)
#' identified. If both \code{barcodes} and \code{sample.sheet} are supplied, only \code{sample.sheet} is used in loading
#' methylation data. The value of \code{barcodes} is tested for validity but it is not used as a filter.
#'
#' @seealso \code{\link{methylumIDAT}} in package \pkg{methylumi}
#'
#' @author Pavlo Lutsik
#' @export
read.idat.files <- function(base.dir, 
		barcodes = NULL, 
		sample.sheet = NULL,
		sep.samples=rnb.getOption("import.table.separator"),
		load.chunk=NULL,
		keep.methylumi=FALSE,
		verbose = TRUE){
	
	if(!(is.character(base.dir) && length(base.dir) == 1 && (!is.na(base.dir)))) {
		stop("invalid value for base.dir")
	}
	if(!is.null(barcodes)) {
		if(!is.character(barcodes) && length(barcodes) != 0 && (!any(is.na(barcodes)))) {
			stop("invalid value for barcodes; expected non-empty character vector or NULL")
		}
	}
	if(!parameter.is.flag(verbose)){
		stop("invalid value for verbose; expected TRUE or FALSE")
	}

	if(verbose){
		rnb.logger.start("Loading Data from IDAT Files")
	}

	if(is.null(barcodes) && is.null(sample.sheet)){
		## Attempt to locate sample annotation file
		sample.sheet <- list.files(base.dir, pattern = "[Ss]ample", full.names = TRUE)
		sample.sheet <- sample.sheet[file.info(sample.sheet)[, "isdir"] == FALSE]
		if(length(sample.sheet) > 1) {
			rnb.error(paste("No unique candidate for sample annotation found in", base.dir))
		}
		if(length(sample.sheet) == 0 || is.na(sample.sheet)){
			rnb.error(paste("No candidate for sample information found in", base.dir))
		}
	}
	if(!is.null(sample.sheet)) {
		if(is.character(sample.sheet) && (!is.matrix(sample.sheet))){
			if(length(sample.sheet) != 1) {
				stop("invalid value for sample.sheet")
			}
			sample.sheet <- read.sample.annotation(sample.sheet, sep=sep.samples)
		} else if(is.matrix(sample.sheet)){
			sample.sheet <- as.data.frame(sample.sheet, check.names = FALSE, stringsAsFactors = FALSE)
		}
		if(!is.data.frame(sample.sheet)){
			stop("invalid value for sample.sheet")
		}
		if(!("barcode" %in% colnames(sample.sheet))){
			sample.sheet <- check.barcode(sample.sheet)
			if(is.null(sample.sheet)){
				rnb.error("sample.sheet does not contain column barcode")
			}
			if(verbose) {
				rnb.info("Added column barcode to the provided sample annotation table")
			}
		}
	}
	if(is.null(barcodes) && is.null(sample.sheet)){
		stop("invalid value for barcodes and/or sample.sheet")
	}
	if(ncol(sample.sheet) < 2){
		rnb.error("The sample annotation table has less than two columns. Check the \"default.table.separator\" option")
	}

	## Check subdirectories
	if(check.idat.subdirs(base.dir)){
		## Found subdirectories that contain .idat files
		if(verbose) {
			rnb.info("Found subdirectores with idat files; creating symbolic links")
		}
		base.dir <- prepare.idat.dir(base.dir)
		if(is.null(base.dir)){
			rnb.error("Could not create temporary directory or create links to idat files")
		}
		fn <- attr(base.dir, "failed")
		if(fn != 0) {
			rnb.warning(paste("Could not create links to", fn, "idat files; these files will be ignored"))
		}
		rm(fn)
	}
	
	## Check whether all files a present
	
	if(!is.null(sample.sheet)){
		fn.base<-sample.sheet[,"barcode"]
	}else{
		fn.base<-barcodes
	}
	
	fn.expected<-paste(fn.base, c("Red.idat", "Grn.idat"), sep="_")
	fn.available<-list.files(base.dir, pattern="idat")
	fn.missing<-setdiff(fn.expected, fn.available)
	if(length(fn.missing)>0){
		rnb.error(sprintf("Some IDAT files are not present in the supplied base directory, for instance %s", 
				paste(fn.missing[1:min(6,length(fn.missing))], collapse=", ")))
	}
	
	## Read the data using the methylumi package
	methylumi.params <- list(idatPath = base.dir, n = TRUE)
	if(is.null(sample.sheet)){
		if(!is.null(load.chunk) && load.chunk>length(barcodes)){
			rnb.warning("The value of the load.chunk parameter is incorrect.")
			load.chunk=NULL
			nsamp<-length(barcodes)
		}
	}else{
		if(!is.null(load.chunk) && load.chunk>nrow(sample.sheet)){
			rnb.warning("The value of the load.chunk parameter is incorrect.")
			load.chunk=NULL
			nsamp<-nrow(sample.sheet)
		}
	}
		
	if(is.null(load.chunk)){
		if(is.null(sample.sheet)){
			methylumi.params$barcodes <- barcodes
		}else{ # is.data.frame(sample.sheet)
			methylumi.params$pdat <- sample.sheet
		}
		mls <- suppressWarnings(suppressMessages(do.call(methylumIDAT, methylumi.params)))
		if(!keep.methylumi){
			mls<-as(mls, "RnBeadRawSet")
		}
	}else{
		# load the data chunk-by-chunk, converting to RnBeadRawSet 
        # and combining them after every step	
		
		mls <- NULL
		if(is.null(sample.sheet)) {
			nsamp<-length(barcodes)
		}else{
			nsamp<-nrow(sample.sheet)
		}
		nch <- nsamp %/% load.chunk
		
		## TODO: remove this strange thing
		# temporarily disable region summarization
		rt<-rnb.getOption("region.types")
		rnb.options(region.types=character())
		####
		
		for(ch in 1:(nch+(nsamp %% load.chunk > 0))){
			if(verbose){
				rnb.info(sprintf("Loading sample chunk %d", ch))
			}
			
			ixxs <- (load.chunk*(ch-1L)+1L):min((load.chunk*ch), nsamp)
			
			if(is.null(sample.sheet)){
				methylumi.params$barcodes <- barcodes[ixxs]
			} else { # is.data.frame(sample.sheet)
				methylumi.params$pdat <- sample.sheet[ixxs,]
			}
			
			mlsch <- suppressWarnings(suppressMessages(do.call(methylumIDAT, methylumi.params)))
			
			rnbch <- as(mlsch, "RnBeadRawSet")
			
			if(is.null(mls)){
				mls<-rnbch
			}else{
				mls.old <- mls
				# a couple of disk space saving hacks
				if(rnb.getOption("enforce.destroy.disk.dumps")){
					if(mls.old@status$disk.dump){
						mls.old@status$discard.ff.matrices<-TRUE
					}
					if(rnbch@status$disk.dump){
						rnbch@status$discard.ff.matrices<-TRUE
					}
				}
				mls<-combine(mls.old, rnbch)
				if(rnb.getOption("enforce.destroy.disk.dumps")){
					if(mls@status$disk.dump){
						mls@status$discard.ff.matrices<-NULL
					}
				}
				rnb.call.destructor(mls.old)
				rnb.call.destructor(rnbch)
			}
			
			rm(rnbch)
			rnb.cleanMem()
		}
		#gc()
		
		rnb.options(region.types=rt)
		
		for (region.type in rnb.region.types.for.analysis("hg19")) {
			mls <- summarize.regions(mls, region.type)
		}
		rnb.cleanMem()
		
	}
	
	if(verbose){
		rnb.logger.completed()
	}
	
	return(mls)

}

########################################################################################################################

#' read.idat.files2
#'
#' Reads a directory of \code{.idat} files and initializes an object of type \code{\linkS4class{MethyLumiSet}}.
#'
#' @param base.dir       Directory that contains the \code{.idat} files to be read; or a character vector of such
#'                       directories.
#' @param barcodes       Optional non-empty \code{character} vector listing the barcodes of the samples that should be
#'                       loaded. If supplied, this vector must not contain \code{NA} among its elements.
#' @param sample.sheet   Optional file name containing a table of sample annotation data, or the table itself in the
#'                       form of a \code{\link{data.frame}} or \code{matrix}. Only (and all) samples defined in this
#'                       table will be loaded. The table is expected to contain a column named \code{"barcode"} that
#'                       lists the samples' Sentrix  barcodes. If such a column is not present, this function searches
#'                       for columns \code{"Sentrix_ID"} and \code{"Sentrix_Position"} (or similar) that build a
#'                       barcode.
#' @param sep.samples	 character used as field separator in the sample sheet file. 
#' 						 Default value is taken by the call to \code{rnb.getOption("import.table.separator")}
#' @param load.chunk	 \code{integer} of size one, giving the number of IDAT files which should be loaded in 
#' 						 one loading cycle or \code{NULL}, in which case an attempt will be made to load all files 
#' 						 in one go. Should be assigned in case the number of IDATs is more than one thousand.
#' @param verbose        Flag indicating ifthe messages to the logger should be sent. Note that the logger
#'                       must be initialized prior to calling this function. Logging is useful for keeping a
#'                       record of the downloaded and processed samples. Also, informative messages are
#'                       stored in case of an error.
#' @return Loaded dataset of HumanMethylation450K samples, encapsulated in an object of type \code{MethyLumiSet}.
#'
#' @details
#' If neither \code{barcodes}, nor \code{sample.sheet} are specified, the function attempts to locate a file in
#' \code{base.dir} containing sample annotation information. It fails ifsuch a file cannot be (unambiguously)
#' identified. If both \code{barcodes} and \code{sample.sheet} are supplied, only \code{sample.sheet} is used in loading
#' methylation data. The value of \code{barcodes} is tested for validity but it is not used as a filter.
#'
#' @seealso \code{\link{methylumIDAT}} in package \pkg{methylumi}
#'
#' @author Pavlo Lutsik
#' @noRd 
read.idat.files2 <- function(base.dir, 
		barcodes = NULL, 
		sample.sheet = NULL,
		sep.samples=rnb.getOption("import.table.separator"),
		load.chunk=NULL,
		verbose = TRUE){
	
	if(!(is.character(base.dir) && length(base.dir) == 1 && (!is.na(base.dir)))) {
		stop("invalid value for base.dir")
	}
	if(!is.null(barcodes)) {
		if(!is.character(barcodes) && length(barcodes) != 0 && (!any(is.na(barcodes)))) {
			stop("invalid value for barcodes; expected non-empty character vector or NULL")
		}
	}
	if(!parameter.is.flag(verbose)){
		stop("invalid value for verbose; expected TRUE or FALSE")
	}
	
	if(verbose){
		rnb.logger.start("Loading Data from IDAT Files")
	}
	
	if(is.null(barcodes) && is.null(sample.sheet)){
		## Attempt to locate sample annotation file
		sample.sheet <- list.files(base.dir, pattern = "[Ss]ample", full.names = TRUE)
		sample.sheet <- sample.sheet[file.info(sample.sheet)[, "isdir"] == FALSE]
		if(length(sample.sheet) > 1) {
			rnb.error(paste("No unique candidate for sample annotation found in", base.dir))
		}
		if(length(sample.sheet) == 0 || is.na(sample.sheet)){
			rnb.error(paste("No candidate for sample information found in", base.dir))
		}
	}
	if(!is.null(sample.sheet)) {
		if(is.character(sample.sheet) && (!is.matrix(sample.sheet))){
			if(length(sample.sheet) != 1) {
				stop("invalid value for sample.sheet")
			}
			sample.sheet <- read.sample.annotation(sample.sheet, sep=sep.samples)
		} else if(is.matrix(sample.sheet)){
			sample.sheet <- as.data.frame(sample.sheet, check.names = FALSE, stringsAsFactors = FALSE)
		}
		if(!is.data.frame(sample.sheet)){
			stop("invalid value for sample.sheet")
		}
		if(!("barcode" %in% colnames(sample.sheet))){
			sample.sheet <- check.barcode(sample.sheet)
			if(is.null(sample.sheet)){
				rnb.error("sample.sheet does not contain column barcode")
			}
			if(verbose) {
				rnb.info("Added column barcode to the provided sample annotation table")
			}
		}
	}
	if(is.null(barcodes) && is.null(sample.sheet)){
		stop("invalid value for barcodes and/or sample.sheet")
	}
	if(ncol(sample.sheet) < 2){
		rnb.error("The sample annotation table has less than two columns. Check the \"default.table.separator\" option")
	}
	
	## Check subdirectories
	if(check.idat.subdirs(base.dir)){
		## Found subdirectories that contain .idat files
		if(verbose) {
			rnb.info("Found subdirectores with idat files; creating symbolic links")
		}
		base.dir <- prepare.idat.dir(base.dir)
		if(is.null(base.dir)){
			rnb.error("Could not create temporary directory or create links to idat files")
		}
		fn <- attr(base.dir, "failed")
		if(fn != 0) {
			rnb.warning(paste("Could not create links to", fn, "idat files; these files will be ignored"))
		}
		rm(fn)
	}
	
	## Check whether all files a present
	
	if(!is.null(sample.sheet)){
		fn.base<-sample.sheet[,"barcode"]
	}else{
		fn.base<-barcodes
	}
	
	fn.expected<-paste(fn.base, c("Red.idat", "Grn.idat"), sep="_")
	fn.available<-list.files(base.dir, pattern="idat")
	fn.missing<-setdiff(fn.expected, fn.available)
	if(length(fn.missing)>0){
		rnb.error(sprintf("Some IDAT files are not present in the supplied base directory, for instance %s", 
						paste(fn.missing[1:min(6,length(fn.missing))], collapse=", ")))
	}
	
	## Read the data using the methylumi package
	methylumi.params <- list(idatPath = base.dir, n = TRUE)
	if(is.null(sample.sheet)){
		if(!is.null(load.chunk) && load.chunk>length(barcodes)){
			rnb.warning("The value of the load.chunk parameter is incorrect.")
			load.chunk=NULL
			nsamp<-length(barcodes)
		}
	}else{
		if(!is.null(load.chunk) && load.chunk>nrow(sample.sheet)){
			rnb.warning("The value of the load.chunk parameter is incorrect.")
			load.chunk=NULL
			nsamp<-nrow(sample.sheet)
		}
	}

	platform<-detect.platform(barcodes)
	
	annot<-rnb.annotation2data.frame(rnb.get.annotation(platform))
	annot.ctrls<-rnb.get.annotation(gsub("probes","controls",platform))
	nprobes<-sum(rnb.annotation.size(platform))
	ncprobes<-nrow(annot.ctrls)
	
	tIred<-which(annot$Color=="Red")
	tIgrn<-which(annot$Color=="Grn")
	tII<-which(annot$Design=="II")
	
	if(useff){
		dpvals<-matrix(NA_real_, nrow=nrpobes, ncol=nsamp)
		M<-matrix(NA_real_, nrow=nrpobes, ncol=nsamp)
		U<-matrix(NA_real_, nrow=nrpobes, ncol=nsamp)
		M0<-matrix(NA_real_, nrow=nrpobes, ncol=nsamp)
		U0<-matrix(NA_real_, nrow=nrpobes, ncol=nsamp)
		beadsM<-matrix(NA_integer_, nrow=nrpobes, ncol=nsamp)
		beadsU<-matrix(NA_integer_, nrow=nrpobes, ncol=nsamp)
	}else{
		dpvals<-ff(NA_real_, dim=c(nprobes, nsamp), dimnames=list(NULL, NULL), vmode="double")
		M<-ff(NA_real_, dim=c(nprobes, nsamp), dimnames=list(NULL, NULL), vmode="double")
		U<-ff(NA_real_, dim=c(nprobes, nsamp), dimnames=list(NULL, NULL), vmode="double")
		M0<-ff(NA_real_, dim=c(nprobes, nsamp), dimnames=list(NULL, NULL), vmode="double")
		U0<-ff(NA_real_, dim=c(nprobes, nsamp), dimnames=list(NULL, NULL), vmode="double")
		beadsM<-ff(NA_integer_, dim=c(nprobes, nsamp), dimnames=list(NULL, NULL))
		beadsU<-ff(NA_integer_, dim=c(nprobes, nsamp), dimnames=list(NULL, NULL))
	}
	
	qc.int<-list()
	qc$Cy3<-matrix(NA_real_, nrow=ncprobes, ncol=nsamp)
	qc$Cy5<-matrix(NA_real_, nrow=ncprobes, ncol=nsamp)
			
			
	INTENSITY.SUMMARIZATION.INFO<-list(
			"typeIred"=list(Design="I", Color="Red", Msource="Red", Usource="Red"),
			"typeIgrn"=list(Design="I", Color="Grn", Msource="Grn", Usource="Grn"),
			"typeII"=list(Design="II", Color="Both", Msource="Grn", Usource="Red")
	)
	
	probe.categories<-INTENSITY.SUMMARIZATION.INFO
	
	get.OOB.channel<-function(channel){
		if(channel=="Red"){
			return("Grn")
		}else if(channel=="Grn"){
			return("Red")
		}else{
			stop("Wrong input for channel")
		}
	}
	
	for(pc in probe.categories){
		pc$Indices<-which(annot$Design==pc$Design, annot$Color==pc$Color)
		pc$Maddress<-annot[pc$Indices,"AddressB"]
		pc$Uaddress<-annot[pc$Indices,"AddressA"]
	}
		
	for(sid in 1:nsamp){
		
		barcode<-sample.sheet$barcodes[sid]
		idatfile.red<-file.path(base.dir, sprintf("%s_Red.idat",barcode))
		idatfile.grn<-file.path(base.dir, sprintf("%s_Red.idat",barcode))
		
		if(!file.exists(idatfile.red)){
			rnb.error(sprintf("IDAT file not found: %s", idatfile.full))
		}
		if(!file.exists(idatfile.grn)){
			rnb.error(sprintf("IDAT file not found: %s", idatfile.grn))
		}
		
		int.files[["Red"]]<-readIDAT(idatfile.red)
		int.files[["Grn"]]<-readIDAT(idatfile.grn)
		
		qc$Cy5<-int.files[["Red"]]@Quant[match(),]
		qc$Cy3<-int.files[["Grn"]]
				
		for(pc in probe.categories){
					
			Mmap<-match(pc$Maddress, int.files[[pc$Msource]]$MidBlock)
			Umap<-match(pc$Uaddress, int.files[[pc$Usource]]$MidBlock)
						
			M[pc$Indices,sid]<-int.files[[pc$Msource]]$Quant[Mmap,1]
			U[pc$Indices,sid]<-int.files[[pc$Usource]]$Quant[Umap,1]
			M0[pc$Indices,sid]<-int.files[[get.OOB.channel(pc$Msource)]]$Quant[Mmap,1]
			U0[pc$Indices,sid]<-int.files[[get.OOB.channel(pc$Usource)]]$Quant[Umap,1]
			
			beadsM[pc$Indices,sid]<-int.files[[pc$Msource]]$Quant[Mmap,2]
			beadsU[pc$Indices,sid]<-int.files[[pc$Msource]]$Quant[Umap,2]
			
			mpval<-M[pc$Indices,sid]
			upval<-U[pc$Indices,sid]
			
			dpval[pc$Indices,sid]<-rowMax(cbind(
							mpval,
							upval
			))
		}				
	}
	
	
	if(is.null(load.chunk)){
		if(is.null(sample.sheet)){
			methylumi.params$barcodes <- barcodes
		}else{ # is.data.frame(sample.sheet)
			methylumi.params$pdat <- sample.sheet
		}
		#mls <- suppressWarnings(suppressMessages(do.call(methylumIDAT, methylumi.params)))
	
				
		#if(!keep.methylumi){
		#	mls<-as(mls, "RnBeadRawSet")
		#}
	}else{
		# load the data chunk-by-chunk, converting to RnBeadRawSet 
		# and combining them after every step	
		
		mls <- NULL
		if(is.null(sample.sheet)) {
			nsamp<-length(barcodes)
		}else{
			nsamp<-nrow(sample.sheet)
		}
		nch <- nsamp %/% load.chunk
		
		## TODO: remove this strange thing
		# temporarily disable region summarization
		rt<-rnb.getOption("region.types")
		rnb.options(region.types=character())
		####
		
		for(ch in 1:(nch+(nsamp %% load.chunk > 0))){
			if(verbose){
				rnb.info(sprintf("Loading sample chunk %d", ch))
			}
			
			ixxs <- (load.chunk*(ch-1L)+1L):min((load.chunk*ch), nsamp)
			
			if(is.null(sample.sheet)){
				methylumi.params$barcodes <- barcodes[ixxs]
			} else { # is.data.frame(sample.sheet)
				methylumi.params$pdat <- sample.sheet[ixxs,]
			}
			
			mlsch <- suppressWarnings(suppressMessages(do.call(methylumIDAT, methylumi.params)))
			
			rnbch <- as(mlsch, "RnBeadRawSet")
			
			if(is.null(mls)){
				mls<-rnbch
			}else{
				mls.old <- mls
				# a couple of disk space saving hacks
				if(rnb.getOption("enforce.destroy.disk.dumps")){
					if(mls.old@status$disk.dump){
						mls.old@status$discard.ff.matrices<-TRUE
					}
					if(rnbch@status$disk.dump){
						rnbch@status$discard.ff.matrices<-TRUE
					}
				}
				mls<-combine(mls.old, rnbch)
				if(rnb.getOption("enforce.destroy.disk.dumps")){
					if(mls@status$disk.dump){
						mls@status$discard.ff.matrices<-NULL
					}
				}
				rnb.call.destructor(mls.old)
				rnb.call.destructor(rnbch)
			}
			
			rm(rnbch)
			rnb.cleanMem()
		}
		#gc()
		
		rnb.options(region.types=rt)
		
		for (region.type in rnb.region.types.for.analysis("hg19")) {
			mls <- summarize.regions(mls, region.type)
		}
		rnb.cleanMem()
		
	}
	
	if(verbose){
		rnb.logger.completed()
	}
	
	return(mls)
	
}


########################################################################################################################

#' read.GS.report
#'
#' Reads in a Genome Studio report, exported as a single file.
#' 
#' @param gsReportFile  location of the GS report file
#' @param pd			alternative sample annotation, if the \code{gsReporFile} is missing the sample section as 
#' 						\code{data.frame} of \code{character} singleton with the file name
#' @param sep       	character used as field separator in the sample sheet file and in the GS report file
#' 						 (should be identical).
#' 						Default value is taken by the call to \code{rnb.getOption("import.table.separator")}
#' @param keep.methylumi a flag indicating whether the a \code{MethyLumiSet}
#'						 object should be returned instead of a \code{RnBeadRawSet}.
#' @param verbose Flag indicating ifthe messages to the logger should be sent. Note that the logger must be initialized prior to calling this function. 
#' 			     Logging is useful for keeping a record of the downloaded and processed samples. Also, informative messages are stored in case of an error.
#' 
#' @return  MethylumiSet object with the data from the report
#'
#' @export
read.GS.report<-function(
		gsReportFile, 
		pd=NULL,
		sep=rnb.getOption("import.table.separator"),
		keep.methylumi=FALSE,
		verbose=TRUE){
	
	if(verbose && logger.isinitialized()){
		rnb.error <- logger.error
	} else {
		rnb.error <- function(msg) { stop(paste(msg)) }
	}
	
	if(verbose){
		rnb.logger.start("Loading Data from a GenomeStudio Report File")
	}
	
	if(!(is.null(pd) || is.character(pd) || is.data.frame(pd))){
		rnb.error("Invalid value for sample descriptions")
	}
	
	if(is.character(pd)){
		pd<-read.sample.annotation(pd, sep=sep)
	}
	
	methylumiSet<-methylumiR(gsReportFile, sampleDescriptions=pd, sep=sep)
	
	if(verbose){
		rnb.logger.completed()
	}
	
	if(keep.methylumi){
		return(methylumiSet)	
	}else{
		return(as(methylumiSet, "RnBeadRawSet"))
	}
	
}
########################################################################################################################

#' read.bed.files
#' 
#' Reads a reduced-representation/whole-genome bisulfite sequencing data set from a set of BED files 
#' 
#' @param base.dir 		Directory with 6+2 formatted BED files contatining processed methylation data
#' @param file.names    Optional non-empty \code{character} vector listing the names of the files that should be
#'                      loaded relative to \code{base.dir}. If supplied, this vector must not contain \code{NA} 
#' 						among its elements.
#' @param sample.sheet  Optional file name containing a table of sample annotation data, or the table itself in the form
#'                      of a \code{\link{data.frame}} or \code{matrix}. Only (and all) samples defined in this table
#'                      will be loaded. The table is expected to contain a column named \code{"barcode"} that lists the
#'                      samples' Sentrix  barcodes. If such a column is not present, this function searches for columns
#'                      \code{"Sentrix_ID"} and \code{"Sentrix_Position"} (or similar) that build a barcode.
#' @param file.names.col Column of the sample sheet which contains the file names (integer singleton). If \code{NA}
#' 						an attempt will be made to find a suiting column automatically. 
#' @param assembly		Genome assembly. Defaults to human (\code{"hg19"})
#' @param region.types	\code{character} vector storing the types of regions for which the methylation information is to
#'                      be summarized. The function \code{\link{rnb.region.types}} provides the list of all supported
#'                      regions. Setting this to \code{NULL} or an empty vector restricts the dataset to site
#'                      methylation only.
#' @param coord.shift.plus	The frame shift between the the CpG annotation (1-based) and the coordinates in the loaded BEDs. 
#' 						If BEDs have 0-based coordinates, \code{coord.shift.plus=1} (default).
#' @param skip.lines	The number of top lines to skip while reading the BED files 
#' @param sep.samples	\code{character} singleton used as field separator in the sample sheet file. 
#' 						Default value is taken by the call to \code{rnb.getOption("import.table.separator")}
#' @param merge.bed.files In case multiple BED files are specified for each sample, the flag indicates whether the 
#' 						methylation calls should be merged after reading
#' @param useff			If \code{TRUE}, functionality provided by the \code{ff} package will be used to read the data efficiently.
#' @param verbose       Flag indicating ifthe messages to the logger should be sent. Note that the logger
#'                      must be initialized prior to calling this function. Logging is useful for keeping a
#'                      record of the downloaded and processed samples. Also, informative messages are
#'                      stored in case of an error.
#' @param ...			Further arguments which are passed to the internal function \code{read.bissnp.bed} and to \code{read.table}
#' 
#' @details 			Each loaded BED file should follow the 6+2 convention: six standard BED columns (chromosome, start, end, feature.name,
#' 						feature value, strand) plus two additional columns -- mean methylation level of a methylation site and the number of reads
#' 						covering it.
#' @author Pavlo Lutsik
#' @export 
read.bed.files<-function(base.dir=NULL, 
		file.names = NULL, 
		sample.sheet = NULL, 
		file.names.col=0, 
		assembly=rnb.getOption("assembly"), 
		region.types=rnb.region.types.for.analysis(rnb.getOption("assembly")),
		coord.shift.plus=1L, 
		skip.lines=1, 
		sep.samples=rnb.getOption("import.table.separator"),
		merge.bed.files=TRUE,
		useff=rnb.getOption("disk.dump.big.matrices"),
		verbose = TRUE,
		...){

	if(!parameter.is.flag(verbose)) {
		stop("Invalid value for verbose; expected TRUE or FALSE")
	}
	if(!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}

	if(verbose){
		rnb.logger.start("Loading Data From BED Files")
	}
	
	if(all(is.null(base.dir), is.null(file.names), is.null(sample.sheet))){
		rnb.error("one of base.dir, file.names, sample.sheet should be supplied")
	}

	if(!is.null(base.dir) && (!is.character(base.dir) || length(base.dir)!=1)){
		rnb.error("Invalide value for parameter base.dir")
	}
	
	if(!is.null(file.names)) {
		if(!is.character(file.names) && length(file.names) != 0 && (!any(is.na(file.names)))) {
			rnb.error("Invalid value for file.names; expected non-empty character vector or NULL")
		}
	}

	if(!(assembly %in% rnb.get.assemblies())) {
		rnb.error("unsupported assembly")
	}
	if(!(is.character(region.types) && (!any(is.na(region.types))))) {
		rnb.error("Invalid value for region.types")
	}
	
	if(is.null(file.names) && is.null(sample.sheet)) {
		## Attempt to locate sample annotation file
		sample.sheet <- list.files(base.dir, pattern = "[Ss]ample", full.names = TRUE)
		sample.sheet <- sample.sheet[file.info(sample.sheet)[, "isdir"] == FALSE]
		if(length(sample.sheet) > 1) {
			rnb.error(c("No unique candidate for sample annotation found in", base.dir))
		}
		if(length(sample.sheet) == 0 || is.na(sample.sheet)) {
			rnb.error(c("No candidate for sample information found in", base.dir))
		}
	}

	if(!is.null(sample.sheet)){
#		if(verbose){
#			rnb.status(c("Supplied a sample annotation, parsing ..."))
#		}
#		
		if(is.character(sample.sheet) && (!is.matrix(sample.sheet))) {
			if(length(sample.sheet) != 1) {
				rnb.error("invalid value for sample.sheet")
			}
			sample.sheet <- read.sample.annotation(sample.sheet, sep=sep.samples)
		} else if(is.matrix(sample.sheet)) {
			sample.sheet <- as.data.frame(sample.sheet, check.names = FALSE, stringsAsFactors = FALSE)
		}
		if(!is.data.frame(sample.sheet)) {
			rnb.error("invalid value for sample.sheet")
		}
		
		if(!is.null(rnb.getOption("identifiers.column"))) 
			sample.names<-sample.sheet[,rnb.getOption("identifiers.column")] else
			sample.names<-sample.sheet[,1]
	}else{
		sample.names<-sub("\\.bed", "", file.names)
	}
	
	if(is.null(file.names)){

		if(is.na(file.names.col)){
			bed.col<-find.bed.column(sample.sheet)
			file.names.col<-bed.col[[1]]
			file.names<-sample.sheet[,bed.col[[1]]]
		}else if(file.names.col>0 && file.names.col<=ncol(sample.sheet)){
			file.names<-sample.sheet[,file.names.col]
		}else{
			file.names.col<-grep("BED|bed", colnames(sample.sheet))
			file.names<-sample.sheet[,file.names.col]
		}
		file.names<-as.character(file.names)
		### manage the multiple file names 
		
		fn.lists<-strsplit(file.names, split=";")
		if(any(sapply(fn.lists, length)>1)){
			multiple.files<-TRUE
			n.files<-sapply(fn.lists, length)
			new.sample.sheet<-do.call("rbind", 
					lapply(1:nrow(sample.sheet), function(index) 
								sample.sheet[rep(index, n.files[index]),]))
						
			if(merge.bed.files){
				new.sample.sheet$AutomaticUniqueSampleID<-do.call("c", lapply(1:length(n.files),function(index) rep(index, n.files[index])))
			}
			file.names<-unlist(fn.lists)
			sample.names<-sub("\\.bed", "", file.names)
			new.sample.sheet[,file.names.col]<-file.names
			sample.sheet<-new.sample.sheet
		}else{
			multiple.files<-FALSE
		}
		
		if(!is.null(base.dir)){
			present.files<-list.files(base.dir, pattern="*", full.names=FALSE)
			if(length(setdiff(file.names, present.files))){
				rnb.error(c("These BED files", setdiff(file.names, present.files),"are missing in", base.dir))
			}
			file.names<-file.path(base.dir, file.names)
		}else{
			existing.files<-sapply(file.names, file.exists )
			if(!all(existing.files)){
				rnb.error(c("Theses BED files do not exist:", file.names[!existing.files]))
			}
		}
	}
	
	#skip.lines <- 1
	data.matrices<-lapply(file.names, read.single.bed, context="cg", 
			...,
#			chr.col=1L,
#			start.col=2L,
#			end.col=NA,
#			c.col=4L,
#			t.col=5L,
#			strand.col=3L,
#			mean.meth.col=NA,
#			coverage.col=NA,
#			coord.shift=0L,
#			nrows=10000L,
			skip=skip.lines, ffread=useff)
	
	if(verbose){
		rnb.status(sprintf("Read %d BED files", length(data.matrices)))
	}

	found.chroms<-character()
	found.strands<-character()
	site.ints<-NULL
	site.indices<-list()
	
	chr.offset<-100e9
	strand.offset<-1e9
	
	for(dmi in 1:length(data.matrices)){
		if(useff){
			open(data.matrices[[dmi]])
		}
		chrs.dm<-as.character(unique(data.matrices[[dmi]][,1L]))
		found.chroms<-c(found.chroms, unique(chrs.dm[!(chrs.dm %in% found.chroms)]))
		chr.match<-match(data.matrices[[dmi]][,1L], found.chroms)
		
		strands.dm<-as.character(unique(data.matrices[[dmi]][,3L]))
		found.strands<-c(found.strands, strands.dm[!(strands.dm %in% found.strands)])
		strand.match<-match(data.matrices[[dmi]][,3L], found.strands)
		site.ints.dm<-chr.match*chr.offset + strand.match*strand.offset + data.matrices[[dmi]][,2L]
		if(is.null(site.ints)){
			site.ints<-site.ints.dm
			site.indices[[dmi]]<-1:nrow(data.matrices[[dmi]])
			next
		}else{
			site.ints<-c(site.ints, unique(site.ints.dm[!(site.ints.dm %in% site.ints)]))
			site.indices[[dmi]]<-match(site.ints.dm, site.ints)
		}
		if(useff){
			close(data.matrices[[dmi]])
		}
	}
	
	dm.subsample<-list()
	if(prod(length(site.ints), length(data.matrices))>.Machine$integer.max){
		sites.allowed <- as.integer(.Machine$integer.max/length(data.matrices))
		sample.site.inds <- sort(sample.int(length(site.ints),sites.allowed))
		msg<-c("Full dataset is too large to be supported by ff. --> downsampling to",sites.allowed,"( of",length(site.ints),") sites")
		rnb.warning(msg)
		all.sites<-site.ints[sample.site.inds]
		for(dmi in 1:length(data.matrices)){
			if(useff){
				open(data.matrices[[dmi]])
			}
			site.ints.dm<-site.ints[site.indices[[dmi]]]
			new.sites<-match(site.ints.dm, all.sites)
			dm.subsample[[dmi]]<-which(!is.na(new.sites))
			site.indices[[dmi]]<-new.sites[!is.na(new.sites)]
			if(useff){
				close(data.matrices[[dmi]])
			}
		}
	}else{
		all.sites<-site.ints
		for(dmi in 1:length(data.matrices)){
			dm.subsample[[dmi]]<-1:dim(data.matrices[[dmi]])[1]
		}
	}
	
	
	covg.present<-sapply(1:length(data.matrices), function(indx){
				if(useff){
					open(data.matrices[[indx]])
				}
				present<-dim(data.matrices[[indx]])[2]==5L
				if(useff){
					close(data.matrices[[indx]])
				}
				return(present)
				})
	
	if(!all(covg.present)){
		rnb.warning(c("Coverage information is not present for the following BED files: ", paste(file.names[which(!covg.present)], sep=", ")))
		rnb.warning(c("Discarded coverage information"))
	}
		
	sites<-matrix(nrow=length(all.sites), ncol=3, dimnames=list(NULL, c("chr", "coord", "strand")))
	if(!useff){
		meth<-matrix(nrow=length(all.sites), ncol=length(data.matrices), dimnames=list(NULL, sample.names))
		if(all(covg.present)){
			covg<-matrix(nrow=length(all.sites), ncol=length(data.matrices), dimnames=list(NULL, sample.names))
		}else{
			covg<-NULL
		}
	}else{
		#sites<-ff(NA_integer_, dim=c(length(all.sites), 3), dimnames=list(NULL, c("chr", "coord", "strand")))
		meth<-ff(NA, dim=c(length(all.sites), length(data.matrices)), dimnames=list(NULL, sample.names), vmode="double")
		if(all(covg.present)){
			covg<-ff(NA_integer_, dim=c(length(all.sites), length(data.matrices)), dimnames=list(NULL, sample.names))
		}else{
			covg<-NULL
		}
	}
	
	chroms<-names(rnb.get.chromosomes(assembly=assembly))
	dump<-sapply(1:length(data.matrices), function(indx){
				if(useff){
					open(data.matrices[[indx]])
				}
				dm<-data.matrices[[indx]]
				#dm.rows2<-intersect(all.sites, rownames(dm))
				#dm.rows<-match(all.sites, rownames(dm))
				#dm<-dm[dm.rows[!is.na(dm.rows)],]
				#rows<-which(all.sites %in% rownames(dm))
				rows<-site.indices[[indx]]
				rows.dm<-dm.subsample[[indx]]
				sites[rows,1]<<-as.integer(match(dm[rows.dm,1L], chroms))
				sites[rows,3]<<-as.integer(factor(dm[rows.dm,3L], levels=c("+","-","*")))
				sites[rows,2]<<-as.integer(dm[rows.dm,2L]) + c(coord.shift.plus, 0L, 0L)[sites[rows,3]]
				meth[rows,indx]<<-as.numeric(dm[rows.dm,4L])/100
				if(!is.null(covg)){
					covg[rows,indx]<<-as.integer(dm[rows.dm,5L])
				}
				if(useff){
					delete(dm)
				}
				return(NULL)		
			})
	rm(dump)
	#clean the memory of big objects
	rm(data.matrices)
	rm(site.indices)
	rm(site.ints)
	rm(dm.subsample)
	rnb.cleanMem()
	
	if(verbose){
		rnb.status(sprintf("Combined a data matrix with %d sites and %d samples", dim(sites)[1], dim(sites)[2]))
	}
	
	#remove sites with illegal chromosomes
	legal.sites <- !is.na(sites[,1])
	num.illegal <- sum(!legal.sites)
	
	sites <- sites[legal.sites,,drop=FALSE]
	#meth <- meth[legal.sites,,drop=FALSE]
	#covg <- covg[legal.sites,,drop=FALSE]

	if(verbose && num.illegal > 0) {
		rnb.info(c("Removed", num.illegal, "sites with unknown chromosomes"))
	}

	if(verbose){
		rnb.status("Processed all BED files")
	}
#	i<-1L
# 	more correct solution, based on the cytosine coordinates (takes too long)
#	annot<-endoapply(rnb.get.annotation(), function(ann){
#				ann[(1:(rnb.annotation.size()[i]/2))*2-1]<-flank(ann[(1:(rnb.annotation.size()[i]/2))*2-1], -1,start=FALSE) 
#				ann[(1:(rnb.annotation.size()[i]/2))*2]<-flank(ann[(1:(rnb.annotation.size()[i]/2))*2], -1)
#				i<<-i+1
#				ann
#			})

	guess.strand <- all(sites[,3L]==3)
	if(verbose && guess.strand) {
		rnb.info("Inferring strand information from annotation enabled")
	}
    # Chromosome ordering fix. 
    #
    #chrdata<-sapply(sort(unique(sites[,1L])), function(chr){
    chrdata<-lapply(sort(unique(sites[,1L])), function(chr){
	#dump<-sapply(unique(sites[,1L]), function(chr){
		chr.subset<-which(sites[,1L]==chr)
		# version with proper ordering
		#
		# chr.subset<-sites[,1L]==chr
		if(guess.strand){
			site.ranges <- GRanges(rep(chroms[chr], length(chr.subset)), 
							IRanges(start=sites[chr.subset,2], 
									width=rep(1,length(chr.subset))), 
									rep("*",length(chr.subset)))
			c.coords.ref <- GenomicRanges::flank(rnb.get.annotation(assembly=assembly)[[chr]],-1)

			match <- GenomicRanges::findOverlaps(c.coords.ref, site.ranges, type="start")
		} else {
			site.ranges<-GRanges(rep(chroms[chr], length(chr.subset)), 
							IRanges(start=sites[chr.subset,2L], 
									width=rep(1,length(chr.subset))), 
									c("+","-","*")[sites[chr.subset,3L]])
						
	#				match<-GenomicRanges::findOverlaps(site.ranges, rnb.get.annotation(assembly=assembly)[[chr]], type="start")
	#				sites[chr.subset[queryHits(match)],2]<<-subjectHits(match)
	#				sites[chr.subset[-queryHits(match)],2]<<-NA
	#				return(NULL)
			
			# version with proper ordering
			#			
			 match<-GenomicRanges::findOverlaps(rnb.get.annotation(assembly=assembly)[[chr]], site.ranges, type="start")
		 }
		 bed.h <- subjectHits(match)
		 ref.h <- queryHits(match)
		 multi.hit.ref <- duplicated(ref.h)
		 if(any(multi.hit.ref)){
		 	rnb.warning(paste("Multiple entries detected in the bed file that match to the same reference C. --> Picking the first one"))
		 	bed.h <- bed.h[!multi.hit.ref]
		 	ref.h <- ref.h[!multi.hit.ref]
		 }
		 multi.hit.bed <- duplicated(bed.h)
		 if(any(multi.hit.bed)){
		 	rnb.warning(paste("Some entries in the bed file match multiple reference Cs. --> Picking the first one"))
		 	bed.h <- bed.h[!multi.hit.bed]
		 	ref.h <- ref.h[!multi.hit.bed]
		 }

		 n.match <- length(match)
		 
		 #chrmeths<-meth[chr.subset,,drop=FALSE][bed.h,,drop=FALSE]
		 #chrcovgs<-covg[chr.subset,,drop=FALSE][bed.h,,drop=FALSE]
		 #indices<-chr.subset[bed.h]
		 chrsites<-cbind(rep(1L,n.match),rep(chr, n.match), ref.h, chr.subset[bed.h])
		 #return(list(chrsites,chrmeths,chrcovgs))
		 chrsites
	
	})
	
	rnb.cleanMem()
	#version with proper ordering
	#sites<-do.call("rbind", chrdata[1,])
	sites<-do.call("rbind", chrdata)
	#meth<-do.call("rbind", chrdata[2,])
    #covg<-do.call("rbind", chrdata[3,])
	valid<-rowSums(is.na(sites))==0
	
	if(verbose){
		rnb.status(c("Matched", sum(valid),"of", length(all.sites),"methylation sites to the annotation"))
	}
	
	if(!is.null(covg)){
		zero.covg.sites <- rowSums(covg[sites[,4L],], na.rm=TRUE)==0
		if(verbose && any(zero.covg.sites)) {
			rnb.status(c("Removed",length(intersect(which(zero.covg.sites), which(valid))),"of",
							sum(valid),"methylation sites because they were not covered in any sample"))
		}
		valid<-which(valid | zero.covg.sites)
		
		if(!useff){
			covg[is.na(covg)]<-0L
		}else{
			for(ci in 1:ncol(covg)){
				covg[is.na(covg[,ci]),ci]<-0L
			}
		}
	}
	
	sites<-sites[valid,]
	if(!useff){
		meth<-meth[which(legal.sites)[sites[,4L]],]
		if(!is.null(covg)){
			covg<-covg[which(legal.sites)[sites[,4L]],]
		}
	}else{
		meth.old<-meth
		meth<-ff(NA, dim=c(nrow(sites), ncol(meth.old)), dimnames=list(NULL, sample.names), vmode="double")
		meth[,]<-meth.old[which(legal.sites)[sites[,4L]],]
		rm(meth.old); rnb.cleanMem()
		if(!is.null(covg)){
			covg.old<-covg
			covg<-ff(NA_integer_, dim=c(nrow(sites), ncol(covg.old)), dimnames=list(NULL, sample.names))
			covg[,]<-covg.old[which(legal.sites)[sites[,4L]],]
			rm(covg.old); rnb.cleanMem()
		}
	}
	sites<-sites[,-4L]
	colnames(sites)<-c("context", "chr", "index")
	
	meth.sample.n.valid<-integer()
	for(ci in 1:ncol(meth)){
		meth.sample.n.valid[ci] <- sum(!is.na(meth[,ci]))
	}
	
	if(any(meth.sample.n.valid<1)){
		rnb.error(c("The following samples have no valid methylation values:",colnames(meth)[meth.sample.n.valid<1]))
		#stop("invalid sample methylation values")
	}

#	sites<-cbind(rep(1, nrow(sites)), sites)
#	sites<-sites[,2:4]
#	colnames(sites)<-c("context", "chr", "index")
	
	if(is.null(sample.sheet)){
		sample.sheet<-data.frame(Sample_Label=sample.names, rownames=sample.names)
	}
	
	#be memory efficient by cleaning up
	# rownames take a lot of space

#	sites <- sites[valid,,drop=FALSE]
#	meth  <- meth[valid,,drop=FALSE]
#	covg  <- covg[valid,,drop=FALSE]
	rnb.cleanMem()
#	if(!is.null(covg)){
#		zero.covg.sites <- rowSums(covg, na.rm=TRUE)==0
#		if(any(zero.covg.sites)){
#			valid <- !zero.covg.sites
#			sites <- sites[valid,,drop=FALSE]
#			meth  <- meth[valid,,drop=FALSE]
#			covg  <- covg[valid,,drop=FALSE]
#			rnb.cleanMem()
#
#			if(verbose) {
#				rnb.status(c("Removed",sum(zero.covg.sites),"of",
#								length(valid),"methylation sites because they were not covered in any sample"))
#			}
#		}
#	}

	if(verbose){
		rnb.logger.start("Creating RnBiseqSet object")
	}
	rnbs<-new("RnBiseqSet", 
			pheno=sample.sheet,
			sites=sites[,],
			meth.sites=meth,
			covg.sites=covg,
			region.types=region.types,
			assembly=assembly,
			verbose=verbose,
			useff=useff)
	
	if(multiple.files && merge.bed.files){
		if(verbose){
			rnb.status("Merging infromation from multiple BED files")
		}
		rnbs<-mergeSamples(rnbs,"AutomaticUniqueSampleID")	
		rnbs@pheno$AutomaticUniqueSampleID<-NULL		
	}
	
	if(verbose){
		rnb.logger.completed()
	}
	if(verbose){
		rnb.logger.completed()
	}
	
	return(rnbs)
}

########################################################################################################################

##	read.bissnp.bed
## 
##	reads a BED file with methylation information
## 
##	@param file											the input BED file
##  @param context										the prefix of the output rownames
##  @param chr.col, start.col,end.col,strand.col		the column with the chromosome information
##	@param mean.meth.col, c.col, t.col					columns with methylation information
##	@param coverage.col									column with coverage information
##  @param is.epp.style									flag for custom Broad Epigenome Pipeline (EPP) bed style: 
##                                                      chrom \t start \t end \t 'methylated_count/total_count' \t meth_score_scaled_0_1000 \t strand
##                                                      OVERWRITES ALL OTHER parameters except \code{file} (also neglects ...)
##  @param coord.shift 									An integer specifying the coordinate adjustment applied to the start and end coordinates.
## 
##  @return a \code{data.frame} object with DNA methylation and coverage information. The row names are formed by the following convension: 
##  		\code{context\\.read.delim(file,...)[,chr.col]\\.read.delim(file,...)[,start.col]\\.read.delim(file,...)[,strand.col]}.
##  @author Pavlo Lutsik
##
read.single.bed<-function(file, 
		context="cg", 
		chr.col=1L, 
		start.col=2L, 
		end.col=3L, 
		strand.col=6L, 
		mean.meth.col=7L, 
		coverage.col=8L, 
		c.col=NA, 
		t.col=NA, 
		is.epp.style=FALSE, 
		coord.shift=0L,
		ffread=FALSE,
		...){
	
	if(all(is.na(mean.meth.col), any(is.na(c(c.col, t.col)))) && !is.epp.style){
		stop("No methylation information columns supplied")
	}
	
	if((is.null(mean.meth.col) || is.na(mean.meth.col)) && all(is.na(coverage.col), any(is.na(c(c.col, t.col)))) && !is.epp.style){
		stop("If the mean methylation column is absent,
			either coverage column and one of C- or T-count column, or both C- and T-count columns should be specified.")
	}
	
	if(any(!is.integer(c(chr.col, start.col, coverage.col, coord.shift, c(end.col, mean.meth.col, c.col,t.col)[!is.na(c(end.col, mean.meth.col, c.col,t.col))])))){
		stop("Invalid parameter values supplied, column indices should be integer")
	}
	
	if(min(na.omit(c(chr.col, start.col, end.col, strand.col, mean.meth.col, coverage.col, c.col, t.col)))<1){
		stop("Invalid parameter values supplied, column indices should not be negative or zero")
	}
	
	if(!is.character(file) || length(file)!=1) 
		stop("Invalid file name")
	
	## read top of the file to determine the column classes
	if(ffread){
		mc<-list(...)
		meth.tab.init <- read.table(file, skip=mc$skip, nrows=1000, header=FALSE, sep="\t", stringsAsFactors=FALSE)
		columnClasses <- get.bed.column.classes(meth.tab.init, 
				chr.col=chr.col, 
				start.col=start.col, 
				end.col=end.col, 
				strand.col=strand.col, 
				mean.meth.col=mean.meth.col, 
				coverage.col=coverage.col, 
				c.col=c.col, 
				t.col=t.col,
				useff=TRUE)
	}
	
	if(is.epp.style){
		
		if(!ffread){
			meth.tab <- read.table(file,...,header=FALSE, sep="\t",stringsAsFactors=FALSE)
			mm <- strsplit(meth.tab[,4],split="/")
			mm <- t(sapply(mm,as.numeric))
			data.set <- data.frame(chrom=meth.tab[,1],start=meth.tab[,2]+1,strand=meth.tab[,6],meth=mm[,1]/mm[,2]*100,cov=mm[,2]) #+1 to account for 0-based coordinates
		}else{
			csf<-getOption("stringsAsFactors")
			options(stringsAsFactors=TRUE)
			data.set<-read.delim.ffdf(x=NULL, file=file, next.rows=10000L, ..., colClasses=columnClasses, header=FALSE, sep="\t")
			options(stringsAsFactors=csf)
			mm <- strsplit(gsub("\\'", "", data.set[,4]),split="/")
			mm <- t(sapply(mm,as.numeric))
			data.set<-data.set[c(1,2,6)]
			colnames(data.set)<-c("chrom", "start", "strand")
			data.set$start<-ff(data.set[,2]+1)
			data.set$meth<-ff(mm[,1]/mm[,2]*100)
			data.set$cov<-ff(mm[,2])
			rm(mm);rnb.cleanMem()
		}
		#rownames(data.set)<-paste(rep(context, dim(data.set)[1]), data.set[,1], data.set[,2], data.set[,3], sep=".")
		if(ffread){	
			close(data.set)
		}
		return(data.set)
	} else {
		
		if(!ffread){
			data.set<-read.delim(file, ..., stringsAsFactors=FALSE, header=FALSE)
		}else{
			csf<-getOption("stringsAsFactors")
			options(stringsAsFactors=TRUE)
			data.set<-read.delim.ffdf(x=NULL,
					file=file, next.rows=10000L, ..., colClasses=columnClasses, header=FALSE)
			options(stringsAsFactors=csf)
		}
		
		#checking the validity
		
		if(max(chr.col, start.col, end.col, strand.col, mean.meth.col, coverage.col, c.col, t.col, na.rm=TRUE)>ncol(data.set)){
			stop(sprintf("The loaded bed file %s has less columns than expected, check the correctness of column assignment (or the BED style option)", file))
		}
		
		#TODO: strand information
		if(is.na(strand.col)){
			if(!ffread){
				data.set$strand <- factor(rep("*",nrow(data.set)))
			}else{
				data.set$strand <- as.ff(factor(rep("*",nrow(data.set))))
			}
			strand.col <- which(colnames(data.set)=="strand")
		} # else {
		# 	possible.strands <- unique(data.set[,strand.col])
		# 	valid.strand.info <- all(possible.strands %in% c("+","-"))
		# 	if(!valid.strand.info){
		# 		stop("Only + and - are allowed in the strand column")
		# 	}
		# }
		data.set[,strand.col]<-sub("\\.", "*", data.set[,strand.col])

		if(coord.shift!=0L){
			data.set[,start.col] <- data.set[,start.col] + coord.shift
			if(!is.na(end.col)){
				data.set[,start.col] <- data.set[,start.col] + coord.shift
			}
		}
#		rn<-paste(rep(context, dim(data.set)[1]), data.set[,chr.col], data.set[,start.col], data.set[,strand.col], sep=".")
#		if(!ffread){
#			rownames(data.set)<-rn
#		}else{
#			row.names(data.set)<-rn
#		}
#		rm(rn)
		rnb.cleanMem()
		
		if(!is.null(mean.meth.col) && !is.na(mean.meth.col)){
			
			columns<-c(chr.col, start.col, strand.col, mean.meth.col)
			if(!is.na(coverage.col)){
				columns<-c(columns, coverage.col)
			}
			if(!ffread){
				return(data.set[,columns])
			}else{
				res<-data.set[columns]
				close(res); close(data.set)
				return(res)
			}
			
		}else if((is.na(coverage.col) && all(!is.na(c(c.col,t.col)))) || (!is.na(coverage.col) && any(!is.na(c(c.col, t.col))))){
			
			if(is.na(coverage.col)){
				cvg <- data.set[,t.col]+data.set[,c.col]
			}else{
				cvg <- data.set[,coverage.col]
			}
			if(!is.na(c.col)){
				if(!ffread){
					return(cbind(data.set[,c(chr.col,start.col,strand.col)],data.set[,c.col]/cvg*100, cvg))
				}else{
					res<-data.set[c(chr.col,start.col,strand.col,c.col,c.col)]
					res[,4]<-data.set[,c.col]/cvg*100
					res[,5]<-cvg
					close(res); close(data.set)
					return(res)
				}
			}else{
				if(!ffread){
					return(cbind(data.set[,c(chr.col,start.col,strand.col)],(1-data.set[,t.col]/cvg)*100, cvg))
				}else{
					res<-data.set[c(chr.col,start.col,strand.col,t.col,t.col)]
					res[,4]<-(1-data.set[,t.col]/cvg)*100
					res[,5]<-cvg
					close(res); close(data.set)
					return(res)
				}
			}
		}else{
			stop("Could not process the BED file: no column with mean methylation specified, and no way to calculate it")
		}
	}
}

########################################################################################################################

##
##	find.bed.column
##
##	Automatically parses the provided bisulfite sequencing annotation file
##	and tries to guess the column with file names
##
##	@author Pavlo Lutsik
##
find.bed.column<-function(annotation, 
		verbose=TRUE){
	
	
	if(verbose){
		rnb.logger.start("Automatically parsing the provided sample annotation file")
	}
	
	# find a column with file names
	
	classes<-sapply(annotation, class)
	
	potential.fnames<-which(classes%in%c("character","factor"))
	
	candidate.fname<-sapply(potential.fnames, function(cix){
				
				extensions<-sapply(strsplit(as.character(annotation[,cix]),"\\."), 
						function(spl) spl[length(spl)] %in% c("bed", "BED", "cov", "COV")
				)
				
				all(extensions)
			})
	
	if(!any(candidate.fname)){
		msg<-"Could not identify a column with file names in the supplied sample annotations file"
		if(verbose){
			rnb.error(msg)
		}
	}else if(length(which(candidate.fname))==1){
		msg<-sprintf("Potential file names found in column %s of the supplied annotation table", 
				potential.fnames[which(candidate.fname)]) 
		
		if(verbose){
			rnb.status(msg)
		}
	}else{
		msg<-paste("Multiple columns containing file names detected (", potential.fnames[which(candidate.fname)],
				"). Please disable automatic parsing and specify the column containing the file names explicitly") 
		if(verbose){
			rnb.error(msg)
		}
	}
	
	absolute<-sapply(potential.fnames[candidate.fname], function(cix){
				paths<-sapply(annotation[,cix], grepl, pattern="/|\\\\" )
				all(paths)
			})
	
	if(verbose){
		rnb.logger.completed()
	}
	
	return(list(potential.fnames[candidate.fname], absolute))
}
########################################################################################################################
get.bed.column.classes<-function(bed.top,
		chr.col=1L, 
		start.col=2L, 
		end.col=3L, 
		strand.col=6L, 
		mean.meth.col=7L, 
		coverage.col=8L, 
		c.col=NA, 
		t.col=NA,
		useff=FALSE){
	
	classes<-sapply(bed.top, class)
	if(useff){
		classes[grep("logical", classes)]<-"numeric"
		classes[grep("character", classes)]<-"factor"
		classes[grep("integer", classes)]<-"numeric"
	}
	if(useff){
		classes[chr.col]<-"factor"
		classes[strand.col]<-"factor"
	}else{
		classes[chr.col]<-"character"
		classes[strand.col]<-"character"
	}
	classes[start.col]<-"integer"
	classes[end.col]<-"integer"
	classes[coverage.col]<-"integer"
	classes[mean.meth.col]<-"numeric"
	if(useff){
		classes[chr.col]<-"factor"
	}else{
		classes[chr.col]<-"character"
	}
	if(!is.na(c.col)){
		classes[c.col]<-"integer"
	}	
	if(!is.na(t.col)){
		classes[t.col]<-"integer"
	}
	names(classes)<-NULL
	return(classes)
}
########################################################################################################################
detect.infinium.platform<-function(barcodes){
	
	inf27k.idats.present<-any(grepl("_[ABCDEFGHIJKL]", barcodes))
	inf450k.idats.present<-any(grepl("_R0[1-6]C0[1-2]", barcodes))
	
	if(inf27k.idats.present && inf450k.idats.present){
		rnb.error("Undefined platform: both 450k and 27k IDATs are present")
	}
	
	if(inf27k.idats.present){
		return("probes27")	
	}
	
	if(inf450k.idats.present){
		return("probes450")
	}
}