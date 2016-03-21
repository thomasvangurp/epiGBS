########################################################################################################################
## dataExport.R
## created: 2013-01-15
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Data exporting routines
########################################################################################################################

## rnb.export.fail
##
## Adds a paragraph describing why exporting data for a specific region type might fail.
##
## @param report Report to contain the description.
## @param txt Text to be prepended to the paragraph.
## @author Yassen Assenov
rnb.export.fail <- function(report, txt = character()) {
	txt <- c(txt, "There are several reasons why a certain output file cannot be (fully) generated. Examples include:")
	rnb.add.paragraph(report, txt)
	txt <- list(
		"The corresponding region type is invalid.",
		"The corresponding region type is not supported by the dataset.",
		"Due to security restrictions, the creation of files in the output directory is not allowed.",
		"A file or directory with the same name exists and cannot be overwritten.",
		"The disk is full or the user quota is exceeded.")
	rnb.add.list(report, txt)
}

########################################################################################################################

#' rnb.RnBSet.to.GRangesList
#'
#' convert an \code{\linkS4class{RnBSet}} object to a \code{GRangesList} object
#' @param rnb.set Object of class \code{\linkS4class{RnBSet}} 
#' @param reg.type region type to be converted
#' @param return.regular.list flag indicating whether a regular \code{list} object should be returned instead
#'             of a \code{GRangesList}. Might improve performance in some cases
#' @return a \code{GRangesList} or \code{list} object with one list element (\code{GRanges}) for each sample in \code{rnb.set}
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' result <- rnb.RnBSet.to.GRangesList(rnb.set.example)
#' }
rnb.RnBSet.to.GRangesList <- function(rnb.set,reg.type="sites",return.regular.list=FALSE){
	if (!inherits(rnb.set, "RnBSet")) {
		stop("Invalid value for rnb.set: Expected RnBSet")
	}
	if (!(reg.type %in% c(rnb.region.types(assembly(rnb.set)),"sites"))){
		stop("Unsupported region type")
	}
	# logger.status("Gathering data") #DEBUG MESSAGE
	anno.sites <- annotation(rnb.set,type=reg.type)
	anno.sites <- anno.sites[,c("Chromosome","Start","End","Strand")]
	granges.coords <- data.frame2GRanges(anno.sites, ids = NULL, chrom.column = "Chromosome", start.column = "Start",
						end.column = "End", strand.column = "Strand", assembly = assembly(rnb.set), sort.result=FALSE)
	# logger.status("...getting data into the correct order") #DEBUG MESSAGE
	oo <- order(as.integer(seqnames(granges.coords)),start(granges.coords), end(granges.coords), as.integer(strand(granges.coords)))
	granges.coords <- granges.coords[oo]
	mm <- meth(rnb.set,type=reg.type)[oo,]
	ccov <- covg(rnb.set,type=reg.type)[oo,]

	# logger.status("Creating list for GRanges objects") #DEBUG MESSAGE
	rnbs.grl <- lapply(samples(rnb.set),FUN=function(ss){
		if (!is.null(ccov)){
			dd <- data.frame(score=mm[,ss],coverage=ccov[,ss])
		} else {
			dd <- data.frame(score=mm[,ss])
		}
		rows.no.nas <- !is.na(dd$score)
		granges.cur <- granges.coords[rows.no.nas]
		elementMetadata(granges.cur) <- dd[rows.no.nas,,drop=FALSE]
		# logger.status(c("...done processing sample",ss)) #DEBUG MESSAGE
		return(granges.cur)
	})
	names(rnbs.grl) <- samples(rnb.set)
	# logger.status("sorting entries")
	# rnbs.grl <- lapply(rnbs.grl , function(y) { y[order(as.integer(seqnames(y)),start(y), end(y), as.integer(strand(y))), ] })
	if (!return.regular.list){
		logger.status("converting to GRangesList") #DEBUG MESSAGE
		rnbs.grl <- GRangesList(rnbs.grl) #this is very memory intense for larger WGBS datasets
	}
	# #alternative: endoapply
	# rnbs.grl <- rep(GRangesList(granges.coords),length(samples(rnb.set)))
	# names(rnbs.grl) <- samples(rnb.set)
	# logger.status("Adding sample information") #DEBUG MESSAGE
	# i <- 1 
	# rnbs.grl <- endoapply(rnbs.grl , function(gr){
	# 	#a bit of a hack getting the sample name using i in endoapply.
	# 	#If you can think of a better way, let me know.
	# 	#endoapply was necessary for performance reasons in large WGBS datasets
	# 	ss <- names(rnbs.grl)[i] 
	# 	i <<- i + 1
	# 	logger.status(c("...processing sample",ss)) #DEBUG MESSAGE
	# 	if (!is.null(ccov)){
	# 		dd <- data.frame(score=mm[,ss],coverage=ccov[,ss])
	# 	} else {
	# 		dd <- data.frame(score=mm[,ss])
	# 	}
	# 	logger.status(c("......setting elementMetadata")) #DEBUG MESSAGE
	# 	elementMetadata(gr) <- dd
	# 	logger.status(c("......retrieving non-NAs")) #DEBUG MESSAGE
	# 	subset(gr,!is.na(elementMetadata(gr)[,"score"]))
	# })

	# rnbs.grl <- endoapply(rnbs.grl , function(y) { y[order(as.integer(seqnames(y)),start(y), end(y), as.integer(strand(y))), ] })
	return(rnbs.grl)
}

#' rnb.RnBSet.to.bed
#'
#' convert an \code{\linkS4class{RnBSet}}  object to \code{*.bed} files.
#' @param rnb.set Object of class \code{\linkS4class{RnBSet}}
#' @param out.dir output directory. If not existing, it will be created. otherwise files in that directory are overwritten.
#' @param reg.type region type to be converted
#' @param names.quant.meth should the names of the bed regions contain information on the methylation level. 
#' 		   If TRUE the following format is applied: meth_percent%[coverage]. Coverage is available only when
#' 			\code{covg(rnb.set)} is not NULL
#' @param verbose More detailed logger output
#' @return (invisibly) a summary list containing information on the conversion step.
#'         elements are \code{filenames} (a table containing information on which sample has been written to what filename)
#'         and \code{assembly} (a string indicating the assembly used by \code{rnb.set}).
#' @details Details on bed can be found in the \href{http://genome.ucsc.edu/FAQ/FAQformat.html}{UCSC Genome Browser
#'          documentation}.  Each methylation site is an entry in the resulting bed file. The Score column corresponds
#'          to a site's methylation value in the interval \code{[0,1]}.
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.RnBSet.to.bed(rnb.set.example,tempdir())
#' }
rnb.RnBSet.to.bed <- function(rnb.set,out.dir,reg.type="sites",names.quant.meth=TRUE,verbose=TRUE){
	if (!file.exists(out.dir)){
		dir.create(out.dir)
	}
	if (!(reg.type %in% c(rnb.region.types(assembly(rnb.set)),"sites"))){
		stop("Unsupported region type")
	}
	#convert RnBSet to GRangesList
	if (verbose) logger.status("Converting to GRangesList")
	rnb.set.grl <- rnb.RnBSet.to.GRangesList(rnb.set,reg.type=reg.type,return.regular.list=TRUE) #might take a while

	#output the text files
	for (i in 1:length(samples(rnb.set))){
		ss <- samples(rnb.set)[i]
		rr <- rnb.set.grl[[ss]]
		if (verbose) logger.status(c("Exporting sample",ss))
		bed.names <- 1:length(rr)
		if (names.quant.meth){
			bed.names <- paste("'",round(elementMetadata(rr)[,"score"]*100),"%",sep="")
			if ("coverage" %in% colnames(elementMetadata(rr))){
				bed.names <- paste(bed.names,"[",elementMetadata(rr)[,"coverage"],"]",sep="")
			}
			bed.names <- paste(bed.names,"'",sep="")
		}
		table.bed <- data.frame(chrom=as.character(seqnames(rr)),
								start=start(rr)-1,#-1 to adjust for 0-based format of bed files
								end=end(rr),
								name=bed.names,
								score=round(elementMetadata(rr)[,"score"]*1000), #scale score to [0,1000]
								strand=as.character(strand(rr)),stringsAsFactors=FALSE)
		table.bed[table.bed[,"strand"]=="*","strand"] <- "."
		write.table(format(table.bed,scientific=FALSE),file=file.path(out.dir,paste("rnbeads_sample",i,".bed",sep="")),
					quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,na=".")
	}
	res <- list()
	res[["filenames"]] <- data.frame(i=1:length(samples(rnb.set)),sample=samples(rnb.set),filename=paste("rnbeads_sample",1:length(samples(rnb.set)),".bed",sep=""))
	res[["assembly"]] <- assembly(rnb.set)
	res[["contains.overlapping.regions"]] <- isDisjoint(rnb.set.grl[[1]])
	invisible(res)
}

#' rnb.RnBSet.to.bedGraph
#'
#' convert an \code{\linkS4class{RnBSet}}  object to \code{*.bedGraph} files.
#' @param rnb.set Object of class \code{\linkS4class{RnBSet}} 
#' @param out.dir output directory. If not existing, it will be created. otherwise files in that directory are overwritten.
#' @param reg.type region type to be converted
#' @return (invisibly) a summary list containing information on the conversion step.
#'         elements are \code{filenames} (a table containing information on which sample has been written to what filename)
#'         and \code{assembly} (a string indicating the assembly used by \code{rnb.set}).
#' @details Details on bedGraph can be found \href{http://genome.ucsc.edu/goldenPath/help/bedgraph.html}{here}. 
#'          Each methylation site is an entry in the resulting bedGraph file. The Score column corresponds to a site's
#'          methylation value in the interval \code{[0,1]}.
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.RnBSet.to.bedGraph(rnb.set.example,tempdir())
#' }
rnb.RnBSet.to.bedGraph <- function(rnb.set,out.dir,reg.type="sites"){
	if (!suppressPackageStartupMessages(require(rtracklayer))) {
		stop("missing required package rtracklayer")
	}
	if (!file.exists(out.dir)){
		dir.create(out.dir)
	}
	if (!(reg.type %in% c(rnb.region.types(assembly(rnb.set)),"sites"))){
		stop("Unsupported region type")
	}
	#convert RnBSet to GRangesList
	rnb.set.grl <- rnb.RnBSet.to.GRangesList(rnb.set,reg.type=reg.type,return.regular.list=TRUE) #might take a while
	
	#output the text files
	for (i in 1:length(samples(rnb.set))){
		ss <- samples(rnb.set)[i]
		rtracklayer::export(rnb.set.grl[[ss]],file.path(out.dir,paste("rnbeads_sample",i,".bedGraph",sep="")),"bedGraph")
	}
	res <- list()
	res[["filenames"]] <- data.frame(i=1:length(samples(rnb.set)),sample=samples(rnb.set),filename=paste("rnbeads_sample",1:length(samples(rnb.set)),".bedGraph",sep=""))
	res[["assembly"]] <- assembly(rnb.set)
	res[["contains.overlapping.regions"]] <- !isDisjoint(rnb.set.grl[[1]])
	invisible(res)
}

## rnb.diffmeth.to.EpiExplorer.file
##
## Export the differential methylation information to a file that can be read by EpiExplorer
## @author Fabian Mueller
rnb.diffmeth.to.EpiExplorer.file <- function(rnb.set, diffmeth, fname, comp.name, reg.type, referenceSet) {
	if (!(is.character(fname) && length(fname) == 1 && (!is.na(fname)))) {
		stop("invalid value for fname")
	}
	if (!(reg.type %in% c(rnb.region.types(assembly(rnb.set)),"sites")) || !(reg.type %in% c(names(rnb.set@regions),"sites"))){
		stop("Unsupported region type")
	}
	get.type <- reg.type
	if (reg.type=="sites") {
		get.type <- "CpG"
		dmt <- get.table(diffmeth,comp.name,"sites",return.data.frame=TRUE)
		col.rank <- "combinedRank"
		col.diff <- "mean.diff"
		col.quot <- "mean.quot.log2"
		col.p.val <- "diffmeth.p.adj.fdr"
		index.map <- rnb.set@sites
	} else {
		dmt <- get.table(diffmeth,comp.name,reg.type,return.data.frame=TRUE)
		col.rank <- "combinedRank"
		col.diff <- "mean.mean.diff"
		col.quot <- "mean.mean.quot.log2"
		col.p.val <- "comb.p.adj.fdr"
		index.map <- rnb.set@regions[[reg.type]]
	}
	aa.rnbd <- rnb.get.annotation(get.type,assembly(rnb.set))
	
	#go over each chromosome and get the corresponding diffmeth data as GRanges object
	dm.grl <- lapply(1:length(aa.rnbd),FUN=function(cci){
		res <- aa.rnbd[[cci]]
		cur.chrom.members <- index.map[,2]==cci #get all the indexes belonging to the current chromosome
		inds <- index.map[cur.chrom.members,3] 
		dm.df <- dmt[cur.chrom.members,c(col.rank,col.diff,col.quot,col.p.val)]
		dm.df <- data.frame(inRnBSet=rep(TRUE,nrow(dm.df)),dm.df)
		na.vec <- rep(NA,length(res))
		empty.df <- data.frame(inRnBSet=rep(FALSE,length(res)),col.rank=na.vec,col.diff=na.vec,col.quot=na.vec,col.p.val=na.vec)
		names(empty.df) <- c("inRnBSet",col.rank,col.diff,col.quot,col.p.val)
		if (length(inds > 0)){
			empty.df[inds,] <- dm.df
		}
		elementMetadata(res) <- empty.df
		return(res)
	})
	dm.gr <- do.call(c,dm.grl)
	df2export <- as.data.frame(elementMetadata(dm.gr))
	
	write.table(format(df2export,scientific=FALSE),file=fname,
			quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE,na="")
	return(fname)
}

## create.ucsc.track.hub
##
## creates the directory structure and essential files for a UCSC track hub.
## @param hub.dir directory where the track hub should be created. Must be nonexisting.
## @param rnb.set Object of class \code{RnBSet} for which the hub is to be created.
## @param data.type either "bigBed" or "bigWig"
## @param ana.name a name for the analysis
## @param ana.desc a decription of the analysis
## @param email email address of the creator
## @return nothing of particular interest
## @details For details on UCSC track hubs see \href{http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html}{here}.
## @author Fabian Mueller
create.ucsc.track.hub <- function(hub.dir,rnb.set,data.type="bigBed",ana.name="RnBeads Analysis",ana.desc="Methylation data processed by RnBeads",email="-@-.com"){
	if (file.exists(hub.dir)){
		stop("Could not create UCSC Hub directory as the specified directory already exists")
	}
	if (!(data.type %in% c("bigBed","bigWig"))){
		stop("Unsupported data type")
	}
	hub.name <- "RnBeads Analysis Hub"
	#follow the steps from http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html
	#step1 is handled by rnb.export.to.ucsc and manual conversion to binary format
	#Step2
	dir.create(hub.dir)
	#Step3 has to be done manually
	#Step4
	fileConn<-file(file.path(hub.dir,"hub.txt"))
	writeLines(c(paste("hub",hub.name),
				 paste("shortLabel",ana.name),
				 paste("longLabel",ana.desc),
				 "genomesFile genomes.txt",
				 paste("email",email)
	),fileConn)
	close(fileConn)
	#Step5
	fileConn<-file(file.path(hub.dir,"genomes.txt"))
	writeLines(c(paste("genome",assembly(rnb.set)),
				 paste("trackDb ",assembly(rnb.set),"/trackDb.txt",sep="")
	),fileConn)
	close(fileConn)
	#Step6
	assembly.subdir <- file.path(hub.dir,assembly(rnb.set))
	if (!file.exists(assembly.subdir)){
		dir.create(assembly.subdir)
	}
	#Step7
	type.str <- data.type
	if (data.type == "bigBed"){
		type.str <- "bigBed 5 ."
	}
	lines2write <- c()
	for (i in 1:length(samples(rnb.set))){
		ss <- samples(rnb.set)[i]
		if (data.type=="bigBed"){
			lines2write <- c(lines2write,
				c(paste("track",paste("RnBeadsMethylationTrack",i,sep="")),
						paste("bigDataUrl",paste("rnbeads_sample",i,".",data.type,sep="")),
						paste("shortLabel",paste(ss)),
						paste("longLabel",paste("RnBeads methylation track for sample",ss)),
						paste("visibility","dense"),
						paste("color","0,60,120"), 
						paste("spectrum","on"),
						paste("type ",type.str,"\n",sep=""))
			)
		} else if (data.type=="bigWig"){
			lines2write <- c(lines2write,
				     c(paste("track",paste("RnBeadsMethylationTrack",i,sep="")),
					 paste("bigDataUrl",paste("rnbeads_sample",i,".",data.type,sep="")),
					 paste("shortLabel",paste(ss)),
					 paste("longLabel",paste("RnBeads methylation track for sample",ss)),
					 paste("visibility","full"),
					 paste("viewLimits","0:1"), #bigWig
					 paste("viewLimitsMax","0:1"), #bigWig
					 paste("maxHeightPixels","100:30:10"), #bigWig
					 paste("color","0,60,120"), 
					 paste("spectrum","on"),
	 				 paste("type ",type.str,"\n",sep=""))
			)
		}
	}
	fileConn<-file(file.path(hub.dir,assembly(rnb.set),"trackDb.txt"))
	writeLines(lines2write,fileConn)
	close(fileConn)
	#Step8
	#Maybe later
}

## rnb.convert.bedGraph.to.bigWig
##
## converts \code{*.bedGraph} files to \code{*.bigWig} files
## @param bedGraphFilenames filenames of the \code{*.bedGraph} files (input files).
## @param bigWig filenames of the \code{*.bigWig} files (output files).
## @param assembly assembly to be used. Important for determining chromosome sizes
## @details Supported operating systems are currently Unix and MacOS only. The corresponding
## 			binaries of \code{bedGraphToBigWig} are installed along with \pkg{RnBeads}. So are
##          the chromosome sizes files.
## @author Fabian Mueller
rnb.convert.bedGraph.to.bigWig <- function(bedGraphFilenames,bigWigFilenames,assembly){
	if (length(bedGraphFilenames)!=length(bigWigFilenames) | length(bedGraphFilenames) < 1){
		stop("non-matching input and output filenames")
	}
	if (!Sys.info()["sysname"] %in% c("Linux","Darwin")){
		stop(paste("Unsupported operating system:",Sys.info()["sysname"]))
	}
	if (!assembly %in% rnb.get.assemblies()){
		stop("unsupported assembly")
	}
	get.f <- function(fname) {
		tryCatch(system.file(fname, package = "RnBeads", mustWork = TRUE),
				error = function(e) { stop(paste("Internal error in RnBeads: required file", fname, "not found")) })
	}
	bedGraphTobigWig.exe <- ""
	if (Sys.info()["sysname"] == "Linux"){
		bedGraphTobigWig.exe <- get.f("bin/linux_x86.64/bedGraphToBigWig")
	} else if (Sys.info()["sysname"] == "Darwin"){
		bedGraphTobigWig.exe <- get.f("bin/macOSX.i386/bedGraphToBigWig")
	}
	chromSizes.file <- get.f(paste("extdata/chromSizes/",assembly,".chrom.sizes",sep=""))
	
	n <- length(bedGraphFilenames)
	for (i in 1:n){
		system(paste(bedGraphTobigWig.exe,bedGraphFilenames[i],chromSizes.file,bigWigFilenames[i]),ignore.stdout=TRUE)
	}
}

## rnb.convert.bed.to.bigBed
##
## converts \code{*.bed} files to \code{*.bigBed} files
## @param bedGraphFilenames filenames of the \code{*.bed} files (input files).
## @param bigBed filenames of the \code{*.bigBed} files (output files).
## @param assembly assembly to be used. Important for determining chromosome sizes
## @details Supported operating systems are currently Unix and MacOS only. The corresponding
## 			binaries of \code{bedToBigBed} are installed along with \pkg{RnBeads}. So are
##          the chromosome sizes files.
## @author Fabian Mueller
rnb.convert.bed.to.bigBed <- function(bedFilenames,bigBedFilenames,assembly){
	if (length(bedFilenames)!=length(bigBedFilenames) | length(bedFilenames) < 1){
		stop("non-matching input and output filenames")
	}
	if (!Sys.info()["sysname"] %in% c("Linux","Darwin")){
		stop(paste("Unsupported operating system:",Sys.info()["sysname"]))
	}
	if (!assembly %in% rnb.get.assemblies()){
		stop("unsupported assembly")
	}
	get.f <- function(fname) {
		tryCatch(system.file(fname, package = "RnBeads", mustWork = TRUE),
				error = function(e) { stop(paste("Internal error in RnBeads: required file", fname, "not found")) })
	}
	bedToBigBed.exe <- ""
	if (Sys.info()["sysname"] == "Linux"){
		bedToBigBed.exe <- get.f("bin/linux_x86.64/bedToBigBed")
	} else if (Sys.info()["sysname"] == "Darwin"){
		bedToBigBed.exe <- get.f("bin/macOSX.i386/bedToBigBed")
	}
	chromSizes.file <- get.f(paste("extdata/chromSizes/",assembly,".chrom.sizes",sep=""))
	
	n <- length(bedFilenames)
	for (i in 1:n){
		system(paste(bedToBigBed.exe,bedFilenames[i],chromSizes.file,bigBedFilenames[i]),ignore.stdout=TRUE)
	}
}

#' rnb.export.to.trackhub
#'
#' convert an \code{\linkS4class{RnBSet}}  object to a UCSC-style track hub. 
#' @param rnb.set Object of class \code{\linkS4class{RnBSet}}
#' @param out.dir output directory. If not existing, it will be created. otherwise files in that directory are overwritten.
#' @param reg.type region type to be converted
#' @param data.type either "bigBed" or "bigWig"
#' @param ... parameters passed on to the track hub generating procedure
#' @return a list containing information on the export
#' @details During execution the RnBSet is converted to bed files. If the operating system is supported (currently Unix and MacOS only)
#'          these are automatically converted to bigBed files. If your operating system is not supported, you need to create them manually
#'          (see the \href{http://genome.ucsc.edu/FAQ/FAQformat.html}{UCSC Genome Browser documentation} for details).
#'          For details on UCSC track hubs see the
#'          \href{http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html}{UCSC tracks help page}.
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.export.to.trackhub(rnb.set.example,tempdir())
#' }
rnb.export.to.trackhub <- function(rnb.set,out.dir,reg.type="sites",data.type="bigBed",...){
	if (!inherits(rnb.set, "RnBSet")) {
		stop("Invalid value for rnb.set: Expected RnBSet")
	}
	if (!file.exists(out.dir)){
		dir.create(out.dir)
	}
	if (!(reg.type %in% c(rnb.region.types(assembly(rnb.set)),"sites"))){
		stop("Unsupported region type")
	}
	if (!(data.type %in% c("bigBed","bigWig"))){
		stop("Unsupported data type")
	}
	res <- list(assembly=assembly(rnb.set))
	if (data.type=="bigBed"){
		track.hub.dir <- file.path(out.dir,"trackHub_bigBed")
		bed.dir <- file.path(out.dir,"bed")
		
		#convert RnBSet to bed
		rnb.logger.start("Conversion to BED")
		bed.conv <- rnb.RnBSet.to.bed(rnb.set,bed.dir,reg.type=reg.type,names.quant.meth=TRUE)
		rnb.logger.completed()
		#convert to binary
		in.file.list <- file.path(bed.dir,bed.conv$filenames$filename)
		out.file.list <- file.path(track.hub.dir,assembly(rnb.set),paste("rnbeads_sample",bed.conv$filenames$i,".bigBed",sep=""))
		
		res[["filenames.bed"]] <- bed.conv$filenames
		
		rnb.logger.start("Creating Track Hub")
		#write the UCSC track hub structure
		create.ucsc.track.hub(track.hub.dir,rnb.set,data.type="bigBed",...)
		if (Sys.info()["sysname"] %in% c("Linux","Darwin")){
			rnb.convert.bed.to.bigBed(in.file.list,out.file.list,assembly=assembly(rnb.set))
			res[["converted.bigBed"]] <- TRUE
		}
		else {
			rnb.info("Skipped conversion of bed to bigBed as no conversion binary is currently available for your operating system")
			res[["converted.bigBed"]] <- FALSE
		}
		rnb.logger.completed()
	} else if (data.type=="bigWig"){
		track.hub.dir <- file.path(out.dir,"trackHub_bigWig")
		bedgraph.dir <- file.path(out.dir,"bedgraph")
		
		#convert RnBSet to bedGraph
		rnb.logger.start("Conversion to bedGraph")
		bedGraph.conv <- rnb.RnBSet.to.bedGraph(rnb.set,bedgraph.dir,reg.type=reg.type)
		rnb.logger.completed()
		#convert to binary
		in.file.list <- file.path(bedgraph.dir,bedGraph.conv$filenames$filename)
		out.file.list <- file.path(track.hub.dir,assembly(rnb.set),paste("rnbeads_sample",bedGraph.conv$filenames$i,".bigWig",sep=""))
		
		res[["contains.overlapping.regions"]] <- bedGraph.conv[["contains.overlapping.regions"]]
		res[["filenames.bedGraph"]] <- bedGraph.conv$filenames
		if (res[["contains.overlapping.regions"]]){
			rnb.warning(c("Skipping conversion to track hub because bedGraph files for region type",reg.type,
					"contain overlapping fragments"))
			res[["converted.bigWig"]] <- FALSE
			unlink(bedgraph.dir, recursive = TRUE, force = TRUE)
			return(res)
		}
		
		rnb.logger.start("Creating Track Hub")
		#write the UCSC track hub structure
		create.ucsc.track.hub(track.hub.dir,rnb.set,data.type="bigWig",...)
		if (Sys.info()["sysname"] %in% c("Linux","Darwin")){
			rnb.convert.bedGraph.to.bigWig(in.file.list,out.file.list,assembly=assembly(rnb.set))
			res[["converted.bigWig"]] <- TRUE
			#delete bedGraph files after successful conversion
			unlink(bedgraph.dir, recursive = TRUE, force = TRUE)
		}
		else {
			rnb.info("Skipped conversion of bedGraph to bigWig as no conversion binary is currently available for your operating system")
			res[["converted.bigWig"]] <- FALSE
		}
		rnb.logger.completed()
	}
	
	invisible(res)
}

#' rnb.export.to.ewasher
#'
#' Data exported to a format compatible with the FaST-LMM-EWASher tool for cell-mixture adjustment.
#' see \href{http://www.nature.com/nmeth/journal/v11/n3/full/nmeth.2815.html}{Zou, J., et al., Nature Methods, 2014}
#' for further details on the tool.
#' @param rnb.set Object of class \code{\linkS4class{RnBSet}} 
#' @param out.dir output directory. If not existing, it will be created and all exported files will be placed here.
#' 				  If existing, this functions results in an error.
#' @param reg.type region type to be exported
#' @param ... passed on to \code{get.comparison.info}
#' @return a list containing information on the export
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.export.to.ewasher(rnb.set.example,tempfile(pattern="forEwasher"))
#' }
rnb.export.to.ewasher <- function(rnb.set, out.dir, reg.type="sites", ...){
	#example: rnb.export.to.ewasher(rnb.set,out.dir="~/tmp/rnb4ewasher",pheno.cols=rnb.getOption("differential.comparison.columns"))
	res <- list()
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set; expected RnBSet")
	}
	if (file.exists(out.dir)){
		stop("Could not create output directory as the specified directory already exists")
	}
	dir.create(out.dir)

	data.dir <- file.path(out.dir,"data_file")
	dir.create(data.dir)
	pheno.dir <- file.path(out.dir,"phenotype_file")
	dir.create(pheno.dir)
	covar.dir <- file.path(out.dir,"covar_file")
	dir.create(covar.dir)
	map.dir <- file.path(out.dir,"map_file")
	dir.create(map.dir)
	res[["data.subdir"]] <- "data_file"
	res[["pheno.subdir"]] <- "phenotype_file"
	res[["covar.subdir"]] <- "covar_file"
	res[["map.subdir"]] <- "map_file"


	cmp.info <- get.comparison.info(rnb.set,adjust.celltype=FALSE,...)
	
	mm <- meth(rnb.set,type=reg.type)
	aa <- annotation(rnb.set,type=reg.type)
	site.ids <- rownames(aa)
	sample.ids <- samples(rnb.set)
	# check if the sample ids start with a letter. If not prepend "S_"
	# because EWASher has problems with sample names starting with non-letter characters
	illegal.sample.id <- !grepl("^[[:alpha:]]", sample.ids)
	sample.ids[illegal.sample.id] <- paste0("S_",sample.ids[illegal.sample.id])

	cmp.ids   <- paste0("cmp",1:length(cmp.info))
	cmp.desc <- sapply(cmp.info,FUN=function(x){x$comparison})
	tab.summary <- data.frame(
		comparison_id=cmp.ids,
		comparison_desc=cmp.desc,
		stringsAsFactors=FALSE)
	res[["comparison.table"]] <- tab.summary

	#write the methylation data
	logger.start("Exporting methylation data")
	tab.data <- data.frame(
		ID=site.ids,
		mm
	)
	colnames(tab.data) <- c("ID",sample.ids)
	fn <- file.path(data.dir,"methylation.tsv")
	write.table(tab.data, file=fn,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
	res[["data.files"]] <- rep("methylation.tsv",length(cmp.info))
	logger.completed()

	#write the location mapping data
	logger.start("Exporting genome mapping data")
	map.data <- data.frame(
		chrom=gsub("chr", "", aa$Chromosome),
		ID=site.ids,
		chromStart=aa$Start,
		chromEnd=aa$End-1
	)
	fn <- file.path(map.dir,"mapping.map")
	write.table(map.data, file=fn,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	res[["map.files"]] <- rep("mapping.map",length(cmp.info))
	logger.completed()

	res[["pheno.files"]] <- rep(NA,length(cmp.info))
	res[["covar.files"]] <- rep(NA,length(cmp.info))
	for (i in 1:length(cmp.info)){
		logger.start(c("Exporting data for comparison:",cmp.info[[i]]$comparison))
		#write the phenotype data
		logger.status("Exporting phenotype table...")
		grps <- rep(NA,length(sample.ids))
		grps[cmp.info[[i]]$group.inds$group1] <- 1
		grps[cmp.info[[i]]$group.inds$group2] <- 0
		tab.pheno <- data.frame(
			sample=sample.ids,
			pheno=grps
		)
		fnn <- paste0("pheno_",cmp.ids[i],".tsv")
		fn <- file.path(pheno.dir,fnn)
		write.table(tab.pheno, file=fn,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
		res[["pheno.files"]][i] <- fnn

		#covariates
		pairing.info <- NULL
		covariates <- cmp.info[[i]]$adjustment.table
		if (cmp.info[[i]]$paired){
			pairing.info <- rep(NA,length(sample.ids))
			#the sample indices are assumed to be in pairing order and the number of samples in each group are equal
			num.pairs <- length(cmp.info[[i]]$group.inds$group1)
			pairing.info[cmp.info[[i]]$group.inds$group1] <- paste0("pair",1:num.pairs)
			pairing.info[cmp.info[[i]]$group.inds$group2] <- paste0("pair",1:num.pairs)
			pairing.info <- as.factor(pairing.info)
			if (is.null(covariates)){
				covariates <- data.frame(matrix(NA,nrow=length(sample.ids),ncol=0))
			}
			covariates <- data.frame(covariates,pairs=pairing.info)
		}
		if (!is.null(covariates)){
			logger.status("Exporting covariate table...")
			opts.na.orig <- options("na.action")[[1]] #for converting the dataframe and keeping NA's temporarily change the NA action
			options(na.action='na.pass')
			covariates <- model.matrix(~0+.,data=covariates)
			options(na.action=opts.na.orig)
			tab.covar <- data.frame(
				intercept=1,
				sample=sample.ids,
				covariates
			)
			fnn <- paste0("covar_",cmp.ids[i],".tsv")
			fn <- file.path(covar.dir,fnn)
			write.table(tab.covar, file=fn,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
			res[["covar.files"]][i] <- fnn
		}
		logger.completed()
	}

	fn <- file.path(out.dir,"summary.tsv")
	write.table(tab.summary, file=fn,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
	return(res)
}


### rnb.section.export.ct.adj
###
### add information to the report for export
### @author Fabian Mueller
### @aliases rnb.section.export
### @param ewasher.obj Ewasher export results object. See \code{\link{rnb.export.to.ewasher}} for details.
### @param rnbSet RnBSet object
### @param report report object to which the content is added
### @return the updated report object
rnb.section.export.ct.adj <- function(ewasher.obj,rnbSet,report){
	res.ewasher <- ewasher.obj
	ewasher.dir <- file.path(rnb.get.directory(report, "data"),"sites","ewasher")
	sectionText <- c("Site-level methylation and differential methylation data has been exported to formats compatible with tools for ",
		"cell-mixture adjustment.")
	report <- rnb.add.section(report, 'Export to Other Cell Type Adjustment Tools', sectionText)
	sectionText <- c("In order to analyze the exported data with FaST-LMM-EWASher, first obtain either the ",
		"<a href=\"http://research.microsoft.com/en-us/downloads/472fe637-7cb9-47d4-a0df-37118760ccd1/\">R-version</a> ",
		"or the <a href=http://research.microsoft.com/en-us/downloads/8af2ab12-a205-4bbb-809c-a333ecaffa40/>Python-version</a> ",
		"of the tool and install it on your machine. ",
		"Note that the FaST-LMM-EWASher software is currently in development. Hence, installation (particularly on non-Windows-computers) ",
		"is a non-trivial task. We recommend to contact the FaST-LMM-EWASher developers in case of technical difficulties. ",
		"RnBeads exported files that can serve as an input to FaST-LMM-EWASher. ",
		"These input files can be found in the ",paste("<a href=\"",ewasher.dir,"\">data directory</a> ",sep=""),
		"accompanying this report. The following table contains the filenames for the input data for each differential methylation comparison:")
	report <- rnb.add.section(report, 'FaST-LMM-EWASher', sectionText, level=2)
	data.file.dir <- file.path(rnb.get.directory(report, "data"),"sites","ewasher",res.ewasher$data.subdir)
	data.file.links <- paste("<a href=\"",file.path(data.file.dir,res.ewasher$data.files),"\">",res.ewasher$data.files,"</a>",sep="")
	map.file.dir <- file.path(rnb.get.directory(report, "data"),"sites","ewasher",res.ewasher$map.subdir)
	map.file.links <- paste("<a href=\"",file.path(map.file.dir,res.ewasher$map.files),"\">",res.ewasher$map.files,"</a>",sep="")
	pheno.file.dir <- file.path(rnb.get.directory(report, "data"),"sites","ewasher",res.ewasher$pheno.subdir)
	pheno.file.links <- paste("<a href=\"",file.path(pheno.file.dir,res.ewasher$pheno.files),"\">",res.ewasher$pheno.files,"</a>",sep="")
	covar.file.dir <- file.path(rnb.get.directory(report, "data"),"sites","ewasher",res.ewasher$covar.subdir)
	covar.file.links <- ifelse(!is.na(res.ewasher$covar.files),
		paste("<a href=\"",file.path(covar.file.dir,res.ewasher$covar.files),"\">",res.ewasher$covar.files,"</a>",sep=""),
		"N/A")
	tt.ewasher <- data.frame(
		comparison=res.ewasher$comparison.table$comparison_desc,
		comparison_id=res.ewasher$comparison.table$comparison_id,
		data_file=data.file.links,
		pheno_file=pheno.file.links,
		map_file=map.file.links,
		covar_file=covar.file.links,
		stringsAsFactors=FALSE
	)
	rnb.add.table(report,tt.ewasher)

	sectionText <- c("To run the Python-version of FaST-LMM-EWASher with the output file use the following syntax for your shell or console ",
		"commands:")
	rnb.add.paragraph(report, sectionText)
	cmdSyntax <- paste0("<code>","python fastlmm-ewasher.py data_file pheno_file [-covar covar_file] [-map map_file]","</code>")
	rnb.add.paragraph(report,cmdSyntax)
	sectionText <- c("To run the R-version of FaST-LMM-EWASher with the output file use the following syntax for your shell or console ",
		"commands:")
	rnb.add.paragraph(report, sectionText)
	cmdSyntax <- paste0("<code>","Rscript FLE-script.r data_file pheno_file [-covar covar_file] [-map map_file]","</code>")
	rnb.add.paragraph(report,cmdSyntax)
	return(report)
}

### rnb.section.tnt
###
### add information to the report for export
### @author Fabian Mueller
### @aliases rnb.section.export
### @param res.exp Exporting results object. See \code{\link{rnb.execute.tnt}} for details.
### @param rnbSet RnBSet object
### @param report report object to which the content is added
### @return the updated report object
rnb.section.tnt <- function(res.exp,rnbSet,report){
	if (is.null(res.exp)){
		sectionText <- "Nothing was exported"
		report <- rnb.add.section(report, 'Nothing Exported', sectionText)
	} else {
		is.biseq <- "RnBiseqSet" %in% class(rnbSet)
		is.beads <- "RnBeadSet" %in% class(rnbSet)
		export.ucsc.types <- unique(unlist(lapply(res.exp,FUN=function(x){x$export.trackhub})))
		any.export.bed <- any(unlist(lapply(res.exp,FUN=function(x){x$export.bed})))
		export.ewasher <- is.element("sites",names(res.exp))
		if (export.ewasher) export.ewasher <- is.element("export.ewasher",names(res.exp$sites))
		sectionText <- "Methylation Data was exported to the following formats:<ul>\n"
		if (any.export.bed){
			sectionText <- paste(sectionText,"<li><a href=\"http://genome.ucsc.edu/FAQ/FAQformat.html#format1\">bed</a></li>")
		}
		if ("bigBed" %in% export.ucsc.types){
			sectionText <- paste(sectionText,"<li><a href=\"http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html\">Track hub</a>: bigBed</li>")
		}
		if ("bigWig" %in% export.ucsc.types){
			sectionText <- paste(sectionText,"<li><a href=\"http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html\">Track hub</a>: bigWig</li>")
		}
		if (export.ewasher){
			refText <- c("Zou, J., Lippert, C., Heckerman, D., Aryee, M., Listgarten, J. (2014). ",
				"Epigenome-wide association studies without the need for cell-type composition. ",
				"<i>Nature Methods.</i>, 11(3), 309-311")
			report <- rnb.add.reference(report, refText)
			sectionText <- paste(sectionText, "<li>Input for the FaST-LMM-EWASher tool for cell-mixture adjustment ",
				rnb.get.reference(report, refText), "</li>")
		}
		sectionText <- paste(sectionText,"</ul>")
		report <- rnb.add.section(report, 'Methylation Data', sectionText)
		summary.tab <- rnb.sample.summary.table(rnbSet)
		fname <- paste("sample_methylation_summary.csv",sep="")
		fname <- rnb.write.table(summary.tab,fname,fpath=rnb.get.directory(report, "data", absolute = TRUE),format="csv",gz=FALSE,row.names = FALSE,quote=FALSE)

		sectionText <- paste("This section contains summary metrics for each sample as well as links to exported files. The summary table below contains the following metrics. It is also available as ",
			"<a href=\"", rnb.get.directory(report, "data"), "/", fname,"\">csv file</a>.",sep="")
		column.description <- c("<ul>",
									"<li>sampleName: Name of the sample</li>",
									"<li>*_num (* can be 'sites' or a region type): Number of sites or regions with coverage in the sample</li>")
		if (is.biseq){
			column.description <- c(column.description,
									"<li>*_covgMean: Mean coverage of sites or regions in the sample</li>",
									"<li>*_covgMedian: Median coverage of sites or regions in the sample</li>",
									"<li>*_covgPerc25,75: 25 and 75 percentile of coverage of sites or regions in the sample</li>",
									"<li>*_numCovg5,10,30,60: Number of sites or regions with coverage greater or equal to 5,10,30,60</li>")
		}
		if (is.beads){
			if (!is.null(dpval(rnbSet))){
				column.description <- c(column.description,
										"<li>sites_numDPval5em2,1em2,1em3: Number of sites with a detection p-value smaller than 0.05,0.01,0.001</li>")
			}
		}
		column.description <- c(column.description,
									"<li>**_numSitesMean (** is any region type): Mean number of sites in a region</li>",
									"<li>**_numSitesMedian (** is any region type): Median number of sites in a region</li>",
									"<li>**_numSites2,5,10,20: Number of regions with at least 2,5,10,20 sites with valid methylation measurements</li>",
							   "</ul>")

		sectionText <- paste(sectionText,paste(column.description,collapse="\n"))
		if (any.export.bed){
			sectionText <- c(sectionText," Bed files for each sample contain the locations of methylation sites and regions the methylation level
						   (in the score column). The following table links to the files. The files can be quite large in size, so you might want to
							consider using the \"save as\" option instead of opening the links directly.")
			tt <- data.frame(sample=as.character(res.exp[[1]]$filenames.bed[,"sample"]),stringsAsFactors=FALSE)
			for (i in 1:length(res.exp)){
				rr <- names(res.exp)[i]
				ttt <- res.exp[[i]]$filenames.bed[,c("sample","filename")]
				ttt$filename <- paste("<a href=\"",file.path(rnb.get.directory(report, "data"),rr,"bed",ttt$filename),"\">","bed","</a>",sep="")
				tt <- data.frame(tt,ttt$filename,stringsAsFactors=FALSE)
			}
			colnames(tt)[2:ncol(tt)] <- paste0(names(res.exp),"_bedFile")
			rownames(tt) <- tt$sample
			bed.df <- data.frame(tt[summary.tab$sampleName,2:ncol(tt)],stringsAsFactors=FALSE)
			colnames(bed.df) <- colnames(tt)[2:ncol(tt)] 
			summary.tab <- data.frame(sampleName=summary.tab$sampleName,bed.df,summary.tab[,2:ncol(summary.tab)],stringsAsFactors=FALSE)
		}
		if ("bigBed" %in% export.ucsc.types){
			tt <- data.frame(sample=as.character(res.exp[[1]]$filenames.bed[,"sample"]),stringsAsFactors=FALSE)
			for (i in 1:length(res.exp)){
				rr <- names(res.exp)[i]
				ttt <- res.exp[[i]]$filenames.bed[,c("sample","filename")]
				ttt$filename <- gsub(".bed$",".bigBed",ttt$filename) #replace bed with bigBed
				ttt$filename <- paste("<a href=\"",file.path(rnb.get.directory(report, "data"),rr,"trackHub_bigBed",assembly(rnbSet),ttt$filename),"\">","bigBed","</a>",sep="")
				tt <- data.frame(tt,ttt$filename,stringsAsFactors=FALSE)
			}
			colnames(tt)[2:ncol(tt)] <- paste0(names(res.exp),"_bigBedFile")
			rownames(tt) <- tt$sample
			bigbed.df <- data.frame(tt[summary.tab$sampleName,2:ncol(tt)],stringsAsFactors=FALSE)
			colnames(bigbed.df) <- colnames(tt)[2:ncol(tt)] 
			summary.tab <- data.frame(sampleName=summary.tab$sampleName,bigbed.df,summary.tab[,2:ncol(summary.tab)],stringsAsFactors=FALSE)
		}

		report <- rnb.add.section(report, "Sample summary", sectionText, level=2)
		rnb.add.table(report,summary.tab)
		
		if (length(export.ucsc.types)>0){
			
			sectionText <- paste("Track Hub data was generated for export to various genome browsers. Note that you need a server that is capable of
							serving the tracks to the genome browser via URL. Below, instructions are provided to view the tracks in the UCSC genome browser. Of course the files can also be viewed in other browsers such as the
							<a href=\"http://www.ensembl.org/\">Ensembl genome browser</a>.
							<ol>\n")
			anything.unconverted <- FALSE
			if ("bigBed" %in% export.ucsc.types) {
				if (!res.exp[[1]]$converted.bigBed) {
					anything.unconverted <- TRUE
				}
			}
			if ("bigWig" %in% export.ucsc.types) {
				if (!res.exp[[1]]$converted.bigWig) {
					anything.unconverted <- TRUE
				}
			}
			if (anything.unconverted) {
				sectionText <- paste(sectionText,"<li>Convert the bed/bedGraph files contained in the bed/bedGraph directories (see table below) attached to this report.",
									 "You can use the UCSC tool <code>bedTobigBed</code>/<code>bedGraphTobigWig</code> for this task. More information (e.g. where to obtain the tool)
									 on how to convert can be found <a href=\"http://genome.ucsc.edu/goldenPath/help/bigBed.html\">here</a>.",
									 "Make sure the resulting bigBed/bigWig files are moved to the corresponding UCSC track hub directory (see table(s) below).</li>")
			}
			sectionText <- paste(sectionText,"<li>Make sure you are viewing this report via the internet and not locally (simply look at the URL in the address line
							of your browser). If this is not the case, copy the files belonging to this report
							to a directory which is accessible via the internet and enter the corresponding URL in your browser.
							If you don't want to copy the entire report, it suffices to copy the track hub directories (see table(s) below).</li>")
			sectionText <- paste(sectionText,"<li>Make sure the",res.exp[[1]][["assembly"]],"subdirectories",
								"in the track hub directory (for its location be referred to the table(s) below) contain all the required bigBed/bigWig files</li>")
			sectionText <- paste(sectionText,"<li>Go to the <a href=\"http://genome.ucsc.edu/cgi-bin/hgHubConnect\">UCSC hub connection website</a> 
							to the \"My Hubs\" tab and add the web-URL of the <code>hub.txt</code> file in the directory you just copied.
							Hint: You can use the copy link functionality of your browser (right mouse click) on the track hub txt column in the table(s) below.
							Don't forget to click on \"Load Selected Hubs\".</li>")
			sectionText <- paste(sectionText,paste("<li>Make sure you are looking at the correct genome (",res.exp[[1]][["assembly"]],").
							The loaded hubs should now appear below the browsing window where you can modify their display settings.</li>",sep=""))
			sectionText <- paste(sectionText,"</ol>\nMore information on UCSC track hubs can be found 
						   <a href=\"http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html\">here</a>. Below, you find table(s) containing the directory names
						   for the bed/bedGraph and track hub directories for all exported data, site and region types.\n")
			report <- rnb.add.section(report, 'Tracks Hubs', sectionText, level=2)
			
			rrs <- names(res.exp)
			if ("bigBed" %in% export.ucsc.types){
				track.hub.dirs <- file.path(rnb.get.directory(report, "data"),rrs,"trackHub_bigBed")
				track.hub.urls <- file.path(rnb.get.directory(report, "data"),rrs,"trackHub_bigBed","hub.txt")
				bed.dirs <- file.path(rnb.get.directory(report, "data"),rrs,"bed")
				track.hub.dirs.links <- paste("<a href=\"",track.hub.dirs,"\">",track.hub.dirs,"</a>",sep="")
				track.hub.urls.links <- paste("<a href=\"",track.hub.urls,"\">",track.hub.urls,"</a>",sep="")
				bed.dirs.links <- paste("<a href=\"",bed.dirs,"\">",bed.dirs,"</a>",sep="")
				tt.bb <- data.frame("bed directory"=bed.dirs.links,"track hub directory"=track.hub.dirs.links,"track hub txt"=track.hub.urls.links)
				rownames(tt.bb) <- rrs
				report <- rnb.add.section(report, "bigBed", NULL, level=3)
				rnb.add.table(report,tt.bb)
			}
			if ("bigWig" %in% export.ucsc.types){
				regions.contain.overlap <- unlist(lapply(res.exp,FUN=function(x){
					base::ifelse(is.null(x[["contains.overlapping.regions"]]),FALSE,x[["contains.overlapping.regions"]])
				}))
				
				if (sum(!regions.contain.overlap)>0){
					track.hub.dirs <- file.path(rnb.get.directory(report, "data"),rrs[!regions.contain.overlap],"trackHub_bigWig")
					track.hub.urls <- file.path(rnb.get.directory(report, "data"),rrs[!regions.contain.overlap],"trackHub_bigWig","hub.txt")
					track.hub.dirs.links <- paste("<a href=\"",track.hub.dirs,"\">",track.hub.dirs,"</a>",sep="")
					track.hub.urls.links <- paste("<a href=\"",track.hub.urls,"\">",track.hub.urls,"</a>",sep="")
					tt.bw <- data.frame("track hub directory"=track.hub.dirs.links,"track hub txt"=track.hub.urls.links)
					if (!res.exp[[1]]$converted.bigWig){
						bed.dirs <- file.path(rnb.get.directory(report, "data"),rrs,"bedgraph")
						bed.dirs.links <- paste("<a href=\"",bed.dirs,"\">",bed.dirs,"</a>",sep="")
						tt.bw <- data.frame("bedgraph directory"=bed.dirs.links,tt.bw)
					}
					rownames(tt.bw) <- rrs[!regions.contain.overlap]
					report <- rnb.add.section(report, "bigWig", NULL, level=3)
					rnb.add.table(report,tt.bw)
				}
			}
		}
	}
	return(report)
}
#report <- rnb.section.export(res,report)

#' rnb.execute.tnt
#'
#' export RnBSet to various output data formats
#' @author Fabian Mueller
#' @aliases rnb.execute.tnt
#' @param rnb.set \code{\linkS4class{RnBSet}} object
#' @param out.dir output directory.
#' @param exp.bed A character vector indicating which data types should be exported to UCSC. Possible values in the vector are \code{bigBed} and \code{bigWig}.
#' 				  If \code{NULL}, UCSC export is disabled
#' @param exp.trackhub file types which should be exported to a trackhub structure.
#' @param region.types a character vector indicating region types to be exported
#' @param ... Arguments passed to \code{\link{rnb.export.to.trackhub}}
#' @return a list containing information on the export
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.execute.tnt(rnb.set.example,tempdir())
#' }
rnb.execute.tnt <- function(rnb.set,out.dir,exp.bed=rnb.getOption("export.to.bed"),exp.trackhub=rnb.getOption("export.to.trackhub"),
		                       region.types=rnb.getOption("export.types"),...){
	logger.start("Generating Tracks and Tables")
	reg.types <- intersect(region.types,c(rnb.region.types(assembly(rnb.set)),"sites"))
	if (length(reg.types)<1){
		stop("Region types not found")	
	}
	if (length(reg.types) < length(region.types)){
		logger.warning("Not all region types were found in the annotation")
	}
	ress <- list()
	for (rr in reg.types){
		logger.start(c("Exporting",rr))
		out.dir.cur <- file.path(out.dir,rr)
		dir.create(out.dir.cur)
		res <- list(export.bed=FALSE)
		if ("bigBed" %in% exp.trackhub){
			logger.start("Creating Track Hub -- bigBed")
			res.exp <- rnb.export.to.trackhub(rnb.set,out.dir.cur,reg.type=rr,
							ana.name=paste("RnBeads Analysis -- ",rr,",bigBed",sep=""),
							ana.desc=paste(rr, "methylation (bigBed)"),data.type="bigBed",...)
			logger.completed()
			res <- c(res,res.exp)
			res[["export.bed"]] <- TRUE
			res[["export.trackhub"]] <- c("bigBed")
		} else if (exp.bed){
			logger.start("Creating BED Files")
			bed.dir <- file.path(out.dir.cur,"bed")
			res.exp <- rnb.RnBSet.to.bed(rnb.set,bed.dir,reg.type=rr,names.quant.meth=TRUE)
			logger.completed()
			res[["filenames.bed"]] <- res.exp[["filenames"]]
			res[["export.bed"]] <- TRUE
			res[["export.trackhub"]] <- c()
		}
		if ("bigWig" %in% exp.trackhub){
			logger.start("Creating UCSC Track Hub -- bigWig")
			res.exp <- rnb.export.to.trackhub(rnb.set,out.dir.cur,reg.type=rr,
					ana.name=paste("RnBeads Analysis -- ",rr,",bigWig",sep=""),
					ana.desc=paste(rr, "methylation (bigWig)"),data.type="bigWig",...)
			logger.completed()
			res <- c(res,res.exp)
			if (!res[["contains.overlapping.regions"]]){
				res[["export.trackhub"]] <- c(res[["export.trackhub"]],"bigWig")
			}
			
		}
		ress <- c(ress,list(res))
		logger.completed()
	}
	names(ress) <- reg.types
	logger.completed()
	invisible(ress)
}

########################################################################################################################

## rnb.execute.export.csv.save
##
## Attempts to save the methylation values matrix of the given region type to the specified file.
##
## @param rnb.set         Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
## @param region.type     Site or region type to be exported.
## @param output.location Target file path.
## @param fname           Target file name. If this file already exists, it will be overwritten.
## @param gz              Flag indicating whether the file should be zipped in \code{gz} format.
## @return The (possibly updated) target file name, invisibly. If \code{gz} is \code{TRUE}, the string \code{".gz"} will
##         be appended to \code{fname}.
##
## @seealso \code{\link{rnb.write.table}}
## @author Yassen Assenov
rnb.execute.export.csv.save <- function(rnb.set, region.type, output.location, fname,
	gz = rnb.getOption("gz.large.files")) {
	mm <- meth(rnb.set, region.type, row.names = TRUE)
	reg.annotation <- annotation(rnb.set, region.type)
	accepted.columns <- tolower(colnames(reg.annotation)) %in% c("chrom", "chromosome", "start", "end", "strand")
	mm <- cbind(data.frame(ID = rownames(mm), check.names = FALSE, stringsAsFactors = FALSE),
		reg.annotation[, accepted.columns], mm)
	rnb.write.table(mm, fname, output.location, gz = gz, row.names = FALSE, quote=FALSE)
}

########################################################################################################################

#' rnb.execute.export.csv
#'
#' Exports (selected) methylation tables of the given dataset to comma-separated value files.
#' 
#' @param rnb.set         Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param output.location \code{character} or \code{\linkS4class{Report}} specifying the output directory. If this is a
#'                        report, the output directory is set to be a subdirectory named \code{csv} of the report's data
#'                        directory. Set this parameter to the empty string (\code{""}) or \code{NA} to use the current
#'                        working directory. If the given path does not exist, this function attempts to create it.
#' @param region.types    \code{character} vector indicating region types to be exported.
#' @return \code{character} vector containing the names of the files to which data were exported; prepended by
#'         \code{output.location}. In case a certain region type could not be exported (see the \emph{Details} section),
#'         the corresponding element of this vector is \code{NA}.
#'
#' @details
#' The names of the generated output files are formed by the prefix \code{"betas_"}, followed by a number between
#' 1 and \code{length(region.types)}. The extension is \code{.csv} or \code{.csv.gz}, depending on the value of the
#' \pkg{RnBeads} option \code{"gz.large.files"}. Any such files that already exist in the output directory, are
#' overwritten.
#' 
#' There are several reasons why a certain output file cannot be (fully) generated. Examples for failures are listed
#' below:
#' \itemize{
#'   \item{}{The corresponding region type is invalid.}
#'   \item{}{The corresponding region type is not supported by the dataset. If the type is loaded in \pkg{RnBeads},
#'           use the \code{\link[=summarize.regions,RnBSet-method]{summarize.regions}} method prior to calling this function,
#'           in order to include the support of this region type in the dataset.}
#'   \item{}{Due to security restrictions, the creation of files in the output directory is not allowed.}
#'   \item{}{A file or directory with the same name exists and cannot be overwritten.}
#'   \item{}{The disk is full or the user quota is exceeded.}
#' }
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.execute.export.csv(rnb.set.example, "", summarized.regions(rnb.set.example))
#' }
#' @author Yassen Assenov
#' @export
rnb.execute.export.csv <- function(rnb.set, output.location, region.types = rnb.getOption("export.types")) {
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (inherits(output.location, "Report")) {
		fprefix <- paste0(rnb.get.directory(output.location, "data", absolute = FALSE), "/csv/")
		output.location <- file.path(rnb.get.directory(output.location, "data", absolute = TRUE), "csv")
	} else {
		if (!(length(output.location) == 1 && (is.character(output.location) || is.na(output.location)))) {
			stop("invalid value for output.location; expected Report or character of length 1")
		}
		if (is.na(output.location) || output.location == ".") {
			output.location <- ""
		} else if (output.location != "") {
			output.location <- sub("^(.+)[/\\\\]$", "\\1", output.location) # remove trailing slash or backslash
		}
		fprefix <- ifelse(output.location == "", "", paste0(gsub("\\", "/", output.location, fixed = TRUE), "/"))
	}
	if (is.null(region.types)) {
		region.types <- character()
	} else if (!(is.character(region.types))) {
		stop("invalid value for region.types")
	}
	if (output.location != "") {
		if (!file.exists(output.location)) {
			## Attempt to create output.location
			if (!dir.create(output.location, showWarnings = FALSE, recursive = TRUE)) {
				rnb.error(c("Could not create directory", output.location))
			}
		} else if (!file.info(output.location)[1, "isdir"]) {
			rnb.error(c("Specified path is not a directory:", output.location))
		}
	}

	result <- region.types
	names(result) <- region.types
	if (length(region.types) != 0) {
		for (i in 1:length(region.types)) {
			fname <- paste0("betas_", i, ".csv")
			fname <- tryCatch(
				rnb.execute.export.csv.save(rnb.set, region.types[i], output.location, fname),
				error = function(e) { return(as.character(NA)) })
			if (!is.na(fname)) {
				fname <- paste0(fprefix, fname)
			}
			result[i] <- fname
		}
	}
	return(result)
}

########################################################################################################################

#' rnb.section.export.csv
#'
#' Creates a report section dedicated to exporting (selected) methylation tables to comma-separated value files.
#'
#' @param report         Report on data export to contain the CSV section. This must be an object of type
#'                       \code{\linkS4class{Report}}.
#' @param export.summary Result of the exporting procedure, as returned by \code{\link{rnb.execute.export.csv}}.
#' @return The (possibly) modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.export.csv <- function(report, export.summary) {
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (!(is.character(export.summary) && (!is.null(names(export.summary))))) {
		stop("invalid value for export.summary")
	}
	if (length(export.summary) == 0) {
		return(report)
	}

	txt <- "This section summarizes the result of exporting methylation data to comma-separated value (CSV) files."
	report <- rnb.add.section(report, "Export to Comma-separated Value Files", txt)

	if (all(is.na(export.summary))) {
		txt <- "Methylation matrices for none of the requested region types could be exported. "
		rnb.export.fail(report, txt)
		return(report)
	}
	txt <- c("The table below contains links to all generated CSV files that store exported methylation &beta; values.")
	if (any(is.na(export.summary))) {
		txt <- c(txt, " Note that the data for ", sum(is.na(export.summary)), " region types could not be exported.")
	}
	rnb.add.paragraph(report, txt)
	ftypes <- rep(".csv", length(export.summary))
	ftypes[grep("\\.csv\\.gz$", export.summary)] <- ".csv.gz"
	tsummary <- data.frame(
		"Region type" = names(export.summary),
		"File" = paste("<a href=\"", export.summary, "\">", ftypes, "</a>", sep = ""),
		check.names = FALSE, stringsAsFactors = FALSE)
	tsummary[is.na(export.summary), "File"] <- as.character(NA)
	rnb.add.table(report, tsummary, row.names = FALSE)
	if (any(is.na(export.summary))) {
		rnb.export.fail(report)
	}

	return(report)
}
