########################################################################################################################
## annotations.R
## created: 2012-08-16
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Collection of helper constants and functions related to the management of probe and region annotations.
########################################################################################################################
########################################################################################################################

#' rnb.fix.strand
#'
#' Converts, if necessary, the provided strand information into a factor vector with levels \code{"+"}, \code{"-"} and
#' \code{"*"}.
#'
#' @param values Vector of strand information to be converted.
#' @return The (possibly modified) strand information as a factor vector.
#'
#' @author Yassen Assenov
#' @noRd
rnb.fix.strand <- function(values) {
	if (is.factor(values) && setequal(levels(values), c("+", "-", "*"))) {
		values[is.na(values)] <- "*"
	} else {
		values <- as.character(values)
		values[is.na(values)] <- "*"
		i.positive <- values %in% c("+", "1", "+1", "F", "f", "TOP")
		i.negative <- values %in% c("-", "-1", "R", "r", "BOT")
		values[i.positive] <- "+"
		values[i.negative] <- "-"
		values[!(i.positive | i.negative)] <- "*"
		values <- factor(values, levels = c("+", "-", "*"))
	}
	return(values)
}

########################################################################################################################

#' rnb.sort.regions
#'
#' Sorts the given regions based on start and end position.
#'
#' @param x Genomic regions as an object of type \code{\link{GRanges}} or \code{\link{GRangesList}}.
#' @return Set of the same regions as x, sorted based on start and end positions.
#'
#' @author Yassen Assenov
#' @noRd
rnb.sort.regions <- function(x) {
	if (inherits(x, "GRanges")) {
		return(x[order(start(x), end(x), as.integer(strand(x))), ])
	}
	if (inherits(x, "GRangesList")) {
		return(endoapply(x, function(y) { y[order(start(y), end(y), as.integer(strand(y))), ] }))
	}
	stop("invalid value for x")
}

########################################################################################################################

#' rnb.load.bed
#'
#' Loads a BED file into a \code{data.frame} with fixed column names. The file contents is validated for structure
#' (at least 3 columns), as well as for integer values in columns 2 and 3.
#'
#' @param fname BED file to load.
#' @return \code{data.frame} with at least 3 and at most 6 columns. The column names are: \code{"chromosome"},
#'         \code{"start"}, \code{"end"}, \code{"id"}, \code{"score"} and \code{"strand"}. Columns after the sixth one,
#'         if present, are dropped.
#'
#' @author Yassen Assenov
#' @noRd
rnb.load.bed <- function(fname) {
	BED.COLUMNS <- c("chromosome", "start", "end", "id", "score", "strand")
	tbl <- tryCatch(suppressWarnings(read.delim(fname, header = FALSE, quote = "", comment.char = "#",
				stringsAsFactors = FALSE, na.strings = "")), error = function(e) { e })
	if (inherits(tbl, "error")) {
		if (grepl("cannot open", tbl, fixed = TRUE)) {
			stop("cannot open file")
		}
		stop("invalid file format")
	}
	if (ncol(tbl) < 3) {
		stop("invalid file format; expected at least 3 columns")
	}
	if (!(is.integer(tbl[[2]]) && is.integer(tbl[[3]]))) {
		stop("invalid file format; expected start and end positions in columns 2 and 3, respectively")
	}
	if (length(BED.COLUMNS) < ncol(tbl)) {
		tbl <- tbl[, 1:length(BED.COLUMNS)]
	}
	colnames(tbl) <- BED.COLUMNS[1:ncol(tbl)]
	tbl[[1]] <- as.character(tbl[[1]])
	if (4 <= ncol(tbl)) {
		tbl[[4]] <- as.character(tbl[[4]])
		if (anyDuplicated(tbl[[4]]) == 0) {
			rownames(tbl) <- tbl[[4]]
			tbl <- tbl[, -4]
		}
	}
	return(tbl)
}

########################################################################################################################

#' get.cpg.stats
#'
#' Computes CpG-related statistics for the specified regions.
#'
#' @param chrom.sequence Chromosome sequence, usually obtained from the assembly's genome definition. This must be an
#'                       object of type \code{MaskedDNAString}.
#' @param starts         \code{integer} vector of start positions for the regions of interest. 
#' @param ends           \code{integer} vector of end positions for the regions of interest.
#' @return Table of statistics for the regions in the form of a \code{matrix} with the following columns:
#'         \code{"CpG"} and \code{"GC"}. The columns contain the number of CpG dinucleoties and the number of C and G
#'         bases in each region.
#'
#' @author Yassen Assenov
#' @export
get.cpg.stats <- function(chrom.sequence, starts, ends) {
	#if (!inherits(chrom.sequence, "MaskedDNAString")) {
		#stop("invalid value for chrom.sequence")
	#}
	if (!(is.integer(starts) && !any(is.na(starts)))) {
		stop("invalid value for starts")
	}
	if (!(is.integer(ends) && !any(is.na(ends)))) {
		stop("invalid value for ends")
	}
	chrom.regions <- suppressWarnings(Views(chrom.sequence, start = starts, end = ends))
	cbind(
		"CpG" = dinucleotideFrequency(chrom.regions)[, "CG"],
		"GC" = as.integer(rowSums(letterFrequency(chrom.regions, c("C", "G")))))
}

########################################################################################################################

#' append.cpg.stats
#'
#' Appends additional metadata columns for CpG count and GC density to the specified regions.
#'
#' @param genome.data Genome of interest.
#' @param regionlist  Genomic regions as a list of \code{GRanges} objects (or an obect of type \code{GRangesList}),
#'                    containing one set of regions per chromosome.
#' @return The modified \code{regionlist}. Two columns are appedned to the metadata of each element in this list -
#'         \code{"CpG"} and \code{"GC"}. If the metadata already contains these columns, this function appends columns
#'         with similar names.
#'
#' @author Yassen Assenov
#'
append.cpg.stats <- function(genome.data, regionlist) {
	cpg.stats <- function(chrom) {
		stats <- get.cpg.stats(genome.data[[chrom]], start(regionlist[[chrom]]), end(regionlist[[chrom]]))
		result <- regionlist[[chrom]]
		mcols(result) <- IRanges::cbind(mcols(result), DataFrame(stats))
		result
	}
	regions.enriched <- foreach(chrom = names(regionlist), .packages = "GenomicRanges",
		.export = c("get.cpg.stats", "genome.data", "regionlist")) %dopar% cpg.stats(chrom)
	names(regions.enriched) <- names(regionlist)
	return(GRangesList(regions.enriched))
}

########################################################################################################################

#' data.frame2GRanges
#'
#' Converts a \code{data.frame} that defines genomic regions to object of type \code{GRanges}.
#'
#' @param dframe        Table defining genomic regions.
#' @param ids           Region names (identifiers) as a \code{character} vector, or \code{NULL} if no names are present.
#' @param chrom.column  Column name or index that lists the chromosome names.
#' @param start.column  Column name or index that lists the start positions of the regions.
#' @param end.column    Column name or index that lists the end positions of the regions.
#' @param strand.column Column name or index that lists the strands on which the regions are located. Set this to
#'                      \code{NULL} if this region set is not strand-specific.
#' @param assembly      Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported
#'                      genomes.
#' @param sort.result   Should the resulting table be sorted
#' @return \code{GRanges} object encapsulating all well defined regions on supported chromosomes, contained in
#'         \code{dframe}. Columns other that the ones listed as parameters in this function are included as metadata.
#'
#' @author Yassen Assenov
#' @export
data.frame2GRanges <- function(dframe, ids = rownames(dframe), chrom.column = "chromosome", start.column = "start",
                               end.column = "end", strand.column = NULL, assembly = "hg19", sort.result = TRUE) {
  if (is.character(chrom.column)) { chrom.column <- which(colnames(dframe) == chrom.column) }
  if (is.character(start.column)) { start.column <- which(colnames(dframe) == start.column) }
  if (is.character(end.column)) { end.column <- which(colnames(dframe) == end.column) }
  if (is.character(strand.column)) { strand.column <- which(colnames(dframe) == strand.column) }
  
  ## Process the chromosome names
  chroms <- as.character(dframe[, chrom.column])
  chroms <- paste0("chr", sub("^chr", "", chroms))
  param.list <- list()
  param.list[["seqnames"]] <- chroms
  
  ## Process the regions
  dframe <- dframe[!(is.na(dframe[, start.column]) | is.na(dframe[, end.column])), ]
  param.list[["ranges"]] <- IRanges(start = dframe[, start.column], end = dframe[, end.column], names = ids)
  
  ## Process the strands
  if (!is.null(strand.column)) {
    param.list[["strand"]] <- rnb.fix.strand(dframe[, strand.column])
  }
  
  for (cname in colnames(dframe)[-c(chrom.column, start.column, end.column, strand.column)]) {
    param.list[[cname]] <- dframe[[cname]]
  }
  if (!is.null(assembly)) {
    i.valid <- (chroms %in% names(CHROMOSOMES[[assembly]]))
    param.list <- lapply(param.list, function(x) { x[i.valid] })
  }
  param.list[["check.names"]] <- FALSE
  result <- do.call(GRanges, param.list)
  if (!is.null(assembly)) {
    seqlevels(result) <- names(CHROMOSOMES[[assembly]])
  }
  genome(result) <- assembly
  if (sort.result){
    result <- rnb.sort.regions(result)
  }
  return(result)
}

########################################################################################################################

#' load.annotations
#'
#' Loads the specifieid annotation tables and mappings, if they are not already loaded.
#'
#' @param assembly Genome assembly regions to load.
#' @param sites    Sites and mappings to load. This parameter is considered only if \code{assembly} is not \code{NULL}.
#' @return \code{TRUE} if the required assembly (and sites) were successfully loaded; \code{FALSE} if such annotations
#'         do not exist.
#'
#' @author Yassen Assenov
#' @noRd
load.annotations <- function(assembly = NULL, sites = NULL) {
	load.f <- function(fname, pname = paste0("RnBeads.", assembly)) {
		if (pname != "RnBeads") {
			suppressPackageStartupMessages(library(pname, character.only = TRUE))
		}
		tryCatch(load(system.file(fname, package = pname), envir = .rnb.annotations),
			error = function(e) {
				stop(paste("Internal error in RnBeads: loading required file", fname, "failed"))
			}
		)
	}
	if (length(ls(envir = .rnb.annotations)) == 0) {
		## Load the basic structure, showing supported assemblies, regions and sites
		load.f(file.path("data", "annotations.RData"), pname = "RnBeads")
		for (o.name in ls(envir = .rnb.annotations)) {
			genome.info <- get(o.name, envir = .rnb.annotations, inherits = FALSE)
			a.names <- c(names(genome.info$sites), names(genome.info$regions))
			lengths <- matrix(as.integer(NA), nrow = length(CHROMOSOMES[[o.name]]), ncol = length(a.names),
				dimnames = list(as.character(CHROMOSOMES[[o.name]]), a.names))
			.rnb.annotations[[o.name]][["lengths"]] <- lengths
		}
		rm(o.name, genome.info, a.names, lengths)
	}

	if (!is.null(assembly)) {
		updated <- c()
		if (!exists(assembly, where = .rnb.annotations)) {
			return(FALSE)
		}
		if (is.null(.rnb.annotations[[assembly]][["regions"]][[1]])) {
			## Load assembly regions
			load.f(file.path("data", paste(assembly, ".regions.RData", sep = "")))
			.rnb.annotations[[assembly]][["regions"]] <- get("regions", envir = .rnb.annotations, inherits = FALSE)
			rm("regions", envir = .rnb.annotations)
			updated <- names(.rnb.annotations[[assembly]][["regions"]])
		}
		if (!is.null(sites)) {
			if (sites %in% names(.rnb.annotations[[assembly]][["controls"]])) {
				c.sites <- attr(.rnb.annotations[[assembly]][["controls"]], "sites")
				sites <- c.sites[sites]
				rm(c.sites)
			}
			if (!(sites %in% names(.rnb.annotations[[assembly]][["sites"]]))) {
				return(FALSE)
			}
			if (is.null(.rnb.annotations[[assembly]][["sites"]][[sites]])) {
				## Load the required sites
				load.f(file.path("data", paste(assembly, ".", sites, ".RData", sep = "")))
				sinfo <- get("sites", envir = .rnb.annotations, inherits = FALSE)
				rm("sites", envir = .rnb.annotations)
				.rnb.annotations[[assembly]][["sites"]][[sites]] <- sinfo$sites
				for (cnames in setdiff(names(sinfo), c("sites", "mappings"))) {
					.rnb.annotations[[assembly]][["controls"]][[cnames]] <- sinfo[[cnames]]
				}
				## Assign and/or create mappings
				region.bi <- attr(.rnb.annotations[[assembly]][["regions"]], "builtin")
        if (!is.null(region.bi)){
				for (rname in names(.rnb.annotations[[assembly]][["regions"]])) {
					if (region.bi[rname]) {
						rmapping <- sinfo$mappings[[rname]]
					} else {
						rmapping <- rnb.regions2sites(.rnb.annotations[[assembly]][["regions"]][[rname]], sinfo$sites)
					}
					.rnb.annotations[[assembly]][["mappings"]][[rname]][[sites]] <- rmapping
					rm(rmapping)
				}
				suppressWarnings(rm(sinfo, region.bi, rname))
				updated <- sites
			}}
		}
		if (length(updated) != 0) {
			## Update annotations lengths for the newly loaded regions and/or sites
			genome.info <- get(assembly, envir = .rnb.annotations, inherits = FALSE)
			for (a.name in updated) {
				.rnb.annotations[[assembly]][["lengths"]][, a.name] <- 0L
				if (a.name %in% names(genome.info$regions)) {
					a.lengths <- elementLengths(genome.info$regions[[a.name]])
				} else { # a.name %in% names(genome.info$sites)
					a.lengths <- elementLengths(genome.info$sites[[a.name]])
				}
				.rnb.annotations[[assembly]][["lengths"]][names(a.lengths), a.name] <- a.lengths
			}
		}
	}

	return(TRUE)
}

########################################################################################################################

#' unload.annotations
#'
#' Removes all loaded annotation tables and mappings (if such exist).
#'
#' @author Yassen Assenov
#' @noRd
unload.annotations <- function() {
	rm(list = ls(envir = .rnb.annotations, all.names = TRUE), envir = .rnb.annotations)
}

########################################################################################################################

#' rnb.update.annotation.infos
#'
#' Updates or adds supporting attributes for a new region annotation.
#'
#' @param type     One-element \code{character} vector giving the name of the freshly added/updated region annotation.
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#'
#' @author Yassen Assenov
#' @noRd
rnb.update.annotation.infos <- function(type, assembly) {

	## Set the annotation for built-in regions
	region.builtin <- attr(.rnb.annotations[[assembly]][["regions"]], "builtin")
	if (type %in% names(region.builtin)) {
		region.builtin[type] <- FALSE
	} else {
		region.builtin <- c(region.builtin, FALSE)
		names(region.builtin)[length(region.builtin)] <- type
	}
	attr(.rnb.annotations[[assembly]][["regions"]], "builtin") <- region.builtin

	## Set the information for number of regions per chromosome
	a.lengths.full <- rep.int(0L, length(CHROMOSOMES[[assembly]]))
	names(a.lengths.full) <- names(CHROMOSOMES[[assembly]])
	a.lengths <- elementLengths(.rnb.annotations[[assembly]][["regions"]][[type]])
	a.lengths.full[names(a.lengths)] <- a.lengths
	i.type <- which(colnames(.rnb.annotations[[assembly]][["lengths"]]) == type)
	if (length(i.type) == 0) {
		lengths <- cbind(.rnb.annotations[[assembly]][["lengths"]], a.lengths.full)
		colnames(lengths)[ncol(lengths)] <- type
	} else {
		lengths <- .rnb.annotations[[assembly]][["lengths"]]
		lengths[, i.type] <- a.lengths.full
	}
	.rnb.annotations[[assembly]][["lengths"]] <- lengths
}

########################################################################################################################

#' rnb.annotation.size
#'
#' Gets the size, in number of genomic elements, of the specified annotation.
#'
#' @param type     Name of annotation. Control probe annotations are not accepted.
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#' @return \code{integer} vector showing the number of elements the specified annotation contains per chromosome. The
#'         names of the vector are the names of \code{\link{rnb.get.chromosomes}} for the given genome assembly.
#'         Chromosomes that are not covered by the annotation have their respective value set to \code{0} (zero).
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' rnb.annotation.size("probes450")
#' }
#' @seealso \code{\link{rnb.region.types}} for a list of supported region annotations
#'
#' @author Yassen Assenov
#' @export
rnb.annotation.size <- function(type = "CpG", assembly = "hg19") {
	if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
		stop("invalid value for type")
	}
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	load.annotations(assembly, type)
	if (!(type %in% colnames(.rnb.annotations[[assembly]][["lengths"]]))) {
		stop("unsupported annotation type")
	}
	.rnb.annotations[[assembly]][["lengths"]][, type]
}

########################################################################################################################

#' rnb.annotation2data.frame
#'
#' Transform the specified site, probe or region annotation to \code{data.frame}.
#'
#' @param annotation.table Annotation in the form of non-empty \code{GRangesList} object, as returned by
#'                         \code{\link{rnb.get.annotation}}.
#' @param add.names        Flag indicating if element names should be extracted and returned also as a column named
#'                         \code{"ID"} in the resulting \code{data.frame}. Note that element names, if present, are
#'                         set to be the row names of the table.
#' @return Annotation in the form of a single \code{data.frame}. The columns in this table include, among other,
#'         \code{"Chromosome"}, \code{"Start"} and \code{"End"}.
#'
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' head(rnb.annotation2data.frame(rnb.get.annotation("probes450")))
#' }
#' @author Yassen Assenov
#' @export
rnb.annotation2data.frame <- function(annotation.table, add.names = TRUE) {
	if (!(inherits(annotation.table, "GRangesList") && length(annotation.table) != 0)) {
		stop("invalid value for annotation.table")
	}
	if (!parameter.is.flag(add.names)) {
		stop("invalid value for add.names")
	}

	## Convert the GRangesList object to a data.frame
	el.names <- do.call("c", lapply(annotation.table, names))
	if (is.null(el.names) || anyDuplicated(el.names) != 0) {
		result <- GenomicRanges::as.data.frame(annotation.table, row.names = NULL)
		el.names <- NULL
	} else {
		result <- GenomicRanges::as.data.frame(annotation.table)
		rownames(result) <- el.names
	}
	if (ncol(result) > 2 && identical(colnames(result)[1:2], c("group", "group_name"))) {
		result <- result[, c(-1, -2), drop = FALSE]
	}

	## Adjust column and row names
	c.names <- colnames(mcols(annotation.table[[1]]))
	c.names <- c(colnames(result)[1:(ncol(result) - length(c.names))], c.names)
	c.names[c.names == "CGI.Relation"] <- "CGI Relation"
	c.names[c.names == "strand"] <- "Strand"
	c.names[c.names == "start"] <- "Start"
	c.names[c.names == "end"] <- "End"
	c.names[c.names == "seqnames"] <- "Chromosome"
	c.names[c.names == "element"] <- "ID"
	colnames(result) <- c.names
	if (add.names && (!is.null(el.names))) {
		result[["ID"]] <- el.names
		i <- which(c.names == "width")
	} else {
		i <- which(c.names %in% c("ID", "width"))
	}
	if (length(i) != 0) {
		result <- result[, -i]
	}
	return(result)
}

########################################################################################################################

#' rnb.get.assemblies
#'
#' Gets the supported genome assemblies.
#'
#' @return All supported genome assemblies in the form of a \code{character} vector. These are \code{"hg19"},
#'         \code{"mm10"}, \code{"mm9"} and \code{"rn5"}.
#'
#' @examples
#' \dontrun{
#' "hg19" %in% rnb.get.assemblies()
#' }
#' @author Yassen Assenov
#' @export
rnb.get.assemblies <- function() {
	load.annotations()
	return(ls(envir = .rnb.annotations))
}

########################################################################################################################

#' rnb.get.annotation
#'
#' Extracts the requested annotation for the given genome.
#' 
#' @param type     Name of annotation.
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#' @return Probe, site or region annotation table. If the specified type refers to control probes, the returned value is
#'         a \code{data.frame} listing all respective control probes. Otherwise, this function returns an object of type
#'         \code{\link{GRangesList}} - a list of consistent \code{\link{GRanges}} objects, one per chromosome.
#' 
#' @details
#' When the returned value is of type \code{GRangesList}, it defines the genomic positions of the requested sites,
#' probes or regions. Identifiers, if present, can be obtained using the \code{names} method. Strand information is also
#' included when applicable. Any additional annotation is stored as metadata in the respective \code{GRanges} objects.
#'
#' @seealso \code{\link{rnb.set.annotation}} for adding annotation;
#'   \code{\link{rnb.region.types}} for all loaded region types in a genome assembly
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' rnb.get.annotation("promoters")
#' }
rnb.get.annotation <- function(type = "CpG", assembly = "hg19") {
	if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
		stop("invalid value for type")
	}
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	load.annotations(assembly, type)
	if (type %in% names(.rnb.annotations[[assembly]][["controls"]])) {
		result <- .rnb.annotations[[assembly]][["controls"]][[type]]
	} else if (type %in% names(.rnb.annotations[[assembly]][["sites"]])) {
		result <- .rnb.annotations[[assembly]][["sites"]][[type]]
	} else if (type %in% names(.rnb.annotations[[assembly]][["regions"]])) {
		result <- .rnb.annotations[[assembly]][["regions"]][[type]]
	} else {
		stop("unsupported annotation type")
	}
	return(result)
}

########################################################################################################################

#' rnb.set.annotation
#'
#' Adds or replaces a region annotation table.
#'
#' @param type        One-element \code{character} vector giving the name of the annotation. If this region type is
#'                    already available, it will be overwritten for the current session. The type cannot be one of
#'                    \code{"CpG"}, \code{"probes450"} or \code{"controls450"}, because these names are reserved for the
#'                    annotation tables of CpG dinucleotides, and Infinium methylation and control probes, respectively.
#' @param regions     BED file defining regions (see \emph{Details}). Alternatively, the value of this parameter can be
#'                    a table of genomic regions in the form of a \code{\link{data.frame}}, containing at least the
#'                    following three columns - \code{"chromosome"}, \code{"start"} and \code{"end"} (notice the lower
#'                    case). The \code{"chromosome"} column must be a \code{character} or \code{factor} vector that
#'                    lists chromosome names. The \code{"start"} and \code{"end"} columns are expected to contain
#'                    genomic positions as \code{integer}s. The row names of this \code{data.frame} are used as region
#'                    identifiers.
#' @param description Optional; short description in the form of a non-empty \code{character} vector. The elements in
#'                    this vector are concatenated without a separator to form the description of the annotation.
#' @param assembly    Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported
#'                    genomes.
#'
#' @details
#' In case the parameter \code{regions} specifies an existing BED file, regions are loaded from this file. The number of
#' columns defined must be at least 3. Columns after the sixth one, if present, are dropped. The columns are given the
#' following names: \code{"chromosome"}, \code{"start"}, \code{"end"}, \code{"id"}, \code{"score"} and \code{"strand"}.
#'
#' The annotation tables in \pkg{RnBeads} focus on chromosomes \code{"chr1"}, \code{"chr2"}, ..., \code{"chr22"},
#' \code{"chrX"} and \code{"chrY"}. Regions on other chromosomes are ignored. This function also recognizes the
#' convention of chromosome names such as \code{"1"}, adopted, for example, by \href{http://www.ensembl.org/}{Ensembl}.
#' Apart from this, the region definition table is not examined in details by this function; therefore, regions located
#' on unsupported chromosomes or having invalid (e.g., negative) genomic coordinates are simply not mapped to any sites
#' or probes.
#'
#' @examples
#' \dontrun{
#' my.regions <- data.frame(
#'     chromosome = c("chr1", "chr1"),
#'     start = c(49242278L, 49242372L),
#'     end = c(49242590L, 49242810L),
#'     rownames = c("BEND5E1", "CpG:38"))
#' txt <- "First exon of the BEND5 gene and an overlapping CpG island."
#' rnb.set.annotation("my regions", my.regions, txt)
#' }
#' @seealso \code{\link{rnb.get.annotation}} for extracting annotation;
#'   \code{\link{rnb.region.types}} for all loaded region types in a genome assembly
#' @author Yassen Assenov
#' @export
rnb.set.annotation <- function(type, regions, description = NULL, assembly = "hg19") {
	if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
		stop("invalid value for type")
	}
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}
	if (!is.data.frame(regions)) {
		if (!(is.character(regions) && length(regions) == 1 && (!is.na(regions)))) {
			stop("invalid value for regions")
		}
		if (!file.exists(regions)) {
			stop("invalid value for regions; file does not exist")
		}
		regions <- rnb.load.bed(regions)
	}
	if (!(ncol(regions) >= 3 && all(c("chromosome", "start", "end") %in% colnames(regions)) &&
			(is.character(regions[["chromosome"]]) || is.factor(regions[["chromosome"]])) &&
			is.integer(regions[["start"]]) && is.integer(regions[["end"]]))) {
		stop("invalid format for regions")
	}
	if (!is.null(description)) {
		if (!(is.character(description) && length(description) != 0)) {
			stop("invalid value for description")
		}
		description <- paste(description, collapse = "")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	s.names <- c(names(.rnb.annotations[[assembly]][["controls"]]), names(.rnb.annotations[[assembly]][["sites"]]))
	if (type %in% s.names) {
		stop("invalid value for type; specified name is not allowed")
	}

	## Construct list of region annotation tables with sorted coordinates, one per chromosome
	strand.column <- which(colnames(regions) == "strand")
	if (length(strand.column) == 0) {
		strand.column <- NULL
	}
	regions <- data.frame2GRanges(regions, strand.column = strand.column, assembly = assembly)
	regions <- GenomicRanges::split(regions, seqnames(regions))
	attr(regions, "description") <- description

	## Create mappings from the regions to all defined sites
	sites <- .rnb.annotations[[assembly]][["sites"]]
	sites <- sites[!sapply(sites, is.null)]
	mappings <- lapply(sites, function(x) { rnb.regions2sites(regions, sites = x) })

	## Add the region annotation table to the environment
	.rnb.annotations[[assembly]][["regions"]][[type]] <- regions
	## Add mappings to the environment
	.rnb.annotations[[assembly]][["mappings"]][[type]] <- mappings
	## Set annotation attributes and sizes (lengths)
	rnb.update.annotation.infos(type, assembly)
}

########################################################################################################################

#' rnb.set.annotation.and.cpg.stats
#'
#' wrapper for \code{\link{rnb.set.annotation}} to accept the region format as output by \code{annotation(rnb.set)}.
#' Additionally, CpG statistics are added to the annotation.
#' @param regions a data.frame handled similarly as by \code{\link{rnb.set.annotation}} with the exception that
#'                 the genomic location columns should be specified using upper case first letters
#' @param type,description,assembly Parameters handled exactly as in \code{\link{rnb.set.annotation}}
#' @seealso \code{\link{rnb.set.annotation}}
#' @author Fabian Mueller
#' @export
rnb.set.annotation.and.cpg.stats <- function(type, regions, description = NULL, assembly = "hg19"){
	## TODO: Incorporate this as a parameter in rnb.set.annotation 
	genome.data <- get.genome.data(assembly)
	regs.gr <- data.frame2GRanges(regions, chrom.column = "Chromosome", start.column = "Start",
		end.column = "End", strand.column = "Strand", assembly = assembly)
	regs.grl <- GenomicRanges::split(regs.gr, seqnames(regs.gr))

	regs.grl <- append.cpg.stats(genome.data, regs.grl)
	regs.gr <- GenomicRanges::unlist(regs.grl,use.names=FALSE)

	regs.df <- GenomicRanges::as.data.frame(regs.gr)
	colnames(regs.df)[colnames(regs.df)=="seqnames"] <- "chromosome"
	regs.df$width <- NULL

	rnb.set.annotation(type=type, regions=regs.df, description=description, assembly = assembly)
}

########################################################################################################################

#' rnb.save.annotation
#'
#' Saves the specified region annotation table and its accompanying data structures to a binary file.
#'
#' @param fname    One-element \code{character} vector giving the name of the file to contain the annotation data. If
#'                 this file already exists, it will be overwritten.
#' @param type     One-element \code{character} vector giving the name of the region annotation.
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#'
#' @details
#' This function is used in combination with \code{\link{rnb.load.annotation}} to enable fast reloading of custom region
#' annotations. If can also be used to save a build-in region annotation (e.g. before overwriting it) but not site or
#' control probe annotations.
#'
#' @seealso \code{\link{rnb.load.annotation}} for loading a saved annotation
#' @author Yassen Assenov
#' @export
rnb.save.annotation <- function(fname, type, assembly = "hg19") {
	if (!(is.character(fname) && length(fname) == 1 && (!is.na(fname)))) {
		stop("invalid value for fname")
	}
	if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
		stop("invalid value for type")
	}
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	if (!(type %in% names(.rnb.annotations[[assembly]][["regions"]]))) {
		stop("unsupported annotation type")
	}

	regions <- .rnb.annotations[[assembly]][["regions"]][[type]]
	mappings <- .rnb.annotations[[assembly]][["mappings"]][[type]]
	tryCatch(save(assembly, regions, mappings, file = fname, compress = "xz"),
		error = function(er) { stop(paste("unable to save objects to", fname)) })
}

########################################################################################################################

#' rnb.export.annotation
#'
#' Export the annotation to a defined format (currently only bed is supported
#'
#' @param fname    One-element \code{character} vector giving the name of the file to contain the annotation data. If
#'                 this file already exists, it will be overwritten.
#' @param type     One-element \code{character} vector giving the name of the region annotation.
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#' @param format   Output format. currently only \code{"bed"} is supported.
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' rnb.export.annotation(tempfile(pattern="promoters",fileext=".bed"),"promoters")
#' }
rnb.export.annotation <- function(fname, type, assembly = "hg19", format = "bed") {
	if (!(is.character(fname) && length(fname) == 1 && (!is.na(fname)))) {
		stop("invalid value for fname")
	}
	if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
		stop("invalid value for type")
	}
	if (!(is.character(format) && length(format) == 1 && (!is.na(format)) && format=="bed")) {
		stop("invalid value for format")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	if (!(type %in% c("CpG",names(.rnb.annotations[[assembly]][["regions"]])))) {
		stop("unsupported annotation type")
	}
	aa <- rnb.get.annotation(type,assembly)
	names(aa) <- NULL
	aa <- GenomicRanges::unlist(aa)
	names(aa) <- NULL
	aa.df <- GenomicRanges::as.data.frame(aa)
	colnames(aa.df)[colnames(aa.df)=="seqnames"] <- "chrom"
	aa.df$width <- NULL
	ids <- rownames(aa.df)
	if (length(ids)!=nrow(aa.df)){
		ids <- 1:nrow(aa.df)
	}
	if (format=="bed"){
		tab.2.export <- data.frame(chrom=aa.df$chrom,start=aa.df$start,end=aa.df$end,strand=aa.df$strand,id=ids)
		write.table(format(tab.2.export,scientific=FALSE),file=fname,
				quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE,na=".")
	}
	return(TRUE)
}

########################################################################################################################

#' rnb.export.all.annotation
#'
#' Wrapper for exporting all annotation sets
#'
#' @param out.dir  The directory to write the files to
#' @param types    One-element \code{character} vector giving the name of the region annotation.
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#' @param format   output format. currently only \code{"bed"} is supported.
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' logger.start(fname=NA)
#' rnb.export.all.annotation(tempdir(),c("genes","promoters"))
#' }
rnb.export.all.annotation <- function(out.dir, types=c("CpG",rnb.region.types(assembly)), assembly = "hg19", format = "bed") {
	if (!(is.character(out.dir) && length(out.dir) == 1 && (!is.na(out.dir)))) {
		stop("invalid value for fname")
	}
	if (!(is.character(types) && (!any(is.na(types))))) {
		stop("invalid value for type")
	}
	if (!(is.character(format) && length(format) == 1 && (!is.na(format)) && format=="bed")) {
		stop("invalid value for format")
	}
	for (tt in types){
		if (logger.isinitialized()) {
			logger.start(c("Exporting Region Type",tt))
		}
		rnb.export.annotation(file.path(out.dir,paste(assembly,tt,format,sep=".")), tt, assembly = assembly, format = "bed")
		if (logger.isinitialized()) {
			logger.completed()
		}
	}
	return(TRUE)
}

########################################################################################################################

#' rnb.load.annotation
#'
#' Loads a previously saved custom region annotation from a binary (RData) file.
#' 
#' @param fname One-element \code{character} vector giving the name of the file that contains the annotation data.
#' @param type  One-element \code{character} vector giving the name of the region annotation. If this annotation
#'              is already available, it will be overwritten for the current session.
#' @return Invisibly, \code{TRUE} if the annotation was loaded successfully; an error message if the objects in the
#'         given file do not encode an annotation.
#'
#' @details
#' If the region annotation cannot be loaded from the specified location, this function exits with an error message in
#' the form \code{"unable to load object from ..."}. This could happen, for example, when \code{fname} does not refer to
#' a valid RData file, or the file cannot be accessed due to security restrictions.
#' 
#' If the file is loaded in the current session, but no annotation was added, the function returns invisibly one of the
#' following short failure messages:
#' \describe{
#'   \item{\code{"invalid format"}}{The RData file does not store exactly the following three objects - \code{assembly},
#'        \code{regions}, and \code{mapping}, or they are not of the expected type.}
#'   \item{\code{"unsupported assembly"}}{The specified assembly is unknown.}
#'   \item{\code{"invalid format of regions"}}{The specified region annotation table is invalid.}
#'   \item{\code{"invalid format of mappings"}}{The specified region mapping tables are invalid.}
#' }
#'
#' @seealso \code{\link{rnb.save.annotation}} for saving annotation to a binary file; \code{\link{rnb.set.annotation}}
#'          for loading an annotation from a BED file.
#' @author Yassen Assenov
#' @export
rnb.load.annotation <- function(fname, type) {
	if (!(is.character(fname) && length(fname) == 1 && (!is.na(fname)))) {
		stop("invalid value for fname")
	}
	if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
		stop("invalid value for type")
	}
	result <- tryCatch(load(fname), error = function(err) { return(NULL) })
	if (is.null(result)) {
		stop(paste("unable to load object from", fname))
	}
	if (!setequal(c("assembly", "regions", "mappings"), result)) {
		return(invisible("invalid format"))
	}
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		return(invisible("invalid format"))
	}
	if (!load.annotations(assembly)) {
		return(invisible("unsupported assembly"))
	}
	## TODO: validate regions and validate mappings

	## Compare currently loaded sites with the ones avaialable in the mappings
	site.names.all <- names(.rnb.annotations[[assembly]][["sites"]])
	site.names.mapped <- names(mappings)
	if (!all(site.names.mapped %in% site.names.all)) {
		return(invisible("invalid format of mappings"))
	}
	## Load all unloaded sites
	## TODO: Would it be better to simply ignore the unloaded sites?
	sites.loaded <- .rnb.annotations[[assembly]][["sites"]]
	sites.loaded <- sites.loaded[!sapply(sites.loaded, is.null)]
	site.names.loaded <- names(sites.loaded)
	for (sname in setdiff(site.names.mapped, site.names.loaded)) {
		load.annotations(assembly, sname)
	}
	## Construct additional mappings if necessary
	for (sname in setdiff(site.names.loaded, site.names.mapped)) {
		mappings[[sname]] <- rnb.regions2sites(regions, sites = sites.loaded[[sname]])
	}
	## Put the objects in the .rnb.annotation environment
	sites.all <- .rnb.annotations[[assembly]][["sites"]]
	site.names.loaded <- names(sites.all[!sapply(sites.all, is.null)])
	.rnb.annotations[[assembly]][["regions"]][[type]] <- regions
	.rnb.annotations[[assembly]][["mappings"]][[type]] <- mappings[site.names.loaded]
	rnb.update.annotation.infos(type, assembly)
	invisible(TRUE)
}

########################################################################################################################

#' rnb.remove.annotation
#'
#' Deletes a region annotation table. Use this function with caution; its operation cannot be undone.
#'
#' @param type     One-element \code{character} vector giving the name of the region annotation.
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#' @return Invisibly, \code{TRUE} if the annotation has been successfully deleted, or \code{FALSE} if the specified
#'         region type is not supported.
#'
#' @examples
#' \dontrun{
#' t.regions <- rnb.get.annotation("tiling")
#' rnb.remove.annotation("tiling")
#' }
#'
#' @seealso \code{\link{rnb.get.annotation}}, \code{\link{rnb.region.types}}
#' @author Fabian Mueller
#' @export
rnb.remove.annotation <- function(type, assembly = "hg19") {
	if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
		stop("invalid value for type")
	}
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	if (!(type %in% names(.rnb.annotations[[assembly]][["regions"]]))) {
		warning("unsupported annotation type")
		return(invisible(FALSE))
	}
	.rnb.annotations[[assembly]][["regions"]][[type]]  <- NULL
	.rnb.annotations[[assembly]][["mappings"]][[type]] <- NULL
	i.type <- which(colnames(.rnb.annotations[[assembly]][["lengths"]]) == type)
	.rnb.annotations[[assembly]][["lengths"]] <- .rnb.annotations[[assembly]][["lengths"]][, -i.type]
	invisible(TRUE)
}

########################################################################################################################

#' rnb.region.types
#'
#' Gets the supported region annotations for a given genome assembly.
#'
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#' @return Region types supported by \pkg{RnBeads} in the form of a \code{character} vector. The built-in ones are
#'         \code{"cpgislands"}, \code{"genes"}, \code{"promoters"} and \code{"tiling"}. The names of all custom region
#'         definitions are also included in the returned vector.
#'
#' @examples
#' \dontrun{
#' "promoters" %in% rnb.region.types() # TRUE
#' }
#'
#' @seealso \code{\link{rnb.get.annotation}}, \code{\link{rnb.set.annotation}}
#' @author Yassen Assenov
#' @export
rnb.region.types <- function(assembly = "hg19") {
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	names(.rnb.annotations[[assembly]][["regions"]])
}

########################################################################################################################

#' rnb.get.mapping
#'
#' Gets the mapping information used for a region type. These are structures used to map regions to the genomic loci (or
#' Infinium probes) that target them.
#'
#' @param region.type Region type. The built-in types are \code{"cpgislands"}, \code{"genes"}, \code{"promoters"} and
#'                    \code{"tiling"}.
#' @param target.type Target type for sites.
#' @param assembly    Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported
#'                    genomes.
#' @return \code{list} of mapping structures, one per chromosome. Every mapping structure is an object of type
#'         \code{\link{IRanges}} and stores the range of indices of all sites contained in the respective region.
#'         Regions that do not contain sites are left out of the mapping.
#' @examples
#' \dontrun{
#' promoters2probes <- rnb.get.mapping("promoters", "probes450")
#' promoters2probes[["chr21]]
#' }
#' @author Yassen Assenov
#' @export
rnb.get.mapping <- function(region.type, target.type, assembly = "hg19") {
	if (!(is.character(region.type) && length(region.type) == 1 && (!is.na(region.type)))) {
		stop("invalid value for region.type")
	}
	if (!(is.character(target.type) && length(target.type) == 1 && (!is.na(target.type)))) {
		stop("invalid value for target.type")
	}
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	if (!load.annotations(assembly, target.type)) {
		stop("unsupported target.type")
	}
	if (!(target.type %in% names(.rnb.annotations[[assembly]][["sites"]]))) {
		stop("unsupported target.type")
	}
	if (is.null(.rnb.annotations[[assembly]][["regions"]][[region.type]])) {
		stop("unsupported region.type")
	}
	.rnb.annotations[[assembly]][["mappings"]][[region.type]][[target.type]]
}

########################################################################################################################

#' rnb.get.chromosomes
#' 
#' Gets the chromosome names supported for the specified assembly.
#'
#' @param assembly Genome assembly of interest. See \code{\link{rnb.get.assemblies}} for the list of supported genomes.
#' @return \code{character} vector of supported chromosomes for the specified genome assembly. The elements of the
#'         vector follow the \href{http://www.ensembl.org/}{Ensembl} convention (\code{"1"}, \code{"2"}, ...), and the
#'         names of this vector - the convention of the \href{http://genome.ucsc.edu/}{UCSC Genome Browser}
#'         (\code{"chr1"}, \code{"chr2"}, ...).
#'
#' @examples
#' \dontrun{
#' "chrX" %in% names(rnb.get.chromosomes())
#' }
#' @author Pavlo Lutsik
#' @export
rnb.get.chromosomes <- function(assembly = "hg19") {
	if (!(is.character(assembly) && length(assembly) == 1 && (!is.na(assembly)))) {
		stop("invalid value for assembly")
	}
	assembly <- tolower(assembly)
	if (!load.annotations(assembly)) {
		stop("unsupported assembly")
	}
	return(CHROMOSOMES.L2S[[assembly]])
}

########################################################################################################################

#' rnb.infinium.control.targets
#'
#' Extracts all control probe types in the HumanMethylation450 assay.
#' 
#' @param target  A singleton of type \code{character}, specifying the microarray platform. 
#' 				  \code{"probes450"} and \code{"probes27"} correspond to HumanMethylation450 
#' 				 respectively HumanMethylation27 microarrays
#' 
#' @return \code{character} vector of control targets.
#'
#' @examples
#' \dontrun{
#' "NEGATIVE" %in% rnb.infinium.control.targets()
#' }
#' @author Pavlo Lutsik
#' @export 
rnb.infinium.control.targets <- function(target="probes450") {
	if(target=="probes450"){
		return(HM450.CONTROL.TARGETS)
	}else if(target=="probes27"){
		return(HM27.CONTROL.TARGETS)
	}else{
		stop("Unsupported platform")
	}
}
