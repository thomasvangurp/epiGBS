# Runs the RnBeads analysis with the given parameters. 
# If a annotation file is given: calculate the the CpG of the given chromosomes too. 
.libPaths(c('%(lib_path)s',.libPaths()))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RnBeads))

suppressPackageStartupMessages(require(%(package)s))
genome.data <- %(species)s

annotation.names <- c(%(annotation_names)s)
annotation.files <- c(%(annotation_files)s)
index <- 0

# If there are annotation files:
if (!is.null(annotation.files)){
  # Sets the RnBeads options with the annotation names 
  rnb.options(region.types=annotation.names)
  # Iterates through the annotation files (could be more than one).
  for (file.name in annotation.files){
    # Puts the file data in a R table. 
    chr.list <- read.table(file.name)
    index <- index + 1
    chromosomes <- character()
    startings <- integer()
    ends <- integer()
    names <- character()
    strand <- character()
    added.length <- 0
    
    for (i in 1:length(chr.list[[1]])){
      if (length(chr.list) > 1){
        # If there is specific position based annotation, e.g. chr1 50  100 \n chr1 150 200 (values separated by tab). 
        cpg.pos.list <- matchPattern("CG", genome.data[[toString(chr.list[[1]][i])]][chr.list[[2]][i]:chr.list[[3]][i]])
      } else {
        # If the whole chromosome is 1 annotation, e.g. chr1 \n chr2 
        cpg.pos.list <- matchPattern("CG", genome.data[[toString(chr.list[[1]][i])]])
      }
      # If there are CpG sites in the annotation:
      if (length(start(cpg.pos.list) >= 1)){
        # Puts every position/chromosome in an variable. 
        for (j in 1:length(cpg.pos.list)){
          chromosomes[[j+added.length]] <- toString(chr.list[[1]][i])
          startings[[j+added.length]] <- start(cpg.pos.list)[[j]]
          ends[[j+added.length]] <- start(cpg.pos.list)[[j]]+1L
          names[[j+added.length]] <- paste("annotation_",toString(chr.list[[1]][j]),sep = "")
          strand[[j+added.length]] <- "+"
        }
        added.length <- added.length + length(cpg.pos.list)
      } 
    }
    # Appends every CpG position to a matrix.
    gene.regions <- data.frame(Chromosome=chromosomes, Start=startings, 
                               End=ends,rownames=names, strand=strand)
    # Appends the RnBeads package with the matrix and calculates each CpG site.
    rnb.set.annotation.and.cpg.stats(annotation.names[[index]], gene.regions, assembly="%(assembly)s") 
}
# If there is no annotation file. 
} else {
  rnb.options(region.types="")
}
# TODO: read options from file that can be submitted by user.
rnb.options(assembly="%(assembly)s", import.bed.style="EPP",
            preprocessing = TRUE, filtering.greedycut = FALSE, filtering.high.coverage.outliers = TRUE,
            filtering.snp = "no", filtering.context.removal = NULL, exploratory = TRUE, 
            exploratory.intersample = FALSE, exploratory.region.profiles = character())

report.dir <- "%(directory)s/analysis/%(time)s"
data.source <- c("%(directory)s/bs.bed/","%(directory)s/sample.csv")

# Sets the number of cores to use for the analysis.
logger.start(fname=NA)
num.cores <- %(cores)s
parallel.setup(num.cores)

# Main analysis function of RnBeads. 
rnb.run.analysis(dir.reports=report.dir, data.source=data.source,
                 data.type="bs.bed.dir")

# Quits the program after execution
q(save = "no")