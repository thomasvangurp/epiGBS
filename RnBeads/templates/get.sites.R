## Calculates each CpG position in the given genomes and saves it in an .RData variable in the assembly template folder. 

# If custom installation folder: .lib_paths will be appended with the given installation folder.  
.libPaths(c('%(lib_path)s',.libPaths()))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RnBeads))

# Writes the chromosomes to the assembly code variable. 
# The chromosomes file is located in the /assemblies/assembly_code/bs.bed/chromosomes.txt 
%(assembly)s.chr <- read.table("%(chromosomes)s")
CHROMOSOMES.L2S <- list("%(assembly)s" = %(assembly)s.chr[[1]])
CHROMOSOMES.S2L <- lapply(CHROMOSOMES.L2S, function(x) { paste0("chr", x) })
CHROMOSOMES <- CHROMOSOMES.S2L
for (assembly.name in names(CHROMOSOMES)) {
  names(CHROMOSOMES.S2L[[assembly.name]]) <- CHROMOSOMES.L2S[[assembly.name]]
  names(CHROMOSOMES[[assembly.name]]) <- names(CHROMOSOMES.L2S[[assembly.name]]) <- CHROMOSOMES[[assembly.name]]
}
rm(assembly.name)

## Dinucleotide patterns to be annotated
NUCLEOTIDE.PATTERNS <- c("CG")
names(NUCLEOTIDE.PATTERNS) <- sapply(strsplit(NUCLEOTIDE.PATTERNS, ""), paste, collapse = "p")

# Loads all the genomic information of the the given genome file. 
suppressPackageStartupMessages(require(%(package)s))
genome.data <- %(species)s

processed.chrom <- 0L
count.chrom <- 0L

# Function: Calculates every CpG position on the given chromosomes and writes it in a GRangesList object. 
rnb.update.sites <- function(assembly = "%(assembly)s") {
  chrom.lengths <- seqlengths(genome.data)[CHROMOSOMES[[assembly]]]
  sites <- list()
  pp.dnas <- DNAStringSet(NUCLEOTIDE.PATTERNS)
  for (i in names(NUCLEOTIDE.PATTERNS)){
    pp.p <- pp.dnas[[i]]
    pp.m <- reverseComplement(pp.p)
    curSites <- lapply(CHROMOSOMES[[assembly]], function(chrom) {
      matches.st <- lapply(list(pp.p, pp.m), function(x) {
        ranges(matchPattern(x, genome.data[[chrom]], fixed=FALSE))
      })
      count.chrom <- count.chrom + 1L
      # For every 50 processed chromosomes, it will print the number of chromosomes that already have been analysed.  
      if (count.chrom == 50){
        processed.chrom <- processed.chrom + count.chrom
        print(paste("Processed:", toString(processed.chrom), "chromosomes."))
        count.chrom <- 0L
      }
      length.neighborhood = length(genome.data[[chrom]])
      cp.starts <- start(matches.st[[1]])
      cp.ends <- cp.starts + 1L 
      
      dna_seq <- genome.data[[chrom]]
      cpg.stats <- suppressWarnings(get.cpg.stats(dna_seq, cp.starts, cp.ends))
      matches.gr <- mapply(GRanges, seqnames = list(chrom), ranges = matches.st, strand = list("+", "-"),
                           "CpG" = list(cpg.stats[, "CpG"]), "GC" = list(cpg.stats[, "GC"]))
      matches.gr <- RnBeads:::rnb.sort.regions(do.call(c, matches.gr))
      seqlevels(matches.gr) <- names(CHROMOSOMES[[assembly]])
      seqlengths(matches.gr) <- chrom.lengths
      matches.gr
    })
    sites[[i]] <- GRangesList(curSites)
  }
  return(sites)
}

# Writes every CpG position of every chromosomes to the sites_raw variable.
sites_raw <- rnb.update.sites()
# Makes/writes the CpG positions to the sites variable with the right structure RnBeads can accept. 
sites <- vector(mode="list", length=1)
names(sites) <- "sites"
sites[["sites"]] <- sites_raw[[1]]
# Saves the sites variable to a RData file and saves it to the /assembly/data/ folder. 
save(sites, file="%(sites_output)s")

# Makes a regions .RData file just because it it required for RnBeads but it will be empty.  
regions <- vector(mode="list", length=1)
names(regions) <- "genes"
regions[["genes"]] <- GRangesList()

# The only the thing the regions .RData file gets is some attribute information that is required for the analysis.
# The actual (annotation) regions will be appended with the run.analysis.R script. 
builtin <- logical(length = 4)
names(builtin) <-  c("tiling", "genes", "promoters", "cpgislands")
builtin[[2]] <- TRUE
attr(regions, "builtin") <- builtin
save(regions, file="%(regions_output)s")

# Quit the script after execution. 
q(save = "no")
