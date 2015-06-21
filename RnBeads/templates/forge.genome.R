## Forges a genome directory with a DCF file and creates the directory at a specific location (destdir).  
## The directory is then made with the genome (.2bit) file in it and is ready to be build (R CMD BUILD)
## and installed (R CMD INSTALL).
## DCF = Description file which got the information where the 2bit file is located and what the 
## genomic package name should be, species name, etc.. 
## tmp = the given tmp folder destination. 
suppressPackageStartupMessages(library(BSgenome))
forgeBSgenomeDataPkg("%(DCF)s", destdir="%(tmp)s")
q(save = "no")
