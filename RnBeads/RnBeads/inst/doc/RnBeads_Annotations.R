### R code from vignette source 'RnBeads_Annotations.Rnw'

###################################################
### code chunk number 1: RnBeads_Annotations.Rnw:152-155
###################################################
suppressPackageStartupMessages(library(RnBeads))
#for knitr: avoid tidy because it messes up line breaks
# opts_chunk$set(tidy=FALSE)


###################################################
### code chunk number 2: RnBeads_Annotations.Rnw:162-163
###################################################
rnb.get.assemblies()


###################################################
### code chunk number 3: RnBeads_Annotations.Rnw:195-197 (eval = FALSE)
###################################################
## control.annotation <- rnb.get.annotation("controls450")
## head(control.annotation)


###################################################
### code chunk number 4: RnBeads_Annotations.Rnw:228-230 (eval = FALSE)
###################################################
## probe.annotation <- rnb.get.annotation("probes450")
## probe.annotation


###################################################
### code chunk number 5: RnBeads_Annotations.Rnw:252-255 (eval = FALSE)
###################################################
## gene.annotation <- rnb.get.annotation("genes")
## attr(gene.annotation, "version")
## gene.annotation


###################################################
### code chunk number 6: RnBeads_Annotations.Rnw:262-264 (eval = FALSE)
###################################################
## promoters <- rnb.annotation2data.frame(rnb.get.annotation("promoters"))
## write.csv(promoters, file = "promoters.csv", row.names = FALSE)


