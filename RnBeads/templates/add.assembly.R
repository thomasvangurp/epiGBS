## Appends the given assembly code the RnBeads`s assemblies annotations.RData file. 
## The file`s location is: RnBeads/data/annotation.RData. 

# Removes all current variables in the R .env() 
rm(list=ls(all=TRUE)) 
# Load the .annotation.RData file.
load(%(file)s)
#Coppies the mm9 data to the new assembly code. (mm9 just as a tempalte)
%(assembly)s <- mm9
save.image(file="%(tmp)sannotations.RData")
q(save = "no")