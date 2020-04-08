#!/usr/bin/env Rscript
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

d <- read.table(input, sep="\t",header=TRUE, stringsAsFactors = FALSE)
# stringsAsFactors=F to avoid problems when melting due to different levels
d$N <- NULL
longd <- melt(d, id.vars=c("Genealogy.ID"))

#correct!
#> length(unique(longd$value))
#[1] 1232
#> length(longd[longd$value!="", "value"])
#[1] 1231
# vs:
#$ bawk '$1 != ""' have_metadata | wc -l
#1231

longd$Genealogy.ID <- NULL
longd <- longd[longd$value != "",]
write.table(longd, file=output, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

