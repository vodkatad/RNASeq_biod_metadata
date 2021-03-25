library(NOISeq)
setwd('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK')
counts <- 'H_readcount.tsv.gz'
meta <- 'selected_metadata_annot_final_nolinfo_nooutlier'
counts <- read.table(counts, sep="\t", header=T, row.names=1)
meta <- read.table(meta, sep="\t", header=T, row.names=1)
library(NOISeq)
head(colnames(counts))
colnames(counts)[grepl('.',(colnames(counts)), fixed=TRUE)]
colnames(counts) <- gsub('.','-', colnames(counts), fixed=TRUE)
counts <- counts[, colnames(counts) %in% rownames(meta)]
dim(meta)
dim(counts)
noi <- readData(counts, meta)

library("EnsDb.Hsapiens.v86")
ens <- EnsDb.Hsapiens.v86
txtypes <- genes(ens, columns=c("symbol", "gene_biotype"))
tx <- data.frame(entrez=mcols(txtypes)$symbol, biotype=mcols(txtypes)$gene_biotype)

txfu <- unique(tx)
biot <- data.frame(row.names=unique(as.character(txfu$entrez)))
d <- data.frame(row.names=make.unique(as.character(txfu$entrez)), biot=txfu$biotype)
biotm <- merge(biot, d, by="row.names")

noi <- addData(noi, biotype=biotm)

myexplodata <- dat(noi, type = "biodetection")

explo.plot(myexplodata, plottype = "persample")


myexplodata <- dat(noi, type = "biodetection", k=2, type="biodetection", factor="type")

for (i in seq(1, length(levels(meta$type)))) {
  png(paste0('noiseq',i, '.png'))
  explo.plot(myexplodata, plottype = "persample", samples=i)
  dev.off()
}

mycountsbio = dat(noi, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")

mycountsbio = dat(noi, factor = "type", type = "countsbio",k=2)
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")

mysaturation = dat(noi, k = 2, ndepth = 7, type = "saturation")

explo.plot(mysaturation, toplot = "protein_coding", samples = 657:658, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = "protein_coding", samples = 266:659, yleftlim = NULL, yrightlim = NULL)
