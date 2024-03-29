library(ggplot2)
reads <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/reads_info.tsv.gz'), sep="\t", header=T, row.names=1)
pc_lympho <- read.table('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_values', sep="\t", header=F, row.names=1)

summary(reads)
reads$type <- substr(rownames(reads), 8, 10)

ggplot(data=reads, aes(frac_h, ..count.., fill=type))+geom_histogram(position='dodge', bins=10)
ggplot(data=reads, aes(frac_h, fill=type))+geom_histogram(position='dodge', bins=10)
# plots and info on 210215_MAbd-growth_EGFPaneth_rnaqc

m <- merge(reads, pc_lympho, by="row.names")
colnames(m)[9] <- 'leuco_score'
colnames(m)[10] <- 'PC2'
colnames(m)[11] <- 'type_ext'
ggplot(data=m, aes(leuco_score, fill=type))+geom_histogram(position='dodge', bins=10)

ggplot(data=m, aes(x=leuco_score, y=PC2, color=type))+geom_point()
ggplot(data=m, aes(x=leuco_score, y=frac_h, color=type))+geom_point()

m$model <- substr(m$Row.names, 0, 7)


table(m[m$leuco_score > 50 & m$PC2 > 30,"model"])
length(unique(m[m$leuco_score > 50 & m$PC2 > 30,"model"]))

# TODO curiosity other scores?

all_scores <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK/fpkm_lymphoma_scores.tsv.gz'), sep='\t', header=T, row.names=1)
m <- merge(reads, all_scores, by="row.names")


ggplot(data=m, aes(Leucocyte, fill=type))+geom_histogram(position='dodge', bins=10)
ggplot(data=m, aes(Endothelial, fill=type))+geom_histogram(position='dodge', bins=10)
ggplot(data=m, aes(CAF, fill=type))+geom_histogram(position='dodge', bins=10)
ggplot(data=m, aes(x=CAF, y=Leucocyte, color=type))+geom_point()


#

reads <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/reads_info.tsv.gz'), sep="\t", header=T, row.names=1)
mscores <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_M/M_lymphoma_scores.tsv.gz'), sep="\t", header=T, row.names=1)
hscores <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK/fpkm_H_lymphoma_scores.tsv.gz'), sep="\t", header=T, row.names=1)
colnames(mscores) <- paste0("M_", colnames(mscores))
colnames(hscores) <- paste0("H_", colnames(hscores))
m <- merge(reads, mscores, by="row.names")
rownames(m) <- m$Row.names
m$Row.names <- NULL
m <- merge(m, hscores, by="row.names")
rownames(m) <- m$Row.names
m$Row.names <- NULL


ggplot(data=m, aes(M_Leucocyte, fill=type))+geom_histogram(position='dodge', bins=10)
ggplot(data=m, aes(M_CAF, fill=type))+geom_histogram(position='dodge', bins=10)
ggplot(data=m, aes(M_Endothelial, fill=type))+geom_histogram(position='dodge', bins=10)
ggplot(data=m, aes(x=M_Leucocyte, y=m_tot, color=type))+geom_point()
ggplot(data=m, aes(x=M_CAF, y=m_tot, color=type))+geom_point()
ggplot(data=m, aes(x=M_Endothelial, y=m_tot, color=type))+geom_point()



pc_lympho <- read.table('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_values', sep="\t", header=F, row.names=1)
meta <- pc_lympho[,c(4,5)]
colnames(meta) <- c('batch', 'type')
m <- merge(m, meta, by="row.names")
rownames(m) <- m$Row.names
m$Row.names <- NULL
