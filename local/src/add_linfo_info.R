#!/usr/bin/env Rscript

set.seed(42)


annot <- snakemake@input[['annot']]
fra <- snakemake@input[['fra']]
marco <- snakemake@input[['marco']]
output <- snakemake@output[[1]]

annot <- read.table(annot, sep="\t", quote="", header=FALSE)
fra <- read.table(fra, sep="\t", quote="", header=TRUE) # new file of 2020
for ( i in 1:length(fra$linfoma) ) {
    if ( fra$linfoma[i]=='noinfo' ) {
        fra$linfoma[i] <- 'no'
    }
}
marco <- read.table(marco, sep="\t", quote="", header=TRUE)

fra$branch <- substr(fra[,2], 1, 12)
fra$type <- NULL
annot$branch <- substr(annot[,1], 1, 12)
marco$branch <- substr(marco[,2], 1, 12)
marco$type <- NULL

colnames(annot) <- c("sample_id_R", "RNA_marker", "RNA_PC", "batch", "type", "w3_cetuxi", "w3_irino", "branch")
# colnames(fra) <- c("sample_id_FF", "which_linfo", "branch")
save.image("pippo.RData")
df <- merge(annot, fra, by="branch", all.x = TRUE)
df$FRA_L <- ifelse(df$linfoma == "no", "pass", ifelse(df$linfoma=='NT', "NotTested", "FRA_L"))

df <- merge(df, marco, by="branch", all.x = TRUE)
df$METHYL_L <- ifelse(df$PC2>500, "METHYL_L", ifelse(is.na(df$cluster), "NA", "pass"))

# df[,c("sample_id_FF","which_linfo","crc","mouse","branch","cluster")] <- NULL
df <- df[,c("sample_id_R","RNA_marker","RNA_PC","METHYL_L","FRA_L","batch","type","w3_cetuxi","w3_irino")]
df <- df[order(df$batch),]

write.table(df, output, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
