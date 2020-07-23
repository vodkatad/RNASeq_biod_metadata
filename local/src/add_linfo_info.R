annot <- snakemake@input[['annot']]
fra <- snakemake@input[['fra']]
marco <- snakemake@input[['marco']]
output <- snakemake@output[[1]]

annot <- read.table(annot, sep="\t", quote="", header=FALSE)
fra <- read.table(fra, sep="\t", quote="", header=FALSE)
marco <- read.table(marco, sep="\t", quote="", header=TRUE)

fra$branch <- substr(fra$V1, 1, 12)
annot$branch <- substr(annot$V1, 1, 12)

colnames(annot) <- c("sample_id_R", "RNA_marker", "RNA_PC", "batch", "type", "w3_cetuxi", "w3_irino", "branch")
colnames(fra) <- c("sample_id_FF", "which_linfo", "branch")

df <- merge(annot, fra, by="branch", all.x = TRUE)
df$FRA_L <- ifelse(df$which_linfo == 'r', "FRA_L", ifelse(df$which_linfo=='b', "NotTested", ifelse(is.na(df$which_linfo), "NA", "pass")))

df <- merge(df, marco, by="branch", all.x = TRUE)                         
df$METHYL_L <- ifelse(df$cluster==1, "METHYL_L", ifelse(is.na(df$cluster), "NA", "pass"))

df[,c("sample_id_FF","which_linfo","crc","mouse","branch","cluster")] <- NULL
df <- df[,c("sample_id_R","RNA_marker","RNA_PC","METHYL_L","FRA_L","batch","type","w3_cetuxi","w3_irino")]
df <- df[order(df$batch),]

write.table(df, output, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)