#!/usr/bin/env Rscript

set.seed(42)

score <- snakemake@input[[1]]
output <- snakemake@output[[1]]
differ <- snakemake@output[[2]]


score <- read.table(gzfile(score), sep='\t', quote="", header=TRUE, row.names=1)

info <- data.frame(row.names=rownames(info), stringsAsFactors=FALSE)


for ( i in 1:length(type) ) {
  tmp <- df[df$type==type[i],]
  if ( nrow(tmp) == 1 ) {
    info[i,"Endothelial"] <- round(tmp$Endothelial, digits=8)
    info[i,"CAF"] <- round(tmp$CAF, digits=8)
    info[i,"Leucocyte"] <- round(tmp$Leucocyte, digits=8)
  } else {
    keep <- FALSE
    j <- 1
    while (!keep & j <= nrow(tmp) ) {
      if ( tmp[j,"mouse_score"] ==  tmp[j,"mouse_info"] ) {
        info[i,"Endothelial"] <- round(tmp[j,"Endothelial"], digits=8)
        info[i,"CAF"] <- round(tmp[j,"CAF"], digits=8)
        info[i,"Leucocyte"] <- round(tmp[j,"Leucocyte"], digits=8)
        # info[i,"same_mouse"] <- "yes"
        keep <- TRUE
      }
    j <- j+1
    }
    if ( !keep ) {
      j <- 1
      info[i,"Endothelial"] <- round(mean(tmp$Endothelial), digits=8)
      info[i,"CAF"] <- round(mean(tmp$CAF), digits=8)
      info[i,"Leucocyte"] <- round(mean(tmp$Leucocyte), digits=8)
    #   info[i,"same_mouse"] <- "no"
    #   diff[i,"Endothelial"] <- round(abs(tmp[j,"Endothelial"] - tmp[j+1,"Endothelial"]), digits=8)
    #   diff[i,"CAF"] <- round(abs(tmp[j,"CAF"] - tmp[j+1,"CAF"]), digits=8)
    #   diff[i,"Leucocyte"] <- round(abs(tmp[j,"Leucocyte"] - tmp[j+1,"Leucocyte"]), digits=8)
    }
  }
}

info[is.na(info)] <- ""
# diff <- diff[!is.na(diff$Endothelial),]
names(info)[c(1,2)] <- c("type","mouse")


write.table(info, output, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
# write.table(diff, differ, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
