library(reshape2)
library(ComplexHeatmap)
setwd("/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/april2020")
data <- read.table("metadata_batch_w3vivocetuxi_w3vivoirino", sep="\t", header=F)
colnames(data) <- c('lgene','type','batch','cetuxi','irino')
data$tissue <- substr(data$lgene,8,9)
data$case <- substr(data$lgene,0,7)
data$case_ext <- substr(data$lgene,0,9)

data$batch <- NULL
#data$lgene <- NULL
data <- unique(data)
data$typenonum <- unlist(lapply(data$type, function(x) { 
  y <- strsplit(as.character(x), '.', fixed=TRUE)[[1]][1]
  # y2 <- strsplit(y, '_', fixed=TRUE)[[1]]
  # last <- length(y2)
  # if (!is.na(as.integer(y2[last]))) {
  #   y2 <- y2[-last]
  # }
  # return(paste(y2, collapse="_"))
  return(y)
}))

widetype <- dcast(data, case_ext ~ typenonum)

treat <- unique(data[,c("cetuxi","irino","tissue","case_ext")])
#batch <- data[,c("batch","case_ext")]
#batchu <- data.frame(case_ext=unique(data$case_ext))
#batchu$batch <- sapply(unique(data$case_ext), function(x) { min(batch[batch$case_ext==x,"batch"])})

whole <- merge(widetype, treat, by="case_ext")
rownames(whole) <- whole$case_ext
whole$case_ext <- NULL

dcetux <- whole[!is.na(whole$cetuxi),]
dcetux <- dcetux[order(dcetux$cetuxi),]
column_order_cetux <- rownames(dcetux)

dirino <- whole[!is.na(whole$irino),]
dirino <- dirino[order(dirino$irino),]
column_order_irino <- rownames(dirino)

#whole <- whole[order(whole$tissue, rowSums(whole)),]

whole <- whole[order(whole$tissue),]
whole$tissue <- NULL

cetuxi <- whole$cetuxi
irino <- whole$irino
whole$cetuxi <- NULL
whole$irino <- NULL


col = c("one" = "blue", "two" = "red", "three"= "#e86f0a", "four" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  one = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["one"], col = NA))
  },
  # bug red
  two = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["two"], col = NA))
  },
  three = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["three"], col = NA))
  },
  # small green
  # MISSENSE = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["MISSENSE"], col = NA))
  # },
  # STOP = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["STOP"], col = NA))
  # },
  # INDEL = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["INDEL"], col = NA))
  # },
  four = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["four"], col = NA))
  }
)

column_title = "OncoPrint remastered"

heatmap_legend_param = list(title = "Replicates", at = c("one", "two", "three", "four"))

det <- t(whole)

column_order_all <- colnames(det)
det[det==0] <- ""
det[det==1] <- "one"
det[det==2] <- "two"
det[det==3] <- "three"
det[det==4] <- "four"


oncoPrint(det,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
          row_names_gp = gpar(fontsize=8), column_order=column_order_all)
          #top_annotation = NULL)


det_irino <- det[, colnames(det) %in% rownames(dirino)]
det_irino <- det_irino[ ,match(column_order_irino, colnames(det_irino))]
det_cetuxi <- det[, colnames(det) %in% rownames(dcetux)]
det_cetuxi <- det_cetuxi[ ,match(column_order_cetux, colnames(det_cetuxi))]

oncoPrint(det_cetuxi,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
          row_names_gp = gpar(fontsize=8), column_order=column_order_cetux,
          top_annotation = HeatmapAnnotation(W3=anno_barplot(dcetux$cetuxi)))

oncoPrint(det_irino,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
          row_names_gp = gpar(fontsize=8), column_order=column_order_irino,
          top_annotation = HeatmapAnnotation(W3=anno_barplot(dirino$irino)))

