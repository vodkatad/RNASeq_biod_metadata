#!/usr/bin/env Rscript

library(ggplot2)
library(getopt)


opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'input_in', 'i', 1, 'character',
  'y_in', 'y', 1, 'character',
  'color_in', 'c', 1, 'character',
  'shape_in', 's', 2, 'character',
  'output', 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)


df <- read.table(gzfile(opt$input_in), sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
y <- opt$y
color <- opt$color_in
names(df) <- c("sample","Leucocyte","PC1","PC2","batch","type","w3_cetuxi","w3_irino")
if ( opt$shape_in == "not") {
    shape <- NULL
} else {
    shape <- opt$shape_in
    if ( color == shape ) {
    stop("Please use different color and shape")
    }
}


if ( color == "Leucocyte" | color == "CAF" | color == "Endothelial") {
p <- ggplot(df, aes_string(x="PC1", y=y, color=color, shape=shape)) +
  labs(title="Scatter plot PC") +
  scale_color_distiller(type="seq", palette="OrRd", direction=1) +
  scale_shape_manual(values=c(0,1,2,3,7,16,17,8)) +
  geom_point() +
  scale_y_continuous(breaks=seq(-50, 70, 10)) +
  geom_hline(yintercept=0, color="black") +
  geom_vline(xintercept=0, color="black") +
  theme_bw() #+
  # theme(axis.line=element_line(size=2, colour = "black"))
} else {
p <- ggplot(df, aes_string(x="PC1", y=y, color=color, shape=shape)) +
  labs(title="Scatter plot PC") +
  scale_color_brewer(type='qual', palette="Dark2") +
  scale_shape_manual(values=c(0,1,2,3,7,16,17,8)) +
  geom_point() + 
  scale_y_continuous(breaks=seq(-50, 70, 10)) +
  geom_hline(yintercept=0, color="black") +
  geom_vline(xintercept=0, color="black") +
  theme_bw() #+
  # theme(axis.line=element_line(size=2, colour = "black"))
}

png(opt$output, width=8, height=5, units="in", type="cairo", res=300)
p
dev.off()


# png(opt$output, width=8, height=5, units="in", type="cairo", res=300)
# ggplot(df, aes_string(x="PC1", y=y, color=color, shape=shape)) +
#   labs(title="Scatter plot PC") +
#   scale_color_brewer(type='qual', palette="Dark2") +
#   # geom_density_2d() +
#   geom_point() + 
#   coord_cartesian(ylim=c(-1000,1000), xlim=c(-1000,1000)) +
#   theme_minimal()
# dev.off()
