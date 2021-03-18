#!/usr/bin/env Rscript

# Script to perform stepwise (chr by chr) Roar analysis. 
# Requires a gtf with _PRE and _POST gene_ids and bam files from the two 
# conditions to be compared.
# Count post randomizations and expr/roar correlations.

meanAcrossAssays <- function(assays, wantedColumns) {
   # Is the conversion to dataframe slow? Is there a more efficient way without
   # changing the SE/Assays structure?
   wantedCols <- lapply(assays, function(x) { x[,wantedColumns] } )
   return(rowMeans(as.data.frame(wantedCols)))
}

getmM <- function(alignments, gtf, stranded) {   
   summOv <- function(x) {
      summarizeOverlaps(features=gtf, reads=x, ignore.strand=!stranded, mc.cores=rds@cores)
   }
   summOvPost <- function(x) {
      summarizeOverlaps(features=postCoords, reads=x, ignore.strand=!stranded, mc.cores=rds@cores)
   } 
   
   # Now we need to keep means and totals of counts over PRE/POST for the two lists.
   # In the simpler case with a single alignment for both conditions we just keep the counts.
   preElems <- grep("_PRE$", mcols(gtf)$gene_id)
   postElems <- grep("_POST$", mcols(gtf)$gene_id)
   if (length(preElems)+length(postElems) != length(mcols(gtf)$gene_id)) {
      stop("The prePostCoords given for this RoarDataset are wrong: some of the gene_id
           does not end in _PRE/_POST.")
   }
   if (length(preElems) != length(postElems)) {
      stop("The prePostCoords given for this RoarDataset are wrong: the number of PRE is different
           from the number of POST.")
   }
   # The number of elements has to be checked because some horrible special case with recycling
   # of vector elements in the comparison are possible.
   if (!all(sub("_PRE", "", mcols(gtf)$gene_id[preElems]) ==
               sub("_POST", "", mcols(gtf)$gene_id[postElems]))) {
      stop("The prePostCoords given for this RoarDataset are wrong: not all prefixes of PRE-POST
           correspond.")
   }
   # Check uniqueness of gene_id.
   geneIds <- sub("_PRE", "", mcols(gtf)$gene_id[preElems])
   if (!all(geneIds==make.unique(geneIds))) {
      stop("The prePostCoords given for this RoarDataset are wrong: gene_ids (prefixes of PRE-POST)
           are not unique.")
   }
   preCoords <- gtf[preElems,]
   postCoords <- gtf[postElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=length(gtf)/2, ncol=2),
                              rowRanges=preCoords, 
                              colData=DataFrame(row.names=c("pre","post"))
   )
   if (length(alignments) == 1) {
      SEpre <- summOv(alignments[[1]])
      SEpost <- summOvPost(alignments[[1]])
      assay(se,1)[,"pre"] <- assays(SEpre)$counts[preElems,]
      assay(se,1)[,"post"] <- assays(SEpost)$counts 
   } else {
      len <- length(alignments)
      counts <- SummarizedExperiment(assays = matrix(nrow=length(gtf)/2, ncol=2),
                                     rowRanges=preCoords, 
                                     colData=DataFrame(row.names=c("pre","post"))
      )
      
      for (i in 1:length(alignments)) {
         SEpre <- summOv(alignments[[i]])
         SEpost <- summOvPost(alignments[[i]])
         assay(counts,i) <- matrix(nrow=length(gtf)/2, ncol=2)
         assay(counts,i)[,"pre"] <- assays(SEpre)$counts[preElems,]
         assay(counts,i)[,"post"] <- assays(SEpost)$counts 
      }
   }
   preLen <- end(preCoords) - start(preCoords) + 1
   postLen <- end(postCoords) - start(postCoords) + 1
   # I had to add "+1" as long as the end coords are not inclusive.
   # Then the bam lengths to correct our lengths, ie: postLen+ReadLength-1
   if (length(alignments) > 1) {
      assay(se,1)[,"pre"] <- meanAcrossAssays(assays(counts), "pre") # here peak of memory usage?
      assay(se,1)[,"post"] <- meanAcrossAssays(assays(counts),"post")
      # Also the length correction should consider all the samples!
      len <- unlist(lapply(alignments, qwidth))
      corr <- mean(len)
   } else {
      corr <- mean(qwidth(alignments[[1]]))
                   # qwidth(x): Returns an integer vector of length length(x) containing the length 
                   # of the query *after* hard clipping (i.e. the length of the query sequence 
                   # that is stored in the corresponding SAM/BAM record).
   }
   # Ok, now if we had a single sample for both conditions we had the data charged in
   # countPrePost, otherwise we have the means (in the same SE/RDS object).
   # If there is a single sample for one condition and more than one for the other there is
   # a little (I hope) unuseful overload to get the mean for the single sample. 
   # The countsTreatment/control slot are still kept as long as we will need them in computePvals.
   postLen <- postLen + corr - 1
   mM <- (assay(se,1)[,"pre"]*postLen)/(assay(se,1)[,"post"]*preLen)-1
   res <- data.frame(row.names=sub("^\\s+","",sub("_POST","",mcols(postCoords)$gene_id)), 
                     mM=mM,
                     counts=assay(se,1)[,"pre"])
}

checkReadable <- function(filename) {
   res <- file.access(names=filename, mode=4) == 0
   if (!res) {
      warning(paste(filename, "is not readable", sep=" "))
   }
   res
}

arguments <- matrix(c(
   'help', 'h', 0, "logical",
   'debug', 'd', 1, "character",
   'gtf' , 'a', 1, "character",
   'stranded'  , 's', 1, "character",
   'bams'  , 'b', 1, "character"
), ncol=4, byrow=T)

library(getopt)
opt <- getopt(arguments)

if (!is.null(opt$help)) {
   stop(getopt(arguments, command=get_Rscript_filename(), usage=TRUE))
}

if (is.null(opt$gtf)) {
   stop("Missing gtf [-a filename] annotation option\n")
}

if (is.null(opt$bams)) {
   stop("Missing -b param")
}

if (is.null(opt$stranded) | (opt$stranded != "yes" & opt$stranded != "no")) {
   stop("The stranded option is needed and it should be 'yes' or 'no'")
}
stranded <- FALSE
if (opt$stranded == "yes") {
	stranded <- TRUE
}


library(roar)
library(rtracklayer)
library(Rsamtools)
library(GenomicAlignments)

bams <- as.vector(unlist(strsplit(opt$bams, ",")))

if (!all(sapply(c(bams, opt$gtf), checkReadable))) {
   stop("One of the given files does not exist or is not readable")  
}

# @gilo
# FIXME gtfGRanges<- import(opt$gtf, asRangedData=FALSE)
# Error in .local(con, format, text, ...) :
#       unused argument (asRangedData = FALSE)
gtfGRanges <- import(opt$gtf)
chrs <- seqlevels(gtfGRanges)

orderBam <- function(bam) {
   tmp <- tempfile()
   ordered <- sortBam(bam, tmp, byQname=FALSE, maxMemory=512)
   garbage <- indexBam(ordered)
   return(ordered)
}

orderedBams <- lapply(bams, orderBam)

workOnChr <- function(chr) {
   write(paste("Working on", chr), stderr())
   reduced <- keepSeqlevels(gtfGRanges, chr, pruning.mode="tidy") 
   coords <- c(start(reduced), end(reduced))  # To keep strandness in consideration! Is this needed?
   begin <- min(coords)
   end <- max(coords)
   spanChr <- GRanges(seqnames=chr,ranges=IRanges(start=begin,width=end-begin+1))
   loadBam <- function(bam) {
      param <- ScanBamParam(which=spanChr)
      res <- readGAlignments(file=bam, param = param)
      return(res)
   } 

   bamsGenomicAlignments <- lapply(orderedBams, loadBam)
   
   res <- getmM(bamsGenomicAlignments, reduced, stranded)
   return(res)
}

allRes <- lapply(chrs, workOnChr)
meltedRes <- do.call("rbind", allRes)
preElems <- grep("_PRE$", mcols(gtfGRanges)$gene_id)
pre <- gtfGRanges[preElems,]
preLen <- end(pre) - start(pre) + 1
names <- sub("^\\s+","",sub("_PRE", "",mcols(pre)$gene_id))
meltedRes <- meltedRes[match(names, rownames(meltedRes)),]
sumPre <- sum(meltedRes[,"counts"])
meltedRes$Fpkm <- (meltedRes[,"counts"]*1000000000)/(preLen*sumPre)
write.table(meltedRes, sep="\t", quote=FALSE)

unlink(orderedBams)
bai <- sub(".bam", ".bai", orderedBams)
unlink(bai)

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}


