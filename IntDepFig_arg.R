#!/usr/bin/R

library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(ggbio)

setwd("~/*")

args<-commandArgs(TRUE)
Int <-toString(args)

# Int <- "2369CxCML341C_ZmPR1.bed"
Comp <- substr(Int,1,nchar(Int)-4)
Comp2 <- (paste(Comp,"C.bed",sep=""))

fig_title <- substr(Int,1,nchar(Int)-4)
marker <- (paste(substr(fig_title,1,nchar(fig_title)-6),"_convert.bed",sep=""))
marker_file <- (paste("/Users/Todd/Dropbox/GBSThesis/NILAS_Final/BED/convert/chr_bed/",marker,sep=""))
group <- lapply(strsplit(fig_title, "_"), `[`, 1)
QTL <- lapply(strsplit(fig_title, "_"), `[`, 2)
rec <- lapply(strsplit(fig_title, "x"), `[`, 1)
donor <- lapply(strsplit(fig_title, "x"), `[`, 2)
donor2 <- substr(donor,1,nchar(donor)-6)
fig_title2 <- (paste(donor2,"x",rec,QTL,"Subgroup"))


avs.file <- marker_file
 avs.data <- read.table(avs.file,header=F,sep="\t",stringsAsFactors=F)
 avs.granges <- GRanges(seqnames=avs.data[,1],
                         ranges=IRanges(start=avs.data[,2],
                                       end=avs.data[,3]),
                              strand="*")
                              
avs3.file <- "ZMPR.bed"
avs3.data <- read.table(avs3.file,header=F,sep="\t",stringsAsFactors=F)
avs3.granges <- GRanges(seqnames=avs3.data[,1],
                        ranges=IRanges(start=avs3.data[,2],
                                      end=avs3.data[,3]),
                              strand="*")
                              
avs4.file <- Int
avs4.data <- read.table(avs4.file,header=F,sep="\t",stringsAsFactors=F)
avs4.granges <- GRanges(seqnames=avs4.data[,1],
                        ranges=IRanges(start=avs4.data[,2],
                                      end=avs4.data[,3]),
                              strand="*")
                              
avs5.file <- Comp2
avs5.data <- read.table(avs5.file,header=F,sep="\t",stringsAsFactors=F)
avs5.granges <- GRanges(seqnames=avs5.data[,1],
                        ranges=IRanges(start=avs5.data[,2],
                                      end=avs5.data[,3]),
                              strand="*")
                              
Combine <- rbind(avs4.data,avs5.data)
Max_Int <- max(Combine[,4], na.rm = TRUE)
                              
mark <- autoplot(avs.granges, layout = "karyogram", cytoband = FALSE, alpha=.5)

title = "Introgression\nDepth"

Z <- mark + layout_karyogram(avs4.granges, geom = "rect", ylim = c(11, 21), aes(fill=avs4.data[,4]),color = "black") + layout_karyogram(avs5.granges, geom = "rect", ylim = c(22, 32), aes(fill=avs5.data[,4]),color = "black") + layout_karyogram(avs3.granges,geom = "rect", ylim = c(11, 21), color = "green", fill=NA) + scale_fill_gradient(title,low = "yellow", high = "red",limits=c(0,Max_Int)) + ggtitle(fig_title2)

ggsave(filename=(paste(fig_title2,".png",sep="")))
