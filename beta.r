#!/usr/bin/env Rscript
suppressWarnings(library('picante'))

in_path = '/PATH/TO/THE/ASVs/TABLE' #Edit your PATH here!
###NOTE:
###The table should be in the format: Features(rows)*Samples(columns)
###All values should be numerical

otu <- read.delim(in_path, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu<-t(otu)
sample.dist<-vegdist(numOtudf,method = "bray")
####################################### NMDS #################################
# 计算nmds
otu.nmds<-metaMDS(sample.dist)
otu.stress<-otu.nmds$stress
otu.nmds.points<-otu.nmds$points
otu.nmds.points<-data.frame(Samples=rownames(otu.nmds.points),otu.nmds.points)
