#!/usr/bin/env Rscript
suppressWarnings(library('picante'))

in_path = '/PATH/TO/THE/ASVs/TABLE' #Edit your PATH here!
###NOTE:
###The table should be in the format: Features(rows)*Samples(columns)
###All values should be numerical
######## functions #########

alpha <- function(x, tree = NULL, base = exp(1)) {
        require('picante')
        est <- estimateR(x)
        Richness <- est[1, ]
        Chao1 <- est[2, ]
        ACE <- est[4, ]
        Shannon <- diversity(x, index = 'shannon', base = base)
        Simpson <- diversity(x, index = 'simpson')
        Shannoneven <- Shannon / log(Richness, base)
        goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
        result <- data.frame(Richness, Shannon, Shannoneven, Simpson, Chao1, ACE, goods_coverage)
        return(result)
}

otu <- read.delim(in_path, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu<-t(otu)
alpha_all <- alpha(otu, treefile, base = 2)
