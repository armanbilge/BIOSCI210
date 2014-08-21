#!/usr/bin/env Rscript

library("HWEintrinsic")

b0 <- c(33, 24, 5)
b1 <- c(32, 24, 5)
g0 <- c(35, 5, 21)
g1 <- c(23, 29, 9)
loci.genotype.counts <- list(b0, b1, g0, g1)
loci.ids <- c("blue gen. 0", "blue gen. 1", "green gen. 0", "green gen. 1")

HWE.expected <- function(genotypes) {
    N = sum(genotypes)
    p = (2 * genotypes[1] + genotypes[2]) / (2 * N)
    q = 1 - p
    return(c(p^2, 2 * p * q, q^2) * N)
}

for (i in 1:length(loci.ids)) {
    print(loci.ids[i])
    dataset <- new("HWEdata", data = loci.genotype.counts[[i]])
    print(HWE.expected(dataset@data.vec))
    print(hwe.bf(dataset))
    print(hwe.ibf(dataset, sum(dataset@data.vec)))
    print("")
}
