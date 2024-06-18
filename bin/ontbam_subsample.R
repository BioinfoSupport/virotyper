#!/usr/bin/env Rscript



#-#-#-#-#-#-#-#-#-#-#-#-#
# Argument parsing
#-#-#-#-#-#-#-#-#-#-#-#-#
library(optparse)
option_list <- list( 
  make_option("--target-cov",help="Subsample reads to the given target coverage.",type="numeric",default = 10000),
  make_option("--out",help="Name of the bam file to generate [required]",type="character"),
  make_option("--random-seed",help="Random seed to use when subsampling",type="integer",default = 12345)
)
opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "usage: %prog --out <output-bam-file> <input-bam-file>",
  description = "Randomly subsample reads in a BAM file to a given target coverage along each chromosome"
),positional_arguments = 1)
if (is.null(opt$options$out)) stop("--out argument is required")

# Testing parameters
#opt <- list(args = c("data/tests/bam/240222-HE.HHV1.bam"),options = list("min-len" = 300, "min-qs" = 20, "target-cov"=10000))


#-#-#-#-#-#-#-#-#-#-#-#-#
# Methods definitions
#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages({
  #BiocParallel::register(BiocParallel::MulticoreParam(workers=4))
  BiocParallel::register(BiocParallel::SerialParam(progressbar = TRUE))
  library(GenomicAlignments)
})





set.seed(opt$options$"random-seed")
bf <- BamFile(opt$args,asMates = FALSE,yieldSize = 1e5)
sbparam <- ScanBamParam(what=c("rname"),flag = scanBamFlag(isSecondaryAlignment = FALSE,isUnmappedQuery = FALSE))


# Compute coverage per contig
cov <- mean(coverage(bf,param=sbparam))
cov_pval <- opt$options$"target-cov"/cov

# Define filtering rules
rules <- FilterRules(list(
  maxcov = function(x) {runif(length(x$rname),0,1) <= unname(cov_pval)[x$rname]} # Generate a random number and 
))

# Apply the filters
filterBam(bf,filter = rules,destination = opt$options$out,param = sbparam)






