#!/usr/bin/env Rscript



#-#-#-#-#-#-#-#-#-#-#-#-#
# Argument parsing
#-#-#-#-#-#-#-#-#-#-#-#-#
library(optparse)
option_list <- list( 
  make_option("--min-qs",help="Minimum average quality scores of reads.",type="numeric",default = 20),
  make_option("--min-len",help="Minimum read length.",type="integer",default = 300),
  make_option("--out",help="Name of the bam file to generate [required]",type="character")
)
opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "usage: %prog --out <output-bam-file> <input-bam-file>",
  description = "Filter out reads in a bam file that are too short or where average quality score is too low"
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

# Define filtering rules
rules <- FilterRules(list(
  minqs = function(x) {mean(as(x$qual,"IntegerList")) >= opt$options$"min-qs"},
  minlen = function(x) {nchar(x$qual) >= opt$options$"min-len"}
))

# Apply the filters
filterBam(
  BamFile(opt$args,asMates = FALSE,yieldSize = 1e5),
  filter = rules,destination = opt$options$out,
  param = ScanBamParam(what="qual")
)






