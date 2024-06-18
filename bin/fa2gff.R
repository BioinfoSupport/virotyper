#!/usr/bin/env Rscript


#-#-#-#-#-#-#-#-#-#-#-#-#
# Argument parsing
#-#-#-#-#-#-#-#-#-#-#-#-#
library(optparse)
option_list <- list( 
  make_option("--out",help="Path to the output GFF file [required]",type="character")
)
opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "usage: %prog --out <output-gff-file> <input-fasta-file>",
  description = "Generate a GFF compatible with 'bcftools csq' from a FASTA file of genes sequences."
),positional_arguments = 1)
if (is.null(opt$options$out)) stop("--out argument is required")



#-#-#-#-#-#-#-#-#-#-#-#-#
# Main
#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages({
  library(Rsamtools)
  library(Biostrings)
  library(rtracklayer)
})


# Generate FAI
indexFa(FaFile(opt$args))

# Load the FASTA containing dna sequences of genes
fa <- readDNAStringSet(opt$args)

# Generate features on each sequence
gene <- GRanges(seqinfo(fa),type="gene",strand="+",biotype="protein_coding")
gene$Name <- as.character(seqnames(gene))
gene$ID <- paste0("gene:",gene$Name)

tx <- GRanges(seqinfo(fa),type="transcript",strand="+",biotype="protein_coding")
tx$Name <- as.character(seqnames(tx))
tx$ID <- paste0("transcript:",tx$Name)
tx$Parent <- paste0("gene:",tx$Name)

exon <- GRanges(seqinfo(fa),type="exon",strand="+",biotype="protein_coding")
exon$Name <- as.character(seqnames(exon))
exon$ID <- paste0("exon:",exon$Name)
exon$Parent <- paste0("transcript:",exon$Name)


cds <- GRanges(seqinfo(fa),type="CDS",strand="+",biotype="protein_coding")
cds$Name <- as.character(seqnames(cds))
cds$ID <- paste0("CDS:",cds$Name)
cds$Parent <- paste0("transcript:",cds$Name)
cds$phase <- 0


# Concatenate features and export the GFF
gff <- c(gene,tx,exon,cds)
stopifnot(all(!duplicated(gff$ID)))

seqinfo(gff) <- seqinfo(fa)
rtracklayer::export.gff3(gff,opt$options$out)


