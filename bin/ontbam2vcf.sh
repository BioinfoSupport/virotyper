#!/usr/bin/env bash

function usage() {
  cat <<-EOF
    USAGE: $0 [OPTIONS]

    DESCRIPTION
      Wrapper arround BCFTOOLS to perform variant calling on a BAM file 
      containing an alignment of ONT reads on a reference sequence.
      The BAM file may be generated with command:
      minimap2 -cx ont -z1000000 -p 0.1 -Y --cs -a <ref.fasta> <assembly.fasta> | samtools sort -o <alignment.bam> -

    OPTIONS
      --ref-fasta   Path to FASTA containing the reference sequence
      --ref-gff     Path to GFF containing annotations compatible with "bcftools csq"
      --bam         Path to BAM containing aligned assembly
      --out         Path to the generated .vcf.gz file
EOF
}


# Parse arguments
while getopts "h-:" opt; do
  case ${opt} in
    h) usage; exit 2;;
    -) case "${OPTARG}" in
        ref-fasta=*) REF_FASTA=${OPTARG#*=};;
        ref-gff=*) REF_GFF=${OPTARG#*=};;
        bam=*) BAM=${OPTARG#*=};;
        out=*) OUT=${OPTARG#*=};;
        *) usage; exit 1;;
       esac;;
   esac
done


if ! [ -f "$BAM" ] ; then echo "--bam missing or file doesn't exists"; exit 1; fi
if ! [ -f "$REF_FASTA" ] ; then echo "--ref-fasta missing or file doesn't exists"; exit 1; fi
if ! [ -f "$REF_GFF" ] ; then echo "--ref-gff missing or file doesn't exists"; exit 1; fi
if [ -z "$OUT" ] ; then echo "--out is missing"; exit 1; fi


# Run command
bcftools mpileup -Ou --max-depth 10000 --no-BAQ -Q5 --max-BQ 30 --indel-size 1000 \
  	--ff=UNMAP --ambig-reads=drop --gap-frac=0 --min-ireads=0 --per-sample-mF \
    -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR \
  	--fasta-ref="$REF_FASTA" "$BAM" \
	| bcftools filter -Ou -i 'DP>=20 && (I16[2] > DP*0.1) && (I16[3] > DP*0.1)' \
	| bcftools norm -Ou -m- -f "$REF_FASTA" \
  | bcftools filter -Ou -i '(FORMAT/ADF[:1] >= DP*0.1) && (FORMAT/ADR[:1] >= DP*0.1)' \
  | bcftools +fill-tags -Ou -- -t 'INFO/AF:1=Float(FORMAT/AD[0:1]/DP)' \
  | bcftools csq -Oz --local-csq --write-index --force -o "$OUT" \
      --fasta-ref="$REF_FASTA" \
      --gff="$REF_GFF"






