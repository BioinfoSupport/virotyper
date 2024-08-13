#!/usr/bin/env bash

function usage() {
  cat <<-EOF
    USAGE: $0 [OPTIONS]

    DESCRIPTION
      Wrapper arround MINIMAP2 + BCFTOOLS to perform variant calling between
      an assembly and a reference sequence. The given assembly sequences are 
      mapped on the reference with minimap2 and the alignment is stored into 
      <out-prefix.bam>. Then variant calling is performed with bcftools and stored
      into <out-prefix.vcf.gz>.

    OPTIONS
      --ref-fasta   Path to FASTA containing the reference sequence
      --ref-gff     Path to GFF containing annotations compatible with "bcftools csq"
      --asm-fasta   Path to FASTA containing the assembly sequence
      --out-prefix  Path to the generated .vcf.gz file
EOF
}


# Parse arguments
while getopts "h-:" opt; do
  case ${opt} in
    h) usage; exit 2;;
    -) case "${OPTARG}" in
        ref-fasta=*) REF_FASTA=${OPTARG#*=};;
        ref-gff=*) REF_GFF=${OPTARG#*=};;
        asm-fasta=*) ASM_FASTA=${OPTARG#*=};;
        out-prefix=*) OUT_PREFIX=${OPTARG#*=};;
        *) usage; exit 1;;
       esac;;
   esac
done


if ! [ -f "$ASM_FASTA" ] ; then echo "--asm-fasta missing or file doesn't exists"; exit 1; fi
if ! [ -f "$REF_FASTA" ] ; then echo "--ref-fasta missing or file doesn't exists"; exit 1; fi
if ! [ -f "$REF_GFF" ] ; then echo "--ref-gff missing or file doesn't exists"; exit 1; fi
if [ -z "$OUT_PREFIX" ] ; then echo "--out-prefix is missing"; exit 1; fi


# Run alignment command
minimap2 -cx asm5 \
  -z1000000 -p 0.1 -Y \
  --cs -a "${REF_FASTA}" "${ASM_FASTA}" \
  | samtools sort -o "${OUT_PREFIX}.bam" -

# Index the bam file
samtools index "${OUT_PREFIX}.bam"

# Run variant calling command
bcftools mpileup -Ou --max-depth 10000 --no-BAQ --indel-size 1000 \
  	--ff=UNMAP --ambig-reads=drop --gap-frac=0 --min-ireads=0 --per-sample-mF \
    -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR \
  	--fasta-ref="${REF_FASTA}" "${OUT_PREFIX}.bam" \
	| bcftools norm -Ou -m- -f "${REF_FASTA}" \
  | bcftools filter -Ou -i '(FORMAT/AD[:1] >= 1)' \
  | bcftools +fill-tags -Ou -- -t 'INFO/AF:1=Float(FORMAT/AD[0:1]/DP)' \
  | bcftools csq -Oz --local-csq --write-index --force -o "${OUT_PREFIX}.vcf.gz" \
      --fasta-ref="${REF_FASTA}" \
      --gff="${REF_GFF}"





