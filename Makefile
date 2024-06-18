

DB_DIR ?= ./data/HHV1/db
MINIMAP2_FLAGS = -z1000000 -p 0.1

# Main rule
%.report:%.bam
	$(MAKE) \
	  $*.vcf.gz \
	  $*.report.html

.PRECIOUS: %.asm.bam %.asm.bcf %.filt.bam %.10k.bam %.bam.bai %.bcf %.bcf.txt

# BAM indexing and ONT read filtering rules
%.bam.bai:%.bam; samtools index "$<"
%.filt.bam:%.bam; ./bin/ontbam_filter.R --out "$@" "$<"
%.10k.bam:%.bam; ./bin/ontbam_subsample.R --out "$@" "$<"



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 1a - Alignment and variant calling from a FASTA containing an assembly
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
%.asm.bam:%.fasta; minimap2 -cx asm5 $(MINIMAP2_FLAGS) -Y --cs -a "$(DB_DIR)/ref.fasta" "$<" | samtools sort -o "$@" -

%.asm.vcf.gz %.asm.vcf.gz.csi:%.asm.bam %.asm.bam.bai
	./bin/asm2vcf.sh --ref-gff="$(DB_DIR)/ref.gff" --ref-fasta="$(DB_DIR)/ref.fasta" --bam="$<" --out="$@"

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 1b - Variant calling from a BAM generated by minimap2 and containing aligned reads
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
%.vcf.gz %.vcf.gz.csi:%.bam %.bam.bai
	./bin/ontbam2vcf.sh --ref-gff="$(DB_DIR)/ref.gff" --ref-fasta="$(DB_DIR)/ref.fasta" --bam="$<" --out="$@"

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 4 - Generate reports
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
%.report.html:$(DB_DIR)/resistances.chk.xlsx %.bam %.bam.bai %.vcf.gz
	Rscript -e 'rmarkdown::render("scripts/02_report.Rmd",params=list(input_bam_file="$*.bam",input_vcf_file="$*.vcf.gz",input_db_dir="$(DB_DIR)",output_docx_report="$*.report.docx"),knit_root_dir="$(PWD)",output_dir="$(@D)",output_file="$(@F)")'


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Clean rule
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
%/clean:
	rm -f \
	  '$*'/*.filt.bam \
	  '$*'/*.10k.bam \
	  '$*'/*.asm.bam \
	  '$*'/*.bai \
    '$*'/*.filt.10k.vcf.gz \
    '$*'/*.filt.10k.vcf.gz.csi \
    '$*'/*.asm.vcf.gz \
    '$*'/*.asm.vcf.gz.csi \
	  '$*'/*.report.html \
	  '$*'/*.report.docx


