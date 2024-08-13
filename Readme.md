

# Disclaimer

```
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
```

# Usage

## In the cloud

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/BioinfoSupport/virotyper.git/HEAD?urlpath=rstudio)

[![Binder](http://mybinder.org/badge_logo.svg)](http://mybinder.org/v2/gh/BioinfoSupport/virotyper.git/HEAD?urlpath=shiny/app)


# Locally

```bash
repo2docker --Repo2Docker.platform=linux/amd64 ./
```

```bash
# Run a bash in the container
docker run --rm -it -v $(pwd):/data --workdir /data unigebsp/virotyper bash


# Clean reporting files in the given directory
make data/HHV1/bam/clean
make data/HHV1/fasta/clean
make data/HHV2/bam/clean
make data/HHV2/fasta/clean

# Run variant calling from a single existing BAM file
make data/HHV1/bam/240222-HE.HHV1.filt.10k.report
make DB_DIR=data/HHV2/db data/HHV2/bam/240222-RO.HHV2.filt.10k.report

# Run variant calling for all BAM files in a directory
ls data/HHV1/bam/*.bam | grep -v filt | sed s/.bam$/.filt.10k.report/ | xargs make -kj6
ls data/HHV2/bam/*.bam | grep -v filt | sed s/.bam$/.filt.10k.report/ | xargs make -j6 DB_DIR=data/HHV2/db

# Run variant calling for all FASTA files in a directory
ls data/HHV1/fasta/*.fasta | sed s/$/.all/ | xargs make -kj6
ls data/HHV2/fasta/*.fasta | sed s/.fasta$/.asm.report/ | xargs make -kj6 DB_DIR=data/HHV2/db

# Generate reports for all VCF files in a directory
ls data/HHV1/vcf/*.vcf.gz data/HHV1/vcf/*.bam | sed s/$/.all/ | xargs make -kj6
```




# Build the image
```bash
docker build -t unigebsp/virotyper ./
```


# Parking



```bash
# Query mutations of a given strain in the database
bcftools view --private --samples 'HSV1-UL30:L702H:1:T>A' /repository/pipeline/ref/hsv1.db.vcf.bgz
bcftools view --private --samples 'HSV1-UL30:A910V:1:C>T' /repository/pipeline/ref/hsv1.db.vcf.bgz

# Generate some consensus sequence from the DB to test the pipeline
bcftools consensus --sample 'HSV1-UL30:L702H:1:T>A' --fasta /repository/pipeline/ref/hsv1.fasta /repository/pipeline/ref/hsv1.db.vcf.bgz > out/tests/consensus01.fasta


# Merge multiple VCF into a single one
bcftools merge -Oz --force-samples out/tests/{consensus01,consensus02}.hsv1.mm2.csq.ann.vcf.bgz > out/ab.hsv1.mm2.csq.ann.vcf.bgz
bcftools index out/ab.hsv1.mm2.csq.ann.vcf.bgz
bcftools view out/ab.hsv1.mm2.csq.ann.vcf.bgz
make -f /repository/pipeline/Makefile out/ab.hsv1.mm2.csq.ann.csq.ann.vcf.bgz
bcftools query -f'[%CHROM\t%POS\t%SAMPLE\t%TBCSQ\n]' out/ab.hsv1.mm2.csq.ann.csq.ann.vcf.bgz
```


## IGV

Example IGV script to generate screenshots of given locus
```bash
./IGV_2.16.2/igv.sh --genome ./pipeline/ref/hsv1.fasta --batch /dev/stdin <<EOF
load ./pipeline/ref/hsv1.gff
load data/HHV1_UL23_UL30_231108.hsv1.mm2.bam
colorBy READ_STRAND
goto JN555585.1-UL30:2105 JN555585.1-UL30:2729
expand hsv1.gff
snapshot out/snapshot.png
EOF
```



