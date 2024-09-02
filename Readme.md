

# Disclaimer

```
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
```

# Usage

- prerequisite: a resistance database must exist in folder `data/<DBID>/db`, see section below for the format.


```bash
# Build the docker image
docker build -t virotyper ./

# Run a bash in the container
docker run --rm -it -v $(pwd):/data --workdir /data virotyper bash

# Generate a report for a VCF/BAM/FASTA file
make DB_ID=hsv2 <filename>.vcf.all
make DB_ID=hsv2 <filename>.bam.all
make DB_ID=hsv1 <filename>.fasta.all
```



# Resistance database

A resistance database is a subdirectory located in `data/db/<DB_ID>/` with the following content:
```
data/db/hsv1/
  - ref.fasta         reference sequence(s)
  - ref.gff           reference annotation (optional: generated from ref.fasta if missing, assuming each sequence is a gene)
  - resistances.xlsx  an excel containing known resitances   
  - template.docx     a template word document used to generate the resistance report
```


`resistances.xlsx` must be an excel file containing a sheet named `resistances` with at least columns:

 - `strain_id`: format is `UL23-HHV2:A19V` to encode a mutation on the 19th amino-acid of gene `UL23-HHV2` which is `A` in the reference genome, and became a `V` in the mutated strain. `UL23-HHV2:t42a` encode a DNA mutation at nucleotide 42 which is `t` in the reference and become a `a` in the mutated strain. To encode the insertion of a `g` at position 366 on gene `UL23-HHV2` the format is `UL23-HHV2:t365tg` (Note: the mutation refer to nucleotide before). Similarly to encode the deletion of sequence `ac` at position 862 on gene `UL23-HHV2`, the format is: `UL23-HHV2:cac861c`.
 - `mutation_type`: Type of mutation (Either: Natural, Resistant, Susceptible)
 - `drug`: Name of the drug the mutation as an impact on (missing if mutation type is Natural).
 
The same mutation can appear several time in the file if it confere a resistance to multiple drugs.





# Parking

```bash
# Clean reporting files in the given directory
make data/HHV1/bam/clean
make data/HHV1/fasta/clean
make data/HHV2/bam/clean
make data/HHV2/fasta/clean

# Run variant calling from a single existing BAM file
make data/HHV1/bam/240222-HE.HHV1.filt.10k.bam.all
make DB_ID=hsv2 data/HHV2/bam/240222-RO.HHV2.filt.10k.bam.all

# Run variant calling for all BAM files in a directory
ls data/HHV1/bam/*.bam | grep -v filt | sed s/.bam$/.filt.10k.bam.all/ | xargs make
ls data/HHV1/bam/*.filt.10k.bam | sed s/.bam$/.vcf.gz.all/ | xargs make
ls data/HHV2/bam/*.bam | grep -v filt | sed s/.bam$/.filt.10k.bam.all/ | xargs make DB_DIR=data/HHV2/db

# Run variant calling for all FASTA files in a directory
ls data/HHV1/fasta/*.fasta | sed s/$/.all/ | xargs make
ls data/HHV2/fasta/*.fasta | sed s/$/.all/ | xargs make DB_DIR=data/HHV2/db

# Generate reports for all VCF files in a directory
ls data/HHV1/vcf/*.vcf.gz data/HHV1/vcf/*.bam | sed s/$/.all/ | xargs make -kj6
```



