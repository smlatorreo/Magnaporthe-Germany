# Magnaporthe-Germany

Methods for the manuscript: "The blast fungus *Magnaporthe* (Syn. *Pyricularia*) occurs in Germany on wild grasses and has the capacity to infect cereal crops. AC. Barragán, SM. Latorre, A. Harant, J. Win, P. Mock, A. Malmgren, HA. Burbano, S. Kamoun, and T. Langner"

DOI: XX.XX.XX

## 1. Mapping, variant calling and genetic clustering of blast fungus isolates
Requirements:

Program                  | Location
------------------------ | ----------------------------
*AdapterRemoval v2*      | (https://github.com/mikkelschubert/adapterremoval)
*Bwa-mem2 v.2.1*         | (https://github.com/bwa-mem2/bwa-mem2)
*samtools v.1.11*        | (https://github.com/samtools/samtools)
*sambamba v0.8.0*        | (https://github.com/biod/sambamba)
*GATK v4.2*              | (https://github.com/broadinstitute/gatk/releases)
*bcftools v.1.11*        | (https://github.com/samtools/bcftools)

### 1.1. Alignment of short reads to the rice-infecting *M. oryzae* reference genome

Raw .fastq sequences were trimmed with *AdapterRemoval2*
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --gzip --basename $sample.trimmed
```

We used the rice-infecting *Magnaporthe oryzae* 70-15 assembly as the reference genome and indexed this genome using *Bwa-mem2*.
```bash
bwa index 70-15.fa
```

*BWA mem2* was used to map the trimmed reads to the reference genome, *samtools* to discard non-mapped reads, and *sambamba* to sort and mark PCR optical duplicates.
```bash
bwa-mem2 mem -R "@RG\tID:$sample\tSM:$sample" 70-15.fa $sample1.trimmed.R1.fastq.gz $sample1.trimmed.R2.fastq.gz > sample1.sam
samtools view -SbhF 4 > sample1_mapped.bam
sambamba sort -o sample1_mapped_sorted.bam sample1_mapped.bam
sambamba markdup sample1_mapped_sorted.bam sample1_mapped_sorted.dd.bam
```

### 1.2. Variant calling
We used the *HaplotypeCaller* from *GATK* to generate genomic haplotype calls per individual using the duplicate-marked BAM file as input.
```bash
gatk HaplotypeCaller -R 70-15.fa -I sample1_mapped_sorted.dd.bam -O sample1.g.vcf.gz
```

We used *CombineGVCFs*, *GenotypeGVCFs* and *SelectVariants* from *GATK* to combine the individual genomic VCFs, call genotypes and filter SNPs, respectively.
```bash
gatk CombineGVCFs -R 70-15.fa -V sample1.g.vcf.gz -V sample2.g.vcf.gz -V sampleN.g.vcf.gz -O blast.g.vcf.gz
gatk GenotypeGVCFs -R 70-15.fa -ploidy 1 -V blast.g.vcf.gz -O blast.raw.vcf.gz
gatk SelectVariants -select-type SNP -V blast.raw.vcf.gz -O blast.raw.snps.vcf.gz
```

We extracted all Quality-by-Depth (QD) values
```bash
bcftools view -H blast.raw.snps.vcf.gz | cut -f8 |
awk -F "QD=" '{print $2}' | cut -f1 -d ";" | gzip >  blast.raw.snps.QD.gz
```

Based on the distribution of Quality-by-Depth values, we set filters of one standard deviation around the median value (See [Latorre et al, 2022](https://doi.org/10.1101/2022.03.06.482794)).
```python
# Python
import pandas as pd
QD = pd.read_csv('blast.raw.snps.QD.gz', header = None, compression = 'gzip')
med = QD.median()
lower = med - QD.std()
upper = med + QD.std()
print(lower, upper)
```

Finally, using the above-mentioned scheme, we filtered SNPs using *GATK VariantFiltration* and created a new VCF file.
```bash
gatk VariantFiltration --filter-name "QD" \
--filter-expression "QD <= $lower || QD >= $upper" \
-V blast.raw.snps.QD.gz \
-O blast.snps.filter.vcf.gz
```
## 2. Reference-free clustering of blast fungus isolates based on K-mer sharing

Requirements:
Program    | Location
---------- | --------
Mash v.3.2 | https://github.com/marbl/Mash/releases/tag/v2.3
ape (R)    | http://ape-package.ird.fr/


Per-sample kmer sketches were created with Mash:
```bash
# Options
# -r : Assummes the input as a read set
# -g 43M : Assumes a genome size of 43 Megabases to calculate p-values
# -s 1000000 : Sketch size
# -m 3 : Minimum 3 copies per kmer are required
# -k 21 : Kmer size
# -I $SAMPLE : Sample ID
# -o $SAMPLE : Sample output prefix
# $SAMPLE.R1.fastq.gz : R1 pairwise raw read
# $SAMPLE.R2.fastq.gz : R2 pairwise raw read

mash sketch -r -g 43M -s 1000000 -m 3 -k 21 -I $SAMPLE -o $SAMPLE $SAMPLE.R1.fastq.gz $SAMPLE.R2.fastq.gz
```
The output files with the Mash sketches can be located [here](/data/Mash_kmer_sketches/)

We then computed pariwise mash distances
```bash
# Options
mash triangle -E *.msh > distances_mash.out
```
The pairwise distance files can be found in the [data folder](/data/)

Finally, we used NJ to cluster the samples based on Mash distances

```r
# Load pariwise distances (Mash output)
dat <- read.table('distances.mash.out', header = F)
isolates <- sort(unique(c(dat[,1], dat[,2])))
N <- length(isolates)
# Create a symmetric matrix
m <- matrix(ncol = N, nrow = N)
co <- combn(N, 2, simplify = T)
for(i in 1:ncol(co)){
  indA <- isolates[co[1,i]]
  indB <- isolates[co[2,i]]
  d <- dat[((dat[,1] == indA) & (dat[,2] == indB)) | ((dat[,1] == indB) & (dat[,2] == indA)), 3]
  m[co[1,i], co[2,i]] <- d
  m[co[2,i], co[1,i]] <- d
}
colnames(m) <- isolates
rownames(m) <- isolates

# Load ape package and compute a NJ tree
library(ape)
tree <- nj((m))
plot(tree)
```
The trees can be found in the [data folder](/data/)
