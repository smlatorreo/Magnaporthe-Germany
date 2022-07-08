# Magnaporthe-Germany

Methods for the manuscript: "The blast fungus *Magnaporthe* (Syn. *Pyricularia*) occurs in Germany on wild grasses and has the capacity to infect cereal crops. AC. Barragán, SM. Latorre, A. Harant, J. Win, P. Mock, A. Malmgren, HA. Burbano, S. Kamoun, and T. Langner"

DOI: XX.XX.XX

## 1. Mapping, variant calling and genetic clustering of blast fungus isolates
Requirements:
Program    | Location
---------- | --------
bwa


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
