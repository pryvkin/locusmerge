locusmerge
==========

An algorithm for clustering multi-mapping RNA-seq reads into representative loci

## Installing dependencies
#### samtools
```
git clone https://github.com/samtools/samtools -b standalone
cd samtools
make
cd ..
```
#### bedtools
```
git clone https://github.com/arq5x/bedtools2
cd bedtools2
make 
sudo make install
```
#### R packages:
In R:
```
install.packages(c("igraph", "fastcluster"))
```

## Building
```
make
```
## Usage
### Building genome annotations
First, you need a genome annotation:
```
./locusmerge_build_annot_hsa
./locusmerge_build_annot_ath
```

### Running locusmerge
Example for H. sapiens, output being placed in data/out*:
```
./locusmerge in.bam hsa data/out
```

Example for A. thaliana with parallelized sorting:
```
./locusmerge --ncpu=3 in.bam ath data/out
```

To see more options:
```
./locusmerge
```

## Output

Reads that map to clustered loci:

*.bamfile.bam.clustered

Clustered loci and information about each cluster:

*.bed.clustered
*.clusters

