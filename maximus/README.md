# maximus

snakemake-powered commandline to map ONT reads to a reference pan-genome and assemble 

```
pip install -e .
maximus --help
maximus run --help
```

* maximus requires an input csv file with 4 columns
* Each row is a sample
* Column 1 is the sample name, column 2 is the long read ONT fastq file format, columns 3 and 4 are the short read (MGI files)
* No headers

e.g.

sample1,sample1_long_read.fastq.gz, sample1_short_read_R1.fastq.gz, sample1_short_read_R2.fastq.gz
sample2,sample2_long_read.fastq.gz, sample2_short_read_R1.fastq.gz, sample2_short_read_R2.fastq.gz

If you specify ```--short_reads```, the required input sheet is 3 columns:

* maximus requires an input csv file with 3 columns
* Each row is a sample
* Column 1 is the sample name,  columns 2 and 3 are the short read (MGI files)
* No headers

e.g.

sample1, sample1_short_read_R1.fastq.gz, sample1_short_read_R2.fastq.gz
sample2, sample2_short_read_R1.fastq.gz, sample2_short_read_R2.fastq.gz

# Usage

Requires a reference (pan) genome to map the reads to passed as --reference. Otherwise just the input, output and threads


```
 maximus run --input test_pseudomonas.csv --output test_pseudomonas --reference Genome_Refs/PA01.fasta --threads 8
```


For short read only, specify --short_reads

```
 maximus run --input test_pseudomonas.csv --output test_pseudomonas --reference Genome_Refs/PA01.fasta --threads 8 --short_reads
```
