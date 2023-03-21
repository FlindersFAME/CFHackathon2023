# maximus

snakemake-powered commandline to map ONT reads to a reference pan-genome and assemble 

```
pip install -e .
maximus --help
maximus run --help
```

* maximus requires an input csv file with 4 columns
* Each row is a sample
* Column 1 is the sample name, column 2 is the long read ONT fastq file in *fastq.gz format, columns 3 and 4 are the short read (MGI files)
* No headers

e.g.

sample1,sample1_long_read.fastq.gz, sample1_short_read_R1.fastq.gz, sample1_short_read_R2.fastq.gz
sample2,sample2_long_read.fastq.gz, sample2_short_read_R1.fastq.gz, sample2_short_read_R2.fastq.gz