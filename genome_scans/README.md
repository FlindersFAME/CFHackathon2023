# Compare genome scans

This is some straightforward ways to test the presence of a genome in the metagenomes. There are probably better ways to do this (see [singlem](https://wwood.github.io/singlem/)) but this works quite effectively.

1. Generate a set of fasta files you want to compare to.

For my tests, I had a set of genomes for _Pseudomonas_, _Mycobacterium abscessus_ complex, _Staphylococcus aureus_, and _Streptococcus salivarus_. I also have a fasta file of 16S sequences. I used a somewhat old one from the [RDP](http://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz), mostly because it is in fasta format, while GreenGenes and NCBI data was not. The 16S hasn't changed that much in the last couple of years.

2. Use `minimap2` to map all the reads. I used this command


```bash
SAMPLE=$1
minimap2 -t 8 -a --sam-hit-only -xsr  Genomes/StrepSalivarus.fna.gz fastq/${SAMPLE}_R1.fastq.gz fastq/${SAMPLE}_R2.fastq.gz   | samtools view -@ 8 -bh | samtools sort -o bam/${SAMPLE}.Strep.bam -
minimap2 -t 8 -a --sam-hit-only -xsr Genomes/Staph_aureus.fna.gz fastq/${SAMPLE}_R1.fastq.gz fastq/${SAMPLE}_R2.fastq.gz   | samtools view -@ 8 -bh | samtools sort -o bam/${SAMPLE}.Staph.bam -
minimap2 -t 8 -a --sam-hit-only -xsr Mycobacteria/myco100.fna.gz fastq/${SAMPLE}_R1.fastq.gz fastq/${SAMPLE}_R2.fastq.gz   | samtools view -@ 8 -bh | samtools sort -o bam/${SAMPLE}.Myco100.bam -
minimap2 -t 8 -a --sam-hit-only -xsr Genomes/Pseudomonas_aeruginosa_noNs.fasta.gz fastq/${SAMPLE}_R1.fastq.gz fastq/${SAMPLE}_R2.fastq.gz   | samtools view -@ 8 -bh | samtools sort -o bam/${SAMPLE}.P_aeruginosa.bam -
minimap2 -t 8 -a --sam-hit-only -xsr rdp/current_Bacteria_unaligned.fa.gz fastq/${SAMPLE}_R1.fastq.gz fastq/${SAMPLE}_R2.fastq.gz   | samtools view -@ 8 -bh | samtools sort -o bam/${SAMPLE}.16S.bam -
```

3. Extract read counts or sequences that map in fasta format.

Now you can provide `CFHackathon2023/genome_scans/bam_to_reads.py` with a file to either get the reads from (`-b`) or one or more bam files to ignore reads in this file (`-n`). In other words, if you extract all the reads that match just _Staphylococcus_ (for example), you can ignore the reads that match the 16S, transposons, and so on.

```
python CFHackathon2023/genome_scans/bam_to_reads.py -p -b bam/1651490_20171010_S.Myco100.bam -n bam/1651490_20171010_S.Staph.bam -n bam/1651490_20171010_S.Strep.bam -n bam/1651490_20171010_S.P_aeruginosa.bam -n bam/1651490_20171010_S.16S.bam
```
