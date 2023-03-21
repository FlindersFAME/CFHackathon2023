rule flye:
    input:
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_long.fastq")
    output:
        directory(os.path.join(FLYE,"{sample}")),
        os.path.join(FLYE,"{sample}", "assembly.fasta"),
        os.path.join(FLYE,"{sample}", "assembly_info.txt")
    conda:
        os.path.join('..', 'envs','assemble.yaml')
    resources:
        mem_mb=BigJobMem,
        time=BigTime
    params:
        FLYE_MODEL
    threads:
        BigJobCpu
    shell:
        """
        flye {params[0]} {input[0]} -t {threads}  --out-dir {output[0]}
        """

rule unicycler_hybrid:
    input:
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_long.fastq"),
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R1.fastq"),
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R2.fastq")
    output:
        directory(os.path.join(UNICYCLER_HYBRID,"{sample}")),
        os.path.join(UNICYCLER_HYBRID,"{sample}", "assembly.fasta")
    conda:
        os.path.join('..', 'envs','unicycler.yaml')
    resources:
        mem_mb=BigJobMem,
        time=BigTime
    threads:
        BigJobCpu
    shell:
        """
        unicycler -1 {input[1]}  -2 {input[2]}  -l {input[0]}  -o {output[0]} 
        """


rule unicycler_short:
    input:
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R1.fastq"),
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R2.fastq")
    output:
        directory(os.path.join(UNICYCLER_SHORT,"{sample}")),
        os.path.join(UNICYCLER_SHORT,"{sample}", "assembly.fasta")
    conda:
        os.path.join('..', 'envs','unicycler.yaml')
    resources:
        mem_mb=BigJobMem,
        time=BigTime
    params:
        FLYE_MODEL
    threads:
        BigJobCpu
    shell:
        """
        unicycler -1 {input[0]}  -2 {input[1]}    -o {output[0]} 
        """

rule aggr_assemble:
    input:
        expand(os.path.join(FLYE,"{sample}", "assembly.fasta"), sample = SAMPLES),
        expand(os.path.join(FLYE,"{sample}", "assembly_info.txt"), sample = SAMPLES),
        expand(os.path.join(UNICYCLER_SHORT,"{sample}", "assembly.fasta"), sample = SAMPLES),
        expand(os.path.join(UNICYCLER_HYBRID,"{sample}", "assembly.fasta"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_assemble.txt")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    shell:
        """
        touch {output[0]}
        """



