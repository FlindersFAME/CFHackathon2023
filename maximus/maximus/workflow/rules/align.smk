rule minimap_ont:
    input:
        get_input_lr_fastqs
    output:
        os.path.join(ALIGN,"{sample}_long.bam")
    threads:
        BigJobCpu
    resources:
        mem_mb=MassiveJobMem,
        time=BigTime
    params:
        Reference
    conda:
        os.path.join('..', 'envs','align.yaml')
    shell:
        '''
        minimap2 -ax map-ont -t {threads} {params[0]} {input[0]} | samtools view -@ {threads} -S -b > {output[0]}
        '''

rule bam_sort_ont:
    input:
        os.path.join(ALIGN,"{sample}_long.bam")
    output:
        os.path.join(ALIGN,"{sample}_sorted_long.bam")
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    conda:
        os.path.join('..', 'envs','align.yaml')
    shell:
        '''
        samtools sort -@ {threads} {input[0]} -o {output[0]}
        '''

rule bam_stats:
    input:
        os.path.join(ALIGN,"{sample}_sorted_long.bam")
    output:
        os.path.join(BAM_STATS,"{sample}_bam_long.stats")
    threads:
        BigJobCpu
    resources:
        mem_mb=SmallJobMem,
        time=MediumTime
    conda:
        os.path.join('..', 'envs','align.yaml')
    shell:
        '''
        samtools stats -@ {threads} {input[0]} | grep ^SN | cut -f 2- > {output[0]} 
        '''

rule bam_map_fastq_long:
    """converted mapped reads to fastq"""
    input:
        os.path.join(ALIGN,"{sample}_sorted_long.bam")
    output:
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_long.fastq")
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools view -F 4 -F 2048 -@ {threads} -h -b {input[0]} | samtools fastq - > {output[0]}
        """

rule minimap_short:
    input:
        get_input_r1_fastqs,
        get_input_r2_fastqs
    output:
        os.path.join(ALIGN,"{sample}_short.bam")
    threads:
        BigJobCpu
    resources:
        mem_mb=MassiveJobMem,
        time=BigTime
    params:
        Reference
    conda:
        os.path.join('..', 'envs','align.yaml')
    shell:
        '''
        minimap2 -ax sr -t {threads} {params[0]} {input[0]} {input[1]}  | samtools view -@ {threads} -S -b > {output[0]}
        '''

rule bam_sort_short:
    input:
        os.path.join(ALIGN,"{sample}_short.bam")
    output:
        os.path.join(ALIGN,"{sample}_sorted_short.bam")
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    conda:
        os.path.join('..', 'envs','align.yaml')
    shell:
        '''
        samtools sort -@ {threads} {input[0]} -o {output[0]}
        '''

rule bam_stats_short:
    input:
        os.path.join(ALIGN,"{sample}_sorted_short.bam")
    output:
        os.path.join(BAM_STATS,"{sample}_bam_short.stats")
    threads:
        BigJobCpu
    resources:
        mem_mb=SmallJobMem,
        time=MediumTime
    conda:
        os.path.join('..', 'envs','align.yaml')
    shell:
        '''
        samtools stats -@ {threads} {input[0]} | grep ^SN | cut -f 2- > {output[0]} 
        '''

rule bam_map_fastq_short:
    """converted mapped reads to fastq"""
    input:
        os.path.join(ALIGN,"{sample}_sorted_short.bam")
    output:
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R1.fastq"),
        os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R2.fastq")
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools fastq -f 3 -F 4 -1 {output[0]} -2 {output[1]} {input[0]}
        """

#### aggregation rule
rule align_aggr:
    """aggregate qc"""
    input:
        expand(os.path.join(ALIGNED_FASTQ,"{sample}_mapped_long.fastq"), sample = SAMPLES),
        expand(os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R1.fastq"), sample = SAMPLES),
        expand(os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R2.fastq"), sample = SAMPLES),
        expand(os.path.join(BAM_STATS,"{sample}_bam_short.stats"), sample = SAMPLES),
        expand(os.path.join(BAM_STATS,"{sample}_bam_long.stats"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "align.txt")
    threads:
        1
    shell:
        """
        touch {output[0]}
        """

#### aggregation rule
rule align_aggr_short:
    """aggregate qc"""
    input:
        expand(os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R1.fastq"), sample = SAMPLES),
        expand(os.path.join(ALIGNED_FASTQ,"{sample}_mapped_short_R2.fastq"), sample = SAMPLES),
        expand(os.path.join(BAM_STATS,"{sample}_bam_short.stats"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "align_short.txt")
    threads:
        1
    shell:
        """
        touch {output[0]}
        """


