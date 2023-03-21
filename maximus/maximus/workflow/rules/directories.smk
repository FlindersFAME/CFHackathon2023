"""
Database and output locations for the pipeline
"""

### reference directory
#ReferenceDir = config["referenceDir"]

### OUTPUT DIRs
FLAGS = os.path.join(OUTPUT, 'FLAGS')
PROCESSING = os.path.join(OUTPUT, 'PROCESSING')
RESULTS = os.path.join(OUTPUT, 'RESULTS')
ASSEMBLIES = os.path.join(OUTPUT, 'ASSEMBLIES')
FLYE = os.path.join(ASSEMBLIES, 'FLYE')
UNICYCLER_HYBRID = os.path.join(ASSEMBLIES, 'UNICYCLER_HYBRID')
UNICYCLER_SHORT = os.path.join(ASSEMBLIES, 'UNICYCLER_SHORT')




#align.smk 
ALIGN = os.path.join(PROCESSING, 'ALIGN')
BAM_STATS = os.path.join(OUTPUT, 'BAM_STATS')
ALIGNED_FASTQ =  os.path.join(OUTPUT, 'ALIGNED_FASTQ')




