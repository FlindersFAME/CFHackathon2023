"""
All target output files are declared here
"""

if SHORT_READS == True:
    TargetFiles = [
    os.path.join(FLAGS, "align_short.txt"),
    os.path.join(FLAGS, "aggr_assemble_short.txt")
    ]
else:
    TargetFiles = [
    os.path.join(FLAGS, "align.txt"),
    os.path.join(FLAGS, "aggr_assemble.txt")
    ]
