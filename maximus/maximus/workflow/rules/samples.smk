"""
Function for parsing the 'Reads' config and identifying samples and read files
"""


def samplesFromCsv(csvFile):
    """
    Read samples and files from CSV 
    4 cols
    1 = sample 
    2 = Long read fastq file
    3 = Short read fastq R1
    4 = Short read fastq R2
    """
    outDict = {}
    with open(csvFile,'r') as csv:
        for line in csv:
            l = line.strip().split(',')
            if len(l) == 4:
                outDict[l[0]] = {}
                if os.path.isfile(l[1]) and os.path.isfile(l[2]) and os.path.isfile(l[3]) :
                    outDict[l[0]]['LR'] = l[1]
                    outDict[l[0]]['R1'] = l[2]
                    outDict[l[0]]['R2'] = l[3]
                else:
                    sys.stderr.write("\n"
                                     f"    FATAL: Error parsing {csvFile}. {l[1]} or {l[2]} or {l[3]} does not exist. \n"
                                     "    Check formatting, and that \n" 
                                     "    file names and file paths are correct.\n"
                                     "\n")
                    sys.exit(1)
    return outDict


def samplesFromCsvShort(csvFile):
    """
    Read samples and files from CSV 
    3 cols
    1 = sample 
    2 = Short R1
    3 = Short R2
    """
    outDict = {}
    with open(csvFile,'r') as csv:
        for line in csv:
            l = line.strip().split(',')
            if len(l) == 3:
                outDict[l[0]] = {}
                if os.path.isfile(l[1]) and os.path.isfile(l[2])  :
                    outDict[l[0]]['R1'] = l[1]
                    outDict[l[0]]['R2'] = l[2]
                else:
                    sys.stderr.write("\n"
                                     f"    FATAL: Error parsing {csvFile}. {l[1]} or  {l[2]}  does not exist. \n"
                                     "    Check formatting, and that \n" 
                                     "    file names and file paths are correct.\n"
                                     "\n")
                    sys.exit(1)
    return outDict


def parseSamples(csvfile):
    # for reading from directory
    #if os.path.isdir(readFileDir):
    #   sampleDict = samplesFromDirectory(readFileDir)
    if os.path.isfile(csvfile):
        sampleDict = samplesFromCsv(csvfile)
    else:
        sys.stderr.write("\n"
                         f"    FATAL: {csvfile} is neither a file nor directory.\n"
                         "\n")
        sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any samples at all.\n"
                         "\n")
        sys.exit(1)
    return sampleDict

def parseSamplesShort(csvfile):
    # for reading from directory
    #if os.path.isdir(readFileDir):
    #   sampleDict = samplesFromDirectory(readFileDir)
    if os.path.isfile(csvfile):
        sampleDict = samplesFromCsvShort(csvfile)
    else:
        sys.stderr.write("\n"
                         f"    FATAL: {csvfile} is neither a file nor directory.\n"
                         "\n")
        sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any samples at all.\n"
                         "\n")
        sys.exit(1)
    return sampleDict