import csv
import os

def findFileName(start, target_dir):
    for filename in os.listdir(target_dir):
        if filename.startswith(start):
            return filename

def parseFiles():
    f = open("/home/aliu/Projects/CausalAssociation/data/rawData/SRR_USE_LIST.txt", "r")
    f = f.read().splitlines()
    lines = []
    for line in f:
      lines.append('/proj/omics4tb2/aliu/projects/causalAssociation/data/rawData/' + line + '.fastq')

    print(",".join(lines))

def combineCSVFiles():
    OMICS4TB2 = "/home/aliu/omics4tb2" #"/proj/omics4tb2"
    HOME = "/home"
    SCC_USE_LIST_FILE = HOME + "/aliu/Projects/causalAssociation/data/rawData/SRR_USE_LIST.txt"
    GCM_USE_LIST_FILE = HOME + "/aliu/Projects/causalAssociation/data/rawData/GCM_USE_LIST"

    SCC_IDENTIFIERS = []
    GCM_IDENTIFIERS = []
    with open(SCC_USE_LIST_FILE, "r") as file:
        lines = file.read().splitlines()
        SCC_IDENTIFIERS = lines
    with open(GCM_USE_LIST_FILE) as gsmfile:
        gsmfile = csv.reader(gsmfile, delimiter='\t')
        for line in gsmfile:
            GCM_IDENTIFIERS.append(line)

    combined = []

    for x in range(len(GCM_IDENTIFIERS)):
        row = GCM_IDENTIFIERS[x]
        row.append(SCC_IDENTIFIERS[x])
        combined.append(row)

    with open(HOME + '/aliu/Projects/causalAssociation/data/rawData/COMBINED.csv', 'w') as writeOutFile:
        writer = csv.writer(writeOutFile, dialect='excel-tab')
        writer.writerows(combined)

if __name__ == "__main__":
    combineCSVFiles()
