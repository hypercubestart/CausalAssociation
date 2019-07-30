def parseFiles():
    f = open("/home/aliu/Projects/CausalAssociation/data/rawData/SRR_USE_LIST.txt", "r")
    f = f.read().splitlines()
    lines = []
    for line in f:
      lines.append('/proj/omics4tb2/aliu/projects/causalAssociation/data/rawData/' + line + '.fastq')

    print(",".join(lines))

if __name__ == "__main__":
    parseFiles()
