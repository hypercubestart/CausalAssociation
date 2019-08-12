# Quantify data using kallisto
# currently not being used in pipeline
import os, datetime

WORKDIR = "/Volumes/omics4tb2/aliu/projects/causalAssociation/"
RAW_DATA_DIR = WORKDIR + "data/rawData/"
NUM_THREADS = 8

def printt(message):
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

programStart = datetime.datetime.now()

files = open(RAW_DATA_DIR + "filenames.txt")
files = files.read().splitlines()

for file in files:
    if not os.path.isdir(WORKDIR + "results/kallisto/" + file):
        start = datetime.datetime.now()
        forward_read = RAW_DATA_DIR + file + "1_001.fastq"
        backward_read = RAW_DATA_DIR + file + "2_001.fastq"

        command = 'kallisto quant -i {} -o {} --threads={} \
            {} {}' \
            .format(WORKDIR + "data/transcriptomeIndex/Mus_musculus.GRCm38.cdna.all.idx", WORKDIR + "results/kallisto/" + file, NUM_THREADS, forward_read, backward_read)

        printt(command)
        os.system(command)
        printt("completed in {}".format(datetime.datetime.now() - start))

printt("finished mapping/quantifying {} reads in {}".format(len(files), datetime.datetime.now() - programStart))