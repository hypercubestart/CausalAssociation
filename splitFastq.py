# split fastq file by single cells

import os, multiprocessing
import datetime
import csv
import time

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

# def splitFastQ(fastq_file):
#     # too slow!!!
#     if not os.path.isdir("/proj/omics4tb2/aliu/projects/causalAssociation/data/scData/" + fastq_file):
#         os.mkdir("/proj/omics4tb2/aliu/projects/causalAssociation/data/scData/" + fastq_file)
#
#     umi_set = set()
#     library_barcode_counts = {}
#
#     fastq_file_path = "/proj/omics4tb2/aliu/projects/causalAssociation/data/rawData/" + fastq_file + ".fastq"
#     i = 0
#     with open(fastq_file_path, "r") as file:
#         write_file = ""
#         for line in file:
#             if i % 4 == 0:
#                 split = line[:-1].split(" ")
#                 split = split[-1].split("_")
#                 write_file = "/proj/omics4tb2/aliu/projects/causalAssociation/data/scData/" + fastq_file + "/" + split[-2] + '.fastq'
#
#                 if split[-3] not in library_barcode_counts:
#                     library_barcode_counts[split[-3]] = 0
#                 library_barcode_counts[split[-3]] += 1
#                 umi_set.add(split[-1])
#                 if i % 20000 == 0:
#                     printt('reached line {} in fastq file: {}'.format(i, fastq_file))
#
#             appendLineToFile(line, write_file)
#             i += 1
#     metadata_path = "/proj/omics4tb2/aliu/projects/causalAssociation/data/scData/" + fastq_file + ".metadata.txt"
#     for key in library_barcode_counts:
#         appendLineToFile('Library Barcode: {}, Counts: {}\n'.format(key, library_barcode_counts[key]), metadata_path)
#     appendLineToFile('\n', metadata_path)
#     appendLineToFile('unique UMIs: {}\n'.format(len(umi_set)), metadata_path)
#     appendLineToFile('number of lines: {}'.format(i), metadata_path)
#     return None

def appendLineToFile(line, file):
    if not os.path.isfile(file):
        open(file, "w")
    with open(file, "a") as file:
        file.write(line)

def filterFastQ(arguments):
    fastq_file = arguments[0]
    cell_barcodes = arguments[1]

    OUTPUT_DATA_DIR = HOME + "/aliu/Projects/causalAssociation/data/scData/" + fastq_file # OMICS4TB2/projects FOR EAGER, HOME/Projects for LOCAL
    if not os.path.isdir(OUTPUT_DATA_DIR):
        os.mkdir(OUTPUT_DATA_DIR)

    fastq_file_path = OMICS4TB2 + "/aliu/projects/causalAssociation/data/rawData/" + fastq_file + ".fastq"

    for cell_barcode in cell_barcodes:
        start_time = time.time()
        printt('starting to filter cell {} from fastq file {}'.format(cell_barcode, fastq_file))
        cell_fastq_path = OUTPUT_DATA_DIR + "/" + cell_barcode + '.fastq'
        with open(cell_fastq_path, "a+") as sc_fastq_file:
            with open(fastq_file_path, "r") as file:
                i = 4
                append = 0
                for line in file:
                    if i == 4:
                        split = line[:-1].split(" ")
                        split = split[-1].split("_")
                        if split[-2] == cell_barcode:
                            append+= 4
                        i = 0
                    if append > 0:
                        sc_fastq_file.write(line)
                        append-=1
                    i += 1
        printt('done with cell {} from fastq file {} in {:.2f} minutes'.format(cell_barcode, fastq_file, (time.time() - start_time) / 60.))

def getCellBarcodes(dem_file_path):
    with open(dem_file_path) as csvfile:
        in_txt = csv.reader(csvfile, delimiter = '\t')
        for line in in_txt:
            cell_barcodes = line[1: ]
            cell_barcodes = [x.split('_')[1] for x in cell_barcodes]
            if len(cell_barcodes) != len(set(cell_barcodes)):
                printt('Cell Barcode Overlap in dem file: {}'.format(dem_file_path))
            return set(cell_barcodes)

def findFileName(start, target_dir):
    for filename in os.listdir(target_dir):
        if filename.startswith(start):
            return filename

if __name__ == "__main__":
    OMICS4TB2 = "/home/aliu/omics4tb2" #"/proj/omics4tb2"
    HOME = "/home"
    SCC_USE_LIST_FILE = OMICS4TB2 + "/aliu/projects/causalAssociation/data/rawData/SRR_USE_LIST.txt"
    GCM_USE_LIST_FILE = OMICS4TB2 + "/aliu/projects/causalAssociation/data/rawData/GCM_USE_LIST"
    
    SCC_IDENTIFIERS = []
    GCM_IDENTIFIERS = []
    with open(SCC_USE_LIST_FILE, "r") as file: 
        lines = file.read().splitlines()
        SCC_IDENTIFIERS = lines
    with open(GCM_USE_LIST_FILE) as gsmfile:
        gsmfile = csv.reader(gsmfile, delimiter='\t')
        for line in gsmfile:
            gsm_file_name = findFileName(line[0], OMICS4TB2 + '/aliu/projects/causalAssociation/results/expected/')
            GCM_IDENTIFIERS.append(OMICS4TB2 + '/aliu/projects/causalAssociation/results/expected/' + gsm_file_name)

    file_names = list(zip(SCC_IDENTIFIERS, GCM_IDENTIFIERS))

    file_names = file_names[20:] #[:20] for eager, [20:] for local

    args = [[x[0], getCellBarcodes(x[1])] for x in file_names]

    # for filterFastQArg in args:
    #     filterFastQ(filterFastQArg)

    num_cores = 4 # multiprocessing.cpu_count() # 2 for eager
    pool = multiprocessing.Pool(num_cores)
    pool.map(filterFastQ, args)
    pool.close()

    # with open(SCC_USE_LIST_FILE, "r") as file:
    #     lines = file.read().splitlines()
    #
    #     num_cores = multiprocessing.cpu_count()
    #     pool = multiprocessing.Pool(num_cores)
    #     pool.map(splitFastQ, lines)
    #     pool.close()
