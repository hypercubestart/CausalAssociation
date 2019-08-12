# map each single cell fastq using STAR

import os
import multiprocessing
import datetime
import time

def printt(message):
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

def STAR_map_dir(argument):
    """Map fastq file using STAR
    :param argument: [fastq_file_path, STAR output folder]
    """
    start_time = time.time()
    fastq_file_path = argument[0]
    fastq_output_path = argument[1]

    if not os.path.isdir(fastq_output_path):
        os.mkdir(fastq_output_path)

    os.system('STAR --runThreadN 2 --genomeDir /proj/omics4tb2/aliu/projects/causalAssociation/data/genomeIndex '
              '--readFilesIn {} --outFileNamePrefix {} --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate'
              .format(fastq_file_path, fastq_output_path))

    printt('finished {} in {:.2f} minutes'.format(fastq_file_path, ((time.time() - start_time)/60.)))

if __name__ == "__main__":
    OMICS4TB2 = "/proj/omics4tb2"
    # SCC_USE_LIST_FILE = OMICS4TB2 + "/aliu/projects/causalAssociation/data/rawData/SRR_USE_LIST.txt"
    FASTQ_DIR = OMICS4TB2 + "/aliu/projects/causalAssociation/data/scData/SRR8521710/"
    ARGS = []
    for filename in os.listdir(FASTQ_DIR):
        file_path = FASTQ_DIR + filename
        output_path = OMICS4TB2 + "/aliu/projects/causalAssociation/results/STAR/SRR8521710/" + filename.split('.')[0] + '/'
        ARGS.append([file_path, output_path])

    num_cores = 2
    pool = multiprocessing.Pool(num_cores)
    pool.map(STAR_map_dir, ARGS)
    pool.close()

    # with open(SCC_USE_LIST_FILE, "r") as file:
    #     lines = file.read().splitlines()
    #
    #     for srr_identifier in lines:
    #         AML_SAMPLE_FILE_LOCATION = OMICS4TB2 + "/aliu/projects/causalAssociation/data/scData/" + srr_identifier + "/"
    #
    #         for filename in os.listdir(AML_SAMPLE_FILE_LOCATION):
    #             if filename.endswith('fastq'):

