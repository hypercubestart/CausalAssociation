import datetime
import multiprocessing
import os
import time

def cell_analyzer(arguments):
    start_time_timer = time.time()
    CELL_BARCODE = arguments

    # generate vcf file using freebayes
    target_vcf_file = RES_FOLDER + "/" + CELL_BARCODE + '.vcf'
    bam_file_path = SAMPLE_FOLDER + "/" + CELL_BARCODE + "/" + "Aligned.sortedByCoord.out.bam"
    command='{} --fasta-reference {} {} --min-alternate-count 2 --min-alternate-fraction 0.66 --skip-coverage 200 --min-mapping-quality 30 --min-base-quality 20 > {}'.format(FREEBAYES, FASTA_REFERENCE_FILE,bam_file_path,target_vcf_file)
    printt(command)
    os.system(command)
    printt('time elapsed: {:.2f} minutes\n'.format((time.time() - start_time_timer)/60))

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

if __name__ == "__main__":
    OMICS4TB2 = "/proj/omics4tb2"
    FREEBAYES = OMICS4TB2 + "/aliu/tools/freebayes-v1.3.1"

    RES_FOLDER = OMICS4TB2 + "/aliu/projects/causalAssociation/results/freebayes/SRR8521710"

    if not os.path.isdir(RES_FOLDER):
        os.mkdir(RES_FOLDER)

    FASTA_REFERENCE_FILE = OMICS4TB2 + "/aliu/projects/causalAssociation/data/referenceGenome/GCF_000001405.39_GRCh38.p13_genomic.fna"

    SAMPLE_FOLDER = OMICS4TB2 + "/aliu/projects/causalAssociation/results/STAR/SRR8521710"
    CELL_FOLDERS = []
    for CELL_FOLDER_NAME in os.listdir(SAMPLE_FOLDER):
        CELL_FOLDERS.append(CELL_FOLDER_NAME)

    NUM_CORES = 8
    pool = multiprocessing.Pool(NUM_CORES)
    pool.map(cell_analyzer, CELL_FOLDERS)
    pool.close()
