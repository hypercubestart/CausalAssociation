import os,datetime,sys

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

def trimmomatic_caller(path,label,lane):

    detected_file=os.listdir(path)[0]

    executable='time java -jar {}trimmomatic-0.39.jar SE -threads {} -phred33 '.format(trimmomatic_path,number_threads)

    input_file=path+'/'+detected_file+' '
    output_file=fastq_clean_dir+'{}.clean.fastq'.format(lane.split('-')[0])
    
    options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)

    command=executable+input_file+output_file+options

    printt('about to clean {}'.format(label))
    print('')
    print(command)
    print('') 
    os.system(command)
    print('')

    return None

# 0. user defined variables
fastq_dir='/Volumes/omics4tb2/alomana/projects/i18/data/bulk_rna/FASTQ_Generation_2019-06-26_12_35_52Z-188844669/'
fastq_clean_dir='/Volumes/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/'
trimmomatic_path='/Users/alomana/software/Trimmomatic-0.39/'
adapter_file=trimmomatic_path+'adapters/TruSeq3-SE.fa'
number_threads=8

# 1. recover labels
printt('recover reads')

# 1.1. locating input files
all_files=os.listdir(fastq_dir)
elements=list(set([element.split('_L')[0] for element in all_files if element[0] != '.']))
labels=sorted(elements)

sample_correspondance={}
for label in labels:
    sample_correspondance[label]=[]
    for element in all_files:
        if label in element:
            sample_correspondance[label].append(element)
    sample_correspondance[label].sort()

# 2. iterate Trimmomatic command
for label in labels:
    for lane in sample_correspondance[label]:
        path=fastq_dir+lane
        trimmomatic_caller(path,label,lane)
