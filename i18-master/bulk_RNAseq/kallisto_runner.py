import sys,numpy,os,datetime

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

def kallisto_caller(label,all_files):

    printt('about to quantify {}'.format(label))

    for index in indices:
        
        index_output_dir=results_dir+index
        if os.path.exists(index_output_dir) == False:
            os.mkdir(index_output_dir)
            
        output_dir=index_output_dir+'/'+label
        executable='time kallisto quant'
        options=' -i {}{}.idx -o {} --bias -t 8 -b {} --single -l 180 -s 20 --rf-stranded '.format(transcriptome_index_dir,index,output_dir,boots)
        files=' '.join([clean_fastq_dir+element for element in all_files if label in element])
        command=executable+options+files

        print('')
        print(command)
        os.system(command)

    print('')

    return None

### MAIN

# 0. user defined variables
clean_fastq_dir='/Volumes/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/clean_fastq/'
boots=160
results_dir='/Volumes/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/kallisto.{}/'.format(boots)
transcriptome_index_dir='/Volumes/omics4tb2/alomana/projects/i18/data/annotation/'
indices=['host','pathogen','host_pathogen']

if os.path.exists(results_dir) == False:
    os.mkdir(results_dir)

# 1. recover labels
printt('recover labels...')

all_files=os.listdir(clean_fastq_dir)
labels=sorted(list(set(([element.split('_L')[0] for element in all_files]))))

# 2. call kallisto quant
for label in labels:
    kallisto_caller(label,all_files)
