### 
### this pipeline calls variants of single cells and all aggregated cells, per genotype specified.
### 

import pickle,datetime,os,sys,numpy
import multiprocessing,multiprocessing.pool

def cell_analyzer(argument):

    i=argument[0]
    cell_indexes=argument[1]
    output_dir=argument[2]

    printt('working with {}/{} cells'.format(i+1,len(cells)))
    working_cells=cells[:i+1]
    target_dir=output_dir+'subset.{:03d}.cells'.format(i+1)
    if os.path.exists(target_dir) == False:
        os.mkdir(target_dir)

    # f.3. define list of barcodes
    one_cell_barcode_file=target_dir+'/one_cell_barcode.txt'
    with open(one_cell_barcode_file,'w') as f:
        f.write('{}\n'.format(working_cells[-1]))
    f.close()

    if i == cell_indexes[-1]:
        aggregated_barcode_file=target_dir+'/aggregated_barcodes.txt'
        with open(aggregated_barcode_file,'w') as f:
            for barcode in working_cells:
                f.write('{}\n'.format(barcode))
        f.close()

    # f.4. call subset, check number of mapped reads and run freebayes for single cells and aggregated bam files
    one_cell_bam_out_file=target_dir+'/one_cell.bam'
    command='subset-bam --bam {} --cell-barcodes {} --out-bam {}'.format(path_to_bamfile,one_cell_barcode_file,one_cell_bam_out_file)
    printt(command)
    os.system(command)
    if i == cell_indexes[-1]:
        aggregated_bam_out_file=target_dir+'/aggregated.bam'
        command='subset-bam --bam {} --cell-barcodes {} --out-bam {}'.format(path_to_bamfile,aggregated_barcode_file,aggregated_bam_out_file)
        printt(command)
        os.system(command)

    mapped_reads_file=target_dir+'/one_cell_mapped_reads.txt'
    command='samtools view -c -F 260 {} > {}'.format(one_cell_bam_out_file,mapped_reads_file)
    printt(command)
    os.system(command)
    if i == cell_indexes[-1]:
        mapped_reads_file=target_dir+'/aggregated_mapped_reads.txt'
        command='samtools view -c -F 260 {} > {}'.format(aggregated_bam_out_file,mapped_reads_file)
        printt(command)
        os.system(command)

    target_vcf_file=target_dir+'/one_cell_variants.vcf'
    command='freebayes --fasta-reference {} {} --min-alternate-count 2 --min-alternate-fraction 0.66 --skip-coverage 200 --min-mapping-quality 30 --min-base-quality 20 > {}'.format(fasta_reference_file,one_cell_bam_out_file,target_vcf_file)
    printt(command)
    os.system(command)
    if i == cell_indexes[-1]:
        target_vcf_file=target_dir+'/aggregated_variants.vcf'
        command='freebayes --fasta-reference {} {} --min-alternate-count 2 --min-mapping-quality 30 --min-base-quality 20 > {}'.format(fasta_reference_file,aggregated_bam_out_file,target_vcf_file)
        printt(command)
        os.system(command)

    return None

def job_distributor(cells,label,path_to_bamfile):

    printt('working with {} cells'.format(len(cells)))

    output_dir=subset_dir+'{}/'.format(label)
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    cell_indexes=numpy.arange(len(cells))

    arguments=[[cell_index,cell_indexes,output_dir] for cell_index in cell_indexes]

    hydra=multiprocessing.pool.Pool(threads)
    empty=hydra.map(cell_analyzer,arguments)
    hydra.close()

    return None        
    
def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

# 0. user defined variables
selected_cellIDs_file='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/scanpy/species.cellIDs.run2.001.pickle'
fasta_reference_file='/Volumes/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0/fasta/genome.fa'
subset_dir='/home/aliu/omics4tb2/aliu/projects/causalAssociation/results/freebayes/'
threads=8

# 1. recover mouse cellIDs
printt('recovering cellIDs')
f=open(selected_cellIDs_file,'rb')
[mouse_cellIDs,human_cellIDs,chimeric_cellIDs]=pickle.load(f)
f.close()
printt('recovered {} mouse cells'.format(len(mouse_cellIDs)))

printt('checking non-overlapping IDs')
trimmed_cellIDs=[element.split('-')[0] for element in mouse_cellIDs]
cellIDs=list(set(trimmed_cellIDs))
if len(mouse_cellIDs) != len(cellIDs):
    printt('error because of overlapping cellIDs')

# 2. build recursive bam files with sequential number of cells
printt('building BAM files for WT')
cells=[element for element in mouse_cellIDs if '-1' in element]
label='wt'
path_to_bamfile='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/cell_ranger/mouse/d2_dmso_mouse/outs/possorted_genome_bam.bam'
job_distributor(cells,label,path_to_bamfile)

printt('building BAM files for MAVS')
cells=[element.replace('-2','-1') for element in mouse_cellIDs if '-2' in element]
label='mavs'
path_to_bamfile='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/cell_ranger/mouse/d4_dmso_mouse/outs/possorted_genome_bam.bam'
job_distributor(cells,label,path_to_bamfile)

printt('building BAM files for NLRP3')
cells=[element.replace('-3','-1') for element in mouse_cellIDs if '-3' in element]
label='nlrp3'
path_to_bamfile='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/cell_ranger/mouse/d4_pstat_r2_mouse/outs/possorted_genome_bam.bam'
job_distributor(cells,label,path_to_bamfile)
