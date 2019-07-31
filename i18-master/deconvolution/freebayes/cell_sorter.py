###
### this script sorts cells based on variants called previously cell-by-cell and aggregated
### across different genotypes, it considers both independently ratios of reads and qualities for reference and alternative calls
###

import sys,datetime,os,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def cell_variants_evaluator(cell_variants,real):

    supporting_read_ratios=numpy.array([0 for genotype in genotypes])
    supporting_quality_ratios=numpy.array([0 for genotype in genotypes])
    
    for cell_variant in cell_variants:
        printt('working with variant {}/{}. {} {}'.format(cell_variants.index(cell_variant)+1,len(cell_variants),cell_variant[0],cell_variant[1]))
        found_in_bulk=0
        
        for genotype in genotypes:
            read_ratios=[0 for genotype in genotypes]
            quality_ratios=[0 for genotype in genotypes]
            variant_label=cell_variant[0]
            if variant_label in population_variants[genotype]:
                found_in_bulk=found_in_bulk+1
                
                printt('originally found in {} as {} {}'.format(genotype,variant_label,population_variants[genotype][variant_label]))

                RO=population_variants[genotype][variant_label][0]
                AO=population_variants[genotype][variant_label][1]
                QR=population_variants[genotype][variant_label][2]
                QA=population_variants[genotype][variant_label][3]

                observation_ratio=AO/(RO+AO)
                quality_ratio=QA/(QR+QA)

                read_ratios[genotypes.index(genotype)]=read_ratios[genotypes.index(genotype)]+observation_ratio
                quality_ratios[genotypes.index(genotype)]=quality_ratios[genotypes.index(genotype)]+quality_ratio
                
                supporting_read_ratios=supporting_read_ratios+numpy.array(read_ratios)
                supporting_quality_ratios=supporting_quality_ratios+numpy.array(quality_ratios)

                printt('updated genotype support: reads {}; qualities {}'.format(supporting_read_ratios,supporting_quality_ratios))
                                
        if found_in_bulk == 0:
            printt('WARNING: variant {} in cell {} {} not found in population variants'.format(cell_variant,cell,genotype))
                
    read_order=numpy.flip(numpy.argsort(supporting_read_ratios))
    quality_order=numpy.flip(numpy.argsort(supporting_quality_ratios))

    prediction_read=genotypes[read_order[0]]
    prediction_quality=genotypes[quality_order[0]]

    read_support=supporting_read_ratios[read_order[0]]/supporting_read_ratios[read_order[1]]

    if real == prediction_read:
        assessment='SUCCESS'
    else:
        assessment='FAILED'

    printt('label: {}; predictions: {} {}; assessment: {}'.format(real,prediction_read,prediction_quality,assessment))
    print('')

    return read_support

def printt(message):
    
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

def variant_reader(input_file):

    variants=[]
    with open(input_file,'r') as f:
        for line in f:
            if line[0] != '#':
                v=line.split('\t')
                alternatives=v[4].split(',') # alternatives may be multiple
                for i in range(len(alternatives)):
                    label='{}.{}.{}'.format(v[0],v[1],v[4].split(',')[i]) # chr.position.alternative
                    info=v[-1].split(':')
                    RO=int(info[3]); QR=int(info[4]); AO=int(info[5].split(',')[i]); QA=int(info[6].split(',')[i])
                    quantification=[RO,AO,QR,QA] # RO,AO,QR,QA
                    variants.append((label,quantification))

    return variants

###
### MAIN
###

# 0. user defined variables
population_variants_files={}
population_variants_files['wt']='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/freebayes/wt/subset.046.cells/aggregated_variants.vcf'
population_variants_files['mavs']='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/freebayes/mavs/subset.043.cells/aggregated_variants.vcf'
population_variants_files['nlrp3']='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/freebayes/nlrp3/subset.008.cells/aggregated_variants.vcf'
single_cell_dir='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/freebayes/'
genotypes=sorted(list(population_variants_files.keys()))

# 1. read population variants
printt('read population variants')
population_variants={}

for genotype in genotypes:
    population_variants[genotype]={}
    variants=variant_reader(population_variants_files[genotype])
    printt('detected {} population variants for genotype {}'.format(len(variants),genotype))
    for variant in variants:
        population_variants[genotype][variant[0]]=variant[1]
    
# 2. sort cells
printt('sort cells')
x=[]; y=[]; z=[]; t=[]

for genotype in genotypes:
    printt('working with genotype {}'.format(genotype))
    detected_folders=os.listdir(single_cell_dir+genotype)
    detected_cells=sorted([element for element in detected_folders if 'cells' in element])
    for cell in detected_cells:
        printt('working with cell {}'.format(cell))

        mapped_reads_file=single_cell_dir+genotype+'/'+cell+'/one_cell_mapped_reads.txt'
        with open(mapped_reads_file,'r') as f:
            v=f.readline()
            mapped_reads=int(v)
        printt('cell with {} mapped reads'.format(mapped_reads))
        
        variant_file=single_cell_dir+genotype+'/'+cell+'/one_cell_variants.vcf'
        cell_variants=variant_reader(variant_file)
        printt('detected {} variants'.format(len(cell_variants)))

        real=genotype
        read_support=cell_variants_evaluator(cell_variants,real)

        # adding data to build figure
        x.append(mapped_reads); y.append(read_support); z.append(len(cell_variants)); t.append(genotype)

# 3. make figure
unique_genotypes=set(t)

for genotype in unique_genotypes:
    a=[]; b=[]; c=[]

    if genotype == 'mavs':
        marker = 'D'
    if genotype == 'nlrp3':
        marker = 's'
    if genotype == 'wt':
        marker = 'o'

    for i in range(len(x)):
        if t[i] == genotype:
            a.append(x[i])
            b.append(y[i])
            c.append(z[i])

    matplotlib.pyplot.scatter(numpy.log10(a),b,c=numpy.log2(c),marker=marker,cmap='viridis',vmin=min(numpy.log2(z)),vmax=max(numpy.log2(z)),label=genotype)

matplotlib.pyplot.axhline(y=1,ls=':',color='black')

cbar=matplotlib.pyplot.colorbar()
cbar.set_ticks([2,3,4,5,6,7,8])
cbar.set_ticklabels(2**(numpy.array([2,3,4,5,6,7,8])))
cbar.ax.set_title('Detected variants',fontsize=14)

matplotlib.pyplot.xlim([3-0.1,numpy.max(numpy.log10(x))+0.1])

matplotlib.pyplot.yticks([1,2,4,6,8,10],['bkg',2,4,6,8,10])

matplotlib.pyplot.xlabel('log$_{10}$ Mapped reads')
matplotlib.pyplot.ylabel('Support ratio')

leg=matplotlib.pyplot.legend(fontsize=14)
for marker in leg.legendHandles:
    marker.set_color('black')

matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/freebayes/figure.variant_support.pdf')
matplotlib.pyplot.clf()

