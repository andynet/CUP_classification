from gwf import Workflow

gwf = Workflow()

working_directory = '/home/andyb/CUP_classification/faststorage/Andrej'

normal_fastq = f'{working_directory}/normal_fastq_files.txt'
tumor_fastq = f'{working_directory}/tumor_fastq_files.txt'

reference_genome = f'{working_directory}/inputs/hg38.fa'
bowtie2_index = f'{working_directory}/bowtie2_index/hg38'
bed_file = f'{working_directory}/inputs/covered_regions.bed'

out_dir = f'{working_directory}/outputs'
n = 1

with open(normal_fastq) as f:
    normal_fastq_lines = f.readlines()

with open(tumor_fastq) as f:
    tumor_fastq_lines = f.readlines()

for i in range(n):

    # <editor-fold desc="create bam from normal samples">
    nsample_id, nfastq1, nfastq2 = normal_fastq_lines[i].split()
    nsam = f'{out_dir}/{nsample_id}.normal.sam'
    nbam = f'{out_dir}/{nsample_id}.normal.bam'

    gwf.target(
        'create_sam_normal',
        inputs=[f'{nfastq1}', f'{nfastq2}'],
        outputs=[f'{nbam}'], 
        walltime="12:00:00", memory="32g"
    ) << f'''
    bowtie2 -x '{bowtie2_index}' \
            -1 '{nfastq1}'  \
            -2 '{nfastq2}'  \
            -S '{nsam}'
            
    samtools view -Sb {nsam} > {nbam}
    '''
    # </editor-fold>

    # # <editor-fold desc="create vcf from normal samples">
    # nvcf = f'{out_dir}/{nsample_id}.normal.vcf'
    #
    # gwf.target(
    #     'create_vcf_normal',
    #     inputs=[f'{reference_genome}', f'{nbam}'],
    #     outputs=[f'{nvcf}']) \
    # << f"""
    #         samtools mpileup    \
    #             -u -tAD     \
    #             -f {reference_genome}   \
    #             -l {bed_file}    \
    #             {nbam} | bcftools view -v snps -m2
    #         """
    # # </editor-fold>

    # <editor-fold desc="create bam from tumor samples">
    tsample_id, tfastq1, tfastq2 = tumor_fastq_lines[i].split()
    tsam = f'{out_dir}/{tsample_id}.tumor.sam'
    tbam = f'{out_dir}/{tsample_id}.tumor.bam'

    gwf.target(
        'create_sam_tumor',
        inputs=[f'{tfastq1}', f'{tfastq2}'],
        outputs=[f'{tbam}'], 
        walltime="12:00:00", memory="32g"
    ) << f'''
        bowtie2 -x '{bowtie2_index}' \
                -1 '{tfastq1}'  \
                -2 '{tfastq2}'  \
                -S '{tsam}'
                
        samtools view -Sb {tsam} > {tbam}
        '''
    # </editor-fold>

#
# gwf.target(
#     'filter_noncancer_variants',
#     inputs=[f'{samples_vcf[0]}', f'{samples_vcf[1]}'],
#     outputs=[f'{filtered_vcf}']) \
#     << f"""
#     vcftools {samples_vcf[1]} --diff {samples_vcf[0]} --diff-site --out {filtered_vcf}
#     """
