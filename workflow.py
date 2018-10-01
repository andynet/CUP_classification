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

locations = [normal_fastq_lines, tumor_fastq_lines]

for i in range(n):
    for j, stype in enumerate(['normal', 'tumor']):

        # <editor-fold desc="fastq to bam">
        sample_id, fastq1, fastq2 = locations[j][i].split()
        sam = f'{out_dir}/{sample_id}.{stype}.sam'
        sbam = f'{out_dir}/{sample_id}.{stype}.sorted.bam'
        ibam = f'{out_dir}/{sample_id}.{stype}.sorted.bam.bai'

        gwf.target(
            f'{sample_id}.{stype}.sorted.bam',
            inputs=[f'{fastq1}', f'{fastq2}'],
            outputs=[f'{sbam}', f'{ibam}'],
            walltime="12:00:00", memory="32g"
        ) << f'''
        bowtie2 -x '{bowtie2_index}' \
                -1 '{fastq1}'  \
                -2 '{fastq2}'  \
                -S '{sam}'
            
        samtools sort {sam} > {sbam}
        samtools index {sbam}
        '''
        # </editor-fold>

        # <editor-fold desc="bam to vcf">
        vcf = f'{out_dir}/{sample_id}.{stype}.vcf'

        gwf.target(
            f'{sample_id}.{stype}.vcf',
            inputs=[f'{reference_genome}', f'{sbam}'],
            outputs=[f'{vcf}'],
            walltime="12:00:00", memory="32g"
        ) << f"""
                samtools mpileup    \
                    -u -tAD     \
                    -f {reference_genome}   \
                    -l {bed_file}    \
                    {sbam} | bcftools view -v snps -m2 > {vcf}
                """
        # </editor-fold>


#
# gwf.target(
#     'filter_noncancer_variants',
#     inputs=[f'{samples_vcf[0]}', f'{samples_vcf[1]}'],
#     outputs=[f'{filtered_vcf}']) \
#     << f"""
#     vcftools {samples_vcf[1]} --diff {samples_vcf[0]} --diff-site --out {filtered_vcf}
#     """
