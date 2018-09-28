from gwf import Workflow

gwf = Workflow()

working_directory = '/home/andyb/CUP_classification/faststorage/Andrej'

reference_genome = f'{working_directory}/hg38.fa'
bowtie2_index = f'{working_directory}/bowtie2_index/hg38'
normal_fastq = f'{working_directory}/normal_fastq_files.txt'
tumor_fastq = f'{working_directory}/tumor_fastq_files.txt'
bed_file = f'{working_directory}/covered_regions.bed'

out_dir = f'{working_directory}/outputs'
n = 1

with open(normal_fastq) as f:
    normal_fastq_lines = f.readlines()

with open(tumor_fastq) as f:
    tumor_fastq_lines = f.readlines()

for i in range(n):

    sample_id, fastq1, fastq2 = normal_fastq_lines[i].split()
    sam = f'{out_dir}/{sample_id}.normal.sam'

    gwf.target(
        'create_bam',
        inputs=[f'{fastq1}', f'{fastq2}', f'{bowtie2_index}'],
        outputs=[f'{sam}']
    ) << f'''
    bowtie2 -x f'{bowtie2_index}' \
            -1 f'{fastq1}'  \
            -2 f'{fastq2}'  \
            -S f'{sam}'
    '''

#     gwf.target(
#         'create_vcf',
#         inputs=[f'{ref_genome}', f'{samples_bam[i]}'],
#         outputs=[f'{samples_vcf[i]}']) \
#         << f"""
#         samtools mpileup -u -f {ref_genome} {samples_bam[i]} | bcftools view -v -c - > {samples_vcf[i]}
#         samtools mpileup -u -tAD -f hg38.fa -l {bedfile}  {bamfile} | bcftools view -v snps -m2
#         """
#
# gwf.target(
#     'filter_noncancer_variants',
#     inputs=[f'{samples_vcf[0]}', f'{samples_vcf[1]}'],
#     outputs=[f'{filtered_vcf}']) \
#     << f"""
#     vcftools {samples_vcf[1]} --diff {samples_vcf[0]} --diff-site --out {filtered_vcf}
#     """
