from gwf import Workflow

gwf = Workflow()

working_directory = '/home/andyb/CUP_classification/faststorage/Andrej'

normal_fastq = f'{working_directory}/normal_fastq_files.txt'
tumor_fastq = f'{working_directory}/tumor_fastq_files.txt'

reference_genome = f'{working_directory}/inputs/hg38.fa'
bowtie2_index = f'{working_directory}/bowtie2_index/hg38'
bed_file = f'{working_directory}/inputs/covered_regions.bed'

out_dir = f'{working_directory}/outputs'
n = 2

with open(normal_fastq) as f:
    normal_fastq_lines = f.readlines()

with open(tumor_fastq) as f:
    tumor_fastq_lines = f.readlines()

locations = [normal_fastq_lines, tumor_fastq_lines]

for i in range(n):

    sample_ids = []
    vcfs = []

    for j, stype in enumerate(['normal', 'tumor']):

        sample_id, fastq1, fastq2 = locations[j][i].split()
        sam = f'{out_dir}/{sample_id}.{stype}.sam'
        sbam = f'{out_dir}/{sample_id}.{stype}.sorted.bam'
        ibam = f'{out_dir}/{sample_id}.{stype}.sorted.bam.bai'
        vcf = f'{out_dir}/{sample_id}.{stype}.vcf'

        sample_ids.append(sample_id)
        vcfs.append(vcf)

        # <editor-fold desc="fastq to bam">
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

    # <editor-fold desc="filter variants">
    sample_normal_tsv = f'{sample_ids[0]}.vcf.tsv'
    sample_tumor_tsv = f'{sample_ids[1]}.vcf.tsv'
    final_tsv = f'{sample_ids[1]}_tumor_{sample_ids[1]}_normal.tsv'

    gwf.target(
        f'{final_tsv}',
        inputs=[f'{vcfs[0]}', f'{vcfs[1]}'],
        outputs=[f'{final_tsv}'],
        walltime="12:00:00", memory="32g"
    ) << f"""
    less {vcfs[0]} | grep -v '^#' | cut -d$'\t' -f1,2,4,5 | sort > {sample_normal_tsv}
    less {vcfs[1]} | grep -v '^#' | cut -d$'\t' -f1,2,4,5 | sort > {sample_tumor_tsv}
    
    comm -13 {sample_normal_tsv} {sample_tumor_tsv} > {final_tsv}
    """
    # </editor-fold>
