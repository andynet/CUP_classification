from gwf import Workflow
from templates import *

gwf = Workflow(defaults={
    "cores": 16,
    "memory": "32g",
    "walltime": "24:00:00",
    "nodes": None,
    "queue": "normal",
    "account": "CUP_classification",
    "constraint": None,
    "mail_type": None,
    "mail_user": None,
    "qos": None,
})

working_directory = '/home/andyb/CUP_classification/faststorage/Andrej'

normal_fastq = f'{working_directory}/normal_fastq_files.txt'
tumor_fastq = f'{working_directory}/tumor_fastq_files.txt'

reference_genome = f'{working_directory}/inputs/hg38.fa'
bowtie2_index = f'{working_directory}/bowtie2_index/hg38'
bed_file = f'{working_directory}/inputs/covered_regions.bed'

out_dir = f'{working_directory}/outputs'
start = 0
end = 1
threads = 8

with open(normal_fastq) as f:
    normal_fastq_lines = f.readlines()

with open(tumor_fastq) as f:
    tumor_fastq_lines = f.readlines()

locations = [normal_fastq_lines, tumor_fastq_lines]

for i in range(start, end):

    sample_ids = []
    vcfs = []
    tsvs = []

    for j, stype in enumerate(['normal', 'tumor']):

        sample_id, fastq1, fastq2 = locations[j][i].split()
        sam = f'{out_dir}/{sample_id}.{stype}.sam'
        sbam = f'{out_dir}/{sample_id}.{stype}.sorted.bam'
        ibam = f'{out_dir}/{sample_id}.{stype}.sorted.bam.bai'
        vcf = f'{out_dir}/{sample_id}.{stype}.vcf'
        tsv = f'{out_dir}/{sample_id}.{stype}.vcf.tsv'

        sample_ids.append(sample_id)
        vcfs.append(vcf)
        tsvs.append(tsv)

        # <editor-fold desc="fastq to bam">
        gwf.target(
            f'{sample_id}.{stype}.sorted.bam',
            inputs=[f'{fastq1}', f'{fastq2}'],
            outputs=[f'{sbam}'],
            options={'cores': threads},
        ) << f'''
        bowtie2 --threads {threads}                                                                                     \
                -x {bowtie2_index}                                                                                      \
                -1 {fastq1}                                                                                             \
                -2 {fastq2}                                                                                             \
        | samtools sort -@ {threads} - > {sbam}
           
        samtools index {sbam}
        '''
        # </editor-fold>

        tmp = filter_reads(sbam, reference_genome)
        gwf.target_from_template(f'filter_{sbam}', tmp)
        final_bam = tmp[1][0]

        # <editor-fold desc="bam to vcf">
        gwf.target(
            f'{sample_id}.{stype}.vcf',
            inputs=[f'{reference_genome}', f'{final_bam}'],
            outputs=[f'{vcf}'],
        ) << f"""
        samtools index {final_bam}
        
        samtools mpileup                                                                                                \
                -u -tAD                                                                                                 \
                -f {reference_genome}                                                                                   \
                -l {bed_file}                                                                                           \
                {final_bam} | bcftools view -v snps -m2 > {vcf}
        """
        # </editor-fold>

        # <editor-fold desc="vcf to tsv">
        gwf.target(
            f'{sample_id}.{stype}.vcf.tsv',
            inputs=[f'{vcf}'],
            outputs=[f'{tsv}'],
        ) << f"""
        less {vcf} | grep -v '^#' | cut -d$'\t' -f1,2,4,5 | sort > {tsv}
        """
        # </editor-fold>

    # <editor-fold desc="filter variants">
    file_name = f'{sample_ids[1]}_tumor_{sample_ids[1]}_normal.tsv'
    final_tsv = f'{out_dir}/{file_name}'

    gwf.target(
        f'{file_name}',
        inputs=[f'{tsvs[0]}', f'{tsvs[1]}'],
        outputs=[f'{final_tsv}'],
    ) << f"""
    comm -13 {tsvs[0]} {tsvs[1]} > {final_tsv}
    """
    # </editor-fold>

    # <editor-fold desc="create final vcf">
    final_vcf = f'{out_dir}/{file_name}.vcf'
    gwf.target(
        f'{file_name}.vcf',
        inputs=[f'{vcfs[1]}', f'{final_tsv}'],
        outputs=[f'{final_vcf}'],
    ) << f"""
    <{vcfs[1]} grep "^#" > {final_vcf}
    
    join -j1 -t$'\t' -o1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11                                                        \
        <(<{vcfs[1]} grep -v '^#' | awk '{{print $1"-"$2"-"$4"-"$5"\t"$0}}' | sort -k1,1)                               \
        <(<{final_tsv} awk '{{print $1"-"$2"-"$3"-"$4"\t"$0}}' | sort -k1,1)                                            \
        >> {final_vcf}
    """
    # </editor-fold>
