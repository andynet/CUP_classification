from gwf import Workflow
import yaml

# <editor-fold desc="initialization">
with open('config.yaml') as f:
    config = yaml.safe_load(f)

data_dir = config['data_dir']
normal_list = config['normal_list']
tumor_list = config['tumor_list']
reference = config['reference']
bowtie2_index = config['bowtie2_index']
bed_file = config['bed_file']

gwf = Workflow(defaults={
    "cores": 16,
    "memory": "32g",
    "walltime": "24:00:00",
    "queue": "normal",
    "account": "CUP_classification",
})
# </editor-fold>

gwf.target(
    name='CreateSequenceDictionary',
    inputs=[reference],
    outputs=[reference.replace('.fa', '.dict')],
) << f"""
    picard CreateSequenceDictionary                                             \
        REFERENCE={reference}                                                   \
        OUTPUT={reference.replace('.fa', '.dict')}
    """


with open(normal_list) as f:
    normal_fastq_lines = f.readlines()

with open(tumor_list) as f:
    tumor_fastq_lines = f.readlines()

locations = [normal_fastq_lines, tumor_fastq_lines]

start = 0
end = len(normal_fastq_lines)
# end = 1

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
        '''
        # </editor-fold>

        tmp = add_groups(sbam)
        gwf.target_from_template(f'add_groups_{sample_id}_{stype}', tmp)
        bam = tmp[1][0]

        gwf.target_from_template(f'index_{sample_id}_{stype}', index(bam))
        tmp = filter_reads(bam, reference_genome)
        gwf.target_from_template(f'filter_{sample_id}_{stype}', tmp)
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

gwf.target_from_template(f'notify_{i}', notify('andrejbalaz001@gmail.com', final_vcf))

