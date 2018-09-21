from gwf import Workflow

gwf = Workflow()

ref_genome = 'hg38.fa'
samples_bam = ['sample_cancer_filter_sort.bam', 'sample_blood_filter_sort.bam']
samples_vcf = ['sample_cancer.vcf', 'sample_blood.vcf']
filtered_vcf = 'filtered.vcf'

for i in range(len(samples_bam)):
    gwf.target(
        'create_vcf',
        inputs=[f'{ref_genome}', f'{samples_bam[i]}'],
        outputs=[f'{samples_vcf[i]}']) \
        << f"""    
        samtools mpileup -u -f {ref_genome} {samples_bam[i]} | bcftools view -v -c - > {samples_vcf[i]}
        """

gwf.target(
    'filter_noncancer_variants',
    inputs=[f'{samples_vcf[0]}', f'{samples_vcf[1]}'],
    outputs=[f'{filtered_vcf}']) \
    << f"""
    vcftools {samples_vcf[1]} --diff {samples_vcf[0]} --diff-site --out {filtered_vcf}       
    """
