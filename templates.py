def filter_reads(bam, ref, cores=16):

    filtered_bam = bam.replace('.bam', '') + '.filtered.bam'

    inputs = [f'{bam}', f'{ref}']
    outputs = [f'{filtered_bam}']
    options = {'cores': cores}
    spec = f'''
        source /com/extra/java/8/load.sh
        source /com/extra/GATK/3.8/load.sh
        gatk                                            \
            -T PrintReads                               \
            -R {ref}                                    \
            -I {bam}                                    \
            -o {filtered_bam}                           \
            -nct {cores}                                \
            --read_filter BadCigar                      \
            --read_filter DuplicateRead                 \
            --read_filter FailsVendorQualityCheck       \
            --read_filter HCMappingQuality              \
            --read_filter MappingQualityUnavailable     \
            --read_filter NotPrimaryAlignment           \
            --read_filter UnmappedRead                  \
            --filter_bases_not_stored                   \
            --filter_mismatching_base_and_quals
    '''

    return inputs, outputs, options, spec


def create_fasta_dict(fasta):

    _dict = fasta.replace('.fa', '') + '.dict'

    inputs = [f'{fasta}']
    outputs = [f'{_dict}']
    options = {}
    spec = f'''
        picard CreateSequenceDictionary     \
            REFERENCE={fasta}               \
            OUTPUT={_dict}
    '''

    return inputs, outputs, options, spec