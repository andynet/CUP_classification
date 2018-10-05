#!/bin/python3

from gwf import Workflow
import os


def read_lines(file):
    with open(file) as f:
        lines = f.readlines()
    return lines


def create_bowtie2_index(reference_fasta, cores=1, memory='1g'):
    """Template for creating bowtie2 index files"""

    index_directory = os.path.dirname(os.path.abspath(reference_fasta))
    index_basename = os.path.basename(os.path.abspath(reference_fasta))

    inputs = [reference_fasta]
    outputs = [
        f'{index_directory}/bowtie2_index/{index_basename}.1.bt2',
        f'{index_directory}/bowtie2_index/{index_basename}.2.bt2',
        f'{index_directory}/bowtie2_index/{index_basename}.3.bt2',
        f'{index_directory}/bowtie2_index/{index_basename}.4.bt2',
        f'{index_directory}/bowtie2_index/{index_basename}.rev.1.bt2',
        f'{index_directory}/bowtie2_index/{index_basename}.rev.2.bt2',
    ]
    options = {
        'cores': cores,
        'memory': memory,
    }
    spec = f"""
        mkdir -p {index_directory}/bowtie2_index
        cd {index_directory}/bowtie2_index
        
        bowtie2-build --threads {cores} ../{index_basename} {index_basename}
    """

    return inputs, outputs, options, spec


def bowtie2_fq_to_ref(fastq1, fastq2, reference, bam, cores=1, memory='1g'):
    """Template for mapping fastq files to bowtie2 reference."""

    inputs = [fastq1, fastq2]
    outputs = [bam, f'{bam}.bai']
    options = {
        'cores': cores,
        'memory': memory,
    }
    spec = f"""
        bowtie2 --threads {cores}                                                                                     \
                -x {reference}                                                                                      \
                -1 {fastq1}                                                                                             \
                -2 {fastq2}                                                                                             \
        | samtools sort -@ {cores} - > {bam}

        samtools index {bam}
    """

    return inputs, outputs, options, spec


def main():

    gwf = Workflow()

    working_directory = '/home/andyb/CUP_classification/faststorage/Andrej'
    reference_genome = f'{working_directory}/inputs/hg38.fa'

    gwf.target_from_template('create_index', create_bowtie2_index(reference_genome, cores=8, memory='16g'))

    #
    # normal_fastq = f'{working_directory}/normal_fastq_files.txt'
    # tumor_fastq = f'{working_directory}/tumor_fastq_files.txt'
    #
    # bowtie2_index = f'{working_directory}/bowtie2_index/hg38'
    # bed_file = f'{working_directory}/inputs/covered_regions.bed'
    #
    # out_dir = f'{working_directory}/outputs'
    # start = 0
    # end = 1
    # threads = 8
    #
    # normal_fastq_lines = read_lines(normal_fastq)
    # tumor_fastq_lines = read_lines(tumor_fastq)
    #
    # locations = [normal_fastq_lines, tumor_fastq_lines]
    #


if __name__ == '__main__':
    main()
