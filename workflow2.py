#!/bin/python3

from gwf import Workflow


def read_lines(file):
    with open(file) as f:
        lines = f.readlines()
    return lines


def create_bowtie2_index(reference):
    return 0, 0


def map_fastqs_to_reference(workflow, name, fastq1, fastq2, reference, bam, threads):

    main_index, indices = create_bowtie2_index(reference)

    workflow.target(
        name=name,
        inputs=[fastq1, fastq2] + indices,
        outputs=[bam],
        options={'cores': threads},
    ) << f'''
        bowtie2 --threads {threads}                                                                                     \
                -x {main_index}                                                                                      \
                -1 {fastq1}                                                                                             \
                -2 {fastq2}                                                                                             \
        | samtools sort -@ {threads} - > {bam}

        samtools index {bam}
    '''


def main():

    gwf = Workflow()
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

    normal_fastq_lines = read_lines(normal_fastq)
    tumor_fastq_lines = read_lines(tumor_fastq)

    locations = [normal_fastq_lines, tumor_fastq_lines]


if __name__ == '__main__':
    main()
