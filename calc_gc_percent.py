"""
    Simple script that outputs the name, length, GC count, and GC percent of a sequence in a fasta format.

    Example:
        NODE_360_length_156_cov_88.910255       2,019,315       781,621 38.71%
        NODE_361_length_33_cov_92.575760        2,019,380       781,644 38.71%
        NODE_362_length_14453_cov_13.237667     2,033,865       787,505 38.72%
        NODE_363_length_47_cov_445.914886       2,033,944       787,549 38.72%

"""
import argparse


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-fasta', required=True, help='Sequences in Fasta format')
    return parser.parse_args()


def get_counts(fasta_input):
    gc_count = 0
    genome_length = 0

    with open(fasta_input, 'r') as myfile:
        title = None
        for line in myfile:
            if line.startswith('>'):
                line = line.rstrip()

                # save the old record
                if title is not None:
                    contig = {'title': title,
                              'gc_count': '{:,}'.format(gc_count),
                              'length': '{:,}'.format(genome_length),
                              'percent': "{0:.2f}".format(float(gc_count) / genome_length * 100)}
                    yield contig
                # start a new record
                title = line[1:]
                continue

            line = line.rstrip().upper()
            gc_count += line.count('G') + line.count('C')
            genome_length += len(line)


def main():
    options = get_options()
    fasta_input = options.fasta

    for count, record in enumerate(get_counts(fasta_input)):
        print('{title}\t{length}\t{gc_count}\t{percent}%'.format(title=record['title'],
                                                                 length=record['length'],
                                                                 gc_count=record['gc_count'],
                                                                 percent=record['percent']))

if __name__ == '__main__':
    main()
