# coding=utf-8
"""This program take a list of locus tags or protein id, and extracts their nucleotide sequence from a full genbank
file.  This could have been done through NCBI's Batch Entrez but for my purposes I feel this is faster and possible
more flexible.


"""
import os
import argparse
from Bio import SeqIO
from collections import namedtuple


def get_options():
    """ Returns the user's command line arguments in an argparse object """
    parser = argparse.ArgumentParser()

    parser.add_argument('-list', required=True, help='Input List of unique identifiers to extract')
    parser.add_argument('-gb', required=True, help='Genbank file containing the unique identifiers to extract')
    # parser.add_argument('-type', choices=['locus_tag', 'protein_id'], required=True,
    # help='Type of unique identifier you are searching with: protein id, locus tag, etc.')
    return parser.parse_args()


def _next_in_request_list(input_list):
    """ Returns the next line in the input list """
    with open(input_list, 'U') as user_list:
        for line in user_list:
            yield line.rstrip()


def _get_genbank_dict(genbank_file):
    """
        Return a dict containing the genbank features

    :param genbank_file:
    """
    return SeqIO.parse(genbank_file, 'genbank').next().__dict__


def _find_item(gb_dict, item_id):
    # named tuple for storing the SeqFeature and locus id
    """
        Search a genbank dictionary for a feature with our target.  The target can be either a locus tag or protein id.

    :param gb_dict:
    :param item_id:
    :return:
    """
    genome_feature = namedtuple('GenomeFeature', ['locus', 'data'])

    # ---------------------------------------------------------------------
    # Sub-function to cleanup the output by creating a Named Tuple.
    # ---------------------------------------------------------------------
    def _make_target(seq_feature):
        tag = seq_feature.qualifiers['locus_tag'][0]
        return genome_feature(locus=tag, data=seq_feature)

    # ---------------------------------------------------------------------
    # Find our item_id target in the genome by looping through each feature
    # and checking if our item is present in the qualifier records.
    # ---------------------------------------------------------------------
    current_target = {'len': 0, 'feature': None}
    for feature in gb_dict['features']:
        # Skip non-locus features
        if 'locus_tag' not in feature.qualifiers:
            continue
        # search for our item in the qualifiers, skip the feature if we don't find our term
        for k in feature.qualifiers:
            if item_id.split('.')[0] == feature.qualifiers[k][0].split('.')[0]:
                # Found item, break to skip the continue in the else clause
                break
        else:
            # No break called, therefore this feature doesn't contain what we need.  Skip
            continue

        # Only keep the feature with the most data.
        # Gene types have less content than their CDS, rRNA,
        # and ncRNA counter parts
        this_feature_size = len(feature.qualifiers)
        if this_feature_size > current_target['len']:
            current_target['len'] = this_feature_size
            current_target['feature'] = feature

    if current_target['feature'] is not None:
        return _make_target(current_target['feature'])
    else:
        return


def main(genbank_file=None, target_list=None):
    """ Extract nucleotide sequences from a genbank file using a list of targets to grab.
     Output results in fasta format

    """
    if genbank_file is None or target_list is None:
        options = get_options()
        genbank_file = options.gb
        target_list = options.list

    genbank_dict = _get_genbank_dict(genbank_file)

    def _print_fasta(nucl, genome_feature):
        """
            Prints out a nucleotide feature in a fasta format
        :param genome: String sequence of the nucleotides from the genbank file
        :param genome_feature: NamedTuple(locus, SeqFeature)
        """

        # For all
        locus_tag = genome_feature.locus

        # Get Protein Accession ID and GI number
        protein_id = gi_number = product = None
        if 'CDS' in genome_feature.data.type:
            protein_id = genome_feature.data.qualifiers['protein_id'][0]
            # get GI number if it's available
            for db_xref in genome_feature.data.qualifiers['db_xref']:
                if str(db_xref).startswith('GI:'):
                    gi_number = db_xref
                    break
        # Get Product description if available
        if 'product' in genome_feature.data.qualifiers:
            product = genome_feature.data.qualifiers['product'][0]

        start = genome_feature.data.location.start
        end = genome_feature.data.location.end
        # Fixme this script isn't taking into account the strand when cutting out the sequence
        # strand = genome_feature.data.location.strand
        seq = nucl[start:end]

        return ('>{locus_tag} | {protein_id} | {gi_number} | {product}\n'
                '{seq}\n'.format(locus_tag=locus_tag,
                                 protein_id=protein_id,
                                 gi_number=gi_number,
                                 product=product,
                                 seq=seq))

    #
    # Loop through the provided list, and build a unique list of features
    result_basket = {}
    for target in _next_in_request_list(target_list):
        # Extract the target from the item on the list
        feature = _find_item(genbank_dict, target)

        if feature is None:
            continue

        if feature.locus in result_basket:
            # only add if the current SeqFeature is a CDS type since it has more
            # information than a gene type
            if 'CDS' not in feature.data.type:
                continue
        result_basket[feature.locus] = feature

    for k in result_basket:
        yield _print_fasta(genbank_dict['_seq'], result_basket[k])


if __name__ == '__main__':
    for result in main():
        print(result)