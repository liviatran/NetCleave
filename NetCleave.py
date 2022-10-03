import argparse
from predictor.database_functions import peptide_extractor, uniprot_extractor, uniparc_extractor
from predictor.core import all_peptide_uniprot_locator, all_training_data_generator
from predictor.ml_main import run_NN
from predictor.predictions import predict_csv

HELP = ' \
Command:\n \
----------\n \
Run: python3 NetCleave.py --ARG\
'

def parse_args():
    """
    Parse command-line arguments given by the user.
    """
    parser = argparse.ArgumentParser(description = 'Arguments necessary for the main program execution.')
    parser.add_argument('--generate',
                            dest = 'generate',
                            help='Generate training data for the neural network',
                            action='store_true')
    parser.add_argument('--input_csv',
                            dest = 'input_csv',
                            help='Path to the input csv file',
                            action='store',
                            default='input/example_file_NetCleave_score.csv')
    parser.add_argument('--mhc_class',
                            dest = 'mhc_class',
                            help='Major Histocompatibility Complex class',
                            action='store',default='I')
    parser.add_argument('--mhc_family',
                            dest = 'mhc_family',
                            help='Major Histocompatibility Complex allele',
                            action='store',
                            default='HLA-A')
    parser.add_argument('--peptide_data',
                            dest = 'peptide_data',
                            help='Path to peptide data to use for generating data to later train the model',
                            action='store',
                            default='./data/databases/iedb/mhc_ligand_full.csv')
    parser.add_argument('--technique',
                            dest = 'technique',
                            help='',
                            action='store',
                            default='mass spectrometry')
    parser.add_argument('--train',
                            dest = 'train',
                            help='Train the neural network',
                            action='store_true')
    parser.add_argument('--score_csv',
                            dest = 'score_csv',
                            help='Predict a set of cleavage sites from csv',
                            action='store_true')
    return parser.parse_args()


def generating_data(peptide_path, uniprot_path, uniparc_path_headers, uniparc_path_sequence, conditions):
    """
    Generate training data that will be later used to retrain the neural network.
    Returns a dictionary with the selected peptides (key: C-terminal residue, value: peptide list)

    PARAMETERS
    ------------------------------------------------------------
    · peptide_path: path to the peptide data (IEDB or others)
    · uniprot_path: path to the UNIPROT data
    · uniparc_path_headers: path to the UNIPARC headers data
    · uniparc_path_sequence: path to the UNIPARC sequence data
    · conditions: dictionary with the conditions that the IEDB data must fulfill
    ------------------------------------------------------------
    """

    # Get data from files
    if './data/databases/iedb/' in peptide_path:
        peptide_data = peptide_extractor.extract_peptide_data(peptide_path, conditions)
    else:
        peptide_data = peptide_extractor.extract_peptide_data(peptide_path, conditions, iedb=False)
    uniprot_data = uniprot_extractor.extract_uniprot_data(uniprot_path)
    uniparc_data = uniparc_extractor.extract_uniparc_data(uniparc_path_headers, uniparc_path_sequence)
    sequence_data = all_peptide_uniprot_locator.join_data(uniprot_data, uniparc_data)

    # Generate dictionary
    selected_dictionary = all_peptide_uniprot_locator.locate_peptides(peptide_data, sequence_data)

    return selected_dictionary


def main(generate=False, train=False, score_csv=False):
    """
    This is the main function of the program.

    PARAMETERS
    ------------------------------------------------------------
    · generate: generate training data for the neural network
    · train: train the neural network with previously generated data
    · score_csv: predict a set of cleavage sites from a csv file
    ------------------------------------------------------------
    """

    training_data_path = 'data/training_data/{}_{}_{}'.format(mhc_class, technique.replace(' ', '-'), mhc_family)
    models_export_path = 'data/models/{}_{}_{}'.format(mhc_class, technique.replace(' ', '-'), mhc_family)

    if not any([generate, train, score_csv]):
        print('Please, provide an argument. See python3 NetCleave.py -h for more information')

    if generate:
        peptide_path = peptide_data
        uniprot_path = 'data/databases/uniprot/uniprot_sprot.fasta' # download and decompress from https://www.uniprot.org/downloads REVIEWED fasta
        uniparc_path_headers = 'data/databases/uniparc/uniparc-yourlist_M20200416A94466D2655679D1FD8953E075198DA854EB3ES.tab'
        uniparc_path_sequence = 'data/databases/uniparc/uniparc-yourlist_M20200416A94466D2655679D1FD8953E075198DA854EB3ES.fasta'

        conditions = {
                    'Description': None, 'Parent Protein IRI': None,
                    # 'Method/Technique': ('contains', technique),
                    # 'MHC allele class': ('match', mhc_class),
                    # 'Allele Name': ('contains', mhc_family),
                     #'Name': ('contains', 'Homo sapiens'),
                     #'Parent Species': ('contains', 'Homo sapiens')
                     }
        selected_dictionary = generating_data(peptide_path, uniprot_path, uniparc_path_headers, uniparc_path_sequence, conditions)
        all_training_data_generator.prepare_cleavage_data(selected_dictionary, training_data_path)

    if train:
        run_NN.create_models(training_data_path, models_export_path)

    if score_csv:
        predict_csv.score_set(input_csv, models_export_path, 'ABC')


if __name__ == '__main__':

    # Get arguments
    arguments = parse_args()
    generate = arguments.generate
    input_csv = arguments.input_csv
    mhc_class = arguments.mhc_class
    mhc_family = arguments.mhc_family
    peptide_data = arguments.peptide_data
    technique = arguments.technique
    train = arguments.train
    score_csv = arguments.score_csv

    # Call main function
    main(generate, train, score_csv)
