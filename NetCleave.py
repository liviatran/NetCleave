import argparse
from predictor.database_functions import peptide_extractor, uniprot_extractor, uniparc_extractor
from predictor.core import all_peptide_uniprot_locator, all_training_data_generator, cleavage_site_generator
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
    parser.add_argument('--input_fasta',
                            dest = 'input_fasta',
                            help='Path to the input FASTA file',
                            action='store')
    parser.add_argument('--mhc_class',
                            dest = 'mhc_class',
                            help='Major Histocompatibility Complex class',
                            action='store',default='I')
    parser.add_argument('--mhc_family',
                            dest = 'mhc_family',
                            help='Major Histocompatibility Complex allele',
                            action='store',
                            default='HLA')
    parser.add_argument('--peptide_data',
                            dest = 'peptide_data',
                            help='Path to peptide data to use for generating data to later train the model',
                            action='store',
                            default='./data/databases/iedb/mhc_ligand_full.csv')
    parser.add_argument('--peptide_data_additional',
                            dest = 'peptide_data_additional',
                            help='Path to additional peptide data to use (combined with IEDB data) for generating data to later train the model',
                            action='store',
                            default='./data/databases/other/HLA.csv')
    parser.add_argument('--technique',
                            dest = 'technique',
                            help='',
                            action='store',
                            default='mass spectrometry')
    parser.add_argument('--train',
                            dest = 'train',
                            help='Train the neural network',
                            action='store_true')
    parser.add_argument('--type',
                            dest = 'type',
                            help='',
                            action='store',default=1,type=int)
    parser.add_argument('--score_csv',
                            dest = 'score_csv',
                            help='Predict a set of cleavage sites from csv',
                            action='store_true')
    return parser.parse_args()


def generating_data(uniprot_path, uniparc_path_headers, uniparc_path_sequence, type, iedb_path=None, conditions=None, other_path=None):
    """
    Generate training data that will be later used to retrain the neural network.
    Returns a dictionary with the selected peptides (key: C-terminal residue, value: peptide list)

    PARAMETERS
    ------------------------------------------------------------
    · uniprot_path: path to the UNIPROT data
    · uniparc_path_headers: path to the UNIPARC headers data
    · uniparc_path_sequence: path to the UNIPARC sequence data

    · conditions: dictionary with the conditions that the IEDB data must fulfill
    ------------------------------------------------------------
    """

    # Get data from files
    if type==1:
        peptide_data = peptide_extractor.extract_peptide_data(iedb_path,conditions,iedb=True)
    elif type==2:
        peptide_data1 = peptide_extractor.extract_peptide_data(iedb_path,conditions,iedb=True)
        peptide_data2 = peptide_extractor.extract_peptide_data(other_path,iedb=False)
        peptide_data = peptide_extractor.merge_peptide_data(peptide_data1,peptide_data2)
    elif type==3:
        peptide_data = peptide_extractor.extract_peptide_data(other_path,iedb=False)

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
        print('Generating training data, type {}...'.format(type))
        uniprot_path = 'data/databases/uniprot/uniprot_sprot.fasta' # download and decompress from https://www.uniprot.org/downloads REVIEWED fasta
        uniparc_path_headers = 'data/databases/uniparc/uniparc-yourlist_M20200416A94466D2655679D1FD8953E075198DA854EB3ES.tab'
        uniparc_path_sequence = 'data/databases/uniparc/uniparc-yourlist_M20200416A94466D2655679D1FD8953E075198DA854EB3ES.fasta'
        if type==1:
            peptide_path = peptide_data
            iedb_conditions = {
                                'Description': None, 'Parent Protein IRI': None,
                                'Method/Technique': ('contains', technique),
                                'MHC allele class': ('match', mhc_class),
                                'Allele Name': ('contains', mhc_family),
                                 #'Name': ('contains', 'Homo sapiens'),
                                 #'Parent Species': ('contains', 'Homo sapiens')
                                 }
            selected_dictionary = generating_data(uniprot_path,uniparc_path_headers,uniparc_path_sequence,type=1,iedb_path=peptide_path,conditions=iedb_conditions)
        elif type==2:
            peptide_path = peptide_data
            peptide_path2 = peptide_data_additional
            iedb_conditions = {
                                'Description': None, 'Parent Protein IRI': None,
                                'Method/Technique': ('contains', technique),
                                'MHC allele class': ('match', mhc_class),
                                'Allele Name': ('contains', mhc_family),
                                 #'Name': ('contains', 'Homo sapiens'),
                                 #'Parent Species': ('contains', 'Homo sapiens')
                                 }
            selected_dictionary = generating_data(uniprot_path,
                                                      uniparc_path_headers,
                                                      uniparc_path_sequence,
                                                      type=2,
                                                      iedb_path=peptide_path,
                                                      conditions=iedb_conditions,
                                                      other_path=peptide_path2)
        elif type==3:
            peptide_path = peptide_data
            selected_dictionary = generating_data(uniprot_path,
                                                      uniparc_path_headers,
                                                      uniparc_path_sequence,
                                                      type=3,
                                                      other_path=peptide_path)


        all_training_data_generator.prepare_cleavage_data(selected_dictionary, training_data_path)

    if train:
        run_NN.create_models(training_data_path, models_export_path)

    if score_csv:
        # predict_csv.score_set(input_csv, models_export_path, 'ABC')
        outfile = cleavage_site_generator.generateCleavageSites(input_fasta)
        predict_csv.score_set(outfile, models_export_path, 'ABC')


if __name__ == '__main__':

    # Get arguments
    arguments = parse_args()
    generate = arguments.generate
    input_csv = arguments.input_csv
    input_fasta = arguments.input_fasta
    mhc_class = arguments.mhc_class
    mhc_family = arguments.mhc_family
    peptide_data = arguments.peptide_data
    peptide_data_additional = arguments.peptide_data_additional
    technique = arguments.technique
    train = arguments.train
    type = arguments.type
    score_csv = arguments.score_csv

    # Call main function
    main(generate, train, score_csv)
