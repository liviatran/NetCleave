import argparse
import os
import re
import pandas as pd
from Bio import SeqIO
import requests as r
from io import StringIO
from predictor.database_functions import uniprot_extractor

def readFasta(file):
    """
    Read a Fasta file and return the protein identifier and sequence.

    PARAMETERS
    ------------------------------------------------------------
    · file: path to the Fasta file
    ------------------------------------------------------------
    """

    for seq_record in SeqIO.parse(file, "fasta"):
        name = seq_record.id
        sequence = seq_record.seq

    return name,sequence


def generateMERS(file,length):
    """
    Generate short peptides from the whole sequence of a protein.

    PARAMETERS
    ------------------------------------------------------------
    · file: path to the Fasta file
    · length: peptide length
    ------------------------------------------------------------
    """

    name,sequence = readFasta(file)

    peptides_mers = {}
    for index,residue in enumerate(sequence):
        peptide = sequence[index:index+length]
        if len(peptide)<length:
            break
        peptides_mers[index] = peptide

    return peptides_mers


def generateCleavageSites(file,mhc=None,custom_length=None):
    """
    Create a csv file with the cleavage sites of a set of peptides (generated using
    the Fasta file of a protein of interest). The cleavage sites contain the necessary
    sequence for NetCleave to predict the cleavage of a peptide, they consist in
    short sequences of 4+3 amino acids.

    PARAMETERS
    ------------------------------------------------------------
    · file: path to the Fasta file
    · mhc: Major Histocompatibility Complex class. If I, generate 8-11 mers; if II, generate 13-17 mers.
    · custom_length:
    ------------------------------------------------------------
    """

    name,sequence = readFasta(file)

    epitopes_id = []
    epitopes_seq = []
    epitopes_len = []
    cleavage_sites = []

    if mhc=='I': # Generate peptides 8 to 11 amino acids long
        length_range = range(8,12)
    if mhc=='II': # Generate peptides 8 to 11 amino acids long
        length_range = range(13,18)
    if custom_length!=None:
        length_range = range(custom_length,custom_length+1)

    for i in length_range:
        epitope_number = 1
        for key in generateMERS(file,i):
            seq = generateMERS(file,i)[key]
            flanking_region = sequence[key+i:key+i+3]
            cleavage_site = seq[-4:] + flanking_region
            if len(cleavage_site)<7:
                break
            identifier = str(i)+'mer|'+str(epitope_number)+'|'+seq+'|'+name
            epitopes_id.append(identifier)
            epitopes_seq.append(str(seq))
            cleavage_sites.append(str(cleavage_site))
            epitopes_len.append(i)
            epitope_number += 1

    df = pd.DataFrame(list(zip(epitopes_id,epitopes_seq,epitopes_len,cleavage_sites)),columns=['epitope_id','epitope','epitope_length','cleavage_site'])

    if not os.path.exists('./output/'):
        os.mkdir('./output/')

    fasta_name = file.split('/')[-1].split('.')[0]
    outfile = 'output/' + fasta_name + '.csv'
    df.to_csv(outfile, header=True, columns=['epitope_id','epitope','epitope_length','cleavage_site'], index=False)

    return outfile


def cleavageMotif(sequence, peptide, start):
    length = len(peptide) -4
    start+=length
    return sequence[start:start+7]


def retrieveSequenceFromUniprot(protein_code):

    try:
        url="http://www.uniprot.org/uniprotkb/"+protein_code+".fasta"
        response = r.post(url)
        data=''.join(response.text)
        seq=StringIO(data)
        seq_list=list(SeqIO.parse(seq,'fasta'))
        protein_seq = str(seq_list[0].seq)
        return protein_seq
    except IndexError:
        return 0


def generateCleavageSitesUniprot(file,uniprot_data):
    """

    PARAMETERS
    ------------------------------------------------------------
    · file: path to the csv file
    · uniprot_data: Uniprot identifiers and the corresponding sequences, obtained from local database
    ------------------------------------------------------------
    """
    if not os.path.exists('./output/'):
        os.mkdir('./output/')

    df = pd.read_csv(file)
    cols = list(df.columns.values)
    ids = df['uniprot_id'].values
    epitopes = df['epitope'].values

    cleavage_sites = []
    protein_sequences = []
    warnings = []
    epitope = []
    uniprot_ids = []


    uniprot_path = 'data/databases/uniprot/uniprot_sprot.fasta' # download and decompress from https://www.uniprot.org/downloads REVIEWED fasta
    changes_in_db = 0
    for identifier in df['uniprot_id'].unique():
        try:
            uniprot_data[identifier]
        except:
            # Append sequence information for identifier not available in local database
            try:
                sequence = retrieveSequenceFromUniprot(identifier)
                if sequence != 0:
                    with open(uniprot_path,'a') as outfile:
                        print('---> Adding {} to local UniProt database...'.format(identifier))
                        outfile.write('>sp|'+identifier+'|added_by_NetCleave\n')
                        outfile.write(sequence+'\n')
                        changes_in_db+=1

            except:
                # avoid TypeError when UniProt is blank
                pass

    ## Update UniProt data in case new sequences are added to local database
    if changes_in_db>0:
        uniprot_data = uniprot_extractor.extract_uniprot_data(uniprot_path)

    ## Obtain cleavage sites
    for i,e in enumerate(epitopes):
        try:
            fasta_id = ids[i]
            sequence = uniprot_data[fasta_id]
            epitope_match = re.finditer(str(e), str(sequence)) # find peptide in sequence
            indices = [m.start(0) for m in epitope_match] # get start index
            if len(indices)==0:
                epitope.append(e)
                uniprot_ids.append(fasta_id)
                cleavage_sites.append('nan')
                protein_sequences.append(sequence)
                warnings.append('epitope_not_found_in_protein_sequence')
            else:
                for ind in indices:
                    cleavage_site=cleavageMotif(sequence,e,ind)
                    if len(cleavage_site)==7:
                        epitope.append(e)
                        uniprot_ids.append(fasta_id)
                        cleavage_sites.append(cleavage_site)
                        protein_sequences.append(sequence)
                        warnings.append('-')
                    else:
                        epitope.append(e)
                        uniprot_ids.append(fasta_id)
                        cleavage_sites.append('nan')
                        protein_sequences.append(sequence)
                        warnings.append('cannot_generate_cleavage_motif')
        except:
            epitope.append(e)
            uniprot_ids.append(fasta_id)
            cleavage_sites.append('nan')
            protein_sequences.append('nan')
            warnings.append('protein_not_found_in_uniprot')
    ## Append cleavage sites to df
    df = pd.DataFrame(columns=['epitope','cleavage_site', 'protein_sequence', 'warnings'])
    df['epitope'] = epitope
    df['uniprot_id'] = uniprot_ids
    df['cleavage_site'] = cleavage_sites
    df['protein_sequence'] = protein_sequences
    df['warnings'] = warnings
    cols.append('cleavage_site')
    cols.append('warnings')
    file = file.split('/')[-1].split('.')[0]
    outfile = 'output/' + file + '.csv'
    df.to_csv(outfile, header=True, columns=['epitope','uniprot_id','cleavage_site','protein_sequence','warnings'], index=False)
    # df.to_csv(outfile, header=True, columns=cols, index=False)

    return outfile


def generateCleavageSitesSequence(file):
    """

    PARAMETERS
    ------------------------------------------------------------
    · file: path to the csv file
    ------------------------------------------------------------
    """
    if not os.path.exists('./output/'):
        os.mkdir('./output/')

    df = pd.read_csv(file)
    cols = list(df.columns.values)

    ids = df['protein_name'].values
    seqs = df['protein_seq'].values
    epitopes = df['epitope'].values
    cleavage_sites = []
    protein_sequences = []
    warnings = []
    epitope = []
    uniprot_ids = []

    ## Obtain cleavage sites
    for i,e in enumerate(epitopes):
        try:
            sequence = seqs[i]
            fasta_id = ids[i]
            epitope_match = re.finditer(str(e), str(sequence)) # find peptide in sequence
            indices = [m.start(0) for m in epitope_match] # get start index
            if len(indices)==0:
                epitope.append(e)
                uniprot_ids.append(fasta_id)
                cleavage_sites.append('nan')
                protein_sequences.append(sequence)
                warnings.append('epitope_not_found_in_protein_sequence')
            else:
                for ind in indices:
                    cleavage_site=cleavageMotif(sequence,e,ind)
                    if len(cleavage_site)==7:
                        epitope.append(e)
                        uniprot_ids.append(fasta_id)
                        cleavage_sites.append(cleavage_site)
                        protein_sequences.append(sequence)
                        warnings.append('-')
                    else:
                        epitope.append(e)
                        uniprot_ids.append(fasta_id)
                        cleavage_sites.append('nan')
                        protein_sequences.append(sequence)
                        warnings.append('cannot_generate_cleavage_motif')
        except:
            epitope.append(e)
            uniprot_ids.append(fasta_id)
            cleavage_sites.append('nan')
            protein_sequences.append('nan')
            warnings.append('incorrect_protein_sequence')
    # Append cleavage sites to df
    df = pd.DataFrame(columns=['epitope','cleavage_site', 'protein_sequence', 'warnings'])
    df['epitope'] = epitope
    df['uniprot_id'] = uniprot_ids
    df['cleavage_site'] = cleavage_sites
    df['protein_sequence'] = protein_sequences
    df['warnings'] = warnings
    cols.append('cleavage_site')
    cols.append('warnings')
    file = file.split('/')[-1].split('.')[0]
    outfile = 'output/' + file + '.csv'
    df.to_csv(outfile, header=True, columns=['epitope','uniprot_id','cleavage_site','protein_sequence','warnings'], index=False)
    # df.to_csv(outfile, header=True, columns=cols, index=False)

    return outfile
