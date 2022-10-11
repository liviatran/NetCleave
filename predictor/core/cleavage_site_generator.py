import argparse
import os
import pandas as pd
from Bio import SeqIO


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


def generateCleavageSites(file):
    """
    Create a csv file with the cleavage sites of a set of peptides (generated using
    the Fasta file of a protein of interest). The cleavage sites contain the necessary
    sequence for NetCleave to predict the cleavage of a peptide, they consist in
    short sequences of 4+3 amino acids.

    PARAMETERS
    ------------------------------------------------------------
    · file: path to the Fasta file
    ------------------------------------------------------------
    """

    name,sequence = readFasta(file)

    epitopes_id = []
    epitopes_seq = []
    epitopes_len = []
    cleavage_sites = []

    # Generate peptides 8 to 11 amino acids long
    for i in range(8,12):
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
    outfile = 'output/' + fasta_name + '_epitopes.csv'
    df.to_csv(outfile, header=True,index=True)

    return outfile
