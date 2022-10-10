import argparse
import os
import pandas as pd
from Bio import SeqIO


def readFasta(file):

    for seq_record in SeqIO.parse(file, "fasta"):
        sequence = seq_record.seq

    return sequence


def generateMERS(file,length):

    sequence = readFasta(file)

    peptides_mers = {}
    for index,residue in enumerate(sequence):
        peptide = sequence[index:index+length]
        if len(peptide)<length:
            break
        peptides_mers[index] = peptide

    return peptides_mers


def generateCleavageSites(file):

    protein_sequence = readFasta(file)

    peptides = []
    cleavage_sites = []
    peptides_description = []

    for key in generateMERS(file,8):
        peptide = generateMERS(file,8)[key]
        flanking_region = protein_sequence[key+8:key+8+3]
        cleavage_site = peptide[-4:] + flanking_region
        if len(cleavage_site)<7:
            break
        peptides.append(str(peptide))
        cleavage_sites.append(str(cleavage_site))
        peptides_description.append('8mer')

    for key in generateMERS(file,9):
        peptide = generateMERS(file,9)[key]
        flanking_region = protein_sequence[key+9:key+9+3]
        cleavage_site = peptide[-4:] + flanking_region
        if len(cleavage_site)<7:
            break
        peptides.append(str(peptide))
        cleavage_sites.append(str(cleavage_site))
        peptides_description.append('9mer')

    for key in generateMERS(file,10):
        peptide = generateMERS(file,10)[key]
        flanking_region = protein_sequence[key+10:key+10+3]
        cleavage_site = peptide[-4:] + flanking_region
        if len(cleavage_site)<7:
            break
        peptides.append(str(peptide))
        cleavage_sites.append(str(cleavage_site))
        peptides_description.append('10mer')

    for key in generateMERS(file,11):
        peptide = generateMERS(file,11)[key]
        flanking_region = protein_sequence[key+11:key+11+3]
        cleavage_site = peptide[-4:] + flanking_region
        if len(cleavage_site)<7:
            break
        peptides.append(str(peptide))
        cleavage_sites.append(str(cleavage_site))
        peptides_description.append('11mer')

    df = pd.DataFrame(list(zip(peptides,peptides_description,cleavage_sites)),columns=['epitope','epitope_length','cleavage_site'])

    if not os.path.exists('./output/'):
        os.mkdir('./output/')

    fasta_name = file.split('/')[-1].split('.')[0]
    outfile = 'output/' + fasta_name + '_epitopes.csv'
    df.to_csv(outfile, header=True,index=True)

    return outfile
