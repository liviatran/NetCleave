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
    outfile = 'output/' + fasta_name + '_epitopes.csv'
    df.to_csv(outfile, header=True, columns=['epitope_id','epitope','epitope_length','cleavage_site'], index=False)

    return outfile


def generateCleavageSitesUniprot(file,uniprot_data):
    """

    PARAMETERS
    ------------------------------------------------------------
    · file: path to the csv file
    · uniprot_data: Uniprot identifiers and the corresponding sequences,  obtained from local database
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
    accepted_epitopes = []

    for identifier in df['uniprot_id'].unique():

        # Generate FASTA file
        fasta_name = 'output/'+identifier+'.fasta'
        with open(fasta_name,'w') as outfile:
            print('Generating FASTA file: {} ...'.format(fasta_name))
            outfile.write('>'+identifier+'\n')
            outfile.write(uniprot_data[identifier])

    # Obtain cleavage sites
    ne = 0
    for e in epitopes:
        fasta_name = 'output/'+ids[ne]+'.fasta'
        name,sequence = readFasta(fasta_name)
        peptides_dict = generateMERS('output/'+ids[ne]+'.fasta',len(e))
        for key, v in peptides_dict.items():
            if v == e:
                index = key
                flanking_region = sequence[key+len(e):key+len(e)+3]
                cleavage_site = e[-4:] + flanking_region
                cleavage_sites.append(cleavage_site)
                protein_sequences.append(sequence)
        ne+=1
        if e not in accepted_epitopes:
            cleavage_sites.append('N/A')
            protein_sequences.append(sequence)

    # Append cleavage sites to df
    df['cleavage_site'] = cleavage_sites
    df['protein_sequence'] = protein_sequences
    cols.append('cleavage_site')
    file = file.split('/')[-1].split('.')[0]
    outfile = 'output/' + file + '_epitopes.csv'

    df.to_csv(outfile, header=True, columns=cols, index=False)

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
    accepted_epitopes = []
    rejected_epitopes = []

    for i,ps in enumerate(seqs):
        identifier = ids[i]
        # Generate FASTA file
        fasta_name = 'output/'+identifier+'.fasta'
        if not os.path.exists(fasta_name):
            with open(fasta_name,'w') as outfile:
                print('Generating FASTA file: {} ...'.format(fasta_name))
                outfile.write('>'+identifier+'\n')
                outfile.write(ps)

    # Obtain cleavage sites
    ne = 0
    for e in epitopes:
        fasta_name = 'output/'+ids[ne]+'.fasta'
        name,sequence = readFasta(fasta_name)
        peptides_dict = generateMERS('output/'+ids[ne]+'.fasta',len(e))
        for key, v in peptides_dict.items():
            if v == e:
                index = key
                flanking_region = sequence[key+len(e):key+len(e)+3]
                cleavage_site = e[-4:] + flanking_region
                cleavage_sites.append(cleavage_site)
                protein_sequences.append(sequence)
                accepted_epitopes.append(e)
        ne+=1
        if e not in accepted_epitopes:
            rejected_epitopes.append(e)

    # Append cleavage sites to df
    try:
        df['cleavage_site'] = cleavage_sites
    except:
        raise Exception('The following epitope/s are not found in the protein sequence: ',rejected_epitopes)

    cols.append('cleavage_site')
    file = file.split('/')[-1].split('.')[0]
    outfile = 'output/' + file + '_epitopes.csv'

    df.to_csv(outfile, header=True, columns=cols, index=False)

    return outfile
