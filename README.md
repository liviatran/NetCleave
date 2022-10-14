# NetCleave

<p align="justify">
NetCleave is a retrainable method for predicting C-terminal peptide processing of MHC-I and MHC-II pathways.
</p>

<p align="center">
<img src="images/draw_scheme_method.png" width="600">
</p>

<p align="justify">
In brief, NetCleave maps reported IEDB peptides to protein sequences in UniProt/UniParc. After the identification of the C-terminal cleavage site, amino acid sequences are coded using QSAR descriptors, including steric, electrostatic and hydrophobic properties. Finally, a neural network architecture is used to generate the predictive model.
</p>

If you use NetCleave, please cite us:

> <p align="justify"> Amengual-Rigo, P., Guallar, V. NetCleave: an open-source algorithm for predicting C-terminal antigen processing for MHC-I and MHC-II. Sci Rep 11, 13126 (2021). https://doi.org/10.1038/s41598-021-92632-y
</p>

NetCleave has the following dependencies:

- [argparse](https://docs.python.org/3/library/argparse.html)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [pathlib](https://docs.python.org/3/library/pathlib.html)
- [sklearn](https://scikit-learn.org/stable/)
- [keras](https://keras.io/)
- [tensorflow](https://www.tensorflow.org/)
- [biopython](https://biopython.org/)

## How to use NetCleave

<p align="justify">

NetCleave is very easy to use. It has two main functionalities, which are **predict** and **retrain**.

### Predict

The scoring option can be used to predict the C-terminal cleavage of peptides. It allows three different types of input:

1. FASTA file of a single protein, from which epitopes (8 to 11 residue-long) will be generated and scored.

2. CSV file with peptide sequences to predict (column name: *epitope*) and the UniProt identifier of the protein where they come (column name: *uniprot_id*). NetCleave will retrieve the sequence of the whole protein, and create the 4+3 sequences that conform the cleavage site and are necessary for the scoring method.

3. CSV file with peptide sequences to predict (column name: *epitope*), the identifier (column name: *protein_id*) and sequence of the protein where they come (column name: *protein_seq*). Same as type 2, those complete protein sequences will be used to create the 4+3 sequences.

To run the **predict** option, please use the arguments **--predict**, followed by the path to the input file, and **input_type**, followed by the type of input you are using (previous types explained 1,2,3).

```
python NetCleave.py --predict (path to file) --input_type (type of input file)
```

As an example, if you use the file `input/toy1.fasta`, which is a FASTA file, you have to run the following command:

```
python NetCleave.py --predict input/toy1.fasta --input_type 1
```

### Retrain

One of the main advantages of having a neural network model is that it can be retrained, so it is updated and adapt to your specific goal. To do so, you can choose between:

- **Generating your own custom data** and then **training** the model with it.

- **Using pre-trained models** during the **training** of the model.

</p>

> <p align="justify"> User should check the quality of the model. If the loss performance between the training and testing groups differ substantially or if there is any sign of overfitting, the user should modify predictor/ml_main/run_NN.py script. Usually, this phenomena happens because a very small dataset is used (a few peptides), which is not enough for building a high quality model. If this happens, consider to generate a more general model (for instance, instead of HLA-A0201, use HLA-A02 or HLA-A).
</p>
</p>
