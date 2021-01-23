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

> NetCleave: an open-source algorithm for predicting C-terminal antigen processing for MHC-I and MHC-II (manuscript in submission)

NetCleave has the following dependencies:

- [argparse](https://docs.python.org/3/library/argparse.html)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [pathlib](https://docs.python.org/3/library/pathlib.html)
- [sklearn](https://scikit-learn.org/stable/)
- [keras](https://keras.io/)
- [tensorflow](https://www.tensorflow.org/)

## How to use NetCleave

<p align="justify">
NetCleave is very easy to use. It has three main functions:

- generate: gets C-terminal data from IEDB and UniProt/UniParc
- train: runs the neural network and saves weights
- predict_csv: scores C-terminal sites

Users can choose between using NetCleave pre-trained models or easily retraining the models:
</p>

### Using pre-trained models

Several pre-trained models are available, which should cover most of the needs of the scientific community. Most common human models available are: HLA-A, HLA-B, HLA-C, HLA-DP, HLA-DQ and HLA-DP.

### Retraining the method and constructing your own models

NetCleave was specificaly build to be easily retrained: IEDB is continously being updated and algorithms should be periodically updated. In order to retrain NetCleave method, user needs first to download the last versions of IEDB and UniParc/UniProt.

