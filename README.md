# Combined HLA Prediction & Viral Mimicry Analysis
Computational Immunology independent study undertaken by Matt Jogodnik under the guidance of Cliburn Chan, PHD and Annette Jackson, PhD, F(ACHI)

Goal: Elucidating peptide epitopes of interest in the autoimmune condition nephrotic syndrome binding to HLA alleles associated with higher disease risk.

## Description
This program functions as a tool to predict peptide epitopes of interest to autoimmune or other diseases associated with specific HLA Class II alleles. It leverages the top 2 HLA Class II peptide binding prediction algorithms, MixMHC2pred and MARIA, and takes the intersection of their results to increase specificity. It also compares all binding peptides against a separate cohort of viral peptides to assess for potential regions of viral mimicry in the targeted genome.

## Instructions
All functions are contained within engine.py, which is already configured to run with default MARIA and MixMHC2pred thresholds of 95.00% and 2.00% respectively. Inputs are configured within the input folder. Sample files are provided for formatting purposes.

## Environment
The following is a list of all packages/versions that must be installed. I recommend creating a Conda environment as you'll need Python 2.7 as well, which is well-outdated but required for MARIA to run.

- Python 2.7
- Keras 2.0.3 (version specific)
- Pandas
- Numpy
- Scipy
- mkl
- Tensorflow
- Theano
- Python-Levenshtein
- Biopython

## The Legal Stuff
MixMHC2pred and MARIA were created by Racle et al. and Chen et al. respectively, which as specified under their licenses are free to use for research purposes. Links to their papers as well as proper citations can be found below.

All other code was written by Matt Jogodnik and is free to be used for non-commercial research purposes consistent with that of MixMHC2pred and MARIA.

### MixMHC2pred

Racle, J., *et al.* Robust prediction of HLA class II epitopes by deep motif deconvolution of immunopeptidomes. *Nat. Biotechnol* **37**, 1283–1286 (2019).

Link: [MixMHCpred](https://www.nature.com/articles/s41587-019-0289-6)

### MARIA

Chen, B., Khodadoust, M.S., Olsson, N. *et al.* Predicting HLA class II antigen presentation through integrated deep learning. *Nat Biotechnol* **37**, 1332–1343 (2019).

Link: [MARIA](https://www.nature.com/articles/s41587-019-0280-2)