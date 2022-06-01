import os
import csv
from tkinter import W
from Bio import SeqIO
from Bio import Entrez

"""
Args:
    email (string): the identifying email used for accessing Entrez
    commonName (string): the protein's commonly referred to name
    proteinUID (string): the uid of the protein on Entrez's server

Returns:
    N/A, but downloads the fasta file of the identified protein into the proteins folder
"""
def downloadSeq(email, commonName, proteinUID):
    Entrez.email = email
    filename = "input/proteins/" + commonName + ".fasta"
    if not os.path.isfile(filename):
        # Downloading...
        net_handle = Entrez.efetch(
            db="protein", id=proteinUID, rettype="fasta", retmode="text"
        )
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print(commonName + " saved")

"""
Args:
    email (string): the identifying email used for accessing Entrez
    proteinCSV (csv): csv file in the format of two columns: protein name and uid, each row being a protein of interest

Returns:
    N/A, but downloads the fasta files of the identified proteins into the proteins folder
"""
def batchDownload(email, proteinCSV):
    with open(proteinCSV) as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            commonName = row[0]
            proteinUID = row[1]
            downloadSeq(email, commonName, proteinUID)

email = input("User email: ")
    
batchDownload(email, 'input/proteins.csv')
batchDownload(email, 'input/viralproteins.csv')