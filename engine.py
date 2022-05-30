import os, csv
from Bio import SeqIO
from Bio import Entrez
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.PDB.Polypeptide import three_to_one
import pandas as pd

class Epitope:
    def __init__(self, seq, protein, allele, core, mix, maria, viral, original):
        self.seq = seq
        self.protein = protein
        self.allele = allele
        self.core = core
        self.blosum = 0
        self.max_blosum = 0
        self.mix = mix
        self.maria = maria
        self.viral = viral
        self.original = original
        self.mutated = (original != None) and (original.seq != seq)

    def __str__(self):
        return self.seq + " Protein: " + self.protein + " Mix: " + str(self.mix) + " Maria: " + str(self.maria)

    def __repr__(self):
        return self.seq + " Protein: " + self.protein + " Core: " + self.core + " Mix: " + str(self.mix) + " Maria: " + str(self.maria)

"""
Arguments:
    retDict (dict): A dictionary storing all peptides and their pertinent information to be assessed by MixMHC2pred
    peptides (list): A list of all peptides to be tested
    allele (str): A string representing the allele of interest to run MixMHC2pred on
Returns:
    A modified copy of the original dictionary with MixMHC2pred scores for each of the specified peptides
"""
def getMix(retDict, peptides, allele):
    with open('MixMHC2pred-master/temp.txt', 'w') as temp:
        for peptide in peptides:
            temp.write(peptide[0]+"\n")
            if peptide[0] not in retDict[allele].keys():
                retDict[allele][peptide[0]] = {}
            retDict[allele][peptide[0]]["gene"] = peptide[1]
    os.system("cd MixMHC2pred-master; ./MixMHC2pred -i temp.txt -o temp_out.txt -a " + allele)
    with open('MixMHC2pred-master/temp_out.txt', 'r') as out:
        lines = out.readlines()[20:]
        for line in lines:
            line = line.split("\t")
            retVal = float(line[6])
            core_pos = int(line[7])
            peptide = line[0]
            core = peptide[core_pos-1:core_pos+8]
            retDict[allele][peptide]["mix"] = retVal
            retDict[allele][peptide]["core"] = core
    os.remove('MixMHC2pred-master/temp.txt')
    os.remove('MixMHC2pred-master/temp_out.txt')
    return retDict

"""
Arguments:
    retDict (dict): A dictionary storing all peptides and their pertinent information to be assessed by MARIA
    peptides (list): A list of all peptides to be tested
    allele (str): A string representing the allele of interest to run MARIA on
Returns:
    A modified copy of the original dictionary with MARIA scores for each of the specified peptides
"""
def getMaria(retDict, peptides, allele):
    with open('MARIA_local_package/temp.txt', 'w') as temp:
        temp.write('\t'.join(["Allele1", "Allele2 (Same as Allele1 if analyzing a single allele)", "Gene Symbol", "Peptide Sequence", "TPM (Optional)"])+'\n')
        for peptide in peptides:
            if peptide[0] not in retDict[mix_allele(allele)].keys():
                retDict[mix_allele(allele)][peptide[0]] = {}
            retDict[mix_allele(allele)][peptide[0]]["gene"] = peptide[1]
            temp.write(allele + "\t" + allele + "\t" + peptide[1] + "\t" + peptide[0] + "\n")
    os.system("cd MARIA_local_package; python maria.py temp.txt")
    with open('MARIA_local_package/temp.txt.output.txt', 'r') as out:
        for line in out.readlines()[1:]:
            line = line.split("\t")
            peptide = line[3]
            retVal = float(line[7])
            commonName = line[2]
            retDict[mix_allele(allele)][peptide]["maria"] = retVal
    os.remove('MARIA_local_package/temp.txt')
    os.remove('MARIA_local_package/temp.txt.output.txt')
    return retDict

"""
Arguments:
    maria_allele (str): A string representation of the HLA allele of interest in the format of MARIA's list of supported alleles
Returns:
    A string representaion of the HLA allele of interest in the format of MixMHC2pred's list of supported alleles
    Ex: HLA-DRB1*07:01        --> DRB1_07_01
        DQA1*02:01-DQB1*02:02 --> DQA1_02_01_DQB1_02_02
"""
def mix_allele(maria_allele):
    if "HLA" in maria_allele:
        maria_allele = maria_allele[4:]
    return maria_allele.replace("-", "__").replace("*", "_").replace(":", "_")

"""
Arguments:
    allele (str): A string formatted in MixMHC2pred's supported alleles format representing the HLA allele of interest
Returns:
    Boolean representing whether or not MixMHC2pred has support for the provided HLA allele
"""
def mix_supported(allele):
    with open("MixMHC2pred-master/alleles.csv", "r") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=",")
        for row in csv_reader:
            if allele in row[0]:
                return True
        return False

"""
Arguments:
    allele (str): A string formatted in MARIA's supported alleles format representing the HLA allele of interest
Returns:
    Boolean representing whether or not MARIA has support for the provided HLA allele
"""
def maria_supported(allele):
    with open("MARIA_local_package/alleles.csv", "r") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=",")
        for row in csv_reader:
            if allele in row[0]:
                return True
        return False

"""
Arguments:
    seq (str): A string representing the amino acid sequence of a given protein, each amino acid represented by its 1 letter code
    size (int): An int representing the desired reading frame size for slicing
Returns:
    A list of all peptides generated using a reading frame of given size on the provided sequence
"""
def slice_seq(seq, size):
    lst = []
    for i in range(len(seq)-size+1):
        newStr = str(seq[i:i+size])
        if "X" not in newStr:
            lst.append(newStr)
    return lst

"""
Arguments:
    k_prot (str): The amino acid sequence of the self-peptide
    v_prot (str): The amino acid sequence of the viral peptide
Returns:
    A list [a, b], where a = blosum score and b = max blosum score:
    a (int): the sum of computing the BLOSUM score between the self and viral peptide across their 9-amino acid core binding region
    b (int): the maximum possible BLOSUM score assuming the self- and viral peptides have the exact same structure across their 9-amino acid core binding region
"""
def blosum(k_prot, v_prot):
    sum = 0
    max_score = 0
    matrix = matlist.blosum62
    if len(k_prot) != len(v_prot):
        raise Exception("Proteins aren't the same length")
    for i in range(len(k_prot)):
        aa_1 = k_prot[i]
        aa_2 = v_prot[i]
        aa_tup = (aa_1, aa_2)
        if aa_tup not in matrix.keys():
            aa_tup = (aa_2, aa_1)
        score = matrix[aa_tup]
        sum += score
        max_score += matrix[(aa_1, aa_1)]
    return [sum, max_score]

"""
Arguments:
    genes (csv): A csv file containing all the proteins from which to populate peptide candidates (see sample files)
    alleles (csv): A csv file containing all the HLA alleles to be tested
    viral_proteins (csv): A csv file containing all the viral proteins from which to populate viral peptides
    length (int): The desired slicing length, uniform across all peptides
Returns:
    A list of epitope objects generated from running all self proteins, HLA alleles, and viral proteins of interest through MixMHC2pred and MARIA
"""
def generate_epitopes(genes, alleles, viral_proteins, length):
    epitopes = []
    retDict = {}
    viralDict = {}
    with open(genes, 'r') as csvfile, open(alleles, 'r') as csvfile3, open(viral_proteins, 'r') as csvfile4:
        genes_csv = csv.reader(csvfile, delimiter=',')
        alleles_csv = csv.reader(csvfile3, delimiter=",")
        viral_csv = csv.reader(csvfile4, delimiter=",")
        next(genes_csv)
        next(viral_csv)
        peptide_list = []
        viral_list = []
        for row in genes_csv:
            commonName = row[0]
            allFrames = []
            fasta_file = "input/proteins/"+commonName+".fasta"
            record = SeqIO.read(fasta_file, "fasta")
            seq = record.seq
            peptides = slice_seq(seq, length)
            peptides2 = [[k, commonName, None] for k in peptides]
            peptide_list += peptides2
            seq_file = "input/proteins/" + commonName + ".fasta"
            record = SeqIO.read(seq_file, "fasta")
        for row in viral_csv:
            commonName = row[0]
            allFrames = []
            fasta_file = "input/proteins/"+commonName+".fasta"
            record = SeqIO.read(fasta_file, "fasta")
            seq = record.seq
            peptides = slice_seq(seq, length)
            peptides2 = [[k, commonName, None] for k in peptides]
            viral_list += peptides2
        for row in alleles_csv:
            print("Ok starting allele: " + row[0])
            allele = row[0]
            mixed_allele = mix_allele(allele)
            retDict[mixed_allele] = {}
            viralDict[mixed_allele] = {}
            if mix_supported(mix_allele(allele)):
                getMix(retDict, peptide_list, mix_allele(allele))
                getMix(viralDict, viral_list, mix_allele(allele))
            if maria_supported(allele):
                getMaria(retDict, peptide_list, allele)
                getMaria(viralDict, viral_list, allele)
            for p in peptide_list:
                peptide = p[0]
                if peptide[0] == "[":
                    continue
                mixVal = -1
                mariaVal = -1
                if "mix" in retDict[mixed_allele][peptide]:
                    mixVal = retDict[mixed_allele][peptide]["mix"]
                if "maria" in retDict[mixed_allele][peptide]:
                    mariaVal = retDict[mixed_allele][peptide]["maria"]
                if p[2] != None:
                    original = [e for e in epitopes if e.allele == allele and e.seq == p[2]][0]
                else:
                    original = None
                epitopes.append(Epitope(peptide, retDict[mixed_allele][peptide]["gene"], allele, retDict[mixed_allele][peptide]["core"], mixVal, mariaVal, False, original))
            for peptide in viralDict[mixed_allele].keys():
                mixVal = -1
                mariaVal = -1
                if "mix" in viralDict[mixed_allele][peptide]:
                    mixVal = viralDict[mixed_allele][peptide]["mix"]
                if "maria" in viralDict[mixed_allele][peptide]:
                    mariaVal = viralDict[mixed_allele][peptide]["maria"]
                original = None
                epitopes.append(Epitope(peptide, viralDict[mixed_allele][peptide]["gene"], allele, viralDict[mixed_allele][peptide]["core"], mixVal, mariaVal, True, original))
    return epitopes

"""
Arguments:
    epitopes (list): A list of Epitope objects containing both self and viral peptides
Returns:
    The same epitopes, but with each self-peptide matched to its closest viral peptide by comparing BLOSUM scores and stored within the Epitope object
"""
def viral_mimicry(epitopes):
    viral_prots = []
    kidney_prots = []
    for e in epitopes:
        if e.viral == True:
            viral_prots.append(e)
        else:
            kidney_prots.append(e)
    retList = []
    for s in kidney_prots:
        best_prot_score = -100
        closest_viral = ""
        prot_max_score = 0
        for v in viral_prots:
            score_pair = blosum(s.core, v.core)
            score = score_pair[0]
            max_score = score_pair[1]
            if score > best_prot_score:
                best_prot_score = score
                prot_max_score = max_score
                closest_viral = v
        s.closest_viral = closest_viral
        s.blosum = best_prot_score
        s.max_blosum = prot_max_score
    return kidney_prots + viral_prots


"""
Arguments:
    alleles (csv): A csv file containing all the HLA alleles to be tested
    proteins (csv): A csv file containing all the proteins from which to populate peptide candidates (see sample files)
    viral_proteins (csv): A csv file containing all the viral proteins from which to populate viral peptides
    mix_threshold (float): The cut-off for assessing binding candidacy from MixMHC2pred, default at 2.00%
    maria_threshold (float): The cut-off for assessing binding candidacy from MARIA, default at 95.00%
Returns:
    The same epitopes, but with each self-peptide matched to its closest viral peptide by comparing BLOSUM scores and stored within the Epitope object
"""
def run(alleles, proteins, viral_proteins, mix_threshold, maria_threshold):
    epitopeList = generate_epitopes(proteins, alleles, viral_proteins, 15)
    with open("input/alleles.csv", 'r') as all_file, open("output/viral_final.csv", "w") as output:
        csv_reader = csv.reader(all_file)
        csv_writer = csv.writer(output)
        csv_writer.writerow(["allele", "peptide", "protein", "core", "virus", "cvp", "cvp_seq", "cvp_core", "blosum", "max_blosum", "mix", "cvp_mix", "maria", "cvp_maria"])
        for row in csv_reader:
            allele = row[0]
            list1 = []
            for epitope in epitopeList:
                if epitope.allele == allele and epitope.mix <= mix_threshold and ((not maria_supported(allele)) or (maria_supported(allele) and epitope.maria >= maria_threshold)):
                        list1.append(epitope)
            lst1 = viral_mimicry(list1)
            for l in lst1:
                if l.viral == False and l.mutated == False:
                    csv_writer.writerow([allele, l.seq, l.protein, l.core, l.closest_viral.protein.split("-")[0], l.closest_viral.protein, l.closest_viral.seq, l.closest_viral.core, l.blosum, l.max_blosum, l.mix, l.closest_viral.mix, l.maria, l.closest_viral.maria])

alleles = "input/alleles.csv"
proteins = "input/old_proteins.csv"
viral_prots = "input/viralproteins.csv"
run(alleles, proteins, viral_prots, 2.00, 95.00)