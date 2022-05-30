#code to run Python backend 01/16/2019
#####################import#################################
#basic
from __future__ import print_function
import numpy as np
import cPickle as pickle
#keras related
import keras
from keras.models import model_from_json
print('Ddeveloped with keras version 2.0.3, and current keras version: ',keras.__version__)
#math
from scipy.stats import percentileofscore
import Levenshtein

####################Get relative path####################################
path_additional = "/".join(__file__.split("/")[:-1])
#path_additional = '/Users/bchen45/Dropbox/Ash/MARIA/Web' #use this line when running in a jupyter notebook
path_model = path_additional+'/models/'
path_dict = path_additional+'/supporting_file/'

#####################Loading dictionary########################
print('Loading data from '+path_dict)

#loading sparse encoding or amino acids
dict_aa = pickle.load(open(path_dict+'aa_21_sparse_encoding.dict','r'))
dict_aa['-'] = [0]*19 + [-1,0]

#loading gene2seq dictionary (Gene ID to seqeunces)
dict_gene2seq = pickle.load(open(path_dict+'gene2seq.dict','r'))

#default cleavage scores set to 0.39
cleavage_default = 'AAAAAAAKP--' #0.39 for probability

#DR encoding
dict_name = 'HLA_DR_encoding.dict'
dict_dr = pickle.load(open(path_dict+dict_name,'r'))

################## Parameters ######################################
len_mhc = 19
print('Each DR allele presented by a '+str(len_mhc)+'-AA long pseudosequence')
MAXLEN = 26 #26 is set to be consistent with training
#MAXLEN_binding = 26+1+19 #length used for training the binding model
print('The maximum length allowed for MARIA MHC-DR is 25')
batch0 = 1024
char_len = 21
len_char = 21
default_linker = '------'


#####################model and weight #################################
#function to import keras NN model
def import_model(path_model, model_name0,weight_name0):
    model_name0 = path_model+ model_name0
    weight_name0 = path_model + weight_name0
    model0 = model_from_json(open(model_name0).read())
    model0.load_weights(weight_name0)
    return model0

################ NN related Encoding functions ################ 
#the following three functions one-hot encode a variable peptide sequence
def encoding_line(str0, max_len, char_len, dict_aa, classn = 2):
    coded0 = np.zeros((max_len,char_len))
    for i,char0 in enumerate(str0):
        coded0[i,:] = dict_aa[char0] 
    return coded0
def encoding(matrix0, input0, len0,dict_aa, char_len):
    for i, sentence in enumerate(input0):
        matrix0[i] = encoding_line(sentence, len0, char_len,dict_aa)
    return matrix0
def encoding_data(list0,MAXLEN=MAXLEN,dict_aa=dict_aa,char_len=char_len ):
    #encoding   
    X_0_m = np.zeros((len(list0), MAXLEN, char_len))
    X_encoded = encoding(X_0_m,list0, MAXLEN, dict_aa, char_len)
    return X_encoded

#encoding a fixed length of peptides
def encoding_fixed(list0,MAXLEN=MAXLEN,dict_aa=dict_aa,char_len=char_len,add_placeh=0):
    X_0_m = np.zeros((len(list0), MAXLEN, char_len))
    X_encoded = encoding(X_0_m,list0,MAXLEN, dict_aa,char_len)
    X_encoded = X_encoded.reshape((X_encoded.shape[0],MAXLEN*char_len))
    if add_placeh > 0:
        add_block = np.ones((X_encoded.shape[0],add_placeh))
        #print(X_encoded.shape)
        X_encoded = np.concatenate((X_encoded,add_block),axis=1)
        #print(X_encoded.shape)
    return np.array(X_encoded)


#combine HLA allele and peptide sequence
#merge model
def encode_apair(data_pair,dict_aa=dict_aa,len_mhc=len_mhc,MAXLEN=MAXLEN,char_len=char_len):
    #print(data_pair)
    x_mhc = encoding_fixed(data_pair[0],len_mhc,dict_aa,char_len,add_placeh=0)
    x_seq = encoding_data(data_pair[1],MAXLEN,dict_aa,char_len)
    return [x_mhc, x_seq]

#old one sequence model
# def encode_apair(data_pair,dict_aa=dict_aa,len_mhc=fix_len,MAXLEN=MAXLEN_binding,char_len=char_len,spacer='-'):
#     data_seq = [mhc0+spacer+seq0 for mhc0,seq0 in zip(data_pair[0],data_pair[1])]
#     x_seq = encoding_data(data_seq,MAXLEN,dict_aa,char_len)
#     return x_seq

#gene expresion normalization
def convert_gn2exp(list_gn,dict_rna,scaling=10,def_val=5):
    list_out = []
    for gn0 in list_gn:
        if gn0 in dict_rna:
            list_out.append(dict_rna[gn0]+0.001)
        else:
            list_out.append(def_val)
    list_out_log = np.log10(list_out)/float(def_val)
    return list_out_log, list_out

#cleavage related functions
#function to match a frament to sequence
def best_match_in_gene(frag0,seq0):
    i_out = -1
    r_best = 0
    for i in range(6,len(seq0)-len(frag0)-5):
        frag_test = seq0[i:i+len(frag0)]
        r0 = Levenshtein.ratio(frag0,frag_test)
        if r0 > r_best:
            i_out = i
            frag_best = frag_test
            r_best = r0
    #for debugging
    #print('original='+frag0+',best_match='+frag_best+' r='+str(r_best))
    return i_out

#cleavage related functions
def get_cleave_ingene(frag0,seq0,in_gene_priority=False,flanking=6,default_linker=default_linker):
    out = []
    seq0 = default_linker+seq0+default_linker
    index0 = seq0.find(frag0)
    #if this sequence can't be exactly matched, sliding over the gene sequences
    if (index0 == -1) & in_gene_priority:
        #find the best match
        index0 = best_match_in_gene(frag0,seq0)
    #if this sequence can be exactly matched
    if index0 > -1:
        part1 = seq0[index0-flanking:index0]
        part2 = seq0[index0+len(frag0):index0+len(frag0)+flanking]
        #print(part1)
        #print(part2)
        out = part1 + part2

    return out

def get_cleave_scores(model_cleavage,pep_list,gn_list,
                      flanking=6, cleavage_default = cleavage_default,char_len=21,
                      dict_aa=dict_aa,in_gene_priority=True,default_linker = default_linker):
    cle_list = []
    list_unfound = []
    out = []
    #get cleavage sites
    for pep0, gn0 in zip(pep_list,gn_list):
        if gn0 in dict_gene2seq:
            seq0 = dict_gene2seq[gn0]
            out = get_cleave_ingene(pep0,seq0,flanking=flanking,in_gene_priority=in_gene_priority,
                                    default_linker=default_linker)
        #if can't find the gene sequence, use default vlaues
        if out == []:
            out = cleavage_default
        cle_list.append(out)
    #run model on cleavage sites
    x_predict = encoding_data(cle_list,MAXLEN=flanking*2,char_len=char_len,
                              dict_aa=dict_aa).reshape(len(cle_list),flanking*2*char_len)
    scores = model_cleavage.predict(x_predict)[:,0]
    return scores

#this version is for generating training too
#in this version, mhc_list is previously translated and has the same size of pep_list
#output encoding sequence matrix too
def run_MARIA_binding_mhclist(model0,mhc_list,pep_list,dict_aa=dict_aa,dict_dr=dict_dr,dict_distr='',
              len_mhc=len_mhc,MAXLEN=MAXLEN,char_len=char_len,batch0=batch0,translate_dr=True,
                      output_ranking=False):
    #convert HLA-DR alleles into pseudosequences
    if translate_dr:
        mhc_list = [dict_dr[allele0] for allele0 in mhc_list]
    #combine mhc list
    data0 = [mhc_list,pep_list]
    #encode for binidng models
    x = encode_apair(data0,dict_aa,len_mhc=len_mhc,char_len=char_len,MAXLEN=MAXLEN)
    #run prediction
    list_predicted = model0.predict(x,verbose=0,batch_size=batch0)[:,0]
    #print(list_predicted)
    if output_ranking:
        list_predicted = get_binding_rank(list_predicted,pep_list,mhc_list,dict_distr)
    return list_predicted


#adding binding, cleavage, and gene expression to a fixed vector
def add_binding_to_dataset(train_list,model_binding,model_cleavage,dict_rna,ranking=True):
    mhc0_list = []
    mhc1_list = []
    seq_list = train_list[1]
    gene_list = []
    
    for fixed0 in train_list[0]:
        allele0 = fixed0[0]
        allele1 = fixed0[1]
        gene_list.append(fixed0[2])
        mhc0_list.append(allele0)
        mhc1_list.append(allele1)
    #print(seq_list)
    score0 = run_MARIA_binding_mhclist(model_binding,mhc0_list,seq_list,dict_distr='',
                 output_ranking=ranking)
    score1 = run_MARIA_binding_mhclist(model_binding,mhc1_list,seq_list,dict_distr='',
                 output_ranking=ranking)
    score_exp,list_tpm = convert_gn2exp(gene_list,dict_rna)
    cleavage_score = get_cleave_scores(model_cleavage,seq_list,gene_list)
    
    #print(len(score_exp))
    score_list = np.array([score0,score1,score_exp,cleavage_score]).T

    #encode peptide sequences
    seq_encoded_list = encoding_data(seq_list,MAXLEN,dict_aa,char_len)

    #get input for MARIA merge model
    train_new = [score_list,seq_encoded_list]

    return train_new, list_tpm

#create a list of HLA in proper formats from HLA-dr input (input a single HLA-DR or two HLA-DR)
def make_mhc_list(mhc_list,pep_list):
    mhc_list = mhc_list.replace(' ','')
    mhc_list = mhc_list.split(',')
    if len(mhc_list) == 1:
        if mhc_list[0] in dict_dr:
            dr_seq = dict_dr[mhc_list[0]]
            out = [[dr_seq,dr_seq]]*len(pep_list)
        else:
            out = 'DR allele not recognized.'
    elif len(mhc_list) == 2:
        if mhc_list[0] in dict_dr and mhc_list[1] in dict_dr:
            dr_seq0 = dict_dr[mhc_list[0]]
            dr_seq1 = dict_dr[mhc_list[1]]
            out = [[dr_seq0,dr_seq1]]*len(pep_list)
        else:
            out = 'DR alleles not recognized.'
    else:
        out = "Please give four digits HLA-DR alleles (e.g. HLA-DRB1*07:01). For single HLA-DR allele, repeat the al"
    return out

#### final scores normalization by comparing to a negative distribution
def get_final_rank(list_predicted,pep_list,dict_distr):
    list_out = []
    for score0, pep0 in zip(list_predicted,pep_list):
        len0 = len(pep0)
        rank0 = percentileofscore(dict_distr[len0],score0)
        list_out.append(rank0)
    return list_out

###core function to predict HLA-DQ from sequences
def predict_dq(model_dq,list_sequence,dict_distr,dict_aa=dict_aa,
                 MAXLEN=25,len_char=len_char,batch0=1024,output_percentile=True ):
    #encoding data
    #MAXLEN is set to 25 to be consisent with the model training
    x = encoding_data(list0=list_sequence,dict_aa=dict_aa,MAXLEN=MAXLEN,char_len=len_char)
    #print(x)
    list_predicted = model_dq.predict(x,verbose=0,batch_size=batch0)[:,1]
    if output_percentile:
        list_predicted_precentile = get_final_rank(list_predicted,list_sequence,dict_distr)
    return list_predicted, list_predicted_precentile

####core function to predict HLA-DR from gene symbol, alleles, and peptide sequences
def predict_with4(pep_list,
                  mhc_list,
                  gn_list,
                  dict_rna,
                  dict_distr,
                  model_merge,
                  model_cleavage,
                  model_binding,
                  binding_ranking=False):
    #make the data into the right format
    data0 = [[],[]]
    #prepare data format
    for mhc0,gn0,pep0 in zip(mhc_list,gn_list,pep_list):
        data0[0].append([mhc0[0],mhc0[1],gn0])
        data0[1].append(pep0)
    #add binding, cleavage, and gene expression
    data0_encoded, list_tpm = add_binding_to_dataset(data0,model_binding,model_cleavage,dict_rna,
        ranking=binding_ranking)
    #run MARIA
    list_scores = model_merge.predict(data0_encoded)[:,1]
    #get percentiles
    list_percentiles = get_final_rank(list_scores,pep_list,dict_distr)
    return list_tpm, list_scores, list_percentiles

#sliding a long sequence and get 15mers
def slide_gene(str0,len0=15):
    list_out = []
    for i in range(0,len(str0)-len0):
        list_out.append(str0[i:i+len0])
    return list_out

#function to run HLA-DR models
#taking three inputs
#1) input 
def run_dr(data_matrix,
           dict_expression,
            tpm_reference):
    ############getting models relavent to HLA-DR
    #binding model
    #path_save = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/model_weight/'
    model1 = 'rnn_merge_all_netmhc_data_n64_sparse_d0.45_dronly.json'
    weight1 = 'rnn_merge_all_netmhc_data_n64_sparse_d0.45_dronly.h5'
    model_binding = import_model(path_model,model1,weight1)
    model_binding.compile(loss='binary_crossentropy', optimizer='adam')
    #cleavage model
    model1 = 'nn_classifier_cleavage_n32_updown6_d0.5.json'
    weight1 = 'nn_classifier_cleavage_n32_updown6_d0.5.weight'
    model_cleavage = import_model(path_model,model1,weight1)
    model_cleavage.compile(loss='binary_crossentropy', optimizer='adam')
    #merge model MARIA
    model1 = 'rnn_withiedb_4feature_merge_d0.4_934auc_norank.json'
    weight1 = 'rnn_withiedb_4feature_merge_d0.4_934auc_norank.h5'

    model_merge = import_model(path_model,model1,weight1)
    model_merge.compile(loss='categorical_crossentropy', optimizer='adam')
    ##loading distributions of scores from random human peptides
    #dict_neg_scores = pickle.load(open(path_dict+'maria_full_neg_score_distr.dict','r'))
    dict_neg_scores = pickle.load(open(path_dict+'neg_score_distr.dict','r'))

    #########getting RNA-Seq dictionary
    #load reference RNA-Seq dictionary if TPM is not specified by the user
    if dict_expression == '':
        #loading gene expression dictionaries
        if tpm_reference == 'MCL':
            dict_rna = pickle.load(open(path_dict+'mcl_patient_tpm_blood_go.dict','r'))
        elif tpm_reference == 'K562':
            dict_rna = pickle.load(open(path_dict+'k562_tpm.dict','r'))
        else:
            dict_rna = pickle.load(open(path_dict+'tcga_cancer_median_tpm.dict','r'))
            dict_rna = dict_rna[tpm_reference]
    else:
        dict_rna = dict_expression

    ########runing prediction
    mhc_list = [[x,y] for x,y in zip(data_matrix[:,0],data_matrix[:,1])]
    list_tpm, list_scores, list_percentiles = predict_with4(pep_list=data_matrix[:,3],
                                                      mhc_list=mhc_list,
                                                      gn_list=data_matrix[:,2],
                                                      dict_rna=dict_rna,
                                                      dict_distr=dict_neg_scores,
                                                      model_merge=model_merge,
                                                      model_cleavage=model_cleavage,
                                                      model_binding=model_binding,
                                                      binding_ranking=False)
    #return values
    return list_tpm, list_scores, list_percentiles

#function to run DQ2.2 model or DQ2.5 model
def run_dq(list_sequence,allele_dq):
    if 'DQ2.2' in allele_dq:
        #load DQ2.2 model and distribution
        model1 = 'lstm_dq2.2_combined_n64_d0.5_rd0.3_auc0.89.json'
        weight1 = 'lstm_dq2.2_combined_n64_d0.5_rd0.3_auc0.89.h5'
        model_dq = import_model(path_model,model1,weight1)
        model_dq.compile(loss='categorical_crossentropy', optimizer='adam')
        ##loading distributions of scores from random human peptides
        dict_neg_scores = pickle.load(open(path_dict+'DQ2.2_dict_distr.dict','r'))
    else:
        #load DQ2.5 model (check 2.2 and 2.5 here)
        model1 = 'lstm_dq2.2_combined_n64_d0.5_rd0.3_auc0.89.json'
        weight1 = 'lstm_dq2.2_combined_n64_d0.5_rd0.3_auc0.89.h5'
        model_dq = import_model(path_model,model1,weight1)
        model_dq.compile(loss='categorical_crossentropy', optimizer='adam')
        ##loading distributions of scores from random human peptides
        dict_neg_scores = pickle.load(open(path_dict+'DQ2.2_dict_distr.dict','r'))

    list_scores, list_percentiles = predict_dq(model_dq=model_dq,
                                                list_sequence=list_sequence,
                                                dict_distr=dict_neg_scores)
    #return values
    return list_scores, list_percentiles 
