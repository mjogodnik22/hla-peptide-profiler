#importing 
import numpy as np
import pandas as pd
import sys
import argparse 

#####dealing with input files####
#function to take in user inputs
#return input file name (including path)
#plus five user defined attributes
def get_attributes():
    #create a parser object
    parser0 = argparse.ArgumentParser(description='MARIA command line tool taking a tab-separated input file')
    #parse mantatory input (input file)
    parser0.add_argument('input_file',default='')
    #parse optional attributes
    parser0.add_argument('-scan',default=0,type=int, help='1 = Analyze each sequence with 15mer sliding windows and use the best score, mandatory for seqeunces > 25 AAs')
    parser0.add_argument('-hladq2',default=0,type=int, help='1 = Analyze HLA-DQ2.2 or HLA-DQ2.5; default 0 = analyze HLA-DR presentation')
    parser0.add_argument('-custom_tpm',default=0,type=int, help='1 = Use user input TPMs for anlaysis; default 0 = use reference RNA-Seq values')
    parser0.add_argument('-tpm_reference',default='MCL',help='Gene expression of tissue/Cancer type to be used as reference gene expression profiles, default MCL. Ignored when -custom_tpm 1')
    parser0.add_argument('-cut_off',default=95.0,type=float,help='Percentile scorec cut_off to indicate whether a peptide will be presented')
    parser0.add_argument('-print_out',default=0,type=int,help='1 = Print the results; default 0 = save result it to a .output.txt file')
    #get attributes from command line
    parser0 = parser0.parse_args()
    #return 'test'
    return parser0.input_file, parser0.scan, parser0.hladq2, parser0.custom_tpm, parser0.tpm_reference, parser0.cut_off, parser0.print_out

#sliding a long sequence and get 15mer sliding windows
def slide_gene(str0,len0=15):
    list_out = []
    if len(str0) < len0:
        return [str0]
    for i in range(0,len(str0)-len0+1):
        list_out.append(str0[i:i+len0])
    return list_out

#input array of [allele0,allele1,gene0,seq0]
#in scan mode, make each sequences into 15mer sliding window (or shorter single window)
#in non_can mode, only convert ultra long sequences (>25AA) into sliding windows
def modify_data_matrix(data_x,scan):
    X = []
    track = {}
    for i,line0 in enumerate(data_x):
        allele0, allele1, gene0, seq0 = line0
        # track which original datapoint each window belongs to
        begin, end = len(X), len(X)
        if scan or len(seq0) > 25:
            for s in slide_gene(seq0):
                X.append([allele0, allele1, gene0, s])
                end += 1
                # data in the range of (begin, end) belongs to datapoint i
        else:
            X.append([allele0, allele1, gene0, seq0])
            end += 1
        track[i] = (begin, end)
    return np.array(X), track

#to recover the best scores based on track 
def select_best(list_tpm,list_maria_scores,list_maria_percentile,list_sequence,track):
    tpm_out,scores_out,percentile_out,sequence_out = [], [], [], []
    #print(outs[0:100])
    for k,vs in track.items():
        best_max = -1000
        #use the 15mer with the highest score
        for j in range(vs[0],vs[1]):
            #get the best segment
            if best_max < list_maria_scores[j]:
                best_max = list_maria_scores[j]
                tpm_selected = list_tpm[j]
                sequence_selected = list_sequence[j]
                percentile_selected = list_maria_percentile[j]
        #update the output
        tpm_out.append(tpm_selected)
        sequence_out.append(sequence_selected)
        percentile_out.append(percentile_selected)
        scores_out.append(best_max) 

    #return the new dataframe
    return tpm_out, scores_out, percentile_out, sequence_out

def print_pd(pd0,sep0='\t'):
    list_columns = list(pd0.columns)
    str_print = sep0.join(list_columns)
    print(str_print)

    for row0 in pd0.iterrows():
        row0 = list(row0[1])
        row0 = [str(x) for x in row0]
        str_print = sep0.join(row0)
        print(str_print)


#get following information based on column orders:
#column 0: allele1
#column 1: allele2 
#column 2: gene symbols
#column 3: sequences
#column 4: custom TPM values (can be empty)
def process_input(pd_input):
    column_index = list(pd_input)
    list_allele1 = list(pd_input[column_index[0]])
    list_allele2 = list(pd_input[column_index[1]])
    list_gene = list(pd_input[column_index[2]])
    list_seq = list(pd_input[column_index[3]])
    list_tpm_custom = list(pd_input[column_index[4]])
    return list_allele1, list_allele2, list_gene, list_seq, list_tpm_custom

#main function of MARIA command line tool
#inputs (see get_attributes for details):
#input_file: input file including path 
#scan: whether scan with 15mer or not (automatically if the peptide is longer than 26)
#hladq2: whether the model is scanning a HLA-DQ2 input or not
#custom_tpm: Whether the model uses user provided tpm values
#tpm_reference: dictionary name to be used as RNA-Seq reference
#cut_off: values to determine positive presenters, default value = 95 (peptides with top 95 percentles or above 
#indicated as positive)
#print_output: true if the result was printed out rather than saved to input_file+'.output.txt'
def main(input_file,scan,hladq2,custom_tpm,tpm_reference,cut_off,print_output):
    
    #get input dataframe
    pd_input = pd.read_csv(input_file,sep='\t')
    #get input information
    list_allele1, list_allele2, list_gene, list_seq, list_tpm_custom = process_input(pd_input)
    #check if it's a DQ allele
    if 'DQ' in list_allele1[0]:
        hladq2 = True
    #check all sequences are longer than 7 AAs
    for seq0 in list_seq:
        if len(seq0) < 8:
            print >>sys.stderr, 'Error: One of peptide sequence lengths is below 8 amino acids.'
            return 'Error: One of peptide sequence lengths is below 8 amino acids.'

    #create a gene-expression dictionary if user specify list_tpm_custom
    if custom_tpm and (not hladq2):
        #test if the value is provided
        list_tpm_custom_new = []
        for x in list_tpm_custom:
            if pd.isnull(x):
                #raise a red flag if the user did not provide TPM values
                print >>sys.stderr, 'Error: Custom RNA-Seq value otpion was chosen, but invalid TPM values were provided. Uncheck Custom option or input valid TPM values.'
                return 'Error: Custom RNA-Seq value option was chosen, but invalid TPM values were provided. Uncheck Custom option or input valid TPM values.'
            else:
                list_tpm_custom_new.append(float(x))

        #build a expression dictionary
        dict_expression = dict(zip(list_gene,list_tpm_custom_new))      
    else:
        dict_expression = ''
        
    #generate input data matrix
    data_matrix = np.array([list_allele1,list_allele2,list_gene,list_seq]).T

    #print(data_matrix)
    #scan the sequences with 15mer sliding windows if necessary
    data_matrix, track = modify_data_matrix(data_matrix,scan)
        
    #main prediction component
    if hladq2:
        #for DQ prediction, only peptide sequences are required and allele should be consistent throughout the input
        print('Running recurrent neural network for HLA-DQ ligand prediction')
        list_score, list_percentile = run_dq(data_matrix[:,3],list_allele1[0])
        list_tpm = [np.nan]*len(list_score)
    else:
        #DR prediction
        print('Running recurrent neural network for HLA-DR ligand prediction')
        list_tpm, list_score, list_percentile = run_dr(data_matrix,
                                                         dict_expression=dict_expression,
                                                         tpm_reference=tpm_reference)
    #get best values based on output and track
    #print(list_tpm,list_score,list_percentile)
    list_tpm,list_maria_score,list_maria_percentile,list_sequence = select_best(list_tpm,list_score,list_percentile,data_matrix[:,3],track)
    #print(list_tpm,list_maria_score,list_maria_percentile)
    if len(list_maria_score) != len(pd_input):
        #print(list_maria_score,len(list_maria_score),len(pd_input))
        print >>sys.stderr,'Error: Score numbers not matched with input query numbers, please check input file format'
        return 'Error: Score numbers not matched with input query numbers, please check input file format'
    
    #indicating positive presenters
    list_pos = []
    for x in list_maria_percentile:
        if x >= cut_off:
            list_pos.append(1)
        else:
            list_pos.append(0)
    
    #make output
    pd_input['TPM estimated'] = np.round(list_tpm,3)
    pd_input['MARIA raw scores'] = np.round(list_maria_score,4)
    pd_input['MARIA percentile scores'] = np.round(list_maria_percentile,3)
    pd_input['15mer core'] = list_sequence
    pd_input['Positive presenters'] = list_pos

    print('MARIA run was successful')
    #output results
    if print_output:
        #print(pd_input.to_string(index=False))
        print_pd(pd_input,sep0='\t')

    else:
        #save
        pd_input.to_csv(input_file+'.output.txt',sep='\t',index=False)
        #indicating finishing running
        print('The output was saved to '+input_file+'.output.txt')

#get input from user command lines
input_file,scan,hladq2,custom_tpm,tpm_reference,cut_off,print_output = get_attributes()
#only run the script if valid input file was provided
if input_file != '':
    from maria_core import run_dr, run_dq
    #run maria and output a text file or print the result
    main(input_file,scan,hladq2,custom_tpm,tpm_reference,cut_off,print_output)

