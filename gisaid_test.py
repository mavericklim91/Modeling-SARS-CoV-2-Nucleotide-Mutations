# -*- coding: utf-8 -*-
import numpy as np

'''
Below functions execute data pre-processing for the gene mutation model
'''

def split_seq(file_name):
    '''
    Function
    ----------
    Opens the 'sequences.fasta' file containing a series of Sars-CoV-2 nucleotide sequences
    Then splits the fasta database into a list of individual sequences
    
    Parameters
    ----------
    Input: file name <str> 
    File name (eg. 'sequences.fasta') downloaded from GISAID database
    
    Output: sequences <list>
    List containing sequences of nucleotides
    Each element is a string of characters 'A,G,T,C' (length ~30,000 chars)
    '''   
    content = open(file_name,'r')              # Opens the fasta database file
    lines = content.readlines()                # Reads the file line by line
    content.close()
    
    headers = []
    sequences = []
    nucleotides = ''
    
    for line in lines:                                           
        if line[0] == ">":                                  # If starts with ">" then it is a header line
            headers.append(line)                            # Append header line to list 'headers'
            if nucleotides != '':                           # If 'nucleotides' is not empty, it is 2nd header onwards
                full_seq = nucleotides.replace('\n','')     # Replace the line break characters to keep only 'A,G,T,C'
                sequences.append(full_seq)                  # Append the full sequence to 'sequences' list
                nucleotides = ''                            # Reset nucleotides to blank to collect next sequence
        else:
            nucleotides += line                             # Start collecting nucleotides into the current sequence
            if line == lines[-1]:                           # Once we have reached the last line
                full_seq = nucleotides.replace('\n','')     # Replace the line break characters to keep only 'A,G,T,C'
                sequences.append(full_seq)                  # Append the final full sequence to list 'sequences'
                nucleotides = ''                            # Reset nucleotides to blank to collect next sequence
    return sequences
    
def convert_seq(sequence):
    '''
    Function
    ----------
    Adds 'N' to make sequence length 30,000 and reshapes list to size (150 x 200)
    Then converts the gene sequence string into numpy array by replacing nucleotides:
    A = 50
    T = 200
    C = 250
    G = 100
    N = 0
    
    Parameters
    ----------
    Input: sequence <str>
    Sequence containing nucleotides A,G,T,C
    
    Output: converted <numpy array>
    Numpy array of size (150 x 200)
    
    '''
    # Add N to make sequence length = 30,000
    addN = 30000 - len(sequence)
    for i in range(addN):
        sequence += 'N'
    
    # Convert sequence from string (A, G, T, C, N) to numpy array (50, 200, 250, 100, 999)
    arr = ''
    for i in sequence:
        if i == 'A':
            arr += '50 '
        if i == 'T':
            arr += '200 '
        if i == 'C':
            arr += '250 '
        if i == 'G':
            arr += '100 '
        if i == 'N':
            arr += '999 '
        if i == '-':
            arr += '999 '
            
    converted = np.fromstring(arr,dtype=np.uint8,sep=' ')     
    converted = np.reshape(converted,(150,200))
    return converted

def mut_finder(seq1,seq2):
    '''
    Function
    ----------
    Compares 2 sequence arrays and finds the differences in each nucleotide index
    If there is a difference, append 1 in the mutation matrix of zeros
    0 = no mutation found
    1 = mutation found

    Parameters
    ----------
    Input: seq1, seq2 <numpy arrays>
    Sequence array of size (150 x 200)
    
    Output: mut <numpy arrays>
    Mutation array of size (150 x 200) --> contains 0s and 1s
    '''
    mut = seq1 - seq2                               # Subtract arrays to find mutations
    
    mut = (mut != 0) & (mut > -500) & (mut < 500)   # Obtain boolean array of mutations
    
    mut = mut.astype(int)                           # Convert boolean array to binary
    return mut
            
def mut_dist(mut_counter):
    '''
    Function
    ----------
    Calculates the inter-occurrence distance of mutations found

    Parameters
    ----------
    Input: mut_counter <numpy array>
    Sequence array of size (150 x 200)
    
    Output: distances <list>
    Nucleotide distances between mutations in the sequence
    '''
    mut_counter = mut_counter.flatten()
    mut_indices = np.where(mut_counter > 1)
    mut_indices = mut_indices[0]

    distances = []
    for i in range(len(mut_indices)-1):
        dist = mut_indices[i+1] - mut_indices[i]
        if dist > 1:                                # Excludes strings of indel mutations
            distances.append(dist)
    return distances
