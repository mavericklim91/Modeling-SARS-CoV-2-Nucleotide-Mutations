# -*- coding: utf-8 -*-
"""
Main program for mutation model analysis
1. Ingest sequences.fasta file containing nucleotide sequences
2. Split individual sequences into a sequence list
3. Converts sequences into numpy array (150x200)
4. Plot the transition matrices as a stacked bar chart
5. Plot the simulated time-series
"""
import matplotlib.pyplot as plt
import numpy as np
from gisaid_test import split_seq, convert_seq

# File names repository
file_name = 'gisaid_test.fasta'
aligned_file_name = 'gisaid_aligned.txt'
ncbi_alignment = 'ncbi_alignment.fasta'
ncbi_ref = 'ncbi_ref_aligned_20.txt'
aligned_193 = 'OM339503.1 and 192 other sequences.aln'

# Select the file to split into list of string sequences
sequences_list = split_seq(aligned_193)

# Create list of numpy arrays; omits consensus (1st sequence in file)
converted_arr = []
for i in range(1,len(sequences_list)):
    converted_arr.append(convert_seq(sequences_list[i]))

# Count frequency of nucleotide substitutions
mut_mat = np.zeros((4,4))

for seq in range(len(sequences_list)-1):            # Iterate over all sequences
    for i in range(len(sequences_list[seq])):       # Iterate over all indices
        a = sequences_list[seq][i]                  # Compare nucleotides
        b = sequences_list[seq+1][i]
        if a == 'T' and b == 'C':
            mut_mat[0][1] += 1
        if a == 'T' and b == 'A':
            mut_mat[0][2] += 1
        if a == 'T' and b == 'G':
            mut_mat[0][3] += 1
        if a == 'C' and b == 'T':
            mut_mat[1][0] += 1
        if a == 'C' and b == 'A':
            mut_mat[1][2] += 1
        if a == 'C' and b == 'G':
            mut_mat[1][3] += 1
        if a == 'A' and b == 'T':
            mut_mat[2][0] += 1
        if a == 'A' and b == 'C':
            mut_mat[2][1] += 1
        if a == 'A' and b == 'G':
            mut_mat[2][3] += 1
        if a == 'G' and b == 'T':
            mut_mat[3][0] += 1
        if a == 'G' and b == 'C':
            mut_mat[3][1] += 1
        if a == 'G' and b == 'A':
            mut_mat[3][2] += 1
 
# Compute Markov transition matrix
print('\nMutation Counter:\n',mut_mat)
total = np.sum(mut_mat)
trans_mat = mut_mat / total
print('\nTransition Rate Matrix:\n',trans_mat)

# Plot stacked bar chart of transition probabilities
for i in range(len(trans_mat)):
    plt.bar(['T','C','A','G'], trans_mat[i], 1)

plt.text(-0.2,0.275,'0.275')
plt.text(0.8,0.275,'0.272')
plt.text(1.8,0.11,'0.103')
plt.text(2.8,0.11,'0.101')

plt.ylim(0, 0.3)
plt.ylabel('Transition Probability')
plt.legend(['T','C','A','G'])

plt.savefig('Mutation_Rate_Stacked')
plt.show()





