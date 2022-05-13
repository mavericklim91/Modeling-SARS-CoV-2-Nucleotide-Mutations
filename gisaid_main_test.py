# -*- coding: utf-8 -*-
"""
Main program for mutation model analysis
1. Ingest sequences.fasta file containing nucleotide sequences
2. Split individual sequences into a sequence list
3. Converts sequences into numpy array (150x200)
4. Plot the sequence arrays in pcolor
5. Plot the heat map of mutations at each index
"""
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import scipy.stats as ss
from matplotlib.patches import Patch
from gisaid_test import split_seq, convert_seq, mut_finder, mut_dist

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
print(converted_arr[0])

# Create heat map of mutations against reference strain
ref_strain = converted_arr[0]
mut_counter = np.zeros((150,200))
for i in range(len(converted_arr)-1):
    mut_counter += mut_finder(converted_arr[0],converted_arr[i+1])

# Plot consensus Nucleotide Map using Matplotlib.pcolor
fig = plt.pcolor(converted_arr[0])
plt.title('Aligned Nucleotide Map', fontweight="bold")
legend_elements = [Patch(facecolor='yellow', label='A'),
                   Patch(facecolor='green', label='G'),
                   Patch(facecolor='blue', label='T'),
                   Patch(facecolor='purple', label='C'),
                   Patch(facecolor='black', label='N')]
plt.legend(handles=legend_elements)
plt.savefig('Nucleotide Map.png', dpi=200)
plt.show()

# Plot Mutation Heatmap using Seaborn
sb.heatmap(mut_counter, cmap='Reds', cbar_kws={'label': 'Mutation Count'})

plt.title("Mutation Heatmap", fontweight="bold")
plt.ylabel('Nucleotide index')

plt.savefig('Mutation Heatmap.png', dpi=200)
plt.show()

# Plot histogram for inter-occurrence distances (spatial model)
distances = mut_dist(mut_counter)
mean = np.mean(distances)

plt.hist(distances, bins=15)
plt.title('Histogram of Inter-occurrence Distances', fontweight='bold')
plt.ylabel('Frequency')
plt.xlabel('Distance')
plt.savefig('Distance_Histograms.png', dpi=200)
plt.show()

# Plot Distances on log-linear scale
n, bins, patches = plt.hist(distances, bins=15, log=True)
plt.show()

X = (bins[:-1] + bins[1:])/2
Y = n
plt.plot(X, Y, marker='o')
plt.yscale('log')
plt.title('Distance Histograms (log-linear scale)', fontweight='bold')
plt.ylabel('log Frequency')
plt.xlabel('Distance')
plt.axvline(x=500, color='b', linestyle='--')
plt.savefig('Distance_Histograms_logscale.png', dpi=200)
plt.show()

# Plot Short Distances (<500) on log-linear scale
dist = np.array(distances)
short = dist[dist < 500]

n, bins, patches = plt.hist(short, bins=15, log=True)
plt.show()

X = (bins[:-1] + bins[1:])/2
Y = n
plt.scatter(X, Y, marker='x')
plt.yscale('log')
plt.title('Short Distance (<500) Histograms (log-linear scale)', fontweight='bold')
plt.ylabel('log Frequency')
plt.xlabel('Distance')
plt.savefig('Short_Distance_Histograms_logscale.png', dpi=200)
plt.show()

# Plot Long Distances (>500) on log-linear scale
dist = np.array(distances)
long = dist[dist > 500]

n, bins, patches = plt.hist(long, bins=15, log=True)
plt.show()

X = (bins[:-1] + bins[1:])/2
Y = n

plt.scatter(X, Y, marker='x')
plt.yscale('log')
plt.title('Long Distance (>500) Histograms (log-linear scale)', fontweight='bold')
plt.ylabel('log Frequency')
plt.xlabel('Distance')
plt.savefig('Long_Distance_Histograms_logscale.png', dpi=200)
plt.show()
