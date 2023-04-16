import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

# Get arguments from command line
parser = argparse.ArgumentParser()
parser.add_argument('--tumor', help='Tumor file')
parser.add_argument('--normal', help='Normal file')
parser.add_argument('--output', help='Output file')
args = parser.parse_args()

# Read txt file with pandas
tumor = pd.read_csv('output/tu.txt', sep='	', header=None)
normal = pd.read_csv('output/wt.txt', sep='	', header=None)

# Rename columns
tumor.columns = ['chrom', 'pos', 'reads']
normal.columns = ['chrom', 'pos', 'reads']

# If position is missing in normal, add it with 0 reads
normal = normal.set_index('pos').reindex(tumor['pos']).reset_index()
tumor = tumor.set_index('pos').reindex(normal['pos']).reset_index()

# Retrieve reads as numpy array
tumor_reads = tumor['reads'].to_numpy()
normal_reads = normal['reads'].to_numpy()

# Retrieve position as numpy array
tumor_pos = tumor['pos'].to_numpy()
normal_pos = normal['pos'].to_numpy()

# Calculate log2
log2 = np.log2(tumor_reads/normal_reads)
log2 = np.nan_to_num(log2, nan=0, posinf=0, neginf=0)

# Plot
plt.figure(figsize=(20,10))
plt.scatter(tumor_pos, log2, color='k', s=1)
plt.xlabel('Position (ChrX)')
plt.ylabel('Log2(Tumor/Normal)')
plt.title('Tumor vs Normal')
plt.savefig('output/plot.png')
