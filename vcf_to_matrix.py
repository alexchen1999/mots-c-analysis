# From NCBI: https://www.ncbi.nlm.nih.gov/gene?cmd=retrieve&dopt=default&list_uids=4549&rn=1
# MT-RNR1 is encoded from position 648-1601 (n=954)

from pysam import VariantFile
from sklearn import decomposition
import numpy as np
import pandas as pd

vcf_filename = "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf"
panel_filename = "integrated_call_samples_v3.20130502.ALL.panel"

genotypes = []
samples = []

with VariantFile(vcf_filename) as vcf_reader:
    counter = 0
    for record in vcf_reader:
        # print('record ', record)
        # print('record keys: ', record)
        # print('record info ', record.info)
        # print('record info keys ', record.info.keys())
        # print('record.pos ', record.pos)
        if (record.pos == 1382):
            print(record)
            counter += 1 
            alleles = [record.samples[x].allele_indices for x in record.samples]
            samples = [sample for sample in record.samples]
            genotypes.append(alleles)
        # print(counter)
        # if counter >= 10000:
        #    break


with open(panel_filename) as panel_file:
    labels = {} # {sample id : population code}
    for line in panel_file:
        line = line.strip().split('\t')
        labels[line[0]] = line[1]

print(samples)
print(genotypes)
print(len(genotypes))
print(len(genotypes[0]))

# print(labels)   

# print(genotypes)
print(len(genotypes))   

genotypes = np.array(genotypes)
print(genotypes.shape)

matrix = np.count_nonzero(genotypes, axis=2)
print(matrix.shape)

matrix = matrix.T
print(matrix.shape)
print(matrix)
'''
pca = decomposition.PCA(n_components=2)
pca.fit(matrix)
print(pca.singular_values_)
to_plot = pca.transform(matrix)
print(to_plot.shape)

df = pd.DataFrame(matrix, index=samples)
df['Population code'] = df.index.map(labels)    
print(df)
df.to_csv("matrix.csv")
'''