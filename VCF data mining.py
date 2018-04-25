
# coding: utf-8

# In[106]:


import pandas as pd
import numpy as np

file = 'G:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\Genomics\M145 Evolved mutants\R2_to_ref_SNPs.csv'

def vcf_clean(file):
    '''
    Input = Exported VCF from Geneious,
    Select columns ['Minimum', 'Maximum', 'Amino Acid Change', 'Change', 'Codon Change', 'locus_tag', 'Polymorphism Type', 'Protein Effect',
                    'Variant Frequency'],
    Remove rows with 'Protein Effect' == NaN or None,
    Remove rows of 'locus_tag' == 'SCO3798' where the barcode is inserted,
    Sort data by 'Variant Frequency' in descending order.
       
    Return cleaned VCF dataframe
    '''
    vcf = pd.read_csv(file, usecols = ['Minimum', 'Maximum', 'Amino Acid Change', 'Change', 'Codon Change', 'locus_tag', 'Polymorphism Type', 'Protein Effect',
       'Variant Frequency'])
    vcf = vcf[~vcf['Protein Effect'].isnull()].loc[vcf['Protein Effect'] != 'None'].loc[vcf['locus_tag'] != 'SCO3798']
    vcf['Variant Frequency']= list(map(lambda x: float(x[0].replace('%', ''))/100, vcf['Variant Frequency'].str.split( '->')))
    vcf.sort_values(by = ['Variant Frequency'], ascending = False, inplace = True)
    
    return vcf


# In[137]:


import glob

def compare_to_wt(wt_vcf, mutant_vcf):
    common = wt_vcf.merge(mutant_vcf, on = ['Minimum', 'Maximum'])
    mutant_to_wt = mutant_vcf[(~mutant_vcf['Minimum'].isin(common['Minimum'])) & (~mutant_vcf['Maximum'].isin(common['Maximum']))].set_index('locus_tag')
    
    return mutant_to_wt

files = glob.glob('G:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\Genomics\M145 Evolved mutants\*.csv')

strain_names = list(map(lambda f: f.split('\\')[6].split('_')[0], files))
vcfs = list(map(lambda f: vcf_clean(f), files))
VCFs = dict(zip(strain_names, vcfs))

wt = VCFs['WT']
mutant_names = strain_names[:-1]
mutants = list(map(lambda key:VCFs[key], mutant_names))
mutants_to_wt = list(map(lambda m: compare_to_wt(wt, m), mutants))
dict_mutants_to_wt = dict(zip(mutant_names, mutants_to_wt))

dict_mutants_to_wt


# In[144]:


# find locus tag with mutations of variant freq > 0.5

def find_locus(dict_mutants_to_wt, variant_freq, mutant_name):
    mutant2wt = dict_mutants_to_wt[mutant_name]
    return mutant2wt[mutant2wt['Variant Frequency'] > 0.5].index.values
    
variant_freq = 0.5
mutation_genes = {}
for name in mutant_names:
    mutation_genes[name] = find_locus(dict_mutants_to_wt, variant_freq, name)
#mutation_genes

df = pd.DataFrame.from_dict(mutation_genes, orient="index")
df.to_csv("G:\\Ye\\Evolution mutants\\mutation_genes.csv")

