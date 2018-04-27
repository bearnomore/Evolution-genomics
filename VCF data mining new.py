
# coding: utf-8

# In[150]:


import pandas as pd
import numpy as np
import glob


# In[180]:


def vcf_clean(file):
    '''
    Input = Exported VCF from Geneious,
    Select columns ['Minimum', 'Maximum', 'Amino Acid Change', 'Change', 'Codon Change', 'locus_tag', 'Polymorphism Type', 'Protein Effect',
                    'Variant Frequency', 'Note'],
    Remove rows with 'Protein Effect' == NaN or None,
    Remove rows of 'locus_tag' == 'SCO3798' where the barcode is inserted,
    Sort data by 'Variant Frequency' in descending order.
       
    Return cleaned VCF dataframe
    '''
    vcf = pd.read_csv(file, usecols = ['Minimum', 'Maximum', 'Change', 'locus_tag', 'Polymorphism Type',
                                       'Protein Effect', 'Variant Frequency', 'Amino Acid Change','Codon Change', 'note'])
    ## vcf = vcf[~vcf['Protein Effect'].isnull()].loc[vcf['Protein Effect'] != 'None'].loc[vcf['locus_tag'] != 'SCO3798']
    # Remove SNPs that caused no effect on protein encoding genes, but maintain all other intra/inter genetic changes
    vcf = vcf.loc[vcf['Protein Effect'] != 'None']
    # Clean the variant frequencies and pick up the floor of the range
    vcf['Variant Frequency']= list(map(lambda x: float(x[0].replace('%', ''))/100, vcf['Variant Frequency'].str.split( '->')))
    vcf.sort_values(by = ['Variant Frequency'], ascending = False, inplace = True)
    
    return vcf


# In[191]:


def clean_WT(wt_vcf, vf_threshold):
    ''' 
    Input = the wt_vcf dataframe from vcf_clean function, the threshold of Variant Frequency
    Remove all rows with 'Variant Frequency' >= vf_threshold (e.g. 0.95) 
    so that the existed variances across the population from time0 in our barcoded strain are excluded from following analysis.
    Output = wt_vcf_new
    '''
    wt_vcf_new = wt_vcf.loc[wt_vcf['Variant Frequency'] < vf_threshold]
    
    return wt_vcf_new


# In[294]:


def compare_to_wt(wt_vcf_new, mutant_vcf):
    '''
    Input = the cleaned wt_vcf_new where variances across the population have been removed, the mutant_vcf cleaned by vcf_clean()
    Output = dataframe common where common variances were found between wt_vcf_new and the mutant_vcf,
             dataframe mutant_to_wt where new variances are found
    '''
    common = wt_vcf_new.merge(mutant_vcf, on = ['Minimum', 'Maximum', 'Change'])
    mutant_to_wt = mutant_vcf[(~mutant_vcf['Minimum'].isin(common['Minimum'])) 
                              & (~mutant_vcf['Maximum'].isin(common['Maximum'])) 
                              & (~mutant_vcf['Change'].isin(common['Change']))]
    mutant_to_wt['locus_tag'].fillna('intergenetic region', inplace = True)
    mutant_to_wt.set_index('locus_tag', inplace = True)
 
    
    # clean up the common df a bit by removing repeated columns
    common = common[['Minimum', 'Maximum', 'Change', 'locus_tag_x',
                    'Polymorphism Type_x', 'Protein Effect_x', 'Variant Frequency_x',
                    'Amino Acid Change_x', 'Codon Change_x', 'note_x', 'Variant Frequency_y']]
    common.rename(index=str, columns={'locus_tag_x':'locus_tag',"Polymorphism Type_x": "Polymorphism Type", "Protein Effect_x": "Protein Effect",
                                     'Variant Frequency_x': 'Variant Frequency WT', 'Amino Acid Change_x': 'Amino Acid Change',
                                     'Codon Change_x': 'Codon Change', 'note_x': 'note', 'Variant Frequency_y': 'Variant Frequency Mutant'}, inplace = True)
    common['locus_tag'].fillna('intergenetic region', inplace = True)
    common.set_index('locus_tag', inplace = True)
    
    return common, mutant_to_wt


# In[295]:


files = glob.glob('G:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\Genomics\M145 Evolved mutants\*.csv')

strain_names = list(map(lambda f: f.split('\\')[6].split('_')[0], files))
vcfs = list(map(lambda f: vcf_clean(f), files))
VCFs = dict(zip(strain_names, vcfs))
wt = VCFs['WT']
vf_threshold = 0.95
wt_new = clean_WT(wt, vf_threshold)

# Compare each mutant df in VCFs to the cleaned new WT vcf, wt_new. 
# This will generate a two layered dictionary where outside keys are the mutant names, 
# and the inside keys are "common" and "variance" that corresponding to variance found 
# in the wt_new and variance not found in wt_new respectively.

mutant_names = strain_names[:-1]
dict_mutants_to_wt = {}
for name in mutant_names:
    [common, variance]= compare_to_wt(wt_new, VCFs[name])
    dict_mutants_to_wt[name] = {}
    dict_mutants_to_wt[name]['common'] = common
    dict_mutants_to_wt[name]['variance'] = variance
wt_new['locus_tag'].fillna('intergenetic region', inplace = True)
wt_new.set_index('locus_tag', inplace = True)


# In[296]:


dict_mutants_to_wt['R1']


# In[298]:


# Compare variance frequencies of mutatations shared by wt and mutants ('the common data set and wt data set')
common_vf_comp = {}
for mutant in mutant_names:
    common_vf_comp[mutant] = dict_mutants_to_wt[mutant]['common'][['Variant Frequency WT', 'Variant Frequency Mutant']]
    
common_vf_comp


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

