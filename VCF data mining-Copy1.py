
# coding: utf-8

# In[138]:


import pandas as pd
import numpy as np
import glob
import matplotlib as mp
from functools import reduce


# In[2]:


def vcf_clean(file):
    '''
    Input = Exported VCF from Geneious,
    Select columns ['Minimum', 'Maximum', 'Amino Acid Change', 'Change', 'Codon Change', 'locus_tag', 'Polymorphism Type', 'Protein Effect',
                    'Variant Frequency', 'Note'],
    Remove rows with 'Protein Effect' == NaN or None,
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


# In[3]:


def list_crosspop_mutation(wt_vcf, vf_threshold):
    ''' 
    Input = the wt_vcf dataframe from vcf_clean function, the threshold of Variant Frequency (e.g. 0.9)
    Keep the existed variances across the population with variant freq above threshold from time0 in our barcoded strain .
    Output = wt_vcf_crosspop
    '''
    wt_vcf_crosspop = wt_vcf.loc[wt_vcf['Variant Frequency'] >= vf_threshold]
    
    return wt_vcf_crosspop


# In[4]:


def clean_WT(wt_vcf, vf_threshold):
    ''' 
    Input = the wt_vcf dataframe from vcf_clean function, the threshold of Variant Frequency
    Remove all rows with 'Variant Frequency' >= vf_threshold (e.g. 0.90) 
    so that the existed variances across the population from time0 in our barcoded strain are excluded from following analysis.
    Output = wt_vcf_new
    '''
    wt_vcf_new = wt_vcf.loc[wt_vcf['Variant Frequency'] < vf_threshold]
    
    return wt_vcf_new


# In[5]:


def compare_to_wt(wt_vcf_new, wt_vcf, mutant_vcf):
    '''
    Input = the cleaned wt_vcf_new where variances across the population have been removed, 
            the wt_vcf including all variances
            the mutant_vcf cleaned by vcf_clean()
    Output = dataframe 'common' where common variances were found in both wt_vcf_new and the mutant_vcf,
             dataframe 'mutant_to_wt' where new variances were found only in mutant_vcf
    '''
    common = wt_vcf_new.merge(mutant_vcf, on = ['Minimum', 'Maximum', 'Change'])
    #mutant_to_wt = mutant_vcf[(~mutant_vcf['Minimum'].isin(common['Minimum'])) 
                              #& (~mutant_vcf['Maximum'].isin(common['Maximum'])) 
                              #& (~mutant_vcf['Change'].isin(common['Change']))]
    mutant_to_wt = mutant_vcf[(~mutant_vcf['Minimum'].isin(wt_vcf['Minimum'])) 
                              & (~mutant_vcf['Maximum'].isin(wt_vcf['Maximum'])) 
                              & (~mutant_vcf['Change'].isin(wt_vcf['Change']))]       
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


# In[6]:


files = glob.glob('G:\Dropbox (Vetsigian lab)\Vetsigian lab Team Folder\Ye\Genomics\M145 Evolved mutants\*.csv')

strain_names = list(map(lambda f: f.split('\\')[6].split('_')[0], files))
vcfs = list(map(lambda f: vcf_clean(f), files))
VCFs = dict(zip(strain_names, vcfs))
wt = VCFs['WT']
vf_threshold = 0.90
wt_new = clean_WT(wt, vf_threshold)

# Compare each mutant df in VCFs to the cleaned new WT vcf, wt_new. 
# This will generate a two layered dictionary where outside keys are the mutant names, 
# and the inside keys are "common" and "variance" that corresponding to the variance shared by wt and mutants 
# and to the variance found only in mutants respectively.

mutant_names = strain_names[:-1]
dict_mutants_to_wt = {}
for name in mutant_names:
    [common, variance]= compare_to_wt(wt_new, wt, VCFs[name])
    dict_mutants_to_wt[name] = {}
    dict_mutants_to_wt[name]['common'] = common
    dict_mutants_to_wt[name]['variance'] = variance
wt_new['locus_tag'].fillna('intergenetic region', inplace = True)
wt_new.set_index('locus_tag', inplace = True)


# In[150]:


dict_mutants_to_wt['W3']['variance']


# In[58]:


# Compare variance frequencies of mutatations shared by wt and mutants ('the common data set and wt data set')
common_vf_comp = {}
for mutant in mutant_names:
    common_vf_comp[mutant] = dict_mutants_to_wt[mutant]['common'][['Variant Frequency WT', 'Variant Frequency Mutant', 'note']]
    

common_var_df = dict_mutants_to_wt['R1']['common']
lab1 = common_var_df.index.values
lab2 = common_var_df['Codon Change']
lab_list = list(zip(lab1, lab2))
tuple(map(lambda tup: tup[0] + ',' + tup[1], lab_list))


# In[62]:


# Build a bar plot to see the common variance change from wt to each mutant
import matplotlib.pyplot as plt

def barplot_common_variance(dict_mutants_to_wt, mutant_name, bar_width, opacity, xtick_rot):
    common_var_df = dict_mutants_to_wt[mutant_name]['common']
    n_group = len(common_var_df.index.values)
    vf_wt = common_var_df['Variant Frequency WT']
    vf_mutant = common_var_df['Variant Frequency Mutant']
    fig, ax = plt.subplots()
    index = np.arange(n_group)
    rects1 = ax.bar(index, vf_wt, bar_width,
                alpha=opacity, color='b', label='WT')
    rects2 = ax.bar(index + bar_width, vf_mutant, bar_width,
                alpha=opacity, color='r', label= mutant_name)
    ax.set_ylabel('Variant Frequency')
    ax.set_title('Frequency of Common Variance shared by WT and ' + mutant_name)
    ax.set_xticks(index + bar_width / 2)
    xticklab1 = common_var_df.index.values 
    xticklab2 = common_var_df['Codon Change']
    lab_list = list(zip(xticklab1, xticklab2))
    xticklab = tuple(map(lambda tup: str(tup[0]) + ',' + str(tup[1]), lab_list))    
    ax.set_xticklabels(xticklab, rotation = xtick_rot)
    ax.legend()
    plt.show()

barplot_common_variance(dict_mutants_to_wt, 'WS1', 0.2, 0.6, 90)   


# In[87]:


# Barplot the variant frequency by mutant group: R (R1 to R5), RH (RH1 to RH5), W(W1 to W5), WS(WS1 and WS2)

# Group the mutants into their phenotypic group
mutant_group = {'R': list(), 'RH':list(), 'W': list(), 'WS': list()}
for key in dict_mutants_to_wt.keys():
    group_name = ''.join([i for i in key if not i.isdigit()])
    if group_name in mutant_group.keys():
        mutant_group[group_name]. append(key)

mutant_names = mutant_group['WS']  
for mut in mutant_names:
    barplot_common_variance(dict_mutants_to_wt, mut, 0.2, 0.6, 90)

mutant_group


# In[152]:


# Investigate variants novel in mutants 

def find_locus(dict_mutants_to_wt, variant_freq, mutant_name):
    '''
    input = dict_mutants_to_wt, variant_frequency and the mutant_name
    output = dataframe of novel mutations only found in the mutants with variant freq > threshold
    '''
    mutant2wt = dict_mutants_to_wt[mutant_name]['variance']
    return mutant2wt[mutant2wt['Variant Frequency'] > variant_freq][['Codon Change', 'Polymorphism Type', 'note', 'Variant Frequency']]

def find_locus4all_mutants(dict_mutants_to_wt, variant_freq):
    '''
    input = dict_mutants_to_wt, variant_frequency 
    output = dictionary of dataframes of novel mutations only found in each mutant with variant freq > threshold
    '''
    mutation_genes = {}
    for name in dict_mutants_to_wt.keys():
        mutation_genes[name] = find_locus(dict_mutants_to_wt, variant_freq, name)
    return mutation_genes

novel_mut_all = find_locus4all_mutants(dict_mutants_to_wt, 0)
#df = pd.DataFrame.from_dict(mutation_genes, orient="index")
#df.to_csv("G:\\Ye\\Evolution mutants\\mutation_genes.csv")
novel_mut_all['W3']


# In[154]:


# Exclude the mutation locus in intergenetic region and SCO3798
def find_mut_in_genes(novel_mut_all, mutant_name, exclude1 = True, exclude2 = True):
    '''
    input = dictionary of df of novel mutations found in each mutant by function find_locus4all_mutants, mutant name
            if exclude1 = True then SCO3798 will be excluded from the final df,
            if exclude2 = True then the intergenetic region will be exclude from the final df
    output = df of mutations in genes other than SCO3798 or intergenetic region with info of locus, change , note and variant freq    
    '''
    
    df = novel_mut_all[mutant_name]
    if exclude1:
        df = df.loc[df.index.values != 'SCO3798'] 
    if exclude2:
        df = df.loc[df.index.values != 'intergenetic region']
    return df

find_mut_in_genes(novel_mut_all, 'R1' , exclude1 = True, exclude2 = True)

def find_mut_in_genes4all_mutants(novel_mut_all, exclude1 = True, exclude2 = True):
    '''
    input = dictionary of df of novel mutations found in each mutant by function find_locus4all_mutants, 
            if exclude1 = True then SCO3798 will be excluded from the final df,
            if exclude2 = True then the intergenetic region will be exclude from the final df
    output = dictionary of dataframes of df of mutations in genes other than SCO3798 or intergenetic region with info of locus, 
             change , note and variant freq    
    '''
    dict_mut_exclude12 = {}
    for name in novel_mut_all.keys():
        dict_mut_exclude12[name] = find_mut_in_genes(novel_mut_all, name, exclude1, exclude2)     
    return dict_mut_exclude12

exclusive_mut_in_mutants = find_mut_in_genes4all_mutants(novel_mut_all, exclude1 = True, exclude2 = True)
exclusive_mut_in_mutants['W3']


# In[157]:



def barplot_mut_var_freq(exclusive_mut_in_mutants, mutant_name, bar_width, opacity, xtick_rot):
    vf = exclusive_mut_in_mutants[mutant_name]['Variant Frequency']
    n_group = len(vf)
    fig, ax = plt.subplots()
    index = np.arange(n_group)
    rects1 = ax.bar(index, vf, bar_width,
                alpha=opacity, color='r', label= mutant_name, align = 'edge')
    ax.set_ylabel('Variant Frequency')
    ax.set_title('Frequency of Exclusive Variance in ' + mutant_name)
    ax.set_xticks(index + bar_width / 2)
    xticklab1 = exclusive_mut_in_mutants[mutant_name].index.values 
    xticklab2 = exclusive_mut_in_mutants[mutant_name]['Codon Change']
    xticklab3 = exclusive_mut_in_mutants[mutant_name]['Polymorphism Type']
    lab_list = list(zip(xticklab1, xticklab2, xticklab3))
    xticklab = tuple(map(lambda tup: str(tup[0]) + ',' + str(tup[1]) + ',' + str(tup[2])[0:3], lab_list))    
    ax.set_xticklabels(xticklab, rotation = xtick_rot)
    ax.legend()
    plt.show()
    
barplot_mut_var_freq(exclusive_mut_in_mutants, 'WS1', 0.4, 0.6, 90)


# In[158]:


mut_in_group = []
for mut_name in mutant_group['R']:
    lab1 = exclusive_mut_in_mutants[mut_name].index.values
    lab2 = exclusive_mut_in_mutants[mut_name]['Codon Change']
    lab3 = exclusive_mut_in_mutants[mut_name]['Polymorphism Type']
    lab = list(zip(lab1, lab2, lab3))
    mut_in_group.extend(tuple(map(lambda tup: str(tup[0]) + ',' + str(tup[1]) + ',' + str(tup[2])[0:3], lab)))
    lab = []
set(mut_in_group)


# In[160]:


# find all non-redundant mutations in each group
def find_non_redund_mut_by_group(exclusive_mut_in_mutants, mutant_group, group_name):
    
    mut_in_group = []
    for mut_name in mutant_group[group_name]:
        lab1 = exclusive_mut_in_mutants[mut_name].index.values
        lab2 = exclusive_mut_in_mutants[mut_name]['Codon Change']
        lab3 = exclusive_mut_in_mutants[mut_name]['Polymorphism Type']
        lab = list(zip(lab1, lab2, lab3))
        mut_in_group.extend(tuple(map(lambda tup: str(tup[0]) + ',' + str(tup[1]) + ',' + str(tup[2])[0:3], lab)))
        lab = [] # set this to empty or multiple run will cause replicated labs
    return set(mut_in_group)

all_labs = find_non_redund_mut_by_group(exclusive_mut_in_mutants, mutant_group, 'WS')
all_labs


# In[165]:


# Merge index (locus_tag) and Change cols in df of exclusive_mut_in_mutants by group

def rearrange_cols_mut(exclusive_mut_in_mutants, mutant_name):
    df = exclusive_mut_in_mutants[mutant_name]
    genes = df.index.values
    changes = df['Codon Change'] 
    pol = df['Polymorphism Type']
    lab = list(zip(genes, changes, pol))
    labs = tuple(map(lambda tup: str(tup[0]) + ',' + str(tup[1]) + ',' + str(tup[2])[0:3], lab))
    df['labels'] = labs
    df = df.set_index('labels').drop(['Codon Change', 'Polymorphism Type'], axis = 1)
    return df

rearrange_cols_mut(exclusive_mut_in_mutants, 'W3')


# In[167]:


def mut_locus_freq_by_group(exclusive_mut_in_mutants, mutant_group, group_name):
    all_labs = find_non_redund_mut_by_group(exclusive_mut_in_mutants, mutant_group, group_name)
    mutant_names = mutant_group[group_name]
    mut_locus_freq = pd.DataFrame(index = all_labs, columns = mutant_names)
    mut_locus_freq['note'] = ''
    for mut_name in mutant_names:
        df = rearrange_cols_mut(exclusive_mut_in_mutants, mut_name)
        for lab in all_labs:
            if lab in df.index.values:
                mut_locus_freq.loc[lab][mut_name] = df.loc[lab]['Variant Frequency']
                mut_locus_freq.loc[lab]['note'] = df.loc[lab]['note']
            else:
                mut_locus_freq.loc[lab][mut_name] = 0
    return mut_locus_freq

mut_locus_freq_by_group(exclusive_mut_in_mutants, mutant_group, 'W')


# In[100]:


def barplot_mut_vf_by_group(exclusive_mut_in_mutants, mutant_group, group_name, bar_width, opacity, xtick_rotate = 90):
    vf = mut_locus_freq_by_group(exclusive_mut_in_mutants, mutant_group, group_name)
    title = 'Variant Frequency of Mutations in Mutant group ' + group_name
    
    ax = vf.plot(kind = 'bar', width = bar_width, title = title)
    ax.set_ylabel('Variant Frequency')
    return ax


# In[101]:



barplot_mut_vf_by_group(exclusive_mut_in_mutants, mutant_group, 'WS', 0.3, 0.2, 90)

