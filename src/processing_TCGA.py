#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Pool, cpu_count
import functools
import os

os.chdir("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/Tools/NBS2_NetPropagation")

# Parse pathway features
entrez2symbol = {}
with open('data/Homo_sapiens.gene_info') as f:
    for line in f.read().splitlines():
        row = line.split('\t')
        entrez2symbol[row[1]]=row[2]

# These are the pathways related to cancer hallmarks, which will be used as edge features.  
# Referance: http://www.cell.com/abstract/S0092-8674(11)00127-9
gene2pathway = {}
pathways = set()
with open('data/hallmarks.txt') as f:
    for line in f.read().splitlines():
        row = line.split('\t')
        if len(row) > 1:
            pathway = row[0].split('|')[1]
            pathways.add(pathway)
            for entrez in row[2:]:
                if entrez not in entrez2symbol:
                    continue
                gene = entrez2symbol[entrez]
                if gene not in gene2pathway:
                    gene2pathway[gene]=set()
                gene2pathway[gene].add(pathway)

with open('data/gene2pathway_genes.txt', 'w') as f:
    f.write("\n".join(gene2pathway.keys()))

with open('data/gene2pathway_pathways.txt', 'w') as f:
    f.write("\n".join([",".join(x) for x in gene2pathway.values()]))

# Parse tumor types and tumor subtypes
def parse_types(fn):
    pat2type = {}
    pat2subtype = {}
    with open(fn) as f:
        for line in f.read().rstrip().splitlines()[1:]:
            row = line.split("\t")
            pat = row[0][:12]
            pat2type[pat] = row[4]
            pat2subtype[pat] = row[3]
        return pat2type, pat2subtype

pat2type, pat2subtype = parse_types('data/Cancer_subtypes_list_updated.tsv')


# Load existing datasets
training_set = pd.read_table('data/somaticAlterations/TCGA_alterations_training.tsv', index_col=0)


# Fisher test against types (Optional)

training_pat2type = training_set.index.to_series().map(pat2type)

gene2fisherp = {}
for gene in training_set.columns:
    gene2fisherp[gene] = 0.
#     gene2fisherp[gene] = {}
    for subtype in training_pat2type.unique():
        tab = pd.crosstab(training_set[gene]>0, training_pat2type==subtype)
        if tab.shape != (2, 2):
            continue
        oddsratio, pvalue = fisher_exact(tab, alternative='greater')
        logp = np.log10(pvalue)
        gene2fisherp[gene] = max(-logp, gene2fisherp[gene])
#         gene2fisherp[gene][subtype] = -logp
        # if pvalue < 0.05/len(training_set.index):
        #     print(gene, subtype, oddsratio, pvalue, tab)


# # Parse PathwayCommons
PathwayCommons = pd.read_table('data/PathwayCommons11.All.hgnc.txt')
PathwayCommons = PathwayCommons.loc[PathwayCommons.loc[:,'INTERACTION_TYPE'].isin(['controls-state-change-of',
                                                                                   'controls-transport-of',
                                                                                   'controls-phosphorylation-of',
                                                                                   'controls-expression-of',
                                                                                   'catalysis-precedes',
                                                                                   'in-complex-with',
                                                                                   'interacts-with', 
                                                                                   'neighbor-of']),:]


def parse_edge_features(mutrates, df):
    edge2features = {}
    for index, row in df.iterrows():
        g1 = row['PARTICIPANT_A']
        g2 = row['PARTICIPANT_B']
        
#         # (Optionally) filter by cancer genes or pathways
#         if not (g1 in cancergenes and g2 in cancergenes):
#             continue

#         # (Optionally) filter by frequently mutated genes
#         if not (g1 in recurrently_mutated_genes and g2 in recurrently_mutated_genes):
#             continue
            
        ty = row['INTERACTION_TYPE']
        ty_d = ty + '_d'
        ty_rev = ty + '_rev'
        sources = row['INTERACTION_DATA_SOURCE'].split(';')
        edge = g1 + '\t' + g2
        edge_rev = g2 + '\t' + g1
        if edge not in edge2features:
            edge2features[edge] = {'gene 1':g1, 'gene 2':g2}
        if edge_rev not in edge2features:
            edge2features[edge_rev] = {'gene 1':g2, 'gene 2':g1}

        # Parse edge type features
        if ty not in ['in-complex-with','interacts-with','neighbor-of']:
            edge2features[edge][ty_d] = 1.
            edge2features[edge_rev][ty_rev] = 1.
        else:
            edge2features[edge][ty] = 1.
            edge2features[edge_rev][ty] = 1.

        # Parse edge source features
        for source in sources:
            edge2features[edge][source] = 1.
            edge2features[edge_rev][source] = 1.

        # Parse pathway features. 
        # If one node is in the pathway, the score is 0.5
        # If both nodes are in the pathway, the score is 1
        if g1 in gene2pathway or g2 in gene2pathway:
            for pathway in pathways:
                edge2features[edge][pathway] = 0.
                edge2features[edge_rev][pathway] = 0.
            for g in [g1, g2]:
                if g in gene2pathway:
                    for pathway in gene2pathway[g]:
                        edge2features[edge][pathway] += 0.5
                        edge2features[edge_rev][pathway] += 0.5

        # Parse mutation features from the training set
        # Calculate mutation rates
        mutrate_g1 = 0
        mutrate_g2 = 0
        if g1 in training_set:
            mutrate_g1 = mutrates.loc[g1]
        if g2 in training_set:
            mutrate_g2 = mutrates.loc[g2]
        edge2features[edge]['mutrate_source'] = mutrate_g1
        edge2features[edge]['mutrate_target'] = mutrate_g2
        edge2features[edge_rev]['mutrate_source'] = mutrate_g2
        edge2features[edge_rev]['mutrate_target'] = mutrate_g1

###### To put later
#         # (optional) add fisher test against subtypes
#         if g1 in gene2fisherp:
#             edge2features[edge]['fisherp_source'] = gene2fisherp[g1]
#             edge2features[edge_rev]['fisherp_target'] = gene2fisherp[g1]
# #             for subtype in gene2fisherp[g1]:
# #                 edge2features[edge]['{}_source'.format(subtype)] = gene2fisherp[g1][subtype]
# #                 edge2features[edge_rev]['{}_target'.format(subtype)] = gene2fisherp[g1][subtype]
#         if g2 in gene2fisherp:
#             edge2features[edge]['fisherp_target'] = gene2fisherp[g2]
#             edge2features[edge_rev]['fisherp_source'] = gene2fisherp[g2]
# #             for subtype in gene2fisherp[g2]:
# #                 edge2features[edge]['{}_target'.format(subtype)] = gene2fisherp[g2][subtype]
# #                 edge2features[edge_rev]['{}_source'.format(subtype)] = gene2fisherp[g2][subtype]

        # Calculate mutual exclusivity / co-occurrence
        ME = 0.
        if g1 in training_set and g2 in training_set:
            if training_set.loc[:,g1].sum() < 4 or training_set.loc[:,g2].sum() < 4:
                continue
            tab = pd.crosstab(training_set.loc[:,g1],training_set.loc[:,g2])
            if tab.shape != (2, 2):
                continue
            if tab.iloc[1,1] * tab.iloc[0,0] >= tab.iloc[0,1] * tab.iloc[1,0]:
                continue
            oddsratio, pvalue = fisher_exact(tab, alternative='less')
            logp = np.log10(pvalue)
            ME = -logp
            if pvalue < 0.05:
                print(g1, g2, oddsratio, pvalue)
        edge2features[edge]['mutual_exclusive'] = ME
        edge2features[edge_rev]['mutual_exclusive'] = ME
            
    return edge2features


training_set_mutrate = training_set.sum() / training_set.shape[0]

n_processes = cpu_count()
pool = Pool(processes=n_processes)

df_split = np.array_split(PathwayCommons, n_processes, axis=0)
parse_edge_features_partial = functools.partial(parse_edge_features, training_set_mutrate)
edge2features_list = pool.map(parse_edge_features_partial, df_split)

pool.close()
pool.join()


from collections.abc import Mapping
def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


edge2features = {}
for i in range(len(edge2features_list)):
    update(edge2features, edge2features_list[i])


edge2features_df = pd.DataFrame.from_dict(edge2features, orient='index')


edge2features_df = edge2features_df.fillna(0).sort_index(1)
cols = list(edge2features_df)
cols.insert(0, cols.pop(cols.index('gene 2')))
cols.insert(0, cols.pop(cols.index('gene 1')))
edge2features_df = edge2features_df[cols]

geneset = set(training_set.sum().sort_values()[-5:].index)


for gene in sorted(geneset):
    edge2features_df[gene+'_source'] = (edge2features_df['gene 1']==gene).astype(int)
    edge2features_df[gene+'_target'] = (edge2features_df['gene 2']==gene).astype(int)

edge2features_df = edge2features_df.fillna(0).sort_index(1)
cols = list(edge2features_df)
cols.insert(0, cols.pop(cols.index('gene 2')))
cols.insert(0, cols.pop(cols.index('gene 1')))
edge2features_df = edge2features_df[cols]


edge2features_df.iloc[:,2:] = edge2features_df.iloc[:,2:].divide(edge2features_df.iloc[:,2:].max())
edge2features_df.dropna(axis=1, inplace=True)


edge2features_df.to_csv('data/TCGA_edge2features.txt', sep='\t', header=False, index=False)



with open('data/TCGA_feature_names.txt', 'w') as f:
    f.write('\n'.join(edge2features_df.columns[2:]))

