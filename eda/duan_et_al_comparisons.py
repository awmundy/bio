import pandas as pd
import numpy as np
pd.options.display.float_format = '{:,.3f}'.format

def safe_rename(df, rename_dict):
    """Rename columns in a dataframe, but only if they exist.
    """
    for old_col, new_col in rename_dict.items():
        assert old_col in df.columns, "Column %s not found in dataframe" % col
        assert new_col not in df.columns, "Column %s already exists in dataframe" % new_col
    df.rename(columns=rename_dict, inplace=True)
    return df


a = pd.read_csv(r'/media/awmundy/Windows/bio/age_related_steatohepatitis/outputs/young_vs_old/all_dge_table.csv')
# supplementary file posted to GEO
b = pd.read_excel(r'/media/awmundy/Windows/bio/age_related_steatohepatitis/supplemental_data/GSE216592_ct-vs-aging.all.annot2.xls.xlsx',
                  engine='openpyxl')

rename_rep = {'ct1_fpkm': 'young_1',
              'ct2_fpkm': 'young_2',
              'ct3_fpkm': 'young_3',
              'aging2_fpkm': 'old_1',
              'aging3_fpkm': 'old_2',
              'aging5_fpkm': 'old_3',
              'Symbol': 'gene_id',
              'FDR': 'adj_p',
              'log2(fc)' : 'logfc'
              }

rename_me = {'adj.P.Val': 'adj_p',
             'logFC': 'logfc',
             }

a = safe_rename(a, rename_me)
b = safe_rename(b, rename_rep)
cols = list(rename_rep.values())
a = a[cols].copy()
b = b[cols].copy()
a['adj_p_sig'] = np.where(a['adj_p'] < 0.05, 1, 0)
b['adj_p_sig'] = np.where(b['adj_p'] < 0.05, 1, 0)
print('my sig genes count:', a['adj_p_sig'].sum()) #1017
print('their sig genes count:', b['adj_p_sig'].sum()) #3186

# drop 24 duplicate records for genes, keeping higher value ones bc sometimes one copy is all 0
b.sort_values(by=['gene_id', 'young_1'], ascending=False, inplace=True)
b.drop_duplicates(subset='gene_id', inplace=True)

m = pd.merge(a, b, on='gene_id', how='outer', validate='1:1', suffixes=('_a', '_b'), indicator=True)
print(m['_merge'].value_counts())

# similar logfc for genes we both call significant
msk =m['_merge'] == 'both'
msk &= m['adj_p_sig_a'] == 1
print('my logfc')
print(m[msk]['logfc_a'].describe()) # mean .340
print('their logfc')
print(m[msk]['logfc_b'].describe()) # mean .307

# read in gene sets, adding flag to gene level data for if it's in at least one of the sets
gene_sets = pd.read_csv('/media/awmundy/Windows/bio/age_related_steatohepatitis/custom_gene_sets/custom_gene_sets.csv')
m['in_gene_sets']= np.where(m['gene_id'].isin(gene_sets['gene_symbol']), 1, 0)

#construct flag cols for if the counts are non 0
m['young_non_zero_b'] = 0
m['old_non_zero_b'] = 0
for col_type in ['young', 'old']:
    count_col = f'{col_type}_non_zero_b'
    cols = [f'{col_type}_{i}_b' for i in range(1, 4)]
    for col in cols:
        msk = ~np.isclose(m[col], 0)
        m.loc[msk, count_col] += 1
m['total_non_zero_b'] = m['young_non_zero_b'] + m['old_non_zero_b']

# 5 or 6 samples with non-zero values for genes in the gene sets
msk = m['in_gene_sets'] == 1
print('distribution of non-zero sample counts across genes in the gene set')
print(m[msk]['total_non_zero_b'].value_counts())

msk = np.isclose(m['adj_p_sig_b'],1)
print('distribution of non-zero sample counts across genes they consider significant')
print(m[msk]['total_non_zero_b'].value_counts())

# generally the right_onlies are less significant, both and left_only seem very similar
print('p values of genes only i have')
print(m[m['_merge'] == 'left_only']['adj_p_a'].describe())
print('my p values for genes we both have')
print(m[m['_merge'] == 'both']['adj_p_a'].describe())
print('their p values for genes we both have')
print(m[m['_merge'] == 'both']['adj_p_b'].describe())

# about 80% are not signficant
msk = m['_merge'] == 'right_only'
print('significance for genes only they have in their data')
m[msk]['adj_p_sig_b'].value_counts(dropna=False)

# same mean logfc for both of us
msk = m['in_gene_sets'] == 1
print('my logfc for genes in the gene sets')
print(m[msk]['logfc_a'].describe())
print('their logfc for genes in the gene sets')
print(m[msk]['logfc_b'].describe())
print('my logfc minus theirs')
print((m[msk]['logfc_a'] - m[msk]['logfc_b']).describe())

# fibrosis and inflammation gene sets have positive logfc which makes sense
g = pd.merge(a, gene_sets, left_on='gene_id', right_on='gene_symbol', how='inner', validate='1:m')
grp = g.groupby('gs_name', as_index=False).agg({'logfc': ['mean']})
print('gene set logfc means, my logfc')
print(grp)


# missing genes in all dge table investigation
gene_sets = pd.read_csv('/media/awmundy/Windows/bio/age_related_steatohepatitis/custom_gene_sets/custom_gene_sets.csv')
gene_sets.sort_values(by='gene_symbol', inplace=True)
# dge = pd.read_csv('/media/awmundy/Windows/bio/age_related_steatohepatitis/outputs/young_vs_old/all_dge_table.csv')
dge = pd.read_csv('/media/awmundy/Windows/bio/age_related_steatohepatitis/outputs/young_vs_old/dge_table.csv')
dge['in_gene_sets']= np.where(dge['gene_id'].isin(gene_sets['gene_symbol']), 1, 0)
# fibrosis and inflammation have almost all genes in the data, but fatty acid is missing half
# RESOLVED, for all but 3 genes in fatty acid and one in fibrosis, issue was whitespace
# fat = gene_sets.loc[gene_sets['gs_name'] == 'Fatty Acid Oxidation Custom'].copy()
gene_sets['in_dge'] = np.where(gene_sets['gene_symbol'].isin(dge['gene_id']), 1, 0)
msk = gene_sets['in_dge'] == 0
print(gene_sets[msk])