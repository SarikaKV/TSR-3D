#!/usr/bin/env python
#generate_samples.py

""" 
"""


import glob,re,random,os,argparse
import urllib, urllib2,requests
import pandas as pd
import numpy as np

__author__ = "Venkata Sarika Kondra"

__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"


parser = argparse.ArgumentParser(description='Adaptive Equal Frequency Binning.')
parser.add_argument('--path', '-path',  metavar='path', 
    default="//media//c00219805//Elements1//Research//Protien_Database//extracted_new_samples//testing//" ,
    help='Directory where the downloaded PDB files are stored.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
	default='sample_ms_serine_protease_kinase_phosphatase_100', help='Name of the sample on which this script should be run.')

args = parser.parse_args()
	
subFolder = args.sample_name
sampleDetailsFile = os.path.join(args.path, subFolder, 'sample_details.csv')
df = pd.read_csv(sampleDetailsFile)
print(df.head(5))
sample_dict = dict(zip(df.protein, df.group))
print(sample_dict)

dali_zscore = os.path.join(args.path, subFolder, 'DALI_zscore_similarity_100proteins.csv')
tm_rmsd = os.path.join(args.path, subFolder, 'TMalign_rmsd_matrix_100prots.csv')
tm_zscore_p1 = os.path.join(args.path, subFolder, 'TMalign_tmscore_matrix_prot1_100prots.csv')
tm_zscore_p2 = os.path.join(args.path, subFolder, 'TMalign_tmscore_matrix_prot2_100prots.csv')
tm_zscore_avg = os.path.join(args.path, subFolder, 'TMalign_tmscore_matrix_avg_100prots.csv')
ce_rmsd = os.path.join(args.path, subFolder, 'CE_rmsd_matrix_100prots_new.csv')
ce_zscore = os.path.join(args.path, subFolder, 'CE_zscore_matrix_100prots_new.csv')
generalised = os.path.join(args.path, subFolder, 'generalised.csv')

df_dali_zscore = pd.read_csv(dali_zscore, sep = '\t', header = None, index_col = 0)
df_dali_zscore.index =  [ sample_dict[x.upper()] +'-'+ x.upper() for x in df_dali_zscore.index]
df_dali_zscore.columns = df_dali_zscore.index
#print(df_dali_zscore)
df_dali_zscore.to_csv(os.path.join(args.path, subFolder, 'dali_zscore.csv'))

df_tm_rmsd = pd.read_csv(tm_rmsd,  header = 0, index_col = 0)
df_tm_rmsd = df_tm_rmsd[df_tm_rmsd.columns[:-1]]
df_tm_rmsd.index =  [sample_dict[x.upper()] +'-'+ x.upper()  for x in df_tm_rmsd.columns]
df_tm_rmsd.columns = df_tm_rmsd.index
#print(df_tm_rmsd.head(5))
df_tm_rmsd.to_csv(os.path.join(args.path, subFolder, 'tm_rmsd.csv'))

df_tm_zscore_p1 = pd.read_csv(tm_zscore_p1,  header = 0, index_col = 0)
new_cols = df_tm_zscore_p1.columns[1:]
df_tm_zscore_p1 = df_tm_zscore_p1[df_tm_zscore_p1.columns[:-1]]
df_tm_zscore_p1.columns = new_cols
df_tm_zscore_p1.index =  [sample_dict[x.upper()] +'-'+ x.upper()  for x in df_tm_zscore_p1.columns]
df_tm_zscore_p1.columns = df_tm_zscore_p1.index
#print(df_tm_zscore_p1.head(5))
df_tm_zscore_p1.to_csv(os.path.join(args.path, subFolder, 'tm_zscore_p1.csv'))

df_tm_zscore_p2 = pd.read_csv(tm_zscore_p2,  header = 0, index_col = 0)
new_cols = df_tm_zscore_p2.columns[1:]
df_tm_zscore_p2 = df_tm_zscore_p2[df_tm_zscore_p2.columns[:-1]]
df_tm_zscore_p2.columns = new_cols
df_tm_zscore_p2.index =  [sample_dict[x.upper()] +'-'+ x.upper()  for x in df_tm_zscore_p2.columns]
df_tm_zscore_p2.columns = df_tm_zscore_p2.index
#print(df_tm_zscore_p2.head(5))
df_tm_zscore_p2.to_csv(os.path.join(args.path, subFolder, 'tm_zscore_p2.csv'))

df_tm_zscore_avg = pd.read_csv(tm_zscore_avg,  header = 0, index_col = 0)
new_cols = df_tm_zscore_avg.columns[1:]
df_tm_zscore_avg = df_tm_zscore_avg[df_tm_zscore_avg.columns[:-1]]
df_tm_zscore_avg.columns = new_cols
df_tm_zscore_avg.index =  [sample_dict[x.upper()] +'-'+ x.upper()  for x in df_tm_zscore_avg.columns]
df_tm_zscore_avg.columns = df_tm_zscore_avg.index
#print(df_tm_zscore_avg.head(5))
df_tm_zscore_avg.to_csv(os.path.join(args.path, subFolder, 'tm_zscore_avg.csv'))

df_ce_rmsd = pd.read_csv(ce_rmsd,  header = 0, index_col = 0)
new_cols = df_ce_rmsd.columns[1:]
df_ce_rmsd = df_ce_rmsd[df_ce_rmsd.columns[:-1]]
df_ce_rmsd.columns = new_cols
df_ce_rmsd.index =  [sample_dict[x[:4].upper()] +'-'+ x[:4].upper()  for x in df_ce_rmsd.columns]
df_ce_rmsd.columns = df_ce_rmsd.index
#print(df_ce_rmsd.head(5))
df_ce_rmsd.to_csv(os.path.join(args.path, subFolder, 'ce_rmsd.csv'))

df_ce_zscore = pd.read_csv(ce_zscore,  header = 0, index_col = 0)
new_cols = df_ce_zscore.columns[1:]
df_ce_zscore = df_ce_zscore[df_ce_zscore.columns[:-1]]
df_ce_zscore.columns = new_cols
df_ce_zscore.index =  [sample_dict[x[:4].upper()] +'-'+ x[:4].upper()  for x in df_ce_zscore.columns]
df_ce_zscore.columns = df_ce_zscore.index
print(df_ce_zscore.head(5))
print(df_ce_zscore.index)
df_ce_zscore.to_csv(os.path.join(args.path, subFolder, 'ce_zscore.csv'))

df_generalised = pd.read_csv(generalised,  header = 0, index_col = 0)
df_generalised.index = df_generalised.columns
df_generalised.to_csv(os.path.join(args.path, subFolder, 'generalised.csv'))

#print(df_generalised)