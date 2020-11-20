#!/usr/bin/env python
#generate_samples.py

""" Generates and downloads new samples of protein PDB IDs.

	Generates new samples of protein PDB IDs from RCSB proteins 
	website and downloads the PDB files to a local directory.

	Attributes:
		path: Directory where the downloaded PDB files are stored.
		all_protein_from_pdb_site_file: A text file containing list of all 
										downloadable PDB IDs at RCSB website.
		no_of_samples: Number of non-overlapping samples required.
		no_of_proteins: Number of protein PDB IDs required per sample.
"""


import re,random,os,argparse
import urllib, urllib2,requests
import pandas as pd
import numpy as np

__author__ = "Venkata Sarika Kondra"

__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"


parser = argparse.ArgumentParser(description='Adaptive Equal Frequency Binning.')
parser.add_argument('--path', '-od',  metavar='path', 
    default="//home//C00219805//Research//Protien_Database//extracted_new_samples//testing//" ,
    help='Directory where the downloaded PDB files are stored.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
	default='t1', help='Name of the sample on which this script should be run.')
parser.add_argument('--all_pdb_file',  metavar='all_pdb_file', 
    default="//home//C00219805//Research//Protien_Database//secondary_structures.csv" ,
    help='A text file containing list of all downloadable PDB IDs at RCSB website.')
parser.add_argument('--all_details_file',  metavar='all_pdb_file', 
    default="/media/c00219805/Elements/Research/Hierarchical_Classification/Database/SCOP/all_details.csv" ,
    help='A text file containing list of all downloadable PDB IDs at RCSB website.')
parser.add_argument('--no_of_samples',  metavar='no_of_samples', 
    default=2 , help='Number of non-overlapping samples required.')
parser.add_argument('--no_of_proteins',  metavar='no_of_proteins', 
    default=199, help='Number of protein PDB IDs required per sample.')
parser.add_argument('--group',  metavar='group', 
    default='all', help='Name of group for which PDB IDs are required.')
parser.add_argument('--prepare_details_file', '-prepare_details_file', \
    action='store_true', default=False, \
    help='Enable this option if details file needed to be created from all_details file.')										


def prepare_details_file(df, no_of_proteins):
	""" Generates new samples of protein PDB IDs.

	Args:
		no_of_proteins: Number of non-overlapping samples required.

	Returns:
		A dataframe with randomly selected proteins in each group.

	"""

	#Remove chain duplicates
	duplicateRowsDF = df[df.duplicated(['protein'])]
	print(len(df), len(duplicateRowsDF))
	df = df[~df.isin(duplicateRowsDF)].dropna()
	print(len(df))
	print(df.head(5))
	print(df["fold_name"].value_counts())
	
	return df


if __name__ == '__main__':
	"""Executable code from command line."""

	args = parser.parse_args()
	non_functional_urls = []

	#Generate random samples.
	#Comment this if you already have the list of PDB ids to be downloaded.
	
	arr = ['1h7w', '1jmx', '1rq6', '2j9u', '1d9c', '2ciw']
	if args.prepare_details_file:
		df = pd.read_csv(args.all_details_file,
			dtype= {"hierarchy": str, "fold_id_short" : str, "class": str, \
			"fold_id_long": str, "fold_name": str, "protein":str, "chain":str, "aa": str},
			usecols = ["hierarchy", "fold_id_short", "class", \
			"fold_id_long", "fold_name", "protein", "chain", "aa"]).fillna("Unknown")
		# if args.group != 'all':
		# 	df = df[df['class'] == args.group]
		df = df[df['protein'].isin(arr)]

		df_details = prepare_details_file(df, args.no_of_proteins)
		df_details.to_csv(os.path.join(args.path, args.sample_name, "sample_details.csv"))

	#Uncomment the next two line, for downloading PDB files of a particular sample when the PDB ids are known.	
	samples_file = pd.read_csv('{}{}//sample_details.csv'.format(args.path,str(args.sample_name)))
	sample = samples_file['protein'].map(str).values
	for id in sample: 
			#Generate the downloadable URL
			url = "http://files.rcsb.org/download/{}.pdb".format(id.replace('+', ''))
			#Check if the URL exists, otherwise skip that PDB ID
			request = requests.get(url)
			if request.status_code == 200:
				outputFilename = "{}{}//{}.pdb".format(args.path,str(args.sample_name),str(id.replace('+', '')))
				response = urllib2.urlopen(url)
				zippedData = response.read()
				# save data to disk
				output = open(outputFilename,'w')
				output.write(zippedData)
				output.close()
				print(url , " extracted to ",outputFilename)
			else:
				non_functional_urls.append(id)
				continue
	print("Un available URLs- ", non_functional_urls)
	print('End sample generation and storing to disk')

