#!/usr/bin/env python
# coding: utf-8

# Filter out viruses
# Split taxonomy into different fields.

import pandas as pd
import argparse
import numpy as np

def argparser():
	"""Arguments"""
	epilog = """Example:
	$  filter_mmseqs_lca_taxonomy.py {SAMPLE}_LCA_tsv.tsv """

	parser = argparse.ArgumentParser(description=__doc__, epilog=epilog,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('filename')
	args = parser.parse_args()
	return args

def read_file_lines(filename):
	expected_num_fields = None
	lines_with_fields = []
	with open(filename, 'r') as file:
		for line in file:
    		# Split the line into fields using the tab separator
			fields = line.rstrip('').split('\t')
			# Some fields are missing the final field; need to ensure all lines have same number of fields
			# Initialize expected number of fields from the first line
			if expected_num_fields is None:
				expected_num_fields = len(fields)
        	# Pad missing fields with empty strings
			while len(fields) < expected_num_fields:
				fields.append('')
			# Select some elements from the fields (replace 0 and 2 with the indices you want)
			selected_elements = [fields[0], fields[1], fields[2], fields[3], fields[8]]
			# Append the selected elements to the list
			lines_with_fields.append(selected_elements)
		df = pd.DataFrame(lines_with_fields, columns=['Contig', 'Taxid', 'Tax_level', 'Best_tax', 'Lineage'])
		dtype_dict = {'Contig':'str', 'Taxid':'int', 'Tax_level':'str', 'Best_tax':'str', 'Lineage': 'str'}
		for k, v in dtype_dict.items():
			df[k] = df[k].astype(v)
	return df

def main():
	"""Main function"""
	args = argparser()
	sample, b, c = args.filename.rsplit('_', maxsplit=2)
	data = read_file_lines(args.filename)

	# Remove trailing whitespace from all cells
	data = data.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
	
	# filter out viruses
	data_filt = data[~data['Lineage'].str.contains('_Viruses')]
	data_filt = data_filt.replace(r'\n',' ', regex=True)

	# Fill in 'Lineage' if missing from unclassified annotation
	# Create condition as mask
	mask = ((data_filt['Lineage'].isna()) | (data_filt['Lineage'] == ' ') | (data_filt['Best_tax'] == 'root') | (data_filt['Tax_level'] == 'no rank') | (data_filt['Lineage'] == '-_cellular organisms') | (data_filt['Tax_level'] == 'cellular root'))
	# Apply mask to add 'unclassified' annotations to 'Lineage' column
	data_filt.loc[mask, 'Lineage'] = '-_Unclassified;d_Unclassified'
	# Add a 'd' for 'domain' to the domain fields (Archaea, Bacteria, and Eukaryota)
	data_filt['Lineage'] = data_filt['Lineage'].replace({';-_Archaea' : ';d_Archaea'}, regex=True)
	data_filt['Lineage'] = data_filt['Lineage'].replace({';-_Bacteria' : ';d_Bacteria'}, regex=True)
	data_filt['Lineage'] = data_filt['Lineage'].replace({';d_Bacteria candidate phyla' : ''}, regex=True)
	data_filt['Lineage'] = data_filt['Lineage'].replace({';d_Bacteria incertae sedis' : ''}, regex=True)
	data_filt['Lineage'] = data_filt['Lineage'].replace({';-_Eukaryota' : ';d_Eukaryota'}, regex=True)

	# Split taxonomic lineage into different fields and use taxonomic rank letters as column headers
	d = pd.DataFrame(data_filt['Lineage'].str.split(';').apply(lambda x:{i.split("_")[0] : i.split("_")[1] for i in x if i}).to_dict()).T
	df = pd.concat([data_filt, d], axis=1)
	df = df.drop(['Lineage', '-', 'Taxid'], axis=1) 

	# Strip trailing whitespace from 'd' column
	df['d'] = df['d'].str.rstrip()

	# Correct class-level taxonomy for haptophytes, chlorophytes, and Oomycota
	df.loc[df['o'] == 'Prymnesiales', 'c'] = 'Prymnesiophyceae'
	df.loc[df['o'] == 'Isochrysidales', 'c'] = 'Prymnesiophyceae'
	df.loc[df['o'] == 'Phaeocystales', 'c'] = 'Prymnesiophyceae'
	df.loc[df['o'] == 'Peronosporales', 'c'] = 'Peronosporales'
	df.loc[df['o'] == 'Saprolegniales', 'c'] = 'Saprolegniales'
	df.loc[df['o'] == 'Lagenidiales', 'c'] = 'Lagenidiales'
	df.loc[df['o'] == 'Pythiales', 'c'] = 'Pythiales'
	df.loc[df['f'] == 'Pycnococcaceae', 'c'] = 'Pseudoscourfieldiales'
	df.loc[df['o'] == 'Pseudoscourfieldiales', 'c'] = 'Pyramimonadophyceae'
	df.loc[df['o'] == 'Pyramimonadales', 'c'] = 'Pyramimonadophyceae'	

	# Export filtered df
	output = sample + "-mmseqs_uniref_lca_filt.tsv"
	df.to_csv(output, sep='\t', index=None)

if __name__ == '__main__':
	main()
