#!/usr/bin/env python
# coding: utf-8

# Extract taxa counts assigned by mmseqs using the UniRef database and LCA
# Filter out uncertain bacterial taxa

import pandas as pd
import argparse
import numpy as np

def argparser():
	"""Arguments"""
	epilog = """Example:
	$  extract_and_count_mmseqs_UniRef_lca.py {SAMPLE}_LCA_tophit_report """

	parser = argparse.ArgumentParser(description=__doc__, epilog=epilog,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('filename')
	args = parser.parse_args()
	return args

def read_file_lines(filename):
	lines_with_fields = []
	with open(filename, 'r') as file:
		for line in file:
    		# Split the line into fields using the tab separator
			fields = line.rstrip().split('\t')
			# Select some elements from the fields (replace 0 and 2 with the indices you want)
			selected_elements = [fields[1], fields[4], fields[7], fields[8]]
			# Append the selected elements to the list
			lines_with_fields.append(selected_elements)
		df = pd.DataFrame(lines_with_fields, columns=['count','avg_seq_id','species','taxonomy'])
		dtype_dict = {'count':'int', 'avg_seq_id':'float', 'species':'str', 'taxonomy':'str'}
		for k, v in dtype_dict.items():
			df[k] = df[k].astype(v)
	return df

def main():
	"""Main function"""
	args = argparser()
	sample, b, c, d = args.filename.rsplit('_', maxsplit=3)
	data = read_file_lines(args.filename)
	# filter out hits with 'root' or 'bacterium' as species
	data = data[data['species'] != 'root']
	data = data[data['species'] != 'bacterium']
	# filter out viruses
	data = data[~data['taxonomy'].str.contains('d_Viruses')]
	# filter out sequences only classified from metagenomes
	data = data[~data['species'].str.contains('metagenome')]
	# filter for average alignment > 50%
	data_filt = data[data['avg_seq_id'] > 0.5]
	# split taxonomy into separate fields
	d = pd.DataFrame(data_filt['taxonomy'].str.split(';').apply(lambda x:{i.split("_")[0] : i.split("_")[1] for i in x if i}).to_dict()).T
	df = pd.concat([data_filt, d], axis=1)
	df = df.drop(['species', 'taxonomy', '-'], axis=1)
	# Export filtered df
	output = sample + "-mmseqs_uniref_lca_counts.tsv"
	output_filt = sample + "-mmseqs_uniref_lca_counts-filt-50percent.tsv"
	#data.to_csv(output, sep='\t', index=None)
	df.to_csv(output_filt, sep='\t', index=None)

if __name__ == '__main__':
	main()
