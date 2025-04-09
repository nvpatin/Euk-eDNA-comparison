#!/usr/bin/env python
# coding: utf-8

# Add sample name to CoverM tpm output files and merge by data set

import os as os
import glob as glob
import pandas as pd
import argparse

def argparser():
	"""Arguments"""
	epilog = """Example:
	$  parse_coverm_tpms_by_dataset.py {TPM_DIRECTORY} {kofams.txt file} """

	parser = argparse.ArgumentParser(description=__doc__, epilog=epilog,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('directory')
	parser.add_argument('kofams')
	args = parser.parse_args()
	return args

def parse_contig_coverages(coverage, kofams):
	tpm = pd.read_csv(coverage, sep='\t')
	sample, x = coverage.split('_O')
	tpm['Sample'] = sample
	tpm = tpm.rename(columns = {'Contig' : 'ORF'})
	tpm = tpm.rename(columns={ tpm.columns[1]: 'TPM' })
	tpm_kofams = tpm.merge(kofams, on='ORF', how='inner')
	return tpm_kofams

def main():
	"""Main function"""
	args = argparser()
	os.chdir(args.directory)

	coverms = []
	kofams = pd.read_csv(args.kofams, sep='\t', names=['ORF', 'KO'])

	for file in glob.glob("*coverm.tsv"):
		df = parse_contig_coverages(file, kofams)
		coverms.append(df)

	coverms_kofams = pd.concat(coverms, ignore_index=True)
	# coverms = coverms.sort_values(by=['Sample', 'ORF'])

	coverms_kofams = coverms_kofams.dropna(subset=['KO'])

	output = args.directory.rsplit("/", maxsplit=4)[1] + "_ORFs_kofams_tpms.tsv"
	coverms_kofams.to_csv(output, sep='\t', index=None)

if __name__ == '__main__':
	main()
