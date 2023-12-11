#!/usr/bin/env python3

## Author: Yuqing Guo
## Version: 2.0
## 2023-10-11
## Last modify:2023-10-16
## Obtain reads SNPs

import os
import argparse
import pysam

import numpy as np
from tqdm import tqdm


class FilestreamVCF():
	"""
	Filestream VCF Object that reads line within a VCF, 
	filters for SNVs only, and returns the given site as a dictionary
	"""
	def __init__(self, vcf, contig = None):
		self.vcf_path = vcf
		self.stats()

	def __len__(self):
		return self.len

	def __iter__(self):
		return self.sites()

	def line_iterator(self):
		with open(self.vcf_path, 'r') as f:
			for line in f:
				if line.startswith("#") or "CHROM" in line:
					continue

				line = line.split("\t")[:6]

				if len(line[3]) != 1 or len(line[4]) != 1: #Filter for only SNPs
					continue

				yield line

	def stats(self):
		unique_chroms = []
		for j, line in enumerate(self.line_iterator()):
			if line[0] not in unique_chroms:
				unique_chroms.append(line[0])

		int_chroms, str_chroms = [], []
		for c in unique_chroms:
			try:
				int_chroms.append(int(c))
			except:
				str_chroms.append(c)

		self.len = j
		self.chromosomes = [str(i) for i in sorted(int_chroms)] + list(sorted(str_chroms))


	def __len__(self):
		return self.len

	def sites(self):
		for line in self.line_iterator():
			try:
				site = {"CHROM" : line[0], "POS" : int(line[1]), "REF" : line[3].upper(), "ALT" : line[4].upper(),"PATMAT" : line[5].upper()}
				yield site
			except (ValueError, IndexError):
				pass


def get_PIR(read, position):
	"""
	Obtains the correct position in read of a given variant position.
	Accounts for soft-clipping and indels.
	"""
	start = read.reference_start
	cigar = read.cigartuples
	if cigar[0][0] == 4:
		start -= cigar[0][1]
	pir_raw = int(position) - start - 1

	if len(cigar) == 1 and cigar[0][0] == 0:
		return pir_raw

	pir = pir_raw
	cigar_counter = 0
	for operator, length in cigar:
		if pir >= cigar_counter + length:
			if operator == 1:
				pir += length
			elif operator == 2:
				pir -= length
		elif pir < cigar_counter:
			return pir
		else:
			if operator == 1:
				pir += length 
			elif operator in [2, 4]:
				return -1
		if operator != 2:
			cigar_counter += length
	return pir

def matches_alt(read, alt):
	"""
	Checks if a read contains the alt at the position specified by the PIR tag

	"""
	return read.query_sequence[read.get_tag("PI")].upper() == alt.upper()


def site_reads(bam, site, filter_dups = False, mapped = False, match_alt = False):
	"""
	Generator that fetchs reads from a BAM that overlap with a given site.
	Applies basic filters (i.e. Is it mapped? Is it a seconday mapping?, etc) and returns passing reads
	If specified, will also filter duplicates and 
	reads with alleles matching that specified as an ALT in the vcf.

	"""
	for read in bam.fetch(site["CHROM"], start = site["POS"]-5000, stop = site["POS"]+5000):
	#	if (read.cigartuples is None) or read.is_unmapped or read.is_secondary or read.is_supplementary or (read.is_duplicate and filter_dups):
	#		continue

		start, end = read.reference_start, read.reference_end
		if start <= site["POS"] and end >= site["POS"]:
			pir = get_PIR(read, site["POS"])
			if pir >= 0:
				read.set_tag("PI", pir)
				if match_alt:
					if matches_alt(read, site["ALT"]):
						yield read
				else:
					yield read

def write_read_features(read, site, tsv):
	"""
	Writes a read with PATMAT features to a TSV.
	"""
	RID = read.query_name

	PIR = read.get_tag("PI")
	BQ_s = read.query_qualities

	if PIR >= len(BQ_s) or PIR < 0:
		return 
	BQ = read.query_qualities[PIR]
	alt = read.query_sequence[PIR]

	frag = read.template_length if read.is_proper_pair else 0
	
	if alt == site["ALT"]:
		PATMAT_info = site["PATMAT"].replace('\n', '')
	elif alt == site["REF"]:
		if site["PATMAT"].replace('\n', '') == "PAT":
			PATMAT_info = "MAT"
		else:
			PATMAT_info = "PAT"
	else:
		PATMAT_info = "Error"	

	MRBQ = np.mean(read.query_qualities)
	MQ = read.mapping_quality
	line = [str(l) for l in [site["CHROM"], site["POS"], site["REF"],site["ALT"],alt,RID,PATMAT_info]]
	line = "\t".join(line) + "\n"
	tsv.write(line)

def check_input_paths(inputs):
	"""
	Checks the inputs to insure they exist.
	"""
	if not os.path.exists(inputs.bam_path):
		raise FileNotFoundError("Input BAM: {} does not exit!".format(inputs.bam_path))

	if not os.path.exists(inputs.vcf_path):
		raise FileNotFoundError("Input VCF: {} does not exit!".format(inputs.vcf_path))

def compare_TSVs(inputs):
	"""
	Benchmarking tool to ensure reads selected match expectations.
	"""
	import pandas as pd
	mine = pd.read_csv(inputs.out_path, sep = '\t', header = None)
	benchmark = pd.read_csv(inputs.benchmark_tsv, sep = '\t', header = None)
	columns = ["chrom", "pos", "ref", "hetSNP","alt", "rid", "PATMAT_info"]
	mine.columns = columns
	benchmark.columns = columns

	mine.sort_values(by = ["rid"], inplace = True), benchmark.sort_values(by = ["rid"], inplace = True)

	mine["idx"] = mine["rid"] + ":" + mine["pos"].apply(str) + mine["pir"].apply(str)
	benchmark["idx"] = benchmark["rid"] + ":" + benchmark["pos"].apply(str) + benchmark["pir"].apply(str)

	added = mine.loc[~mine["idx"].isin(benchmark.idx)]
	missing = benchmark.loc[~benchmark["idx"].isin(mine.idx)]

	print("MISSING: {}".format(len(missing)))
	print(missing)
	print("ADDED: {}".format(len(added)))
	print(added)


def bam_read_filter(inputs):
	"""
	Identifies reads that overlap with a variant site 
	and writes them to a VCF.
	"""
	BAM = pysam.AlignmentFile(inputs.bam_path, "rb")
	VCF = FilestreamVCF(inputs.vcf_path)
	TSV = open(inputs.out_path, "w")
	print()
	for site in tqdm(VCF, desc = "Identifying Variants", total = len(VCF),  unit = " Loci", dynamic_ncols=True):
		for read in site_reads(BAM, site, filter_dups = inputs.filter_dups):
			write_read_features(read, site, TSV)

def main():
	"""
	Main
	"""
	parser = argparse.ArgumentParser(description='PATMAT SNPs Split')
	parser.add_argument("--out", dest="out_path", action = "store", help = "Output Name", required = True)
	parser.add_argument("--bam", dest="bam_path", action = "store", help = "Path to Input BAM", required = True)
	parser.add_argument("--vcf", dest="vcf_path", action = "store", help = "Path to Tumor VCF", required = True)

	parser.add_argument("--filter-dups", dest="filter_dups", action = "store", help = "Filter duplicate reads if specificied", required = False, default = False)
	parser.add_argument("--benchmark-tsv", dest="benchmark_tsv", action = "store", help = "Path to Tumor VCF", required = False)

	inputs = parser.parse_args()
	if not inputs.out_path.endswith(".tsv"):
		inputs.out_path += ".tsv"

	check_input_paths(inputs)

	bam_read_filter(inputs)

	if not inputs.benchmark_tsv is None:
		if not os.path.exists(inputs.benchmark_tsv):
			raise FileNotFoundError("Input Benchmark TSV: {} does not exit!".format(inputs.benchmark_tsv))
		else:
			print("\nComparing identified reads with benchmark TSV\n")
			compare_TSVs(inputs)
	print("Done.")

if __name__ == '__main__':
	main()
