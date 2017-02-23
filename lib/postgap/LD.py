#! /usr/bin/env python

"""

Copyright [1999-2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License")
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

		 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

"""

	Please email comments or questions to the public Ensembl
	developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

	Questions may also be sent to the Ensembl help desk at
	<http://www.ensembl.org/Help/Contact>.

"""
import os
import os.path
import sys
import re
import subprocess
import tempfile
from postgap.DataModel import *
import postgap.Globals

def calculate_window(snp, window_len=500000, population='EUR', cutoff=0.7):
	"""

		Given a SNP id, calculate the pairwise LD between all SNPs within window_size base pairs.

		Args:
		* SNP
		* int, window width
		* string, population name
		* float, r2 cutoff
		Returntype: [ SNP ]

	"""
	### Define the necessary region.
	from_pos = snp.pos - (window_len / 2)
	to_pos = snp.pos + (window_len / 2)

	### Find the relevant 1000 genomes BCF
	chrom_file = os.path.join(postgap.Globals.DATABASES_DIR, '1000Genomes', population, "ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf" % (snp.chrom))
	if not os.path.isfile(chrom_file):
		return [snp]

	### Extract this region out from the 1000 genomes BCF
	region_file, region_file_name = tempfile.mkstemp()
	extract_region_comm = "bcftools view -r %s:%s-%s %s -O v -o %s" % (snp.chrom, from_pos, to_pos, chrom_file, region_file_name)
	sys.stderr.write(extract_region_comm + "\n")
	subprocess.call(extract_region_comm, shell=True)

	### Calculate the pairwise LD using plink2
	ld_file, ld_file_name = tempfile.mkstemp()
	plinkcomm = "plink --vcf %s --r2 --ld-snp %s --ld-window-r2 %f --inter-chr --keep-allele-order --out %s" % (region_file_name, snp.rsID, cutoff, ld_file_name)
	sys.stderr.write(plinkcomm + "\n")
	if subprocess.call(plinkcomm, shell=True) != 0:
		# Catching error reported by plink, generally that the GWAS SNP is not in the VCF as expected
		temp_file_names = [region_file_name] + [ ld_file_name + suffix for suffix in ["", ".ld", ".log", ".nosex"]]
		for temp_file_name in temp_file_names:
			if os.path.isfile(temp_file_name):
				os.remove(temp_file_name)
		return [snp]

	### Read LD file
	ld_snps = []
	for line in open(ld_file_name + '.ld'):
		items = line.split()
		if items[-1] == 'R2':
			continue
		ld = float(items[-1])
		ld_snps.append(SNP(items[-2].split(';')[0], items[-4], int(items[-3])))

	### Clean temp files
	temp_file_names = [region_file_name] + [ ld_file_name + suffix for suffix in ["", ".ld", ".log", ".nosex"]]
	map(os.remove, temp_file_names)
	map(os.close, [region_file, ld_file])

	if len(ld_snps) > 0:
		return ld_snps
	else:
		return [snp]

def get_lds_from_top_gwas(gwas_snp, ld_snps, population='EUR'):
	"""

		For large numbers of SNPs, best to specify SNP region with chrom:to-from, e.g. 1:7654947-8155562

		Args:
		* SNP
		* [ SNP ], SNPs of interest
		* string, population name
		Returntype: dict(SNP => float)

	"""
	### Check inputs
	if gwas_snp.rsID not in [x.rsID for x in ld_snps]:
		ld_snps.append(gwas_snp)
	if len(ld_snps) == 1:
		return dict([(gwas_snp, 1)])

	### Define the region of interest
	assert all(x.chrom == gwas_snp.chrom for x in ld_snps)
	positions = [x.pos for x in ld_snps]
	start = min(positions) - 10
	end = max(positions) + 10

	### Find the relevant BCF file
	chrom_file = os.path.join(postgap.Globals.DATABASES_DIR, '1000Genomes', population, 'ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf' % (gwas_snp.chrom))
	if not os.path.isfile(chrom_file):
		return dict((snp, 1) for snp in ld_snps)

	### Extract the required region from the BCF
	region_file, region_file_name = tempfile.mkstemp()
	extract_region_comm = "bcftools view -r %s:%s-%s %s -O z -o %s" % (gwas_snp.chrom, start, end, chrom_file, region_file_name)
	sys.stderr.write(extract_region_comm + "\n")
	subprocess.call(extract_region_comm, shell=True)

	### Extract the list of SNP ids from this region
	rsID_file, rsID_file_name = tempfile.mkstemp()
	h = open(rsID_file_name, 'w')
	h.write("\n".join(str(snp.rsID) for snp in ld_snps))
	h.close()
	snp_file, snp_file_name = tempfile.mkstemp()
	vcfcomm = "vcftools --gzvcf %s --snps %s --recode --maf 0.01 --min-alleles 2 --max-alleles 2 --out %s" % (region_file_name, rsID_file_name, snp_file_name)
	sys.stderr.write(vcfcomm + "\n")
	subprocess.call(vcfcomm, shell=True)

	### Calculate the pairwise LD using plink2
	ld_file, ld_file_name = tempfile.mkstemp()
	plinkcomm = "plink --vcf %s.recode.vcf --r2 --ld-snp %s --inter-chr --keep-allele-order --out %s" % (snp_file_name, gwas_snp.rsID, ld_file_name)
	sys.stderr.write(plinkcomm + "\n")
	if subprocess.call(plinkcomm, shell=True) != 0:
		temp_file_names = [region_file_name, rsID_file_name] + [snp_file_name + suffix for suffix in ["", ".recode.vcf", ".log"] ] + [ ld_file_name + suffix for suffix in ["", ".ld", ".log", ".nosex"]]
		for temp_file_name in temp_file_names:
			if os.path.isfile(temp_file_name):
				os.remove(temp_file_name)
		return dict([(gwas_snp, 1)])

	### Read LD file
	snp_hash = dict((snp.rsID, snp) for snp in ld_snps)
	r2_dict = {}
	for line in open(ld_file_name + '.ld'):
		items = line.split()
		if items[-1] == 'R2':
			continue
		ld = float(items[-1])
		ld_snp = snp_hash[items[-2].split(';')[0]]
		r2_dict[ld_snp] = ld

	### Clean temp files
	temp_file_names = [region_file_name, rsID_file_name] + [snp_file_name + suffix for suffix in ["", ".recode.vcf", ".log"] ] + [ ld_file_name + suffix for suffix in ["", ".ld", ".log", ".nosex"]]
	for temp_file_name in temp_file_names:
		if os.path.isfile(temp_file_name):
			os.remove(temp_file_name)
	map(os.close, [region_file, rsID_file, snp_file, ld_file])

	return r2_dict
