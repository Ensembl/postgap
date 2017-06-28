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
import tempfile
import numpy
from subprocess import Popen, PIPE
from postgap.DataModel import *
import postgap.Globals
import logging

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

	### use ld_vcf
	ld_comm = [
		"ld_vcf",
		"-f", chrom_file,
		"-r", "%s:%i-%i" % (snp.chrom, from_pos, to_pos),
		"-v", snp.rsID,
		"-w", str(window_len)
	]
	logger = logging.getLogger(__name__)
	logger.debug(" ".join(ld_comm))

	process = Popen(ld_comm, stdout=PIPE)
	(output, err) = process.communicate()
	if process.wait():
		raise Exception(err)

	### Read LD file
	ld_snps = []
	for line in output.split("\n"):

		items = line.split()
		if len(items) == 0:
			continue

		if float(items[6]) >= cutoff:
			ld_pos = 0
			ld_id  = ''

			if items[3] == snp.rsID:
				ld_pos = items[4]
				ld_id  = items[5]
			else:
				ld_pos = items[2]
				ld_id  = items[3]

			ld_snps.append(SNP(ld_id, snp.chrom, int(ld_pos)))

	if any(ld_snp.rsID == snp.rsID for ld_snp in ld_snps):
		return ld_snps
	else:
		return ld_snps + [snp]

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

	### get a list of rsIDs into a file
	rsID_file, rsID_file_name = tempfile.mkstemp()
	h = open(rsID_file_name, 'w')
	h.write("\n".join(str(snp.rsID) for snp in ld_snps))
	h.close()

	### use ld_vcf
	ld_comm = [
		"ld_vcf",
		"-f", chrom_file,
		"-r", "%s:%i-%i" % (gwas_snp.chrom, start, end),
		"-v", gwas_snp.rsID,
		"-n", rsID_file_name,
		"-w", str((end - start) + 1)
	]
	sys.stderr.write(" ".join(ld_comm) + "\n")

	process = Popen(ld_comm, stdout=PIPE)
	(output, err) = process.communicate()
	if process.wait():
		raise Exception(err)

	### Read LD file
	snp_hash = dict((snp.rsID, snp) for snp in ld_snps)
	r2_dict = {}

	for line in output.split("\n"):
		items = line.split()
		if len(items) == 0:
			continue

		ld_id = items[5] if items[3] == gwas_snp.rsID else items[3]
		r2_dict[snp_hash[ld_id]] = float(items[-3])

	os.remove(rsID_file_name)
	os.close(rsID_file)

	r2_dict[gwas_snp] = 1
	return r2_dict

def get_pairwise_ld(ld_snps, population='EUR'):
	"""

		For large numbers of SNPs, best to specify SNP region with chrom:to-from, e.g. 1:7654947-8155562

		Args:
		* [ SNP ], SNPs of interest
		* string, population name
		Returntype: [String (rsID)], Numpy.Array (2D)

	"""
	### Check inputs
	chrom = ld_snps[0].chrom
	assert all(x.chrom == chrom for x in ld_snps)
	positions = [x.pos for x in ld_snps]
	start = min(positions) - 10
	end = max(positions) + 10

	### Find the relevant BCF file
	chrom_file = os.path.join(postgap.Globals.DATABASES_DIR, '1000Genomes', population, 'ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf' % (chrom))
	if not os.path.isfile(chrom_file):
		return [], numpy.zeros((0,0))

	### get a list of rsIDs into a file
	rsID_file, rsID_file_name = tempfile.mkstemp()
	h = open(rsID_file_name, 'w')
	h.write("\n".join(str(snp.rsID) for snp in ld_snps))
	h.close()

	### use ld_vcf
	ld_comm = [
		"ld_vcf",
		"-f", chrom_file,
		"-r", "%s:%i-%i" % (chrom, start, end),
		"-n", rsID_file_name,
		"-w", str((end - start) + 1)
	]
	sys.stderr.write(" ".join(ld_comm) + "\n")

	process = Popen(ld_comm, stdout=PIPE, stderr=PIPE)
	(output, err) = process.communicate()
	if process.returncode:
		print process.returncode
		raise Exception(err)

	### First pass through file: collect rsIDs
	observed_snps = set()
	for line in output.split("\n"):
		items = line.split('\t')
		if len(items) < 9:
			continue
		observed_snps.add(items[3])
		observed_snps.add(items[5])

	### Second pass through file: store LD into matrix
	SNP_ids = [x.rsID for x in ld_snps if x.rsID in observed_snps]
	snp_order = dict((rsID, rank) for rank, rsID in enumerate(SNP_ids))
	r2_array = numpy.zeros((len(SNP_ids), len(SNP_ids)))
	for line in output.split("\n"):
		items = line.split('\t')
		if len(items) < 9:
			continue
		snp_1_rank = snp_order[items[3]]
		snp_2_rank = snp_order[items[5]]
		r2_array[snp_1_rank][snp_2_rank] = float(items[9])
		r2_array[snp_2_rank][snp_1_rank] = float(items[9])

	### Clean up
	os.remove(rsID_file_name)
	os.close(rsID_file)

	return SNP_ids, r2_array
