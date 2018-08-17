#! /usr/bin/env python

"""

Copyright [1999-2018] EMBL-European Bioinformatics Institute

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

import postgap.REST
from postgap.DataModel import *
import postgap.Globals
from postgap.Utils import *

GRCH37_ENSEMBL_REST_SERVER = "http://grch37.rest.ensembl.org"
GRCH38_ENSEMBL_REST_SERVER = "http://rest.ensembl.org"
known_genes = {}
known_snps = {}

def get_gene(gene_name, ENSEMBL_REST_SERVER = GRCH37_ENSEMBL_REST_SERVER):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	key = (gene_name, ENSEMBL_REST_SERVER)
	if key not in known_genes:
		if gene_name[:4] != 'ENSG':
			gene = fetch_gene(gene_name, ENSEMBL_REST_SERVER)
			if gene is None:
				known_genes[key] = None
			else:
				key2 = (gene.id, ENSEMBL_REST_SERVER)
				if key2 in known_genes:
					# Already found the same gene but under a different name (capitalisation changes etc)
					known_genes[key] = known_genes[key2]
				else:
					known_genes[key] = gene
					known_genes[key2] = gene
		else:
			known_genes[key] = fetch_gene_id(gene_name)
	return known_genes[key]

def fetch_gene(gene_name, ENSEMBL_REST_SERVER = GRCH37_ENSEMBL_REST_SERVER):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	server = ENSEMBL_REST_SERVER
	ext = "/lookup/symbol/%s/%s?content-type=application/json" % (postgap.Globals.SPECIES, gene_name)
	try:
		hash = postgap.REST.get(server, ext)
		return Gene(
			name = gene_name,
			id = hash['id'],
			chrom = hash['seq_region_name'],
			tss = int(hash['start']) if hash['strand'] > 0 else int(hash['end']),
			biotype = hash['biotype']
			)
	except:
		return None

def fetch_gene_id(gene_id, ENSEMBL_REST_SERVER = GRCH37_ENSEMBL_REST_SERVER):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	server = ENSEMBL_REST_SERVER
	ext = "/lookup/id/%s?content-type=application/json" % (gene_id)
	try:
		hash = postgap.REST.get(server, ext)
		return Gene(
			name = hash['display_name'],
			id = hash['id'],
			chrom = hash['seq_region_name'],
			tss = int(hash['start']) if hash['strand'] > 0 else int(hash['end']),
			biotype = hash['biotype']
			)
	except:
		return None


def get_ensembl_gene(ensembl_id, ENSEMBL_REST_SERVER = GRCH37_ENSEMBL_REST_SERVER):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	key = (ensembl_id, ENSEMBL_REST_SERVER)
	if key not in known_genes:
		known_genes[key] = fetch_gene_id(ensembl_id, ENSEMBL_REST_SERVER)
	return known_genes[key]

def get_snp_locations(rsIDs, ENSEMBL_REST_SERVER = GRCH37_ENSEMBL_REST_SERVER):
	"""

		Get SNP details from rsID
		* [ string ]
		Returntype: [ SNP ]

	"""
	if len(rsIDs) == 0:
		return []

	unknown_rsIDs = filter(lambda rsID: (rsID, ENSEMBL_REST_SERVER) not in known_snps, rsIDs)

	res = get_snp_locations_simple(unknown_rsIDs, ENSEMBL_REST_SERVER) 

	if len(res) == 0:
		if len(unknown_rsIDs) == 1:
			res = []
		else:
			res = get_snp_locations(unknown_rsIDs[:len(unknown_rsIDs)/2]) + get_snp_locations(unknown_rsIDs[len(unknown_rsIDs)/2:])
	
	return [known_snps[(rsID, ENSEMBL_REST_SERVER)] for rsID in rsIDs if (rsID, ENSEMBL_REST_SERVER) in known_snps]


def get_snp_locations_simple(rsIDs, ENSEMBL_REST_SERVER = GRCH37_ENSEMBL_REST_SERVER):
	"""

		Get SNP details from rsID
		* [ string ]
		Returntype: [ SNP ]

	"""

	assert postgap.Globals.SPECIES is not None, "postgap.Globals.SPECIES must be set or the url will be nonsense"
	
	server = ENSEMBL_REST_SERVER
	ext = "/variation/%s?content-type=application/json" % (postgap.Globals.SPECIES)
	hash = concatenate_hashes(postgap.REST.get(server, ext, data={'ids':chunk}) for chunk in chunks(rsIDs, 199))
	for record in hash.values():
		for synonym in record["synonyms"]:
			hash[synonym] = record

	'''
		Example response:
		{
			"rs56116432": {
				"source": "Variants (including SNPs and indels) imported from dbSNP",
				"mappings": [
					{
						"location": "9:133256042-133256042",
						"assembly_name": "GRCh38",
						"end": 133256042,
						"seq_region_name": "9",
						"strand": 1,
						"coord_system": "chromosome",
						"allele_string": "C/T",
						"start": 133256042
					},
				],
				"name": "rs56116432",
				"MAF": "0.00259585",
				"ambiguity": "Y",
				"var_class": "SNP",
				"synonyms": [],
				"evidence": [
					"Multiple_observations",
					"Frequency",
					"1000Genomes",
					"ESP",
					"ExAC"
				],
				"ancestral_allele": "C",
				"minor_allele": "T",
				"most_severe_consequence": "missense_variant"
			}
		}
	'''
	results = []
	for rsID in rsIDs:
		if rsID in hash:
			for mapping in hash[rsID]['mappings']:
				snp = SNP(
					rsID = rsID,
					chrom = mapping['seq_region_name'],
					pos = (int(mapping['start']) + int(mapping['end'])) / 2,
					approximated_zscore = None
				)
				results.append(snp)
				known_snps[(rsID, ENSEMBL_REST_SERVER)] = snp

	return results

