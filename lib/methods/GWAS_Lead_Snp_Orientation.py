#!/usr/bin/env python

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
import logging
import collections

'''
	Data structure backing the Risk_Allele_Orientation class
'''
Risk_Allele_Orientation_Prototype = collections.namedtuple(
	'Risk_Allele_Orientation_Prototype', 
	[
		# The accession of the SNP
		'rs_id',
		# The base in the allele reported by GWAS
		'base_at_snp_in_gwas',
		# The base in the reference
		'base_at_snp_in_reference',
		# The outcome of: "base_at_snp_in_gwas == base_at_snp_in_reference"
		'risk_allele_present_in_reference',
		# The urls from which the data was obtained.
		'gwas_url',
		'ensembl_url',
	]
)

class Risk_Allele_Orientation(Risk_Allele_Orientation_Prototype):
	'''
		Risk_Allele_Orientation is used to store the orientation of a risk 
		allele from gwas.
		
		The interesting bit it in the field "risk_allele_present_in_reference",
		the other fields are used for generating meaningful messages when
		exceptions are thrown.
	'''
	def __str__(self):
		
		msg = ""
		
		msg +=  "    rs_id = " + self.rs_id + "\n"
		msg +=  "        base_at_snp_in_gwas      = " + str(self.base_at_snp_in_gwas) + "\n"
		msg +=  "        base_at_snp_in_reference = " + str(self.base_at_snp_in_reference) + "\n"
		msg +=  "        \n"
		msg +=  "        Result: risk_allele_present_in_reference = " + str(self.risk_allele_present_in_reference) + "\n"
		msg +=  "        \n"
		msg +=  "        See:\n"
		msg +=  "             gwas_url    = " + self.gwas_url    + "\n"
		msg +=  "             ensembl_url = " + self.ensembl_url + "\n"
		
		return msg

class _exception_with_risk_allele_orientations(Exception):
	
	def __init__(self, risk_allele_orientations, *args):
		self.risk_allele_orientations = risk_allele_orientations

	def description(self):
		return "Overwrite this with a generic description of the exception."
	
	def __str__(self):
		
		risk_allele_orientations = self.risk_allele_orientations
		
		msg = self.description() + "\n"
		msg += "The gwas association has %i risk alleles:\n" % len(risk_allele_orientations)
		
		for risk_allele_orientation in risk_allele_orientations:
			msg +=  str(risk_allele_orientation) + "\n"

		return msg

class some_alleles_present_others_not_exception(_exception_with_risk_allele_orientations):
	def description(self):
		return "Some alleles are present in the reference, others aren't."


class none_of_the_risk_alleles_is_a_substitution_exception(_exception_with_risk_allele_orientations):
	def description(self):
		return "None of the risk alleles is a substitution."


class variant_mapping_is_ambiguous_exception(Exception):
	def description(self):
		return "The variant was mapped to more than one location in the reference. (This never seems to happen.)"

class gwas_data_integrity_exception(Exception):
	pass

def gwas_risk_alleles_present_in_reference(riskAlleles):
	
	risk_allele_orientations = compute_risk_allele_orientations(riskAlleles)
	risk_allele_orientation_consensus = compute_risk_allele_orientation_consensus(risk_allele_orientations)
	
	if risk_allele_orientation_consensus == Risk_Allele_Orientation_Consensus.no_gwas_risk_allele_present_in_reference:
		return False
	
	if risk_allele_orientation_consensus == Risk_Allele_Orientation_Consensus.all_gwas_risk_alleles_present_in_reference:
		return True

	# This should never happen.
	raise Exception

class Risk_Allele_Orientation_Consensus:
	'''
		This class is used like an enum to report the outcome of analysing the
		orientation of the gwas risk alleles belonging to the same association.
		
		It handles the two expected cases that are straightforward to report.
		
		Exceptions to these are handled by raising exceptions.
	'''
	all_gwas_risk_alleles_present_in_reference, \
	no_gwas_risk_allele_present_in_reference    \
	= range(2)

def compute_risk_allele_orientation_consensus(risk_allele_orientations):
	'''
		Figures out, whether the risk alleles reported by gwas are present in 
		the reference.
		
		For any gwas association there may be multiple risk alleles, so this 
		also deals with special cases arising from that.
		
		Return type: A value from the "enum" Risk_Allele_Orientation_Consensus
		
		Returns True, if all risk alleles from gwas are present in the 
		reference, False, if none of the risk alleles are present in the
		reference.
		
		Exceptions:
		
		none_of_the_risk_alleles_is_a_substitution_exception:
		
		if a risk allele in gwas doesn't have a base assigned to it, this risk
		allele is skipped. In that case gwas puts a question mark where the 
		base should be. If none of the gwas associations have a base assigned,
		a "none_of_the_risk_alleles_is_a_substitution_exception" is raised.
		
		some_alleles_present_others_not_exception:
		
		It is possible to have multiple risk alleles in an association. Some 
		may be present in the reference, others not. This is an inconsistent 
		case. If it arises, a "some_alleles_present_others_not_exception" is
		raised.
	'''
	
	# Flags for keeping track, if any of these outcomes it true for all risk 
	# alleles of the gwas association.
	#
	all_gwas_risk_alleles_present_in_reference = True
	no_gwas_risk_allele_present_in_reference   = True
	none_of_the_risk_alleles_is_a_substitution = True

	for risk_allele_orientation in risk_allele_orientations:
		
		no_base_in_gwas    = risk_allele_orientation.base_at_snp_in_gwas == '?'
		no_base_in_ensembl = risk_allele_orientation.base_at_snp_in_reference is None
		
		no_risk_allele_is_a_substitution = no_base_in_gwas and no_base_in_ensembl
		
		none_of_the_risk_alleles_is_a_substitution = none_of_the_risk_alleles_is_a_substitution and no_risk_allele_is_a_substitution
		
		all_gwas_risk_alleles_present_in_reference = all_gwas_risk_alleles_present_in_reference and risk_allele_orientation.risk_allele_present_in_reference
		no_gwas_risk_allele_present_in_reference   = no_gwas_risk_allele_present_in_reference   and not(risk_allele_orientation.risk_allele_present_in_reference)
	
	if none_of_the_risk_alleles_is_a_substitution:
		raise none_of_the_risk_alleles_is_a_substitution_exception(risk_allele_orientations)
	
	if all_gwas_risk_alleles_present_in_reference:
		return Risk_Allele_Orientation_Consensus.all_gwas_risk_alleles_present_in_reference
	
	if no_gwas_risk_allele_present_in_reference:
		return Risk_Allele_Orientation_Consensus.no_gwas_risk_allele_present_in_reference
	
	raise some_alleles_present_others_not_exception(risk_allele_orientations)

def compute_risk_allele_orientations(riskAlleles):
	'''
		For the array of risk alleles returned by gwas, this checks, whether 
		the risk allele is present in the reference or not.
		
		It returns an array of Risk_Allele_Orientation objects.
	'''
	rs_id = assert_risk_alleles_are_from_same_snp(riskAlleles)
	
	base_at_snp_in_reference, ensembl_source_url = fetch_base_at_snp_in_reference(rs_id)
	
	risk_allele_orientations = []

	for riskAllele in riskAlleles:
		
		riskAlleleName = riskAllele["riskAlleleName"]
		(rs_id, base_at_snp_in_gwas) = riskAlleleName.split("-")
		
		current_gwas_risk_allele_present_in_reference = base_at_snp_in_gwas == base_at_snp_in_reference
		
		risk_allele_orientation = Risk_Allele_Orientation(
			
			rs_id = rs_id,
			
			gwas_url    = riskAllele["_links"]["self"]["href"],
			ensembl_url = ensembl_source_url,
			
			risk_allele_present_in_reference = current_gwas_risk_allele_present_in_reference,
			base_at_snp_in_gwas              = base_at_snp_in_gwas,
			base_at_snp_in_reference         = base_at_snp_in_reference,
		)
		
		risk_allele_orientations.append(risk_allele_orientation)
		
	return risk_allele_orientations

def assert_risk_alleles_are_from_same_snp(riskAlleles):
	'''
		This is an integrity check.
		
		All risk alleles from a gwas association should be from the same snp.
	'''
	seen_rs_id = ""
	
	for riskAllele in riskAlleles:
		
		riskAlleleName = riskAllele["riskAlleleName"]
		(raw_rs_id, nucleotide_in_risk_allele) = riskAlleleName.split("-")
		
		# E.g.: On
		# http://wwwdev.ebi.ac.uk/gwas/beta/rest/api/singleNucleotidePolymorphisms/6167/riskAlleles
		#
		# The alleles are from the snp rs4420638. But one of them is: "
		# rs4420638 -?", so the rs_id will appear as  "rs4420638 " with an extra space a the end.
		#
		# The extra whitespace can trip up this assertion.
		#
		rs_id = raw_rs_id.strip()
		
		if seen_rs_id == "":
			seen_rs_id = rs_id
		
		if seen_rs_id != rs_id:
			import json
			raise gwas_data_integrity_exception("More than one snp reported for the same gwas association! One snp is: " + seen_rs_id + " another is: " + rs_id + " \n\n" + json.dumps(riskAlleles))
	
	return seen_rs_id

def fetch_base_at_snp_in_reference(rs_id):
	'''
		Fetches the base for this snp in the reference genome.
		
		If there is more than one, this will raise a 
		variant_mapping_is_ambiguous_exception.
		
		A snp can have more than one base, if it was mapped more than once.
	'''
	bases_at_snp_in_reference, source_url = fetch_bases_at_snp_in_reference(rs_id)
	
	if len(bases_at_snp_in_reference) > 1:
		raise variant_mapping_is_ambiguous_exception()
	
	# It can happen that all of them were insertions.
	if len(bases_at_snp_in_reference) == 0:
		return None, source_url
	
	return bases_at_snp_in_reference[0], source_url

def fetch_bases_at_snp_in_reference(rs_id):
	'''
		Fetches the bases for this snp in the reference genome.
		
		If it the snp was mapped to more than one location and these have 
		a different nucleotide, it will return all nucleotides that occurred 
		once in an array.
		
		If the variant is an insertion, it is skipped.
	'''
	mappings, source_url = fetch_variant_mappings(rs_id)	
	bases_seen_in_reference = dict()
	
	for mapping in mappings:

		allele_string = mapping["allele_string"];
		observed_bases = allele_string.split("/")
		
		# The first component is the reference allele according to:
		#
		# https://genomes-ebi.slack.com/archives/C0JUVJV6W/p1498561584548384
		#
		base_in_reference = observed_bases[0]
		
		# Not using this.
		bases_observed_in_other_populations = observed_bases[1:]
		
		if base_in_reference == "-":
			logging.debug("Skipping variant, because it is not a substitution in the reference, but an insertion: " + base_in_reference);
			continue

		bases_seen_in_reference[base_in_reference] = 1
		
	return bases_seen_in_reference.keys(), source_url

def fetch_variant_mappings(rs_id):
	'''
		Queries the rest server to find all mappings of a given SNP.
		
		Returns the hash returned by the rest server.
	'''
	variant_mappings_rest_call_url = "http://grch37.rest.ensembl.org/variation/homo_sapiens/%s?content-type=application/json" % rs_id
	
	import postgap.REST
	hash = postgap.REST.get(variant_mappings_rest_call_url, ext="")
	'''
		{
		"source": "Variants (including SNPs and indels) imported from dbSNP",
		"mappings": [
			{
				"location": "9:136131429-136131429",
				"assembly_name": "GRCh37",
				"end": 136131429,
				"seq_region_name": "9",
				"strand": 1,
				"coord_system": "chromosome",
				"allele_string": "C/T",
				"start": 136131429
			}
		],
		"name": "rs56116432",
		"MAF": 0.00259585,
		"ambiguity": "Y",
		"var_class": "SNP",
		"synonyms": [ ],
		"evidence": [
			"Frequency",
			"1000Genomes",
			"ESP",
			"ExAC"
		],
		"ancestral_allele": "C",
		"minor_allele": "T",
		"most_severe_consequence": "non_coding_transcript_exon_variant"

		}
	'''
	assert hash["name"] == rs_id
	mappings = hash["mappings"]
	
	return mappings, variant_mappings_rest_call_url
