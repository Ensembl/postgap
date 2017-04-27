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
import re

import postgap.REST
from postgap.Globals import *

import sys
from pprint import pprint


def suggest(term):
	"""

		Find most appropriate EFO term for arbitrary string
		Arg:
		* string
		Returntype: string (EFO ID)

	"""
	server = 'http://www.ebi.ac.uk/spot/zooma/v2/api'
	url_term = re.sub(" ", "%20", term)
	ext = "/summaries?query=%s" % (url_term)
	result = postgap.REST.get(server, ext)
	'''

		Example search result:
		{
			id: "864A7A7335109CD59C2986398637CB519F21DB05",
			semanticTags: [
			"http://www.orpha.net/ORDO/Orphanet_308410"
			],
			annotationURIs: [
					"http://rdf.ebi.ac.uk/resource/zooma/efo/DE1004D39B114BDEEFB7A2EA93D881AD"
					],
			annotationSourceURIs: [
					"http://www.ebi.ac.uk/efo/efo.owl"
					],
			annotationSummaryTypeName: "disease; Orphanet_308410",
			annotatedPropertyType: "disease",
			annotatedPropertyValue: "Autism-epilepsy syndrome due to branched chain ketoacid dehydrogenase kinase deficiency",
			annotatedPropertyUri: "http://rdf.ebi.ac.uk/resource/zooma/FB9E6464F1830B24C08204157E7A088E",
			quality: 72.82359,
			uri: "http://rdf.ebi.ac.uk/resource/zooma/annotation_summary/864A7A7335109CD59C2986398637CB519F21DB05"
		},


	'''

	hits = filter(lambda X: len(X['semanticTags']) == 1, result)
	if len(hits):
		sorted_hits = sorted(hits, key = lambda X: X['quality'])
		selection = sorted_hits[-1]['semanticTags'][0]
		efo = re.sub('.*/', '\1', selection)
		if DEBUG:
			print "Suggested EFO %s" % efo
		return efo

	else:
		if DEBUG:
			print "Suggested EFO N/A"
		return None

def children(efo):
	"""

		Return list of children EFO IDs
		Arg:
		* string (EFO ID)
		Returntype: [ string ] (EFI IDs)

	"""
	server = 'http://www.ebi.ac.uk'
	page = 1
	res = [efo]

	while (True):
		ext = "/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252F%s/descendants?page%i&size=10" % (efo, page)
		hash = postgap.REST.get(server, ext)
		'''

			OWL Output format:
			{
				"_links" : {
					"first" : {
			"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_0001071/descendants?page=0&size=10"
					},
					"prev" : {
			"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_0001071/descendants?page=0&size=10"
					},
					"self" : {
			"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_0001071/descendants?page=1&size=10"
					},
					"last" : {
			"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_0001071/descendants?page=1&size=10"
					}
				},
				"_embedded" : {
					"terms" : [ {
						"iri" : "http://www.ebi.ac.uk/efo/EFO_1000333",
						"label" : "Lung Inflammatory Myofibroblastic Tumor",
						"description" : [ "An intermediate fibroblastic neoplasm arising from the lung. It is characterized by the presence of spindle-shaped fibroblasts and myofibroblasts, and a chronic inflammatory infiltrate composed of eosinophils, lymphocytes and plasma cells." ],
						"annotation" : {
							"NCI_Thesaurus_definition_citation" : [ "NCIt:C39740" ]
						},
						"synonyms" : null,
						"ontology_name" : "efo",
						"ontology_prefix" : "EFO",
						"ontology_iri" : "http://www.ebi.ac.uk/efo",
						"is_obsolete" : false,
						"is_defining_ontology" : true,
						"has_children" : false,
						"is_root" : false,
						"short_form" : "EFO_1000333",
						"obo_id" : "EFO:1000333",
						"_links" : {
							"self" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333"
							},
							"parents" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333/parents"
							},
							"ancestors" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333/ancestors"
							},
							"jstree" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333/jstree"
							},
							"graph" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333/graph"
							}
						}
					} ]
				},
				"page" : {
					"size" : 10,
					"totalElements" : 20,
					"totalPages" : 2,
					"number" : 1
				}
			}

		'''

		if '_embedded' in hash:
			res += [ result['short_form'] for result in hash['_embedded']['terms'] ]

		if page > int(hash['page']['totalPages']):
			break

		page += 1

	if DEBUG:
		print "EFO children: " + "\t".join(res)
	return res

