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


def suggest(term):
	"""

		Find most appropriate EFO term for arbitrary string
		Arg:
		* string
		Returntype: string (EFO ID)

	"""
	server = 'http://www.ebi.ac.uk/spot/zooma/v2/api'
	url_term = re.sub(" ", "%20", term)
	ext = "/services/annotate?propertyValue=%s&filter=required:[none],ontologies:[efo]" % (url_term)
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

[

    {
        "uri": null,
        "annotatedProperty": {
            "uri": null,
            "propertyType": null,
            "propertyValue": "diabetes mellitus"
        },
        "_links": {
            "olslinks": [
                {
                    "href": "http://www.ebi.ac.uk/ols/api/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0000400",
                    "semanticTag": "http://www.ebi.ac.uk/efo/EFO_0000400"
                }
            ]
        },
        "semanticTags": [
            "http://www.ebi.ac.uk/efo/EFO_0000400"
        ],
        "replacedBy": [ ],
        "replaces": [ ],
        "derivedFrom": {
            "uri": "http://rdf.ebi.ac.uk/resource/zooma/annotation_summary/OLS",
            "annotatedProperty": {
                "uri": null,
                "propertyType": null,
                "propertyValue": "diabetes mellitus"
            },
            "_links": {
                "olslinks": [
                    {
                        "href": "http://www.ebi.ac.uk/efo/EFO_0000400",
                        "semanticTag": "http://www.ebi.ac.uk/efo/EFO_0000400"
                    }
                ]
            },
            "semanticTags": [
                "http://www.ebi.ac.uk/efo/EFO_0000400"
            ],
            "replacedBy": [ ],
            "replaces": [ ],
            "annotatedBiologicalEntities": [ ],
            "provenance": {
                "source": {
                    "type": "ONTOLOGY",
                    "name": "http://www.ebi.ac.uk/efo",
                    "uri": "http://www.ebi.ac.uk/efo"
                },
                "evidence": "COMPUTED_FROM_ONTOLOGY",
                "accuracy": null,
                "generator": "http://www.ebi.ac.uk/efo",
                "generatedDate": null,
                "annotator": null,
                "annotationDate": null
            }
        },
        "confidence": "GOOD",
        "annotatedBiologicalEntities": [ ],
        "provenance": {
            "source": {
                "type": "DATABASE",
                "name": "zooma",
                "uri": "http://www.ebi.ac.uk/spot/zooma"
            },
            "evidence": "COMPUTED_FROM_TEXT_MATCH",
            "accuracy": null,
            "generator": "ZOOMA",
            "generatedDate": 1493109533692,
            "annotator": "ZOOMA",
            "annotationDate": 1493109533692
        }
    }

]

	'''

	hits = filter(lambda X: len(X['semanticTags']) == 1, result)
	efo = hits[0]['semanticTags'][0]
	return efo

def children(efo):
	"""

	Return list of children EFO IDs
	Arg:
	* string (EFO ID)
	Returntype: [ string ] (EFI IDs)

	"""

	import logging
	logger = logging.getLogger(__name__)
	
	server = 'http://www.ebi.ac.uk'

	import urllib
	double_quoted_iri = urllib.quote_plus(urllib.quote_plus(efo))

	# E.g.: http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400
	ext = "/ols/api/ontologies/efo/terms/" + double_quoted_iri
	hash = postgap.REST.get(server, ext)
	'''
        {

            "iri": "http://www.ebi.ac.uk/efo/EFO_0000400",
            "label": "diabetes mellitus",
            "description": [
                "A heterogeneous group of disorders characterized by HYPERGLYCEMIA and GLUCOSE INTOLERANCE.",
                "A metabolic disorder characterized by abnormally high blood sugar levels due to diminished production of insulin or insulin resistance/desensitization."
            ],
            "annotation": {
                "DOID_definition_citation": [
                    "DOID:9351"
                ],
                "ICD9_definition_citation": [
                    "ICD9:250"
                ],
                "MSH_definition_citation": [
                    "MSH:D003920"
                ],
                "NCI_Thesaurus_definition_citation": [
                    "NCIt:C2985"
                ],
                "OMIM_definition_citation": [
                    "OMIM:612227"
                ],
                "SNOMEDCT_definition_citation": [
                    "SNOMEDCT:73211009"
                ],
                "bioportal_provenance": [
                    "Diabetes mellitus, NOS[accessedResource: SNOMEDCT:73211009][accessDate: 05-04-2011]",
                    "Diabetes mellitus (disorder)[accessedResource: DOID:9351][accessDate: 05-04-2011]",
                    "DM - Diabetes mellitus[accessedResource: SNOMEDCT:73211009][accessDate: 05-04-2011]",
                    "A heterogeneous group of disorders characterized by HYPERGLYCEMIA and GLUCOSE INTOLERANCE.[accessedResource: MSH:D003920][accessDate: 05-04-2011]",
                    "Diabetes NOS[accessedResource: DOID:9351][accessDate: 05-04-2011]",
                    "Diabetes[accessedResource: DOID:9351][accessDate: 05-04-2011]",
                    "A metabolic disorder characterized by abnormally high blood sugar levels due to diminished production of insulin or insulin resistance/desensitization.[accessedResource: NCIt:C2985][accessDate: 05-04-2011]"
                ],
                "gwas_trait": [
                    "true"
                ],
                "term editor": [
                    "James Malone"
                ]
            },
            "synonyms": [
                "Diabetes mellitus (disorder)",
                "Diabetes",
                "Diabetes mellitus, NOS",
                "DM - Diabetes mellitus",
                "Diabetes NOS"
            ],
            "ontology_name": "efo",
            "ontology_prefix": "EFO",
            "ontology_iri": "http://www.ebi.ac.uk/efo",
            "is_obsolete": false,
            "term_replaced_by": null,
            "is_defining_ontology": true,
            "has_children": true,
            "is_root": false,
            "short_form": "EFO_0000400",
            "obo_id": "EFO:0000400",
            "in_subset": null,
            "obo_definition_citation": null,
            "obo_xref": null,
            "obo_synonym": null,
            "_links": {
                "self": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400"
                },
                "parents": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/parents"
                },
                "ancestors": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/ancestors"
                },
                "hierarchicalParents": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/hierarchicalParents"
                },
                "hierarchicalAncestors": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/hierarchicalAncestors"
                },
                "jstree": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/jstree"
                },
                "children": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/children"
                },
                "descendants": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/descendants"
                },
                "hierarchicalChildren": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/hierarchicalChildren"
                },
                "hierarchicalDescendants": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/hierarchicalDescendants"
                },
                "graph": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/graph"
                }
            }

        }
	'''

	has_children = hash['has_children']
	if (has_children == False):
		return []

	descendants_href = hash['_links']['descendants']['href']
	hash = postgap.REST.get(descendants_href, '')

	terms = hash['_embedded']['terms']
	logger.debug("The Ontology Lookup Service returned the descendants for " + efo)
	for term in terms:
		logger.debug(" - " + term['short_form']);

	result = [ term['iri'] for term in terms ]

	return result

known_terms = dict()

def term(efo):
	"""

	Return term associated to EFO ID
	Arg:
	* string (EFO ID)
	Returntype: string (term)

	"""
	if efo not in known_terms:
		known_terms[efo] = lookup_term(efo)
	return known_terms[efo]

def lookup_term(efo):
	"""

	Return term associated to EFO ID
	Arg:
	* string (EFO ID)
	Returntype: string (term)

	"""
	import logging
	logger = logging.getLogger(__name__)
	
	server = 'http://www.ebi.ac.uk'

	import urllib
	double_quoted_iri = urllib.quote_plus(urllib.quote_plus(query_iris_for_efo_short_form(efo)[0]))

	# E.g.: http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400
	ext = "/ols/api/ontologies/efo/terms/" + double_quoted_iri
	hash = postgap.REST.get(server, ext)
	'''
        {

            "iri": "http://www.ebi.ac.uk/efo/EFO_0000400",
            "label": "diabetes mellitus",
            "description": [
                "A heterogeneous group of disorders characterized by HYPERGLYCEMIA and GLUCOSE INTOLERANCE.",
                "A metabolic disorder characterized by abnormally high blood sugar levels due to diminished production of insulin or insulin resistance/desensitization."
            ],
            "annotation": {
                "DOID_definition_citation": [
                    "DOID:9351"
                ],
                "ICD9_definition_citation": [
                    "ICD9:250"
                ],
                "MSH_definition_citation": [
                    "MSH:D003920"
                ],
                "NCI_Thesaurus_definition_citation": [
                    "NCIt:C2985"
                ],
                "OMIM_definition_citation": [
                    "OMIM:612227"
                ],
                "SNOMEDCT_definition_citation": [
                    "SNOMEDCT:73211009"
                ],
                "bioportal_provenance": [
                    "Diabetes mellitus, NOS[accessedResource: SNOMEDCT:73211009][accessDate: 05-04-2011]",
                    "Diabetes mellitus (disorder)[accessedResource: DOID:9351][accessDate: 05-04-2011]",
                    "DM - Diabetes mellitus[accessedResource: SNOMEDCT:73211009][accessDate: 05-04-2011]",
                    "A heterogeneous group of disorders characterized by HYPERGLYCEMIA and GLUCOSE INTOLERANCE.[accessedResource: MSH:D003920][accessDate: 05-04-2011]",
                    "Diabetes NOS[accessedResource: DOID:9351][accessDate: 05-04-2011]",
                    "Diabetes[accessedResource: DOID:9351][accessDate: 05-04-2011]",
                    "A metabolic disorder characterized by abnormally high blood sugar levels due to diminished production of insulin or insulin resistance/desensitization.[accessedResource: NCIt:C2985][accessDate: 05-04-2011]"
                ],
                "gwas_trait": [
                    "true"
                ],
                "term editor": [
                    "James Malone"
                ]
            },
            "synonyms": [
                "Diabetes mellitus (disorder)",
                "Diabetes",
                "Diabetes mellitus, NOS",
                "DM - Diabetes mellitus",
                "Diabetes NOS"
            ],
            "ontology_name": "efo",
            "ontology_prefix": "EFO",
            "ontology_iri": "http://www.ebi.ac.uk/efo",
            "is_obsolete": false,
            "term_replaced_by": null,
            "is_defining_ontology": true,
            "has_children": true,
            "is_root": false,
            "short_form": "EFO_0000400",
            "obo_id": "EFO:0000400",
            "in_subset": null,
            "obo_definition_citation": null,
            "obo_xref": null,
            "obo_synonym": null,
            "_links": {
                "self": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400"
                },
                "parents": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/parents"
                },
                "ancestors": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/ancestors"
                },
                "hierarchicalParents": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/hierarchicalParents"
                },
                "hierarchicalAncestors": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/hierarchicalAncestors"
                },
                "jstree": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/jstree"
                },
                "children": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/children"
                },
                "descendants": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/descendants"
                },
                "hierarchicalChildren": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/hierarchicalChildren"
                },
                "hierarchicalDescendants": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/hierarchicalDescendants"
                },
                "graph": {
                    "href": "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000400/graph"
                }
            }

        }
	'''
	return hash['label']

def query_iris_for_efo_short_form_list(efo_short_form_list):
	iri_list = []
	for efo_short_form in efo_short_form_list:
		iri_list += query_iris_for_efo_short_form(efo_short_form)
		
	return iri_list

def query_iris_for_efo_short_form(efo_short_form):
	server = 'http://www.ebi.ac.uk'
	ext = "/ols/api/ontologies/efo/terms?short_form=%s" % (efo_short_form)
	result = postgap.REST.get(server, ext)
	terms = result['_embedded']['terms']
	iri_terms = [ term['iri'] for term in terms ]
	return iri_terms
