# -*- coding: utf-8 -*-
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

from postgap.DataModel import *
from pprint import pformat
from postgap.Finemap import *

def summarise(obj, **kwparams):
	
	if type(obj) == list:
		return summarise_list(obj, **kwparams)

	if type(obj) == GWAS_Cluster:
		return summarise_gwas_cluster(obj, **kwparams)

	if type(obj) == SNP:
		return summarise_snp(obj, **kwparams)
	
	if type(obj) == GWAS_SNP:
		return summarise_gwas_snp(obj, **kwparams)

	if type(obj) == GWAS_Association:
		return summarise_generic(obj, **kwparams)
	
	if type(obj) == Cisregulatory_Evidence:
		return summarise_generic(obj, **kwparams)

	if type(obj) == OneDConfigurationSample:
		return str(obj)

	if type(obj) == TwoDConfigurationSample:
		return str(obj)
	
	raise Exception("Don't know how to summarise a " + str(type(obj)))

def summarise_list(list_of_objects, leadstring = "", **kwparams):
	
	string = leadstring + " ┌ List:\n"
	
	for index, current_object in enumerate(list_of_objects):
		
		string += leadstring + " ├─ Item " + str(index) + ":\n"
		string += summarise(current_object, leadstring = leadstring + " │      ") + "\n"
	
	string += leadstring + " └─────"
	
	return string

def summarise_gwas_cluster(gwas_cluster, leadstring = ""):
	
	string =  "GWAS_Cluster\n"
	string += "════════════\n"
	string += "\n"
	
	string += "    Gwas snps:\n"
	string += "    ──────────\n"
	
	string += summarise(gwas_cluster.gwas_snps, leadstring = leadstring + "    ")
	
	string += "\n"
	string += "    LD snps:\n"
	string += "    ────────\n"
	
	string += summarise(gwas_cluster.ld_snps, leadstring = leadstring + "    ")

	return string

def summarise_snp(snp, leadstring = ""):
	return leadstring + pformat( snp )

def summarise_generic(obj, leadstring = ""):
	return leadstring + pformat( obj )

def summarise_gwas_snp(gwas_snp, leadstring = ""):
	
	string  = leadstring + "╔════════════\n"	
	string += leadstring + "║GWAS_SNP\n"	
	string += leadstring + "║ " + "SNP:     " + summarise(gwas_snp.snp, leadstring = "") + "\n"
	string += leadstring + "║ " + "pvalue:  " + str(gwas_snp.pvalue) + "\n"
	string += leadstring + "║ " + "z_score: " + str(gwas_snp.z_score) + "\n"
	string += leadstring + "║ " + "evidence:\n"
	string += summarise(gwas_snp.evidence, leadstring = leadstring + "║ " + "    ") + "\n"
	string += leadstring + "╚════════════"
	return string
