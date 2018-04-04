
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
from flask import Flask, request

import json
import sys
import collections
import re
import inspect
from pprint import pprint

import postgap
import postgap.EFO
import postgap.Globals
import postgap.Integration
from postgap.Utils import *

postgap_server = Flask(__name__)

@postgap_server.route("/")
def hello():
  return "Hello World!"


@postgap_server.route("/query", methods=['GET', 'POST'])
def query():

  options = get_options();

  # query on rsID, rsID+chr+pos or diseases/efos
  if options['rsID'] is not None:
    if len(options['tissues']) == 0:
      options['tissues'] = ["Whole_Blood"]

    if options['chr'] is not None:
      sys.stderr.write("Request type: position\n")
      snp = postgap.DataModel.SNP(rsID = options['rsID'], chrom = options['chr'], pos = int(options['pos']))
      res = postgap.Integration.ld_snps_to_genes([snp], options['tissues'])

    else:
      sys.stderr.write("Request type: rsID\n")
      res = postgap.Integration.rsIDs_to_genes(options['rsID'], options['tissues'])

  else:
    sys.stderr.write("Request type: EFOs\n")
    res = postgap.Integration.diseases_to_genes(options['diseases'], options['efos'], "CEPH", options['tissues'])

  return "\n".join(map(json.dumps, res))


def get_options():

  # set defaults
  options = {
    'tissues': [],
    'diseases': [],
    'efos': [],
    'rsID': None,
    'chr': None,
    'pos': None,
  }

  # get options from URL
  for key in options.keys():
    
    # allow for giving same param twice
    values = request.args.getlist(key)

    # but only for those params allowed to be a list
    if isinstance(options[key], list):
      options[key] += values
    elif len(values) > 0:
      options[key] = values[0]

  pprint(options, stream=sys.stderr)

  # modify
  if len(options['efos']) == 0:
    options['efos'] = filter(lambda X: X is not None, (postgap.EFO.suggest(disease) for disease in options['diseases']))

  # set defaults
  postgap.Globals.DATABASES_DIR = 'databases'
  postgap.Globals.SPECIES = 'Human'

  return options



if __name__ == "__main__":
  postgap_server.run()
