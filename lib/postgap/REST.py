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
import sys
import requests
import json
import time

from postgap.Globals import *

def get(server, ext, data=None):
	"""
		Args:
		* String (server name)
		* String (extension string)
		Return type: JSON object

	"""
	retries = 0

	while True:
		if data is None:
			headers = { "Content-Type" : "application/json" }
			r = requests.get(str(server)+str(ext), headers = headers)
		else:
			headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
			r = requests.post(str(server)+str(ext), headers = headers, data = json.dumps(data))

		if DEBUG:
			sys.stderr.write("REST JSON Query: %s%s\n" % (server, ext))

		if not r.ok:
			sys.stderr.write("Failed to get proper response to query %s%s\n" % (server, ext) )
			sys.stderr.write("With headers:\n" + repr(headers) + "\n")
			if data is not None:
				sys.stderr.write("With data:\n" + repr(data) + "\n")
			if 'Retry-After' in r.headers:
				time.sleep(int(r.headers['Retry-After']))
				retries += 1
				continue
			r.raise_for_status()
			sys.exit()

		try:
			return r.json()
		except:
			sys.stderr.write("Failed to get proper response to query %s%s\n" % (server, ext) )
			raise

	# Failed too many times
	sys.exit()

