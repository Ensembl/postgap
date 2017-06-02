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
import logging
import httplib

from postgap.Globals import *

def get(server, ext, data=None):
	"""
		Args:
		* String (server name)
		* String (extension string)
		Return type: JSON object

	"""
	#logger = logging.getLogger(__name__)
	logger = logging.getLogger(__file__)
	
	for retries in range(3):
		if DEBUG:
			logging.debug("REST JSON Query: %s%s" % (server, ext))
			start_time = time.time()

		try:
			if data is None:
				headers = { "Content-Type" : "application/json" }
				r = requests.get(server.encode('ascii', 'xmlcharrefreplace')+ext.encode('ascii', 'xmlcharrefreplace'), headers = headers, timeout=200)
			else:
				headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
				r = requests.post(server.encode('ascii', 'xmlcharrefreplace')+ext.encode('ascii', 'xmlcharrefreplace'), headers = headers, data = json.dumps(data), timeout=200)
                except requests.exceptions.ReadTimeout:
			continue
		except requests.exceptions.ConnectionError:
			# A timeout can creep up as a connection error, so catching this as well.
			# requests.exceptions.ConnectionError: HTTPConnectionPool(host='grch37.rest.ensembl.org', port=80): Read timed out.
			continue


		if not r.ok:
			
			http_response_code = httplib.responses[r.status_code]
			
			if (http_response_code == ""):
				error_message = "Unknown status code %s" % (r.status_code)
				logging.critical(error_message)
				raise RuntimeError(error_message)
			
			logging.error("Failed to get proper response to query %s%s" % (server, ext) )
			logging.error("With headers:" + repr(headers))
			if data is not None:
				logging.error("With data:" + repr(data))
			
			logging.error("Error code: %s (%s) %s" % (http_response_code, r.status_code, json.dumps(r.json()) ) )

			if retries == 2:
				logging.critical("Giving up.")
				r.raise_for_status()

			if 'Retry-After' in r.headers:
				logging.error("Will try again in %s seconds." % r.headers['Retry-After'])
				time.sleep(int(r.headers['Retry-After']))

			elif r.status_code == requests.codes.forbidden:
				logging.error("Will try again in %s seconds." % 600)
				time.sleep(600) # Sleep 10 minutes while server calms down
			elif r.status_code == 104 or r.status_code == requests.codes.gateway_timeout or r.status_code == requests.codes.request_timeout:
				logging.error("Will try again in %s seconds." % 60)
				time.sleep(60) # Sleep 1 minute while server cools down
			else:
				logging.error("Got status code %s (%s)." % (r.status_code, http_response_code))
				r.raise_for_status()
			continue

		if DEBUG:
			logging.debug("Time: %f" % (time.time() - start_time))

		try:
			return r.json()
		except:
			error_message = "Failed to get proper response to query %s%s" % (server, ext) 
			logging.critical(error_message)
			raise requests.HTTPError(error_message)

	# Failed too many times
	error_message = "Failed too many times to get a proper response for query %s%s !" % (server, ext)
	logging.critical(error_message)
	raise requests.exceptions.ConnectionError(error_message)

