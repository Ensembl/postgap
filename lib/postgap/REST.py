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
import sys
import requests
import json
import time
import logging
import httplib

import signal

class timeout_handler:
	
	url = 'initme'
	
	# Register an handler for the timeout
	def handler_very_short(self, signum, frame):
		logging.info("Url is being a bit sluggish " + self.url)
		signal.signal(signal.SIGALRM, self.handler_long)
		# Set second timeout
		signal.alarm(4)

	# Register an handler for the timeout
	def handler_short(self, signum, frame):
		logging.warning("Waiting for " + self.url)
		signal.signal(signal.SIGALRM, self.handler_long)
		# Set second timeout
		signal.alarm(60)

	def handler_long(self, signum, frame):
		logging.error("Killing request for url: "  + self.url)
		raise requests.exceptions.ReadTimeout("Killed request for url because of timeout: %s" % self.url)

def post_request_with_timeout(url, headers, data, timeout):
	
	th = timeout_handler()
	th.url = url
	
	# Register the signal function handler
	signal.signal(signal.SIGALRM, th.handler_short)

	# Set first timeout
	signal.alarm(1)
	
	r = requests.post(url, headers = headers, data = data, timeout=timeout)
	
	# Cancel alarm
	signal.alarm(0)

	return r

def get_request_with_timeout(url, headers, timeout):
	
	th = timeout_handler()
	th.url = url
	
	# Register the signal function handler
	signal.signal(signal.SIGALRM, th.handler_short)

	# Set first timeout
	signal.alarm(2)
	
	r = requests.get(url, headers = headers, timeout=timeout)
	
	# Cancel alarm
	signal.alarm(0)

	return r

def get(server, ext, data=None):
	"""
		Args:
		* String (server name)
		* String (extension string)
		Return type: JSON object

	"""
	maximum_retries = 10

	
	for retries in range(maximum_retries):
		
		logging.debug("REST JSON Query: %s%s" % (server, ext))
		start_time = time.time()

		try:
			if data is None:
				headers = { "Content-Type" : "application/json" }
				r = get_request_with_timeout(server.encode('ascii', 'xmlcharrefreplace')+ext.encode('ascii', 'xmlcharrefreplace'), headers = headers, timeout=200)
			else:
				headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
				r = post_request_with_timeout(server.encode('ascii', 'xmlcharrefreplace')+ext.encode('ascii', 'xmlcharrefreplace'), headers = headers, data = json.dumps(data), timeout=200)
                except requests.exceptions.ReadTimeout:
			continue
		except requests.exceptions.ConnectionError:
			# A timeout can creep up as a connection error, so catching this as well.
			# requests.exceptions.ConnectionError: HTTPConnectionPool(host='grch37.rest.ensembl.org', port=80): Read timed out.
			logging.error("Got a requests.exceptions.ConnectionError when querying %s%s" % (server, ext) )
			continue
		except requests.exceptions.ChunkedEncodingError:
			# Happens every now and then when the eqtl server feels a bit 
			# stressed.
			logging.error("Got a requests.exceptions.ChunkedEncodingError when querying %s%s" % (server, ext) )
			continue


		if not r.ok:
			
			logging.debug("Something went wrong code: %s" % (r.status_code))

			http_response_code = None
			
			try:
				http_response_code = httplib.responses[r.status_code]
			except KeyError:
				http_response_code = r.status_code
				error_message = "Unknown status code %s" % (r.status_code)
				logging.critical(error_message)
			
			#if (http_response_code == ""):
				#error_message = "Unknown status code %s" % (r.status_code)
				#logging.critical(error_message)
				#raise RuntimeError(error_message)
			
			logging.error("Failed to get proper response to query %s%s" % (server, ext) )
			logging.error("With headers:" + repr(headers))
			if data is not None:
				logging.error("With data:" + repr(data))
			
			response_as_string = None
			
			try:
				response_as_string = json.dumps(r.json())
			except ValueError:
				response_as_string = "<Error when stringifying>" + repr(r) + "</Error when stringifying>"
			
			logging.error("Error code: %s (%s) %s" % (http_response_code, r.status_code, response_as_string ) )

			#if retries == 5:
			#	logging.critical("Giving up.")
			#	r.raise_for_status()

			if 'Retry-After' in r.headers:
				if r.status_code == 429:
					logging.error("Got error 429 'Too Many Requests'" )
				
				logging.error("Will try again in %s seconds." % r.headers['Retry-After'])
				time.sleep(int(r.headers['Retry-After']))

			elif r.status_code == 502:
				logging.error("Got error 502 'Bad Gateway'. Will try again in %s seconds." % 2)
				time.sleep(2) # Sleep while server cools down

			elif r.status_code == requests.codes.forbidden:
				logging.error("Got 'forbidden' error: Will try again in %s seconds." % 600)
				time.sleep(600) # Sleep 10 minutes while server calms down
			elif r.status_code == 104 \
				or r.status_code == requests.codes.gateway_timeout \
				or r.status_code == requests.codes.request_timeout:

				logging.error("Got 'timeout error': Will try again in %s seconds." % 60)
				time.sleep(60) # Sleep 1 minute while server cools down
			elif r.status_code == 400:
				# Check for errors that aren't actually errors
				url = server + ext
				if "/eqtl/" in url:
					logging.info("Error is expected behaviour by the eqtl server and will be passed on.")
					raise EQTL400error(r)

				if "/lookup/symbol" in url or '/lookup/id' in url:
					logging.info("Error is expected behaviour by the Ensembl gene lookup and will be passed on.")
					raise GENE400error(r)
				
				if "/variation/" in url:
					
					response = r.json()
					
					# Happens like this:
					#
					# Failed to get proper response to query http://grch37.rest.ensembl.org/variation/homo_sapiens/rs24449894?content-type=application/json
					# With headers:{'Content-Type': 'application/json'}
					# Error code: Bad Request (400) {"error": "rs24449894 not found for homo_sapiens"}
					#
					if " not found for" in response["error"]:
						logging.info("Error is expected behaviour by the variation endpoint and will be passed on.")
						raise Variation400error(r)
				
				# requests.exceptions.HTTPError: 400 Client Error: Bad Request for url: http://grch37.rest.ensembl.org/overlap/region/Human/5:117435127-119583975?feature=gene;content-type=application/json
				if retries < 5:
					logging.error("Will try again in %s seconds." % 2)
					time.sleep(2)
					continue
				logging.error("Will try again in %s seconds." % 60)
				time.sleep(60) # Sleep 1 minute while server cools down
			else:
				logging.error("Got status code %s (%s)." % (r.status_code, http_response_code))
				r.raise_for_status()
			continue

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

class unhandled_rest_exception(Exception):
    def __init__(self, request):

        # No message to pass, so setting to ""
        super(Exception, self).__init__("")

        # Now for your custom code...
        self.request = request
        self.response = request.json()

class EQTL400error(unhandled_rest_exception):
	pass

class GENE400error(unhandled_rest_exception):
	pass

class Variation400error(unhandled_rest_exception):
	pass




