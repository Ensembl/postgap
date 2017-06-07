#!/usr/bin/env python

import unittest

class Blabladiblupp(unittest.TestCase):
 
	def setUp(self):
		pass
 
	def test_foooooo(self):
		print "Hello world!"
		self.assertEqual(2, 2)
 
if __name__ == '__main__':
	unittest.main()
