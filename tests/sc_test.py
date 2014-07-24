import unittest

import modular_core.libfundamental as lfu
lfu.USING_GUI = False
import stringchemical as sc

import os, sys, pdb

#log = open(os.path.join(os.getcwd(), 'test_ensemble.log'), 'w')
#sys.stdout = log

class dummyTestCase(unittest.TestCase):
	"""Tests for `stringchemical module`."""
	#simple_mcfg = os.path.join(os.getcwd(), 
	#			'stringchemical_dep_mcfgs', 
	#			'MM_kinetics_boring.mcfg')

	def pause(self, *args, **kwargs):
		sys.stdout = sys.__stdout__
		pdb.set_trace()
		sys.stdout = log

	def test_can_make_ensemble(self):
		self.assertFalse(sc.simulate() == None)
		self.assertTrue(lfu.is_module_valid('chemical'))

if __name__ == '__main__':
    unittest.main()


