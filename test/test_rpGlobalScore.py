import unittest
import os
import sys

sys.path.insert(0, '..')

#NOTE: must copy the rpSBML.py file to the root directory
import rpSBML
import rpTool

class TestRPglobalscore(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.rpsbml = rpSBML.rpSBML('test', path=os.path.join('data', 'rpsbml.xml'))

    def test_calculateGlobalScore_rpsbml(self):
        global_score = calculateGlobalScore_rpsbml(self.rpsbml)
        self.assertAlmostEqual(global_score, 0.4208256320472526)
