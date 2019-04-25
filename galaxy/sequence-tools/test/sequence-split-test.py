import unittest
import os
import sys
from split_sequences import test

class TestSplitSequenceMethods(unittest.TestCase):
  def test_test(self):
    self.assertEqual(test(), 42)
