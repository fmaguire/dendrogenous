#!/usr/bin/env python
"""
Base classes for unit testing
"""
import unittest
import os

class BaseTestCase(unittest.TestCase):
    """
    Base test superclass for attributes and methods shared by all
    other unittests
    """
    test_resources = os.path.join('dendrogenous', 'test', 'resources')
    binary_path = os.path.join('dendrogenous', 'dependencies')

    def parse_file(self, fpath):
        """
        Utility function to parse a file to a list
        """
        with open(fpath, 'r') as fh:
            parsed_file = [x.rstrip('\n') for x in fh.readlines()]
        return parsed_file
