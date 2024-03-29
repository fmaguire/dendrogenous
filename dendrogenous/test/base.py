#!/usr/bin/env python
"""
Base classes for unit testing
"""
import unittest
import hashlib
import subprocess
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


    def execute_cmd(self, cmd, input_str=None):
        """
        Execute a command using subprocess using a string as stdin
        and returns the stdout as a string if an input_str is provided
        Otherwise just executes cmd normally
        """
        if input_str is None:
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            (stdout, stderr) = proc.communicate()


        else:
            osnull = open(os.devnull, 'w')
            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                    stderr=osnull)
            (stdout, stderr) = proc.communicate(input_str.encode())

            osnull.close()

       # stdout is returned as a binary string with newline characters
        stdout_list = [x for x in stdout.decode().split(os.linesep) \
                          if len(x) > 0]
        return stdout_list

    def file_hash(self, file_name):
        """
        Calculate md5sum of a given file
        """
        hasher = hashlib.md5()
        with open(file_name, 'rb') as fh:
            buf = fh.read()
            hasher.update(buf)
            return hasher.hexdigest()



    def touch(self, filename):
        """
        UNIX touch facility to make a file
        """
        with open(filename, 'a'):
            os.utime(filename, None)
