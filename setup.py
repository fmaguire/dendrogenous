#!/usr/bin/env python

from distutils.core import setup
import sys
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

setup(name='dendrogenous',
      version='1.0',
      description="Batch create phylogenies with a genome database and list of seeds",
      author="Finlay Maguire",
      author_email="root@finlaymagui.re",
      url="https://github.com/fmaguire/dendrogenous.git",
      tests_require=['pytest'],
      cmdclass = {'test': PyTest},
      packages=["dendrogenous"],)


