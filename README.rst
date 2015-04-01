=============================== 
Dendrogenous
===============================

|Build Status|

|Coverage Status|

|Documentation Status|

Module designed to facilitate the rapid generation of phylogenies using
a MYSQL genome database

-  Free software: BSD license


- Sequences need to be in a MySQL database containing a table "cider" 
   which has sequences and protein_ID fields corresponding to blastdb



Features
--------

Take in a multi-protein fasta file and generates rapid phylogenies for
each sequence based on BLAST homology to genomes in the MYSQL database
specified in SETTINGS.json



.. |Build Status| image:: https://travis-ci.org/fmaguire/dendrogenous.png?branch=develop
   :target: https://travis-ci.org/fmaguire/dendrogenous
.. |Coverage Status| image:: https://coveralls.io/repos/fmaguire/dendrogenous/badge.png?branch=develop
   :target: https://coveralls.io/repos/fmaguire/dendrogenous/branch=develop
.. |Documentation Status| image:: https://readthedocs.org/projects/dendrogenous/badge/?version=latest
   :target: http://dendrogenous.readthedocs.org/en/latest 



