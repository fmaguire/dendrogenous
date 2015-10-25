#!/usr/bin/env python
"""
Settings class that parses user supplied settings.json or command line arguments
"""
import os
import json
import logging
import io
import pymysql
from Bio import SeqIO

class Settings():
    """
    A class to handle all parasing of settings.json as well as checking
    values set and providing default settings.


    run_settings = dg_set.Settings('settings.json')
    """

    def __init__(self, settings_file):
        """
        Parse user settings.json, ensure required fields are provided,
        then add default args and export to class dict
        """


        required_settings = [('input_seqs', str),
                             ('genome_list', list),
                             ('genome_dir', str),
                             ('dbconfig', dict)]

        parsed_settings = self.__parse_settings(settings_file)

        user_input = self.__check_mandatory_fields(parsed_settings, required_settings)

        default_binary_dir = os.path.abspath(os.path.join('dendrogenous',
                                                          'dependencies'))
        default_settings = [('cores', 1),
                            ('output_dir', 'out'),
                            ['minimums', {'min_seqs': 3,
                                          'min_sites': 29}],
                            ['blast_settings', {'num_seqs': 1,
                                                'evalue': 1e-5}],
                            ['binary_paths', {'blastp':   os.path.join(default_binary_dir,
                                                                       'blastp'),
                                              'kalign':   os.path.join(default_binary_dir,
                                                                       'kalign'),
                                              'trimal':   os.path.join(default_binary_dir,
                                                                       'trimal'),
                                              'FastTree': os.path.join(default_binary_dir,
                                                                       'FastTree')}]]

        self.full_settings = self.__collate_default_and_user_settings(default_settings,
                                                                      user_input,
                                                                      required_settings)

        # needs a directory structure to place the logfile
        self.__create_directory_structure()
        self.logger_name = self.__get_logger()

    def __collate_default_and_user_settings(self, default_settings, user_input, required_settings):
        """
        Combine user defined settings with defaults overriding the default where
        appropriate.
        Raises a typeerror if the user defined type doesn't match the default type
        """
        # iterate over default settings populating full_settings
        # with the appropriate user defined setting if it exists or
        # the default value if not
        full_settings = {}
        for field, def_setting in default_settings:
            user_val_or_default = user_input.get(field, def_setting)
            if type(user_val_or_default) != type(def_setting):
                raise TypeError("{0} is incorrectly specified in settings.json should be {1}".format(field, type(def_setting)))
            else:
                full_settings.update({field: user_val_or_default})

        # add default fields from user input to full_settings
        for field, _ in required_settings:
            full_settings.update({field: user_input[field]})

        # make sure all user defined fields are real settings
        for user_field in user_input.keys():
            if user_field not in full_settings.keys():
                raise ValueError("{} setting specified in settings.json does not exist.".format(user_field))
        return full_settings





    def __create_directory_structure(self):
        """
        Generate the required output directory structure in the output_dir
        """
        for folder in self.dir_paths.values():
            if not os.path.exists(folder):
                # exist_ok=True gives mkdir -p type function
                os.makedirs(folder, exist_ok=True)

    def __get_logger(self):
        """
        Initialise and return logger name
        """

        logger_name = "dg_log"
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.DEBUG)

        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch_formatter = logging.Formatter("%(asctime)s - %(message)s")
        ch.setFormatter(ch_formatter)
        logger.addHandler(ch)

        fh = logging.FileHandler(os.path.join(self.dir_paths['run_data'],
                                              "run.log"))
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        fh.setFormatter(fh_formatter)
        logger.addHandler(fh)

        return logger_name


    def __parse_settings(self, settings_file):
        """
        Parse json file or stringIO formatted settings file
        into settings dictionary
        """
        if settings_file.__class__ is io.StringIO:
            settings = json.load(settings_file)

        else:
            if not os.path.exists(settings_file):
                raise ValueError("Settings file does not exist: {}".format(\
                                                                settings_file))
            with open(settings_file, 'r') as settings_fh:
                settings = json.load(settings_fh)

        return settings


    def __check_mandatory_fields(self, parsed_user_input, required_settings):
        """
        Ensure mandatory fields have been provided
        """
        for entry, value_type in required_settings:
            if entry not in parsed_user_input:
                raise ValueError("{} must be defined in settings.json".format(entry))
            elif type(parsed_user_input[entry]) is not value_type:
                raise TypeError("{0} is incorrectly specified in settings.json should be {1}".format(entry, value_type))

        return parsed_user_input

    @property
    def input_seqs(self):
        """
        Validates input seq filename and returns
        """
        putative_file = self.full_settings['input_seqs']
        if os.path.exists(putative_file):
            putative_file = os.path.abspath(putative_file)
        else:
            raise ValueError("specified input_seqs ({}) file does not exist".format(putative_file))
        try:
            with open(putative_file, 'rU') as test_fh:
                sample_seq = next(SeqIO.parse(test_fh, "fasta"))
        except:
            raise ValueError("input_seqs fasta file specified in settings.json {} could not be read".format(putative_file))

        return putative_file

    @property
    def genomes(self):
        """
        Returns list of genome paths and validate it
        """
        genome_list = self.full_settings['genome_list']
        genome_dir = os.path.abspath(self.full_settings['genome_dir'])
        if not os.path.exists(self.full_settings['genome_dir']):
            raise ValueError("{} genome dir path does not exist".format(genome_list))

        genome_paths = [os.path.join(genome_dir, genome_name) for genome_name in genome_list]

        genomes = []
        for genome in genome_paths:
            phr = os.path.exists(genome + '.phr')
            pin = os.path.exists(genome + '.pin')
            psq = os.path.exists(genome + '.psq')

            if not(phr and pin and psq):
                raise ValueError("{} genome blastdb files do not exist".format(genome))
            else:
                genomes.append(genome)

        return genomes

    @property
    def output_dir(self):
        """
        Returns output_dir
        """
        return os.path.abspath(self.full_settings['output_dir'])

    @property
    def dir_paths(self):
        """
        Returns the directory paths for all output
        """
        dir_paths = {"run_data": "0.run_data",
                     "ERROR": "ERROR",
                     "blast_hits": "1.blast_hits",
                     "blast_fail": os.path.join("1.blast_hits", "insufficient_hits"),
                     "alignment": "2.alignment",
                     "mask": "3.mask",
                     "mask_fail": os.path.join("3.mask", "insufficient_sites"),
                     "tree": "4.phylogeny",
                     "named": "5.name"}

        for key in dir_paths.keys():
           dir_paths[key] = os.path.join(self.output_dir, dir_paths[key])

        return dir_paths

    @property
    def dbconfig(self):
        """
        Returns dict of dbconfigs
        """
        self.__test_db_connection(self.full_settings['dbconfig'])
        return self.full_settings['dbconfig']

    def __test_db_connection(self, putative_dbconfig):
        """
        Test it is possible to connect to db
        """

        try:
            con = pymysql.connect(**putative_dbconfig)
            cursor = con.cursor()
            cursor.execute("SELECT VERSION()")
            results = cursor.fetchone()
            if not results:
                raise ValueError("Can't read mysql version in db: {}".format(putative_dbconfig))
        except:
            raise ValueError("Can't connect to mysql db: {}".format(putative_dbconfig))

    @property
    def blast_settings(self):
        """
        Returns blast_settings - a dict containing evalue to use and minimum
        number of sequences to take from each genome
        """
        return self.full_settings['blast_settings']

    @property
    def binary_paths(self):
        """
        Returns binary path dict
        """
        return self.full_settings['binary_paths']

    @property
    def minimums(self):
        """
        Returns minimums - a dict containing minimum number of sequences to align
        and minimum number of sites to build a tree from
        """
        return self.full_settings['minimums']
