#!/usr/bin/env python
"""
Settings class that parses user supplied settings.json or command line arguments
"""
import os
import json
import io
import pymysql

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

        required_settings = ['genome_list', 'genome_dir', 'dbconfig']

        parsed_settings = self.__parse_settings(settings_file)

        user_input = self.__check_mandatory_fields(required_settings,
                                                 parsed_settings)

        default_settings = (('output_dir', 'out'),
                            ('minimums', {'min_seqs': 3,
                                          'min_sites': 30}),
                            ('blast_settings', {'num_seqs': 1,
                                                'evalue': 1e-5})
                            ('binary_paths', {'blastp': 'blastp',
                                              'kalign': 'kalign',
                                              'trimal': 'trimal',
                                              'FastTree': 'FastTree'}),
                            ('log_level', 2))

        self.full_settings = {}
        for field, setting in default_settings:
            self.full_settings[field] = user_input.get(field, setting)

        self.logger = __get_logger(self.log_level)


    def __get_logger(self, log_level):
        """
        Initialise and return logger
        """
        logger = logging.getLogger("dg_log")
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(os.path.join(self.output_dir,
                                                "run_data", "run.log"))
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.ERROR)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        logger.addHandler(fh)
        logger.addHandler(ch)
        return logger


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


    def __check_mandatory_fields(self, required_settings, parsed_user_input):
        """
        Ensure mandatory fields have been provided
        """
        for entry in required_settings:
            if entry not in parsed_user_input:
                raise ValueError("{} must be defined in settings.json".format(entry))

    @property
    def genomes(self):
        """
        Returns list of genome paths and validate it
        """
        genome_list = self.full_settings['genome_list']
        genome_dir = os.path.abspath(self.full_settings['genome_dir'])
        if not os.path.exists(genome_list):
            raise ValueError("{} genome list file does not exist".format(genome_list))

        if not os.path.exists(self.full-settings['genome_dir']):
            raise ValueError("{} genome dir path does not exist".format(genome_list))

        with open(genome_list, 'r') as gl_fh:
            genomes = [os.path.join(genome_dir, genome.strip()) for genome in gl_fh.readlines()]

        for genome in genomes:
            phr = os.path.exists(genome + '.phr')
            pin = os.path.exists(genome + '.pin')
            psq = os.path.exists(genome + '.psq')
            index = os.path.exists(genome + '.index')

            if not(phr and pin and psq and index):
                raise ValueError("{} genome blastdb files do not exist".format(genome))

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
        dir_paths = {"run_data": "run_data",
                     "input_seqs": "0.sequences",
                     "blast_hits": "1.blast_hits",
                     "blast_fail": os.path.join("1.blast_hits", "insufficient_hits"),
                     "alignment": "2.alignment",
                     "mask": "3.mask",
                     "mask_fail": os.path.join("3.mask", "insufficient_sites"),
                     "tree": "4.phylogeny",
                     "name": "5.name"}

        for key in dir_paths.keys():
           dir_paths[key] = os.path.join(self.output_dir, dir_paths[key])

        return dir_paths

    @property
    def dbconfig(self):
        """
        Returns dict of dbconfigs
        """
        self.__test_db_connection(self.full_settings['db_config'])
        return self.full_settings['db_config']

    def __test_db_connection(self, db_config):
        """
        Test it is possible to connect to db
        """
        try:
            config = self.settings['db_config']
            con = pymysql.connect(**config)
            cursor = conn.cursor()
            cursor.execute("SELECT VERSION()")
            results = cursor.fetchone()
            if not results:
                raise ValueError("Can't read mysql version in db: {}".format(config))
        except:
            raise ValueError("Can't connect to mysql db: {}".format(config))

    @property
    def blast_settings(self):
        """
        Returns blast_settings - a dict containing evalue to use and minimum
        number of sequences to take from each genome
        """
        return self.full_settings['blast_settings']

    @property
    def minimums(self):
        """
        Returns minimums - a dict containing minimum number of sequences to align
        and minimum number of sites to build a tree from
        """
        return self.full_settings['minimums']

