import dendrogenous as dg
import dendrogenous.settings

import io
import pickle
import os
import shutil
import pytest

from socket import gethostname

from dendrogenous.test.base import BaseTestCase

class TestSettings(BaseTestCase):
    """
    Unit tests for the core dendrogeous init and reform methods
    """

    def test_clean_init(self):
        """
        test clean initialisation of settings class - pretty much integration test for class
        also checks that outputdir default value "out" is replace with "test_out"
        also checks io string parsing -- bad I know
        """
        needed_settings = io.StringIO('{"input_seqs": "a.fas", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {"a": 4}}')
        test_settings = dg.settings.Settings(needed_settings)
        self.assertEqual(test_settings.minimums, {'min_seqs': 3, 'min_sites': 29})
        self.assertEqual(test_settings.output_dir, os.path.abspath("test_out"))
        self.assertEqual(test_settings.blast_settings, {'num_seqs': 1, 'evalue': 1e-5})
        self.assertEqual(test_settings.binary_paths['blastp'], os.path.join(os.path.abspath(os.path.join('dendrogenous',
                                                                                                         'dependencies')),
                                                                            'blastp'))

    def test_parse_settings_json_file(self):
        """
        Make sure json file is correctly parsed
        """
        needed_settings = os.path.join(self.test_resources, 'test_settings.json')
        test_settings = dg.settings.Settings(needed_settings)
        self.assertEqual(test_settings.minimums, {'min_seqs': 3, 'min_sites': 29})
        self.assertEqual(test_settings.output_dir, os.path.abspath("test_out"))
        self.assertEqual(test_settings.blast_settings, {'num_seqs': 1, 'evalue': 1e-5})
        self.assertEqual(test_settings.binary_paths['blastp'], os.path.join(os.path.abspath(os.path.join('dendrogenous',
                                                                                                         'dependencies')),
                                                                            'blastp'))

    def test__check_mandatory_fields_valueerror_if_default_missing(self):
        """
        ensure Settings raises valueError if mandatory setting absent
        """
        needed_settings = io.StringIO('{}')
        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"output_dir": "blah", "genome_list": ["a", "b", "c"]}')
        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"output_dir": "blah", "genome_dir": "blah"}')
        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"genome_list": ["a", "b", "c"], "genome_dir": "blah"}')
        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)

    def test__check_mandatory_fields_typerror_if_type_wrong(self):
        """
        ensure Settings raises valueError if mandatory setting absent
        """
        needed_settings = io.StringIO('{"input_seqs": "", "genome_list": "", "genome_dir": "foo", "dbconfig": {"a": 4}}')
        with self.assertRaises(TypeError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"input_seqs": "", "genome_list": ["a", "b", "c"], "genome_dir": {}, "dbconfig": {"a": 4}}')
        with self.assertRaises(TypeError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"input_seqs": "", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": []}')
        with self.assertRaises(TypeError):
            test_settings = dg.settings.Settings(needed_settings)


    def test__collate_default_and_user_settings_raises_TypeError_if_new_type_wrong(self):
        """
        Ensure if user defined setting overriding default is of different type to default
        a valuerror is raised
        """
        needed_settings = io.StringIO('{"input_seqs": "", "output_dir": [], "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {"a": 4}}')
        with self.assertRaises(TypeError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"input_seqs": "", "output_dir": "test_out", "blast_settings": [], "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {"a": 4}}')
        with self.assertRaises(TypeError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"input_seqs": "", "output_dir": "test_out", "binary_paths": "", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {"a": 4}}')
        with self.assertRaises(TypeError):
            test_settings = dg.settings.Settings(needed_settings)

    def test__collate_default_and_user_settings_raises_ValueError_if_wrong_user_input_field_identified(self):
        """
        If user defined fields aren't settings names raise ValueError
        """
        needed_settings = io.StringIO('{"input_seqs": "", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {"a": 4}}, "foo": ""')
        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"input_seqs": "", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {"a": 4}}, "bar": {}')
        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)

        needed_settings = io.StringIO('{"input_seqs": "", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {"a": 4}}, "baz": []')
        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)

    def test_genomes_all_exist(self):
        """
        make sure settings correctly find all genomes and return real paths
        """

        genome_list = '["Escherichia_coli_O157_H7_str._Sakai.fas", "Escherichia_coli_IAI39.fas", "Nanoarchaeum_equitans_Kin4-M.fas"]'

        settings_string = '{{"input_seqs": "", "output_dir": "test_out", "genome_list": {0}, "genome_dir": "{1}", "dbconfig": {{"a": 4}}}}'.format(genome_list, self.test_resources)


        needed_settings = io.StringIO(settings_string)

        test_settings = dg.settings.Settings(needed_settings)

        self.assertEqual(len(test_settings.genomes), 3)
        self.assertTrue(os.path.exists(test_settings.genomes[1]))


    def test_genomes_non_existent_valueError(self):
        """
        Make sure self.genomes raises valueerror for non-existent genome
        """

        genome_list = '["foo.fas", "Escherichia_coli_IAI39.fas", "Nanoarchaeum_equitans_Kin4-M.fas"]'

        settings_string = '{{"output_dir": "test_out", "genome_list": {0}, "genome_dir": "{1}", "dbconfig": {{"a": 4}}}}'.format(genome_list, self.test_resources)


        needed_settings = io.StringIO(settings_string)

        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)
            foo = test_settings.genomes

    def test_genomes_dir_non_existent_valueError(self):
        """
        Ensure ValueError is raised when genome directory does not exist
        """
        genome_list = '["Escherichia_coli_O157_H7_str._Sakai.fas", "Escherichia_coli_IAI39.fas", "Nanoarchaeum_equitans_Kin4-M.fas"]'

        settings_string = '{{"output_dir": "test_out", "genome_list": {0}, "genome_dir": "{1}", "dbconfig": {{"a": 4}}}}'.format(genome_list, "foo")

        needed_settings = io.StringIO(settings_string)

        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(needed_settings)
            foo = test_settings.genomes


    def test_dir_paths(self):
        """
        Would be pretty much testing a literal but leaving stub in case later modifications call for it
        """
        pass

    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test_db_config_success(self):
        """
        Ensure dbconfig method works (including test_db_connection) using secret settings
        """
        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret_fh:
            secret_config = pickle.load(secret_fh)

        setting_str = '{{"input_seqs": "", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {0}}}'.format(str(secret_config).replace('\'', '\"'))

        test_settings = dg.settings.Settings(io.StringIO(setting_str))
        foo = test_settings.dbconfig
        self.assertEqual(foo, secret_config)

    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test_db_config_fail(self):
        """
        Ensure dbconfig method fails correctly ValueError (including test_db_connection)
        """
        needed_settings = io.StringIO('{"input_seqs": "", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {}}')
        test_settings = dg.settings.Settings(needed_settings)
        with self.assertRaises(ValueError):
            foo = test_settings.dbconfig

    def test_blast_settings(self):
        """
        Would be pretty much testing a literal but leaving stub in case later modifications call for it
        """
        pass

    def test_minimums(self):
        """
        Would be pretty much testing a literal but leaving stub in case later modifications call for it
        """
        pass

    def test_input_seqs_validation(self):
        """
        Check input seqs validation works properly if file exists and is fasta
        """
        setting_str = '{{"input_seqs": "{}", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {{}}}}'.format(str(os.path.join(self.test_resources, "Nanoarchaeum_equitans_Kin4-M.fas")).replace("\'", ""))

        test_settings = dg.settings.Settings(io.StringIO(setting_str))
        input_seqs = test_settings.input_seqs
        self.assertEqual(os.path.abspath(os.path.join(self.test_resources, "Nanoarchaeum_equitans_Kin4-M.fas")),
                         input_seqs)


    def test_input_seqs_validation_raises_ValueError_if_fasta_wrong(self):
        """
        Check input seqs validation raises a valueerror is input seq file exists but is not fasta
        """
        setting_str = '{{"input_seqs": "{}", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {{}}}}'.format(str(os.path.join(self.test_resources, "broken_blastp_output.xml")).replace("\'", ""))

        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(io.StringIO(setting_str))
            input_seqs = test_settings.input_seqs

    def test_input_seqs_validation_raises_ValueError_if_fasta_does_not_exist(self):
        """
        Check input seqs validation raises a valueerror if input seq file does not exist
        """
        setting_str = '{{"input_seqs": "a.fas", "output_dir": "test_out", "genome_list": ["a", "b", "c"], "genome_dir": "foo", "dbconfig": {{}}}}'
        with self.assertRaises(ValueError):
            test_settings = dg.settings.Settings(io.StringIO(setting_str))
            input_seqs = test_settings.input_seqs


    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')
