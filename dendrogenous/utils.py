#!/usr/bin/env python3

import re
import os
import shutil
import json
import subprocess
import io

def execute_cmd(cmd, input_str=None):
    """
    Execute a command using subprocess using a string as stdin
    and returns the stdout as a string if an input_str is provided
    Otherwise just executes cmd normally
    """
    if input_str is None:
        subprocess.call(cmd.split())

    else:
        proc = subprocess.Popen(cmd,
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        (stdout, stderr) = proc.communicate(input_str.encode())

        # stdout is returned as a binary string with newline characters
        stdout_list = [x for x in stdout.decode().split(os.linesep) \
                            if len(x) > 0]

        return stdout_list


def reformat_accession(seq_record):
    """
    Reformat accessions to be <=20 chars long and not contain any special chars
    """

    if len(seq_record.id) > 20:
        short_id = seq_record.id[:19]
    else:
        short_id = seq_record.id

    reformatted_id = re.sub('[|,/,\,.]', '_', short_id)
    return reformatted_id


def parse_settings(settings_file):
    """
    Parse the settings json and ensure required values are present
    """
    pass


def check_genomes(genome_dir, genome_list):
    """
    Check all the genomes specified in the settings file are in the genome_dir
    Change names and underscores as appropriate
    """
    pass


def create_directory_structure(output_dir):
    """
    Generate the required output directory structure in the output_dir
    """
    pass

