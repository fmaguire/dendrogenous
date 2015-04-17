#!/usr/bin/env python3

import os
import subprocess
from Bio import SeqIO


def execute_cmd(cmd, input_str=None, output_str=False, debug=False):
    """
    Execute a command using subprocess using a string as stdin
    and returns the stdout as a string if an input_str is provided
    Otherwise just executes cmd normally
    """
    if input_str:
        osnull = open(os.devnull, 'w')
        proc = subprocess.Popen(cmd,
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=osnull)
        (stdout, stderr) = proc.communicate(input_str.encode())
        osnull.close()
        return stdout.decode()

    elif output_str and not input_str:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        (stdout, stderr) = proc.communicate()
        return stdout.decode()

    elif debug:
        subprocess.call(cmd, shell=True)

    else:
        with open(os.devnull, 'w') as null:
            subprocess.call(cmd, shell=True,
                             stdout=null, stderr=null)

    return None

def parse_seqs(sequence_file):
    with open(sequence_file, 'rU') as seq_fh:
        sequences = [seq for seq in SeqIO.parse(seq_fh, 'fasta')]
    return sequences

def check_already_run(settings, input_seqs):
    pass
    #settings.output_dir

