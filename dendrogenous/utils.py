#!/usr/bin/env python3

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
        reformatted_id = seq_record.id[:19]
    else:
        reformatted_id = seq_record.id

    reformatted_id.replace('|','_')
    reformatted_id.replace('/','_')
    reformatted_id.replace('\\','_')

    return reformatted_id


