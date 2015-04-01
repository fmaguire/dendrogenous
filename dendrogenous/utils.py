#!/usr/bin/env python3

import os
import subprocess


def execute_cmd(cmd, input_str=None):
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

    else:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        (stdout, stderr) = proc.communicate()

    return stdout.decode()


def create_directory_structure(dir_structure):
    """
    Generate the required output directory structure in the output_dir
    """
    for folder in dir_structure.values():
        if not os.path.exists(folder):
            os.mkdir(folder)
