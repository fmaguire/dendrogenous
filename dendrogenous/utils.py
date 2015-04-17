#!/usr/bin/env python3

import os
import subprocess
import re
from Bio import SeqIO
import glob
import collections

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
    """
    Parse input sequence file while reformatting accessions
    """
    with open(sequence_file, 'rU') as seq_fh:
        sequences = [reformat_accession(seq) for seq in SeqIO.parse(seq_fh, 'fasta')]
    return sequences

def reformat_accession(seq_record):
    """
    Reformat accessions in Seq object to be <=20 chars long and not contain any special chars
    """
    if len(seq_record.id) > 20:
        short_id = seq_record.id[:20]
    else:
        short_id = seq_record.id

    seq_record.id = re.sub('[|,/,\,.,:,\,),(]', '_', short_id)
    return seq_record

def check_already_run(settings, input_seqs):
    """
    Check final output folder and 2 expected fail state folders
    """
    named_phylogenies = glob.glob(os.path.join(settings.dir_paths["named"], '*.named_tre'))
    failed_masks      = glob.glob(os.path.join(settings.dir_paths["mask_fail"], '*.mask_too_short'))
    failed_blasts     = glob.glob(os.path.join(settings.dir_paths["blast_fail"], '*.insufficient_hits'))
    # get rid of file extensions
    remove_extra_path = lambda filename: os.path.basename(os.path.splitext(filename)[0])

    putative_run_seqs = list(map(remove_extra_path, named_phylogenies+failed_masks+failed_blasts))

    run_seqs = set(putative_run_seqs)

    # make sure there aren't multiple terminal output files corresponding to the same
    # input sequence and if there are raise exception with a list of the duplicates
    if len(putative_run_seqs) != len(run_seqs):
        duplicates = [x for x, y in collections.Counter(putative_run_seqs).items() if y > 1]
        error_message = \
                "Duplicate outputs found in [{0}, {1}, {2}]: {3}".format(settings.dir_paths["name"],
                                                                         settings.dir_paths["mask_fail"],
                                                                         settings.dir_paths["blast_fail"],
                                                                         duplicates)

        settings.logger.log.error(error_message)
        raise ValueError(error_message)

    input_seq_ids = set([seq.id for seq in input_seqs])

    # therefore seqs need run is the assymetric set difference
    ids_needing_run = input_seq_ids - run_seqs

    # create list of seqs needing run from this
    seqs_needing_run = [seq for seq in input_seqs if seq.id in ids_needing_run]

    return seqs_needing_run


class PipeError(Exception):
    """
    Custom exception for a raised error in the
    execution of the core dendrogenous class
    Specifically when an output file either isn't
    created or is empty
    """
    def __init__(self, msg):
        self.msg = msg

class MaskFail(Exception):
    """
    Custom exception for the failure of
    masking - specifically not enough sites
    in the mask
    """
    pass

class GetSeqFail(Exception):
    """
    Custom exception for the failure of
    acquiring sequences - specifically not
    enough sequences recovered from blasting genomes
    """
    pass
