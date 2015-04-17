#!/usr/bin/env python

import sys
import dendrogenous as dg
import dendrogenous.settings
import dendrogenous.utils
import dendrogenous.core
import multiprocessing

def main(settings_file):

    settings         = dg.settings.Settings(settings_file)

    input_seqs       = dg.utils.parse_seqs(settings.input_seqs)

    seqs_needing_run = dg.utils.check_already_run(settings, input_seqs)

    processes = [multiprocessing.Process(target=build_phylogeny, args=(seq, settings)) for seq in input_seqs]

    for p in processes:
        p.start()

    for p in processes:
        p.join()

def build_phylogeny(seq, settings):
    seq_job = dg.core.Dendrogenous(seq, settings)
    seq_job.build_named_phylogeny()

if __name__=='__main__':
    main(sys.argv[1])
