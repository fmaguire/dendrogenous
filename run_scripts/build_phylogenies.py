#!/usr/bin/env python
# Automatically generate phylogenies from a settings file
# specifying input fasta and genomes
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

    processes = [multiprocessing.Process(target=build_phylogeny, args=(seq, settings)) for seq in seqs_needing_run]

    for p in processes:
        p.start()

    for p in processes:
        p.join()

def build_phylogeny(seq, settings):
    seq_job = dg.core.Dendrogenous(seq, settings)
    seq_job.build_named_phylogeny()

if __name__=='__main__':
    if len(sys.argv) != 2:
        print("USAGE: build_phylogenies.py settings.json")
        sys.exit(1)

    main(sys.argv[1])