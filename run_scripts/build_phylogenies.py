#!/usr/bin/env python
# Automatically generate phylogenies from a settings file
# specifying input fasta and genomes
import sys
import dendrogenous as dg
import dendrogenous.settings
import dendrogenous.utils
import dendrogenous.core
import joblib
import pickle

def main(settings_file):

    settings = dg.settings.Settings(settings_file)

    input_seqs       = dg.utils.parse_seqs(settings.input_seqs)

    seqs_needing_run = dg.utils.check_already_run(settings, input_seqs)

    r = joblib.Parallel(n_jobs=4, verbose=5)(joblib.delayed(pool_process)\
            (seq, settings) for seq in seqs_needing_run)

def pool_process(seq, settings):
    """
    A hacky and unecessary way to provide a pickle serealisable
    object for multiprocessing to pass off to workers
    - inefficiency in reinstantiating a settings class every time
    """
    seq_job = dg.core.Dendrogenous(seq, settings)
    seq_job.build_named_phylogeny()


if __name__=='__main__':
    if len(sys.argv) != 2:
        print("USAGE: build_phylogenies.py settings.json")
        sys.exit(1)

    main(sys.argv[1])
