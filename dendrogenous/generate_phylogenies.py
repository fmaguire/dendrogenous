#!/usr/bin/env python
"""
Run script for pipeline
"""
import sys
import dendrogenous as dg

if __name__=='__main__':

    # parse CLI input settings or supplied run json
    settings = dg.utils.parse_settings()

    """
    check if directory structure exists:
        make it if it doesn't

      dir structure:
      $run_dir/run_info - input fasta, log data, run settings, sqlite states
               0.input_seqs - single sequence fasta files for each input seq
               1.blast_hits - multifasta of blast_hits for each input seq
               2.alignment - alignment of multifasta for each input seq
               3.mask- mask of multifasta alignment for each input seq
               4.phylogeny - phylogeny for masked alignment for each input seq
               5.named - named phylogeny for each input seq

       Initialise sqlite state database
                sqlite3.connect('run_state.db')
                c = conn.cursor()
                c.execute('''CREATE TABLE sequences
                            (SeqId varchar(255),
                             Blast varchar(255),
                             Align varchar(255),
                             Mask  varchar(255),
                             Tree  varchar(255),
                             Name  varchar(255))''')

                        blast align mask tree name
                seq_id1
                seq_id2
                seq_idN
                with each entry the md5sum of the file or NULL
                conn.commit()
                conn.close()


        if dirs exist:
            if db exists:
                SELECT SeqID FROM sequences ORDER BY desc






    """


    # check what needs run (pre-existing phylogenies)
    # specifically is the timestamp of the single sequence fasta newer than
    # phylogeny (non-existent don't exist)
    #


    #
    if settings.run_seqs == settings.input_seqs:
        dg.utils.generate_directory_structure(settings.run_name)



