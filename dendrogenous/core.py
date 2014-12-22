#!/usr/bin/env python

import dendrogenous.utils as utils

class dendrogenous():
    """
    Tree generation class with methods for all the steps required to generate
    phylogenies
    """

    def __init__(self, seq_record, settings):
        """
        Initialise the class and truncate the accession to 20-chars
        """

        self.seq_name = utils.reformat_accession(seq_record)
        self.seed = ">{0}\n{1}".format(self.seq_name,
                                       seq_record.seq)


        self.output_dir = settings['output_dir']
        # check precomputed state of class by looking for files matching
        # this name in output folders (globs)
        self.state = utils.get_state(self.seq_name, self.output_dir)

    def _blast(self, genome_name):
        """
        Blast seed sequence against db using BLASTP
        """
        blast_settings = self.settings.get('blast_settings', default={})
        num_seqs = blast_settings.get('num_seqs', default=1)
        e_value = blast_settings.get('e_value', default=1e-5)

        blast_cmd = "blastp -db {1} " \
                    " -evalue {3} -max_target_seqs {4}" \
                    " -outfmt 5".format(genome,
                                        evalue,
                                        seqs_to_return_per_genome)

        blast_output = utils.execute_cmd(blast_cmd, stdin_str=self.seed)


    def _parse_blast(self, blast_output):
        """
        Parses a blast_output files and queries the
        MySQLdb to recover the appropriate sequences
        """

        # no get with default vals as this is a mandatory part of the
        # settings.json - use them to connect to mysql server
        db_config = self.settings['db_config']
        con = MySQLdb.connect(host=db_config['host'],
                              user=db_config['user'],
                              passwd=db_config['passwd'],
                              db=db_config['db'])
        cur = con.cursor()

        # Parse the xml output returns an iterator of blast record
        # objects over the xml blastoutput file
        blast_hits = NCBIXML.parse(io.StringIO(blast_output))

        hit_id_set = set()
        for blast_record in blast_hits:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hit_id_set.add(alignment.hit_def)

        hit_records = []
        for hit_id in hit_id_set:
            cur.execute("SELECT sequence FROM cider "
                        "WHERE protein_ID='{0}'".format(hit_id))
            sequence = cur.fetchone()
            if sequence is None:
                with open('protein_IDs_that_dont_exist', 'a+') as fh:
                   fh.write(str(hit_id) + '\n')
                continue
                #raise Exception("Failed to find"
                #                " protein_ID in DB: {0}".format(hit_id))

            sequence_record = SeqRecord(\
                Seq(sequence[0], IUPAC.protein), id=hit_id, description='')
            hit_records.append(sequence_record)
        con.close()


    def seqs(self):
        """
        Get similar sequences to seed by blasting genomes and parsing the output
        """
        self.parsed_blast_hits = os.path.join("1.seqs", self.seq_name)
        # file containing blast hits

        for genome in settings['genome_paths']:
            blast_output = self._blast()
            blast_hits = self._parse_blast(blast_output)
            num_hits += len(hits)
            SeqIO.write(blast_hits, self.parsed_blast_hits, 'fasta')
            os.unlink(temporary_blast_hits)

        if num_seqs < minimum_seqs:
            os.rename(parsed_output_seqs,
                      parsed_output_seqs + '_too_few_seqs_for_alignment')

        # update state correctly


    def align(self):
        """
        Align input seqs using kalign
        """

        if not os.path.exists(ALIGNMENT_FILE):
            self.seqs()

        align_cmd = "kalign -i {0} -o {1}".format(input_seqs,
                                              output_aligned_seqs)
        subprocess.call(align_cmd)

    def mask(self):
        '''Mask input alignment using trimal'''

        if not os.path.exists(alignment_file):
            self.align()

        mask_cmd = "trimal -in {0} -out {1} -nogaps".format(input_aligned_seqs,
                                                            output_masked_seqs)
        accessories.system_call(mask_cmd)

        mask_length = accessories.mask_check(output_masked_seqs)

        if mask_length < cut_off:
            remask_cmd = "trimal -in {0}" \
                         " -out {1} -automated1".format(input_aligned_seqs,
                                                        output_masked_seqs)
            accessories.system_call(remask_cmd)
            new_mask_length = accessories.mask_check(output_masked_seqs)
            if mask_length < cut_off:
                mask_error = "{0} is" \
                             " too short" \
                             " after masking".format(output_masked_seqs)
                os.rename(output_masked_seqs,
                          output_masked_seqs + "_mask_too_short")

                # write to sqlite db


    def phylogeny(self):
        '''Generate phylogeny from masked seqs using FastTree2'''

        if not os.path.exists(MASK_FILE):
            self.mask()

        phylogeny_cmd = "FastTree -bionj -slow" \
                        " -quiet -out {1} {0}".format(input_masked_seqs,
                                                      output_phylogeny)
        subprocess.call(phylogeny_cmd)


    def rename(self):
        """

        """
        if not os.path.exists(PHYLOGENY):
            self.phylogeny()


