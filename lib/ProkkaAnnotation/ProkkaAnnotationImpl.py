# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import subprocess
import uuid
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from BCBio import GFF

from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI
#END_HEADER


class ProkkaAnnotation:
    '''
    Module Name:
    ProkkaAnnotation

    Module Description:
    A KBase module: ProkkaAnnotation
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    def _get_qualifier_value(self, qualifier):
        return qualifier[0] if (qualifier and len(qualifier) > 0) else None
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.scratch = config['scratch']
        #END_CONSTRUCTOR
        pass


    def annotate_contigs(self, ctx, params):
        """
        :param params: instance of type "AnnotateContigsParams" -> structure:
           parameter "assembly_ref" of type "assembly_ref" (Reference to an
           Assembly object in the workspace @id ws
           KBaseGenomeAnnotations.Assembly), parameter "output_workspace" of
           String, parameter "output_genome_name" of String
        :returns: instance of type "AnnotateContigsOutput" -> structure:
           parameter "output_genome_ref" of type "genome_ref" (Reference to
           an Genome object in the workspace @id ws KBaseGenomes.Genome)
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN annotate_contigs
        assembly_ref = params['assembly_ref']
        output_genome_name = params['output_genome_name']
        output_workspace = params['output_workspace']
        au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'], token=ctx['token'])
        orig_fasta_file = au.get_assembly_as_fasta({'ref': assembly_ref})['path']
        renamed_fasta_file = orig_fasta_file + "_renamed.fna"
        records = []
        new_ids_to_old = {}
        contig_counter = 0
        for record in SeqIO.parse(orig_fasta_file, "fasta"):
            contig_counter += 1
            old_id = record.id
            new_id = "contig_" + str(contig_counter)
            sequence = record.seq  # it has type 'Seq'
            record = SeqRecord(sequence, id=new_id, description="(" + old_id + ")")
            records.append(record)
            new_ids_to_old[new_id] = old_id
        SeqIO.write(records, renamed_fasta_file, "fasta")
        output_dir = "/kb/module/work/tmp/temp_" + str(hex(uuid.getnode()))
        subprocess.Popen(["perl", "/kb/prokka/bin/prokka", "--outdir", output_dir, "--prefix",
                          "mygenome", renamed_fasta_file], cwd=self.scratch).wait()
        faa_file = output_dir + "/mygenome.faa"
        cds_to_prot = {}
        for record in SeqIO.parse(faa_file, "fasta"):
            cds_to_prot[record.id] = str(record.seq)
        ffn_file = output_dir + "/mygenome.ffn"
        cds_to_dna = {}
        for record in SeqIO.parse(ffn_file, "fasta"):
            cds_to_dna[record.id] = str(record.seq)
        gff_file = output_dir + "/mygenome.gff"
        if not os.path.isfile(gff_file):
            raise ValueError("PROKKA output GFF file is not found")
        features = []
        with open(gff_file, "r") as f1:
            for rec in GFF.parse(f1):
                contig_id = new_ids_to_old[str(rec.id)]
                print "Contig: id=" + rec.id + ", " + contig_id
                for ft in rec.features:
                    loc = ft.location
                    min_pos = int(loc.start)
                    max_pos = int(loc.end)
                    strand = '+' if loc.strand == 1 else '-'
                    flen = max_pos - min_pos + 1
                    start = min_pos if strand == '+' else max_pos
                    qualifiers = ft.qualifiers
                    generated_id = self._get_qualifier_value(qualifiers.get('ID'))
                    name = self._get_qualifier_value(qualifiers.get('Name'))
                    ec = self._get_qualifier_value(qualifiers.get('eC_number'))
                    gene = self._get_qualifier_value(qualifiers.get('gene'))
                    product = self._get_qualifier_value(qualifiers.get('product'))
                    fid = name if name else generated_id
                    aliases = []
                    if gene:
                        aliases.append(gene)
                    if ec:
                        aliases.append(ec)
                    prot = cds_to_prot.get(generated_id)
                    prot_len = len(prot)
                    dna = cds_to_dna.get(generated_id)
                    dna_len = len(dna)
                    feature = {'id': fid, 'location': [[contig_id, start, strand, flen]], 
                               'type': 'gene', 'aliases': aliases, 'dna_sequence': dna, 
                               'dna_sequence_length': dna_len, 'protein_translation': prot, 
                               'protein_translation_length': prot_len}
                    if product:
                        feature['function'] = product
                    features.append(feature)
        genome = {'id': 'Unknown', 'features': features, 'scientific_name': 'Unknown', 
                  'domain': 'Bacteria', 'genetic_code': 11, 'assembly_ref': assembly_ref}
        ga = GenomeAnnotationAPI(os.environ['SDK_CALLBACK_URL'], token=ctx['token'])
        info = ga.save_one_genome_v1({'workspace': output_workspace, 'name': output_genome_name,
                                      'data': genome})['info']
        genome_ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        returnVal = {'output_genome_ref': genome_ref}
        #END annotate_contigs

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method annotate_contigs return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
