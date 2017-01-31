# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import subprocess
import uuid
import shutil

from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
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
        fasta_file = au.get_assembly_as_fasta({'ref': assembly_ref})['path']
        output_dir = "/kb/module/work/tmp/temp_" + str(hex(uuid.getnode()))
        subprocess.Popen(["perl", "/kb/prokka/bin/prokka", "--outdir", output_dir, 
                          "--prefix", "mygenome", fasta_file], cwd=self.scratch).wait()
        output_file = output_dir + "/mygenome.gbk"
        if not os.path.isfile(output_file):
            raise ValueError("PROKKA output genbank (.gbk) file is not found")
        gfu = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'], token=ctx['token'])
        genome_ref = gfu.genbank_to_genome({'file': {'path': output_file},
                                            'genome_name': output_genome_name,
                                            'workspace_name': output_workspace}
                                           )['genome_ref']
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
