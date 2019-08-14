# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
from pprint import pformat
from ProkkaAnnotation.Util.ProkkaUtils import ProkkaUtils
from installed_clients.WorkspaceClient import Workspace as workspaceService

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
    VERSION = "2.1.1"
    GIT_URL = "https://github.com/slebras/ProkkaAnnotation.git"
    GIT_COMMIT_HASH = "85491b1794f1181f52b33f174990e9a7388ebd1e"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        self.ws_client = workspaceService(config["workspace-url"])
        #END_CONSTRUCTOR
        pass


    def annotate(self, ctx, params):
        """
        :param params: instance of type "AnnotateParams" (Required
           parameters: object_ref - reference to Assembly or Genome object,
           output_workspace - output workspace name, output_genome_name -
           output object name, Optional parameters  (correspond to PROKKA
           command line arguments): --scientific_name  Genome scientific name
           (default 'Unknown') --kingdom [X]      Annotation mode:
           Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria') --genus
           [X]        Genus name (triggers to use --usegenus) --gcode [N]    
           Genetic code / Translation table (set if --kingdom is set)
           (default '11') --rawproduct       Do not clean up /product
           annotation (default OFF) --fast             Fast mode - skip CDS
           /product searching (default OFF) --mincontiglen [N] Minimum contig
           size [NCBI needs 200] (default '1') --evalue [n.n]     Similarity
           e-value cut-off (default '1e-06') --rfam             Enable
           searching for ncRNAs with Infernal+Rfam (SLOW!) (default OFF)
           --norrna           Don't run rRNA search (default OFF) --notrna   
           Don't run tRNA search (default OFF)) -> structure: parameter
           "object_ref" of type "data_obj_ref" (Reference to an Assembly or
           Genome object in the workspace @id ws
           KBaseGenomeAnnotations.Assembly @id ws KBaseGenomes.Genome),
           parameter "output_workspace" of String, parameter
           "output_genome_name" of String, parameter "scientific_name" of
           String, parameter "kingdom" of String, parameter "genus" of
           String, parameter "gcode" of Long, parameter "rawproduct" of type
           "boolean" (A boolean. 0 = false, anything else = true.), parameter
           "fast" of type "boolean" (A boolean. 0 = false, anything else =
           true.), parameter "mincontiglen" of Long, parameter "evalue" of
           String, parameter "rfam" of type "boolean" (A boolean. 0 = false,
           anything else = true.), parameter "norrna" of type "boolean" (A
           boolean. 0 = false, anything else = true.), parameter "notrna" of
           type "boolean" (A boolean. 0 = false, anything else = true.)
        :returns: instance of type "AnnotateOutput" -> structure: parameter
           "output_genome_ref" of type "genome_ref" (Reference to an Genome
           object in the workspace @id ws KBaseGenomes.Genome), parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN annotate
        print(("Input parameters: " + pformat(params)))
        object_ref = params['object_ref']
        object_info = self.ws_client.get_object_info_new({"objects": [{"ref": object_ref}],
                                                           "includeMetadata": 1})[0]
        object_type = object_info[2]

        self.config['ctx'] = ctx
        prokka_runner = ProkkaUtils(self.config)

        if "KBaseGenomeAnnotations.Assembly" in object_type:
            return [prokka_runner.annotate_assembly(params, object_info)]
        elif "KBaseGenomes.Genome" in object_type:
            return [prokka_runner.annotate_genome(params)]
        else:
            raise Exception("Unsupported type" + object_type)
        #END annotate

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method annotate return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def annotate_metagenome(self, ctx, params):
        """
        :param params: instance of type "MetagenomeAnnotateParams" (params:
           optional params:) -> structure: parameter "object_ref" of type
           "data_obj_ref" (Reference to an Assembly or Genome object in the
           workspace @id ws KBaseGenomeAnnotations.Assembly @id ws
           KBaseGenomes.Genome), parameter "output_workspace" of String,
           parameter "output_genome_name" of String, parameter "kingdom" of
           String, parameter "genus" of String, parameter "gcode" of Long,
           parameter "rawproduct" of type "boolean" (A boolean. 0 = false,
           anything else = true.), parameter "fast" of type "boolean" (A
           boolean. 0 = false, anything else = true.), parameter
           "mincontiglen" of Long, parameter "evalue" of String, parameter
           "rfam" of type "boolean" (A boolean. 0 = false, anything else =
           true.), parameter "norrna" of type "boolean" (A boolean. 0 =
           false, anything else = true.), parameter "notrna" of type
           "boolean" (A boolean. 0 = false, anything else = true.)
        :returns: instance of type "MetagenomeAnnotateOutput" -> structure:
           parameter "metagenome_ref" of type "genome_ref" (Reference to an
           Genome object in the workspace @id ws KBaseGenomes.Genome),
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN annotate_metagenome
        print("Input parameters: " + pformat(params))
        if params.get('object_ref'):
          object_ref = params['object_ref']
        else:
          raise ValueError("object_ref argument required.")
        if not params.get("output_workspace"):
          raise ValueError("output_workspace argument required.")
        if not params.get("output_metagenome_name"):
          raise ValueError("output_metagenome_name argument required.")

        object_info = self.ws_client.get_object_info_new({"objects": [{"ref": object_ref}],
                                                           "includeMetadata": 1})[0]
        object_type = object_info[2]
        params['metagenome'] = True
        self.config['ctx'] = ctx
        prokka_runner = ProkkaUtils(self.config)

        if "KBaseGenomeAnnotations.Assembly" in object_type:
            return [prokka_runner.annotate_assembly(params, object_info)]
        elif "KBaseMetagenomes.AnnotatedMetagenomeAssembly" in object_type:
            return [prokka_runner.annotate_metagenome(params)]
        else:
            raise Exception("Unsupported type" + object_type)
        #END annotate_metagenome

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method annotate_metagenome return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {"state": "OK",
                     "message": "",
                     "version": self.VERSION,
                     "git_url": self.GIT_URL,
                     "git_commit_hash": self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
