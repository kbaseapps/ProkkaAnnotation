# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import uuid
from pprint import pformat
from ProkkaAnnotation.Util.ProkkaUtils import ProkkaUtils
from Workspace.WorkspaceClient import Workspace as workspaceService
from KBaseReport.KBaseReportClient import KBaseReport

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
    VERSION = "0.1.0"
    GIT_URL = "https://github.com/bio-boris/ProkkaAnnotation.git"
    GIT_COMMIT_HASH = "faf3ff2ba692cf84681c52e444737b877de8bef4"

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
#        self.output_workspace = None
        self.kbr = KBaseReport(config["SDK_CALLBACK_URL"])
        #END_CONSTRUCTOR
        pass


    def annotate(self, ctx, params):
        """
        :param params: instance of type "AnnotateParams" (Required
           parameters: object_ref - reference to Assembly or Genome object,
           output_workspace - output workspace name, output_genome_name -
           output object name, Optional parameters (correspond to PROKKA
           command line arguments): --scientific_name Genome scientific name
           (default 'Unknown') --kingdom [X]     Annotation mode:
           Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria') --genus
           [X]       Genus name (triggers to use --usegenus) --gcode [N]     
           Genetic code / Translation table (set if --kingdom is set)
           (default '11') --metagenome      Improve gene predictions for
           highly fragmented genomes (default OFF) --rawproduct      Do not
           clean up /product annotation (default OFF) --fast            Fast
           mode - skip CDS /product searching (default OFF) --mincontiglen
           [N] Minimum contig size [NCBI needs 200] (default '1') --evalue
           [n.n]    Similarity e-value cut-off (default '1e-06') --rfam      
           Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default
           OFF) --norrna          Don't run rRNA search (default OFF)
           --notrna          Don't run tRNA search (default OFF)) ->
           structure: parameter "object_ref" of type "data_obj_ref"
           (Reference to an Assembly or Genome object in the workspace @id ws
           KBaseGenomeAnnotations.Assembly @id ws KBaseGenomes.Genome),
           parameter "output_workspace" of String, parameter
           "output_genome_name" of String, parameter "scientific_name" of
           String, parameter "kingdom" of String, parameter "genus" of
           String, parameter "gcode" of Long, parameter "metagenome" of type
           "boolean" (A boolean. 0 = false, anything else = true.), parameter
           "rawproduct" of type "boolean" (A boolean. 0 = false, anything
           else = true.), parameter "fast" of type "boolean" (A boolean. 0 =
           false, anything else = true.), parameter "mincontiglen" of Long,
           parameter "evalue" of String, parameter "rfam" of type "boolean"
           (A boolean. 0 = false, anything else = true.), parameter "norrna"
           of type "boolean" (A boolean. 0 = false, anything else = true.),
           parameter "notrna" of type "boolean" (A boolean. 0 = false,
           anything else = true.)
        :returns: instance of type "AnnotateOutput" -> structure: parameter
           "output_genome_ref" of type "genome_ref" (Reference to an Genome
           object in the workspace @id ws KBaseGenomes.Genome), parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN annotate
        output_workspace = params["output_workspace"]
        print("Input parameters: " + pformat(params))
        object_ref = params['object_ref']
        object_info = self.ws_client.get_object_info_new({"objects": [{"ref": object_ref}],
                                                           "includeMetadata": 1})[0]
        object_type = object_info[2]

        self.config['ctx'] = ctx
        prokka_runner = ProkkaUtils(self.config)
        report_message = ''
        objects_created = []
        
        if "KBaseGenomeAnnotations.Assembly" in object_type:
            #return [prokka_runner.annotate_assembly(params, object_info)]
            ret = prokka_runner.annotate_assembly(params, object_info)
            report_message += ret['report_message']
            objects_created +=  [{"ref": ret['output_genome_ref'], "description": "Annotated genome"}]
            
        elif "KBaseGenomes.Genome" in object_type:
            #return [prokka_runner.annotate_genome(params)]
            ret = prokka_runner.annotate_genome(params)
            report_message += ret['report_message']
            objects_created +=  [{"ref": ret['output_genome_ref'], "description": "Annotated genome"}]
            
        elif "KBaseSets.AssemblySet" in object_type:
            object_set = self.ws_client.get_objects([{"ref": object_ref}])[0]['data']['items']
            object_ref_list = []
            for object_list in object_set:
                object_info = self.ws_client.get_object_info_new({"objects": [{"ref": object_list['ref']}],
                                                           "includeMetadata": 1})[0]
                params['object_ref'] = object_list['ref']
                output_genome_name = object_info[1] + ".prokka"
                params['output_genome_name'] = output_genome_name
                ret = prokka_runner.annotate_assembly(params, object_info)
                report_message += ret['report_message']
                objects_created +=  [{"ref": ret['output_genome_ref'], "description": "Annotated genome"}]
                object_ref_list += [ret['output_genome_ref']]
            new_genomeset_ref = prokka_runner.make_genome_set(output_workspace, object_ref_list,output_genome_name)
            objects_created +=  [{"ref": ret['output_genome_ref'], "description": "Annotated genomeSet"}]
              
        elif "KBaseSearch.GenomeSet" in object_type:
            object_set = self.ws_client.get_objects([{"ref": object_ref}])[0]['data']["elements"]
            object_ref_list = []
            for object_list in object_set:
                object_info = self.ws_client.get_object_info_new({"objects": [{"ref": object_list}],
                                                           "includeMetadata": 1})[0]
                params['object_ref'] = object_list
                output_genome_name = object_info[1] + ".prokka"
                params['output_genome_name'] = output_genome_name
                ret = prokka_runner.annotate_genome(params)
                report_message += ret['report_message']
                objects_created +=  [{"ref": ret['output_genome_ref'], "description": "Annotated genome"}]
                object_ref_list += [ret['output_genome_ref']]
                
            new_genomeset_ref = prokka_runner.make_genome_set(output_workspace, object_ref_list, output_genome_name)
            objects_created +=  [{"ref": ret['output_genome_ref'], "description": "Annotated genomeSet"}]                
        else:
            raise Exception("Unsupported type" + object_type)
        
        report_info = self.kbr.create_extended_report(
            {"message": report_message,
             "objects_created": objects_created,
             "report_object_name": "kb_prokka_report_" + str(uuid.uuid4()),
             "workspace_name": output_workspace
             })
        return [{"report_name": report_info["name"],
               "report_ref": report_info["ref"]}]

        #END annotate

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method annotate return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {"state": "OK",
                     "message": "",
                     "version": self.VERSION,
                     "git_url": self.GIT_URL,
                     "git_commit_hash": self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
