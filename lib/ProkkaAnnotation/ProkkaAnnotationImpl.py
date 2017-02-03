# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import subprocess
import uuid
import shutil
import hashlib
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from BCBio import GFF
from pprint import pformat

from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from GenomeAnnotationAPI.GenomeAnnotationAPIServiceClient import GenomeAnnotationAPI
from KBaseReport.KBaseReportClient import KBaseReport
from biokbase.workspace.client import Workspace as workspaceService  # @UnresolvedImport @IgnorePep8
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
    def _get_input_value(self, params, key):
        if not key in params:
            raise ValueError("Parameter " + key + " should be set in input parameters")
        return params[key]

    def _get_qualifier_value(self, qualifier):
        return qualifier[0] if (qualifier and len(qualifier) > 0) else None
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.scratch = config['scratch']
        self.ws_url = config['workspace-url']
        self.sw_url = config['service-wizard']
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
        print("Input parameters: " + pformat(params))
        assembly_ref = self._get_input_value(params, 'assembly_ref')
        output_genome_name = self._get_input_value(params, 'output_genome_name')
        output_workspace = self._get_input_value(params, 'output_workspace')
        ws_client = workspaceService(self.ws_url, token=ctx['token'])
        
        sso_ret = ws_client.get_objects([{'ref': "KBaseOntology/seed_subsystem_ontology"}])[0]
        sso = sso_ret['data']
        ec_to_sso = {}
        for sso_id in sso['term_hash']:
            sso_name = sso['term_hash'][sso_id]['name']
            if "(EC " in sso_name and sso_name.endswith(")"):
                ec = sso_name[sso_name.index("(EC ") + 4 : -1].strip()
                sso_list = ec_to_sso.get(ec, None)
                if not sso_list:
                    sso_list = []
                    ec_to_sso[ec] = sso_list
                sso_list.append(sso['term_hash'][sso_id])
        print("EC found in SSO: " + str(len(ec_to_sso)))
        print("EC with unique SSO: " + str(len([ec for ec in ec_to_sso if len(ec_to_sso[ec]) == 1])))
        sso_info = sso_ret['info']
        sso_ref = str(sso_info[6]) + '/' + str(sso_info[0]) + '/' + str(sso_info[4])
        with open('/kb/module/work/seed_so.json', 'w') as outfile:
            json.dump(sso, outfile, sort_keys = True, indent = 4)
        
        assembly_info = ws_client.get_object_info_new({'objects': [{'ref': assembly_ref}],
                                                       'includeMetadata': 1})[0]
        assembly_meta = assembly_info[10]
        gc_content = float(assembly_meta.get("GC content"))
        dna_size = int(assembly_meta.get("Size"))
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
        output_dir = "/kb/module/work/tmp/temp_" + str(uuid.uuid4())
        scientific_name = params.get('scientific_name', 'Unknown')
        # --kingdom [X]     Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
        kingdom = str(params.get('kingdom', "Bacteria"))
        prokka_cmd_list = ["perl", "/kb/prokka/bin/prokka", "--outdir", output_dir, "--prefix",
                          "mygenome", "--kingdom", kingdom, "--cpus", '1']
        # --genus [X]       Genus name (triggers to use --usegenus)
        if 'genus' in params and params['genus']:
            prokka_cmd_list.extend(['--genus', str(params['genus']), '--usegenus'])
        # --gcode [N]       Genetic code / Translation table (set if --kingdom is set) (default '0')
        gcode = 11
        if 'gcode' in params and params['gcode']:
            gcode = params['gcode']
        prokka_cmd_list.extend(['--gcode', str(gcode)])
        # --gram [X]        Gram: -/neg +/pos (default '')
        if 'gram' in params and params['gram']:
            prokka_cmd_list.extend(['--gram', str(params['gram'])])
        # --metagenome      Improve gene predictions for highly fragmented genomes (default OFF)
        if 'metagenome' in params and params['metagenome'] == 1:
            prokka_cmd_list.append("--metagenome")
        # --rawproduct      Do not clean up /product annotation (default OFF)
        if 'rawproduct' in params and params['rawproduct'] == 1:
            prokka_cmd_list.append("--rawproduct")
        # --fast            Fast mode - skip CDS /product searching (default OFF)
        if 'fast' in params and params['fast'] == 1:
            prokka_cmd_list.append("--fast")
        # --mincontiglen [N] Minimum contig size [NCBI needs 200] (default '1')
        if 'mincontiglen' in params and params['mincontiglen']:
            prokka_cmd_list.extend(['--mincontiglen', str(params['mincontiglen'])])
        # --evalue [n.n]    Similarity e-value cut-off (default '1e-06')
        if 'evalue' in params and params['evalue']:
            prokka_cmd_list.extend(['--evalue', str(params['evalue'])])
        # --rfam            Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default '0')
        if 'rfam' in params and params['rfam'] == 1:
            prokka_cmd_list.append("--rfam")
        # --norrna          Don't run rRNA search (default OFF)
        if 'norrna' in params and params['norrna'] == 1:
            prokka_cmd_list.append("--norrna")
        # --notrna          Don't run tRNA search (default OFF)
        if 'notrna' in params and params['notrna'] == 1:
            prokka_cmd_list.append("--notrna")
        prokka_cmd_list.append(renamed_fasta_file)
        print(str(prokka_cmd_list))
        subprocess.Popen(prokka_cmd_list, cwd=self.scratch).wait()
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
        cdss = []
        mrnas = []
        features = []
        non_hypothetical = 0
        genes_with_ec = 0
        genes_with_sso = 0
        prot_lengths = []
        with open(gff_file, "r") as f1:
            for rec in GFF.parse(f1):
                contig_id = new_ids_to_old[str(rec.id)]
                for ft in rec.features:
                    loc = ft.location
                    min_pos = int(loc.start)
                    max_pos = int(loc.end)
                    strand = '+' if loc.strand == 1 else '-'
                    flen = max_pos - min_pos + 1
                    start = min_pos if strand == '+' else max_pos
                    location = [[contig_id, start, strand, flen]]
                    qualifiers = ft.qualifiers
                    generated_id = self._get_qualifier_value(qualifiers.get('ID'))
                    if not generated_id:
                        # Skipping feature with no ID (mostly repeat regions)
                        continue
                    dna = cds_to_dna.get(generated_id)
                    if not dna:
                        # Skipping feature with no DNA (mostly repeat regions)
                        continue
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
                        genes_with_ec += 1
                    md5 = hashlib.md5(dna).hexdigest()
                    feature = {'id': fid, 'location': location, 'type': 'gene', 
                               'aliases': aliases, 'md5': md5, 'dna_sequence': dna,
                               'dna_sequence_length': len(dna),
                               }
                    if product:
                        feature['function'] = product
                        if product != "hypothetical protein":
                            non_hypothetical += 1
                    if ec and ec in ec_to_sso:
                        sso_list = ec_to_sso[ec]
                        sso_terms = {}
                        for sso_item in sso_list:
                            sso_terms[sso_item['id']] = {'id': sso_item['id'],
                                                         'evidence': [],
                                                         'term_name': sso_item['name'],
                                                         'ontology_ref': sso_ref,
                                                         'term_lineage': []}
                        feature['ontology_terms'] = {'SSO': sso_terms}
                        genes_with_sso += 1
                    cds = None
                    mrna = None
                    prot = cds_to_prot.get(generated_id)
                    if prot:
                        cds_id = fid + "_CDS"
                        mrna_id = fid + "_mRNA"
                        prot_len = len(prot)
                        prot_lengths.append(prot_len)
                        feature['protein_translation'] = prot
                        feature['protein_translation_length'] = prot_len
                        feature['cdss'] = [cds_id]
                        feature['mrnas'] = [mrna_id]
                        cds = {'id': cds_id, 'location': location, 'md5': md5, 'parent_gene': fid, 
                               'parent_mrna': mrna_id, 'function': (product if product else ''), 
                               'ontology_terms': {}, 'protein_translation': prot, 
                               'protein_translation_length': prot_len, 'aliases': aliases}
                        mrna = {'id': mrna_id, 'location': location, 'md5': md5, 
                                'parent_gene': fid, 'cds': cds_id}
                    features.append(feature)
                    if cds:
                        cdss.append(cds)
                    if mrna:
                        mrnas.append(mrna)
        genome = {'id': 'Unknown', 'features': features, 'scientific_name': scientific_name,
                  'domain': kingdom, 'genetic_code': gcode, 'assembly_ref': assembly_ref,
                  'cdss': cdss, 'mrnas': mrnas, 'source': 'PROKKA annotation pipeline',
                  'gc_content': gc_content,'dna_size': dna_size, 'reference_annotation': 0}
        prov = ctx.provenance()
        ga = GenomeAnnotationAPI(self.sw_url, token=ctx['token'])
        info = ga.save_one_genome_v1({'workspace': output_workspace, 'name': output_genome_name,
                                      'data': genome, 'provenance': prov})['info']
        genome_ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        
        # Prepare report
        report = ''
        report += 'Genome saved to: ' + output_workspace + '/' + output_genome_name + '\n'
        report += 'Number of genes predicted: ' + str(len(features)) + '\n'
        report += 'Number of protein coding genes: ' + str(len(prot_lengths)) + '\n'
        report += 'Number of genes with non-hypothetical function: ' + str(non_hypothetical) + '\n'
        report += 'Number of genes with EC-number: ' + str(genes_with_ec) + '\n'
        report += 'Number of genes with Seed Subsystem Ontology: ' + str(genes_with_sso) + '\n'
        report += 'Average protein length: ' + str(int(sum(prot_lengths) / 
                                                       float(len(prot_lengths)))) + ' aa.\n'
        kbr = KBaseReport(os.environ['SDK_CALLBACK_URL'], token=ctx['token'])
        report_info = kbr.create_extended_report(
            {'message': report,
             'objects_created': [{'ref': genome_ref, 'description': 'Annotated genome'}],
             'report_object_name': 'kb_prokka_report_' + str(uuid.uuid4()),
             'workspace_name': output_workspace
             })

        returnVal = {'output_genome_ref': genome_ref, 'report_name': report_info['name'], 
                     'report_ref': report_info['ref']}
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
