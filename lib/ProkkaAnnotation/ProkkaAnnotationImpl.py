# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import subprocess
import uuid
import hashlib
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
from pprint import pformat,pprint
from subprocess import check_output, CalledProcessError
import sys

from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from KBaseReport.KBaseReportClient import KBaseReport
from biokbase.workspace.client import Workspace as workspaceService  # @UnresolvedImport @IgnorePep8
from DataFileUtil.DataFileUtilClient import DataFileUtil
from collections import namedtuple,defaultdict
#END_HEADER


class ProkkaAnnotation:
    """
    Module Name:
    ProkkaAnnotation

    Module Description:
    A KBase module: ProkkaAnnotation
    """

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.5"
    GIT_URL = "https://github.com/bio-boris/ProkkaAnnotation.git"
    GIT_COMMIT_HASH = "8c6463159e9effb84f719f124e9d82685cb70c9f"

    #BEGIN_CLASS_HEADER

    @staticmethod
    def _get_input_value(params, key):
        """Get value of key after checking for its existence

        :param params: Params dictionary haystack
        :param key: Key to search in Params
        :return: Parameter Value
        :raises ValueError: raises an exception if the key doesn"t exist
        """
        if not key in params:
            raise ValueError("Parameter " + key + " should be set in input parameters")
        return params[key]

    @staticmethod
    def _get_qualifier_value(qualifier):
        """Get first qualifier from the list of qualifiers

        :param qualifier: list contents of the qualifier from BCBio GFF Tools
        :return: first element in the list
        """
        return qualifier[0] if (qualifier and len(qualifier) > 0) else None

    def download_seed_data(self):
        """Download Seed Data Ontology, and set the gene_ontology reference (sso_ref) and
        the create a table from ec numbers to sso (ec_to_sso)

        :return: None
        """
        # Download Seed Reference Data
        sso_ret = self.ws_client.get_objects([{"ref": "KBaseOntology/seed_subsystem_ontology"}])[0]
        sso = sso_ret["data"]
        for sso_id in sso["term_hash"]:
            sso_name = sso["term_hash"][sso_id]["name"]
            if "(EC " in sso_name and sso_name.endswith(")"):
                ec = sso_name[sso_name.index("(EC ") + 4: -1].strip()
                sso_list = self.ec_to_sso.get(ec, None)
                if not sso_list:
                    sso_list = []
                    self.ec_to_sso[ec] = sso_list
                sso_list.append(sso["term_hash"][sso_id])
        print("EC found in SSO: " + str(len(self.ec_to_sso)))
        sso_info = sso_ret["info"]
        sso_ref = str(sso_info[6]) + "/" + str(sso_info[0]) + "/" + str(sso_info[4])
        with open("/kb/module/work/seed_so.json", "w") as outfile:
            json.dump(sso, outfile, sort_keys=True, indent=4)
        self.sso_ref = sso_ref

    def inspect_assembly(self, assembly_meta, assembly_ref):
        """Check to see if assembly has too many contigs and might not be a metagenome or
        non prokaryotic dataset

        :param assembly_meta: information about the assembly reference
        :param assembly_ref: the assembly reference number
        :return: a tuple containing gc_content and dna_size
        """
        gc_content = float(assembly_meta.get("GC content"))
        dna_size = int(assembly_meta.get("Size"))
        n_contigs = 0
        if "N Contigs" in assembly_meta:
            n_contigs = int(assembly_meta.get("N Contigs"))
        else:
            contig = self.ws_client.get_objects([{"ref": assembly_ref}])[0]
            n_contigs = len(contig["data"]["contigs"])
        if n_contigs >= 30000:
            message = """
            Hmmm.  There are over 30,000 contigs in this Assembly. 
            It looks like you are trying to run Prokka on a metagenome or non-prokaryotic data set. 
            If this is a metagenome data set we recommend using an App like MaxBin to first bin the contigs into genome-like bins. 
            These bins can then be individually annotated as a single genome using Prokka. 
            If this data comes from a Eukaryotic sample, KBase does not currently have an annotation app designed for Eukaryotes. 
            Alternatively, you can try reducing the number of contigs using a filter app.")
            raise ValueError("Too many contigs for Prokka.  See logs for details and suggestions
            """
            print(message)
            raise ValueError("Too many contigs for Prokka.  See logs for details and suggestions")

        assembly_info = namedtuple("assembly_info", "gc_content dna_size")
        return assembly_info(gc_content, dna_size)

    @staticmethod
    def create_renamed_assembly(assembly_fasta_filepath):
        """Rename records to be in the format of contig_N and output a new fasta file

        :param assembly_fasta_filepath:
        :return: The path to the fasta file with renamed contigs the number of contigs,
        the mapping from old ids to new ids, and the contigs as SeqRecords
        """
        records = []
        new_ids_to_old = {}
        contig_counter = 0
        for record in SeqIO.parse(assembly_fasta_filepath, "fasta"):
            contig_counter += 1
            old_id = record.id
            new_id = "contig_" + str(contig_counter)
            sequence = record.seq  # it has type "Seq"
            record = SeqRecord(sequence, id=new_id, description="(" + old_id + ")")
            records.append(record)
            new_ids_to_old[new_id] = old_id

        renamed_assembly_fasta_filepath = assembly_fasta_filepath + "_renamed.fna"
        SeqIO.write(records, renamed_assembly_fasta_filepath, "fasta")

        renamed_assembly = namedtuple("renamed_assembly",
                                      "filepath contig_counter new_ids_to_old records")
        return renamed_assembly(renamed_assembly_fasta_filepath, contig_counter, new_ids_to_old,
                                records)

    def run_prokka(self, params, subject_fasta_filepath):
        """Run Prokka

        :param params: Prokka parameters
        :param subject_fasta_filepath: The contigs or genes to run prokka against
        :return: The  directory with all of the prokka output files
        """
        output_dir = "/kb/module/work/tmp/temp_" + str(uuid.uuid4())

        kingdom = str(params.get("kingdom", "Bacteria"))
        # --kingdom [X]     Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default "Bacteria")

        prokka_cmd_list = ["perl", "/kb/prokka/bin/prokka", "--outdir", output_dir, "--prefix",
                           "mygenome", "--kingdom", kingdom]

        # --genus [X]       Genus name (triggers to use --usegenus)
        if "genus" in params and params["genus"]:
            prokka_cmd_list.extend(["--genus", str(params["genus"]), "--usegenus"])
        # --gcode [N]       Genetic code / Translation table (set if --kingdom is set) (default "0")
        if "gcode" in params and params["gcode"]:
            prokka_cmd_list.extend(["--gcode", str(params["gcode"])])
        else:
            prokka_cmd_list.extend(["--gcode", "0"])
        # --gram [X]        Gram: -/neg +/pos (default "")
        if "gram" in params and params["gram"]:
            raise ValueError("gram parameter is not supported in current Prokka installation")
        # --metagenome      Improve gene predictions for highly fragmented genomes (default OFF)
        if "metagenome" in params and params["metagenome"] == 1:
            prokka_cmd_list.append("--metagenome")
        # --rawproduct      Do not clean up /product annotation (default OFF)
        if "rawproduct" in params and params["rawproduct"] == 1:
            prokka_cmd_list.append("--rawproduct")
        # --fast            Fast mode - skip CDS /product searching (default OFF)
        if "fast" in params and params["fast"] == 1:
            prokka_cmd_list.append("--fast")
        # --mincontiglen [N] Minimum contig size [NCBI needs 200] (default "1")
        if "mincontiglen" in params and params["mincontiglen"]:
            prokka_cmd_list.extend(["--mincontiglen", str(params["mincontiglen"])])
        # --evalue [n.n]    Similarity e-value cut-off (default "1e-06")
        if "evalue" in params and params["evalue"]:
            prokka_cmd_list.extend(["--evalue", str(params["evalue"])])
        # --rfam            Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default "0")
        if "rfam" in params and params["rfam"] == 1:
            prokka_cmd_list.append("--rfam")
        # --norrna          Don"t run rRNA search (default OFF)
        if "norrna" in params and params["norrna"] == 1:
            prokka_cmd_list.append("--norrna")
        # --notrna          Don"t run tRNA search (default OFF)
        if "notrna" in params and params["notrna"] == 1:
            prokka_cmd_list.append("--notrna")
        prokka_cmd_list.append(subject_fasta_filepath)
        print("Prokka command line: " + str(prokka_cmd_list))

        try:
            check_output(prokka_cmd_list, cwd=self.scratch)
        except CalledProcessError as e:
            pprint(e)
        return output_dir

    @staticmethod
    def retrieve_prokka_results(output_dir):
        """ Gather up the relevant prokka results, load the records from the results files

        :param output_dir:
        :return: Sequences from the .faa .ffn files and the gff_filepath
        """
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

        prokka_results = namedtuple("prokka_results", "cds_to_prot cds_to_dna gff_filepath")
        return prokka_results(cds_to_prot, cds_to_dna, gff_file)

    def parse_prokka_results(self, **prokka_parse_parameters):
        """ Go through the prokka results from the input contigs and then
        create the features, mrnas and cdss components of the KbaseGenome.Genome object

        :param prokka_parse_parameters: gff_filepath, mappings
        :return: Genome:features Genome:cdss  Genome:mrnas report_message of genes discovered
        """
        gff_filepath = prokka_parse_parameters["gff_filepath"]
        cds_to_dna = prokka_parse_parameters["cds_to_dna"]
        cds_to_prot = prokka_parse_parameters["cds_to_prot"]
        new_ids_to_old = prokka_parse_parameters["new_ids_to_old"]

        cdss = []
        mrnas = []
        features = []
        non_hypothetical = 0
        genes_with_ec = 0
        genes_with_sso = 0
        prot_lengths = []
        with open(gff_filepath, "r") as f1:
            for rec in GFF.parse(f1):
                contig_id = new_ids_to_old[str(rec.id)]
                for ft in rec.features:
                    loc = ft.location
                    min_pos = int(loc.start) + 1
                    max_pos = int(loc.end)
                    strand = "+" if loc.strand == 1 else "-"
                    flen = max_pos - min_pos + 1
                    start = min_pos if strand == "+" else max_pos
                    location = [[contig_id, start, strand, flen]]
                    qualifiers = ft.qualifiers
                    generated_id = self._get_qualifier_value(qualifiers.get("ID"))
                    if not generated_id:
                        # Skipping feature with no ID (mostly repeat regions)
                        continue
                    dna = cds_to_dna.get(generated_id)
                    if not dna:
                        # Skipping feature with no DNA (mostly repeat regions)
                        continue
                    name = self._get_qualifier_value(qualifiers.get("Name"))
                    ec = self._get_qualifier_value(qualifiers.get("eC_number"))
                    gene = self._get_qualifier_value(qualifiers.get("gene"))
                    product = self._get_qualifier_value(qualifiers.get("product"))
                    fid = generated_id
                    aliases = []
                    if name:
                        aliases.append(name)
                    if gene:
                        aliases.append(gene)
                    if ec:
                        aliases.append(ec)
                        genes_with_ec += 1
                    md5 = hashlib.md5(dna).hexdigest()
                    feature = {"id": fid, "location": location, "type": "gene",
                               "aliases": aliases, "md5": md5, "dna_sequence": dna,
                               "dna_sequence_length": len(dna),
                               }
                    if product:
                        feature["function"] = product
                        if product != "hypothetical protein":
                            non_hypothetical += 1
                    if ec and ec in self.ec_to_sso:
                        sso_list = self.ec_to_sso[ec]
                        sso_terms = {}
                        for sso_item in sso_list:
                            sso_terms[sso_item["id"]] = {"id": sso_item["id"],
                                                         "evidence": [],
                                                         "term_name": sso_item["name"],
                                                         "ontology_ref": self.sso_ref,
                                                         "term_lineage": []}
                        feature["ontology_terms"] = {"SSO": sso_terms}
                        genes_with_sso += 1
                    cds = None
                    mrna = None
                    prot = cds_to_prot.get(generated_id)
                    if prot:
                        cds_id = fid + "_CDS"
                        mrna_id = fid + "_mRNA"
                        prot_len = len(prot)
                        prot_lengths.append(prot_len)
                        feature["protein_translation"] = prot
                        feature["protein_translation_length"] = prot_len
                        feature["cdss"] = [cds_id]
                        feature["mrnas"] = [mrna_id]
                        cds = {"id": cds_id, "location": location, "md5": md5, "parent_gene": fid,
                               "parent_mrna": mrna_id, "function": (product if product else ""),
                               "ontology_terms": {}, "protein_translation": prot,
                               "protein_translation_length": prot_len, "aliases": aliases}
                        mrna = {"id": mrna_id, "location": location, "md5": md5,
                                "parent_gene": fid, "cds": cds_id}
                    features.append(feature)
                    if cds:
                        cdss.append(cds)
                    if mrna:
                        mrnas.append(mrna)

        # Prepare report
        report = ""
        report += "Number of genes predicted: " + str(len(features)) + "\n"
        report += "Number of protein coding genes: " + str(len(prot_lengths)) + "\n"
        report += "Number of genes with non-hypothetical function: " + str(non_hypothetical) + "\n"
        report += "Number of genes with EC-number: " + str(genes_with_ec) + "\n"
        report += "Number of genes with Seed Subsystem Ontology: " + str(genes_with_sso) + "\n"
        report += "Average protein length: " + str(int(sum(prot_lengths) /
                                                       float(len(prot_lengths)))) + " aa.\n"

        annotated_assembly = namedtuple("annotated_assembly", "features cdss mrnas report_message")
        return annotated_assembly(features, cdss, mrnas, report)

    def get_new_annotations(self, gff_filepath):
        """

        :param gff_filepath: A dictionary of ids with products and ec numbers
        :return:
        """
        genome = {}
        with open(gff_filepath, "r") as f:
            for rec in GFF.parse(f):
                gid = rec.id
                gene_features = {"id": id}

                for feature in rec.features:
                    qualifiers = feature.qualifiers
                    if "product" in qualifiers:
                        gene_features["function"] = " ".join(qualifiers["product"])

                    if "eC_number" in qualifiers:
                        ec_numbers = qualifiers["eC_number"]
                        sso_terms = dict()
                        for ec in ec_numbers:
                            sso_list = self.ec_to_sso.get(ec, [])
                            for sso_item in sso_list:
                                sso_terms[sso_item["id"]] = {"id": sso_item["id"],
                                                             "evidence": [],
                                                             "term_name": sso_item["name"],
                                                             "ontology_ref": self.sso_ref,
                                                             "term_lineage": []}

                        gene_features["ontology_terms"] = sso_terms
                genome[gid] = gene_features

        return genome

    def write_genome_to_fasta(self, genome_data):
        """

        :param genome_data:
        :return:
        """
        fasta_for_prokka_filepath = os.path.join(self.scratch,
                                                 "features_" + str(uuid.uuid4()) + ".fasta")
        with open(fasta_for_prokka_filepath, "w") as f:
            for item in genome_data["data"]["features"]:
                if "id" not in item and "dna_sequence" not in item:
                    print("COMPLAINING ABOUT THIS BAD SEQUENCE")
                else:
                    f.write(">" + item["id"] + "\n" + item["dna_sequence"] + "\n")
        print("Finished printing to" + fasta_for_prokka_filepath)
        return fasta_for_prokka_filepath

    def annotate_genome_with_new_annotations(self, **annotation_args):
        """

        :param annotation_args: genome_data, new_annotations from prokka, and the output_genome_name
        :type
        :return:
        """
        genome_data = annotation_args["genome_data"]
        new_annotations = annotation_args["new_annotations"]
        output_genome_name = annotation_args["output_genome_name"]

        stats = {"current_functions": len(genome_data["data"]["features"]), "new_functions": 0,
                 "new_ontologies": 0}

        function_report_filepath = os.path.join(self.scratch, "ontology_report")
        ontology_report_filepath = os.path.join(self.scratch, "function_report")

        onto_r = open(function_report_filepath, "w")
        func_r = open(ontology_report_filepath, "w")

        func_r.write("function_id current_function new_function\n")
        onto_r.write("function_id current_ontology new_ontology\n")

        for i, feature in enumerate(genome_data["data"]["features"]):
            fid = feature["id"]
            current_function = feature.get("function", "")
            current_ontology = feature.get("ontology_terms", None)
            new_function = ""
            new_ontology = dict()
            if fid in new_annotations:
                new_function = new_annotations[fid].get("function", "")
                new_ontology = new_annotations[fid].get("ontology_terms", None)

                if new_function is not "" and new_function is not "hypothetical protein":
                    genome_data["data"]["features"][i]["function"] = new_function
                    stats["new_functions"] += 1

                if new_ontology is not None:
                    if "ontology_terms" not in genome_data["data"]["features"][i]:
                        genome_data["data"]["features"][i]["ontology_terms"] = {"SSO": {}}

                    for key in new_ontology.keys():
                        genome_data["data"]["features"][i]["ontology_terms"]["SSO"][key] = \
                            new_ontology[key]
                        stats["new_ontologies"] += 1

            func_r.write(json.dumps([fid, current_function, new_function]) + "\n")
            onto_r.write(json.dumps([fid, current_ontology, new_ontology]) + "\n")

        func_r.close()
        onto_r.close()

        info = self.gfu.save_one_genome({"workspace": self.output_workspace,
                                         "name": output_genome_name,
                                         "data": genome_data["data"],
                                         "provenance": self.ctx.provenance()})["info"]

        genome_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])

        annotated_genome = namedtuple("annotated_genome",
                                      "genome_ref function_report_filepath ontology_report_filepath stats")

        return annotated_genome(genome_ref, function_report_filepath, ontology_report_filepath,
                                stats)

    def upload_file(self, filepath, message="Annotation report generated by kb_prokka"):
        """
        Upload a file to shock
        :param filepath: File to upload
        :param message: Optional Upload Message
        :return:
        """
        output_file_shock_id = self.dfu.file_to_shock({"file_path": filepath})["shock_id"]
        print("Uploaded filepath" + filepath + "to shock and got id" + output_file_shock_id)
        return {"shock_id": output_file_shock_id,
                "name": os.path.basename(filepath),
                "label": os.path.basename(filepath),
                "description": message}

    def report_annotated_genome(self, genome):
        """ Create report output with newly reannotated genome, and some stats

        :param genome: Reannotated Genome Reference, Report Files and Stats
        :return: Reference to Report Object
        """
        genome_ref = genome.genome_ref
        stats = genome.stats

        file_links = [self.upload_file(genome.ontology_report_filepath),
                      self.upload_file(genome.function_report_filepath)]

        report_message = ("Genome Ref:{0}\n"
                          "Number of functions:{1}\n"
                          "Functions found:{2}\n"
                          "Ontology terms found:{3}\n").format(genome_ref,
                                                               stats["current_functions"],
                                                               stats["new_functions"],
                                                               stats["new_ontologies"])

        report_info = self.kbr.create_extended_report(
            {"message": report_message,
             "objects_created": [{"ref": genome_ref, "description": "Annotated genome"}],
             "file_links": file_links,
             "report_object_name": "kb_prokka_report_" + str(uuid.uuid4()),
             "workspace_name": self.output_workspace
             })

        return {"output_genome_ref": genome_ref, "report_name": report_info["name"],
                "report_ref": report_info["ref"]}

    def annotate_genome(self, params):
        """ User input an existing genome to re-annotate.

        :param params: Reference to the genome, Output File Name, UI Parameters
        :return: Report with Reannotated Genome and Stats about it
        """
        genome_ref = self._get_input_value(params, "object_ref")
        output_name = self._get_input_value(params, "output_genome_name")
        genome_data = self.dfu.get_objects({"object_refs": [genome_ref]})["data"][0]
        fasta_for_prokka_filepath = self.write_genome_to_fasta(genome_data)
        output_dir = self.run_prokka(params, fasta_for_prokka_filepath)
        prokka_results = self.retrieve_prokka_results(output_dir)
        new_annotations = self.get_new_annotations(prokka_results.gff_filepath)
        annotated_genome = self.annotate_genome_with_new_annotations(genome_data=genome_data,
                                                                     new_annotations=new_annotations,
                                                                     output_genome_name=output_name)
        return self.report_annotated_genome(annotated_genome)

    def annotate_assembly(self, params, assembly_info):
        """
        Annotate an assembly with Prokka. The steps include to download the assembly as a fasta file,
        rename the contigs, run prokka against the contigs, parse the results, and finally,
        create and upload a genome object.

        :param params: object reference, output_genome_name and output_workspace
        :param assembly_info: Information used to determine if the assembly is too big
        :return: Report with newly annotated assembly as a genome, and stats about it
        """

        assembly_ref = self._get_input_value(params, "object_ref")
        output_genome_name = self._get_input_value(params, "output_genome_name")
        output_workspace = self._get_input_value(params, "output_workspace")
        assembly_info = self.inspect_assembly(assembly_info[10], assembly_ref)
        orig_fasta_file = self.au.get_assembly_as_fasta({"ref": assembly_ref})["path"]

        # Rename Assembly and Keep Track of Old Contigs
        renamed_assembly = self.create_renamed_assembly(orig_fasta_file)
        # Run Prokka with the modified, renamed fasta file
        output_dir = self.run_prokka(params, renamed_assembly.filepath)
        # Prokka_results
        prokka_results = self.retrieve_prokka_results(output_dir)
        # Parse Results
        annotated_assembly = self.parse_prokka_results(gff_filepath=prokka_results.gff_filepath,
                                                       cds_to_dna=prokka_results.cds_to_dna,
                                                       cds_to_prot=prokka_results.cds_to_dna,
                                                       new_ids_to_old=renamed_assembly.new_ids_to_old)
        genome = {"id": "Unknown",
                  "features": annotated_assembly.features,
                  "scientific_name": params.get("scientific_name", "Unknown"),
                  "domain": str(params.get("kingdom", "Bacteria")),
                  "genetic_code": params.get("gcode", 0),
                  "assembly_ref": assembly_ref,
                  "cdss": annotated_assembly.cdss,
                  "mrnas": annotated_assembly.mrnas,
                  "source": "PROKKA annotation pipeline",
                  "gc_content": assembly_info.gc_content,
                  "dna_size": assembly_info.dna_size,
                  "reference_annotation": 0}

        info = self.gfu.save_one_genome({"workspace": output_workspace,
                                         "name": output_genome_name,
                                         "data": genome,
                                         "provenance": self.ctx.provenance()})["info"]

        genome_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])

        report_message = "Genome saved to: " + output_workspace + "/" + \
                         output_genome_name + "\n" + annotated_assembly.report_message
        report_info = self.kbr.create_extended_report(
            {"message": report_message,
             "objects_created": [{"ref": genome_ref, "description": "Annotated genome"}],
             "report_object_name": "kb_prokka_report_" + str(uuid.uuid4()),
             "workspace_name": output_workspace
             })

        return {"output_genome_ref": genome_ref, "report_name": report_info["name"],
                "report_ref": report_info["ref"]}


    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn"t
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.scratch = config["scratch"]
        self.ws_url = config["workspace-url"]
        self.gfu = GenomeFileUtil(os.environ["SDK_CALLBACK_URL"])
        self.ws_client = None
        self.sso_ref = None
        self.au = None
        self.ctx = None
        self.kbr = None
        self.dfu = None
        self.ec_to_sso = {}
        self.output_workspace = None
        #END_CONSTRUCTOR
        pass

    def annotate(self, ctx, params):
        """
        :param params: instance of type "AnnotateParams" (Required
           parameters: object_ref - reference to Assembly or Genome object,
           output_workspace - output workspace name, output_genome_name -
           output object name, Optional parameters (correspond to PROKKA
           command line arguments): --scientific_name Genome scientific name
           (default "Unknown") --kingdom [X]     Annotation mode:
           Archaea|Bacteria|Mitochondria|Viruses (default "Bacteria") --genus
           [X]       Genus name (triggers to use --usegenus) --gcode [N]
           Genetic code / Translation table (set if --kingdom is set)
           (default "11") --metagenome      Improve gene predictions for
           highly fragmented genomes (default OFF) --rawproduct      Do not
           clean up /product annotation (default OFF) --fast            Fast
           mode - skip CDS /product searching (default OFF) --mincontiglen
           [N] Minimum contig size [NCBI needs 200] (default "1") --evalue
           [n.n]    Similarity e-value cut-off (default "1e-06") --rfam
           Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default
           OFF) --norrna          Don"t run rRNA search (default OFF)
           --notrna          Don"t run tRNA search (default OFF)) ->
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
        print("Input parameters: " + pformat(params))
        self.output_workspace = params["output_workspace"]
        object_ref = self._get_input_value(params, "object_ref")
        self.ws_client = workspaceService(self.ws_url, token=ctx["token"])
        self.ctx = ctx
        self.au = AssemblyUtil(os.environ["SDK_CALLBACK_URL"], token=ctx["token"])
        self.kbr = KBaseReport(os.environ["SDK_CALLBACK_URL"], token=ctx["token"])
        self.dfu = DataFileUtil(os.environ["SDK_CALLBACK_URL"])
        self.download_seed_data()

        object_info = self.ws_client.get_object_info_new({"objects": [{"ref": object_ref}],
                                                          "includeMetadata": 1})[0]
        object_type = object_info[2]

        if "KBaseGenomeAnnotations.Assembly" in object_type:
            return [self.annotate_assembly(params, object_info)]
        elif "KBaseGenomes.Genome" in object_type:
            return [self.annotate_genome(params)]
        else:
            raise Exception("Unsupported type" + object_type)
        #END annotate

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {"state": "OK",
                     "message": "",
                     "version": self.VERSION,
                     "git_url": self.GIT_URL,
                     "git_commit_hash": self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
