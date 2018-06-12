# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests
import shutil

from os import environ

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from shutil import copyfile

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from ProkkaAnnotation.ProkkaAnnotationImpl import ProkkaAnnotation
from ProkkaAnnotation.ProkkaAnnotationServer import MethodContext
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from ProkkaAnnotation.authclient import KBaseAuth as _KBaseAuth
from AssemblySequenceAPI.AssemblySequenceAPIClient import AssemblySequenceAPI
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from DataFileUtil.DataFileUtilClient import DataFileUtil


class ProkkaAnnotationTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        config_file = environ.get("KB_DEPLOYMENT_CONFIG", None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items("ProkkaAnnotation"):
            cls.cfg[nameval[0]] = nameval[1]
        # Token validation
        token = environ.get("KB_AUTH_TOKEN", None)
        authServiceUrl = cls.cfg.get("auth-service-url",
                                     "https://kbase.us/services/authorization/Sessions/Login")
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don"t call any logging methods on the context object,
        # it"ll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({"token": token,
                        "user_id": user_id,
                        "provenance": [
                            {"service": "ProkkaAnnotation",
                             "method": "please_never_use_it_in_production",
                             "method_params": []
                             }],
                        "authenticated": 1})
        cls.wsURL = cls.cfg["workspace-url"]
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = ProkkaAnnotation(cls.cfg)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "wsName"):
            cls.wsClient.delete_workspace({"workspace": cls.wsName})
            print("Test workspace was deleted")

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, "wsName"):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_ProkkaAnnotation_" + str(suffix)
        ret = self.getWsClient().create_workspace({"workspace": wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def getBogusAssembly(self):
        # Create a fake assembly with lots of contigs

        assembly_file_name = "bogus.fna"  # "AP009048.fna"
        assembly_temp_file = os.path.join("/kb/module/work/tmp", assembly_file_name)
        with open(assembly_temp_file, "w") as f:
            for i in range(1, 30002):
                f.write("> contig_%d\n" % i)
                f.write("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC\n")

        assembly_name = "Assembly.2"
        au = AssemblyUtil(os.environ["SDK_CALLBACK_URL"], token=self.getContext()["token"])
        assembly_ref = au.save_assembly_from_fasta({"file": {"path": assembly_temp_file},
                                                    "workspace_name": self.getWsName(),
                                                    "assembly_name": assembly_name})
        self.assembly_ref = assembly_ref
        print("Uploaded bogus assembly " + str(assembly_ref))
        return assembly_ref

    def test_validation_integration(self):
        """
        This test does some basic validation tests of the required parameters.
        :return:
        """
        with self.assertRaises(KeyError):
            self.getImpl().annotate(self.getContext(), {})
        with self.assertRaises(KeyError):
            self.getImpl().annotate(self.getContext(), {'object_workspace': '0'})
        with self.assertRaises(Exception):
            self.getImpl().annotate(self.getContext(), {'object_workspace': '0', 'object_ref': 0})
        with self.assertRaises(Exception):
            self.getImpl().annotate(self.getContext(), {'object_workspace': '0', 'object_ref': 0,
                                                        'output_genome_name': 0})

    def test_reannotate_genome(self):
        gfu = GenomeFileUtil(os.environ["SDK_CALLBACK_URL"])

        genome_test_file = os.path.join("/kb/module/test/data/", "rhodobacter_genomic.gbff")
        genome_test_file_scratch = os.path.join("/kb/module/work/tmp", "rhodobacter_genomic.gbff")
        copyfile(genome_test_file, genome_test_file_scratch)

        genome_ref = gfu.genbank_to_genome({'file': {'path': genome_test_file_scratch},
                                            'workspace_name': self.getWsName(),
                                            'genome_name': "rhodobacter_genomic.gbff",
                                            'generate_ids_if_needed': 1})['genome_ref']
        #33281/2/1
        genome_name = "Rhodoannotated_by_prokka"
        print("ABOUT TO ANNOTATE GENOME")
        result = self.getImpl().annotate(self.getContext(),
                                         {"object_ref": genome_ref,
                                          "output_workspace": self.getWsName(),
                                          "output_genome_name": genome_name,
                                          "evalue": None,
                                          "fast": 0,
                                          "gcode": 0,
                                          "genus": "genus",
                                          "kingdom": "Bacteria",
                                          "metagenome": 0,
                                          "mincontiglen": 1,
                                          "norrna": 0,
                                          "notrna": 0,
                                          "rawproduct": 0,
                                          "rfam": 1,
                                          "scientific_name": "RhodoBacter"
                                          })[0]

        genome_ref = self.getWsName() + "/" + genome_name
        re_annotated_genome = self.getWsClient().get_objects([{"ref": genome_ref}])[0]["data"]


    def Xtest_reannotate_genome(self):
        """
        DOESN"T WORK ON CI WITH THIS DATASET, ONLY WITH APPDEV, THIS TEST IS COMMENTED OUT
        This test uploads the genome.json object, replacing the features with a single feature, and runs prokka against this feature.
        The test itself checks the feature to see that it has been updated.
        This test might break if Prokka decides that there is a better function name for this feature.
        :return:
        """
        gfu = GenomeFileUtil(os.environ["SDK_CALLBACK_URL"])

        genome_test_file = os.path.join("/kb/module/test/data/", "RHODO.json")
        genome_test_feature1_file = os.path.join("/kb/module/test/data/", "rsp_0986.json")
        genome_test_feature2_file = os.path.join("/kb/module/test/data/", "rsp_1428.json")
        genome_name = "RhodoBacter2.4"

        assembly_ref = self.getBogusAssembly()

        with open(genome_test_file, "r") as f:
            genome = json.load(f)
            genome['assembly_ref'] = assembly_ref

        # New function found by prokka
        with open(genome_test_feature1_file, "r") as f:
            target_feature_rsp_0986 = json.load(f)

        # Hypothetical Protein as determined by prokka
        with open(genome_test_feature2_file, "r") as f:
            target_feature_rsp_1428 = json.load(f)

        genome["features"] = [target_feature_rsp_0986, target_feature_rsp_1428]

        info = gfu.save_one_genome(
            {"workspace": self.getWsName(), "name": genome_name,
             "data": genome})["info"]

        genome_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])

        result = self.getImpl().annotate(self.getContext(),
                                         {"object_ref": genome_ref,
                                          "output_workspace": self.getWsName(),
                                          "output_genome_name": genome_name,
                                          "evalue": None,
                                          "fast": 0,
                                          "gcode": 0,
                                          "genus": "genus",
                                          "kingdom": "Bacteria",
                                          "metagenome": 0,
                                          "mincontiglen": 1,
                                          "norrna": 0,
                                          "notrna": 0,
                                          "rawproduct": 0,
                                          "rfam": 1,
                                          "scientific_name": "RhodoBacter"
                                          })[0]

        genome_ref = self.getWsName() + "/" + genome_name
        re_annotated_genome = self.getWsClient().get_objects([{"ref": genome_ref}])[0]["data"]

        old_feature = genome["features"][0]
        new_feature = re_annotated_genome["features"][0]

        # TEST NEW PROTEIN FUNCTION
        self.assertNotEquals(old_feature, new_feature)
        self.assertEqual(old_feature["function"], "fructokinase")
        self.assertEqual(new_feature["function"], "Pantothenate kinase")

        old_feature = genome["features"][1]
        new_feature = re_annotated_genome["features"][1]

        # TEST HYPOTHETICAL PROTEIN
        self.assertEqual(old_feature["function"], "putative Pre (Mob) type recombination enzyme")
        self.assertEqual(old_feature["function"], new_feature["function"])

    # NOTE: According to Python unittest naming rules test method names should start from "test". # noqa
    def Xtest_annotate_contigs(self):
        assembly_file_name = "small.fna"  # "AP009048.fna"
        assembly_test_file = os.path.join("/kb/module/test/data/", assembly_file_name)
        assembly_temp_file = os.path.join("/kb/module/work/tmp", assembly_file_name)
        shutil.copy(assembly_test_file, assembly_temp_file)
        assembly_name = "Assembly.1"
        au = AssemblyUtil(os.environ["SDK_CALLBACK_URL"])
        assembly_ref = au.save_assembly_from_fasta({"file": {"path": assembly_temp_file},
                                                    "workspace_name": self.getWsName(),
                                                    "assembly_name": assembly_name})
        # Add a genome to the WS to test ref_paths
        genome_name = "Genome.1"
        genome = {"id": "Unknown", "features": [],
                  "scientific_name": "",
                  "domain": "", "genetic_code": 0,
                  "assembly_ref": assembly_ref,
                  "cdss": [], "mrnas": [],
                  "source": "Magic!",
                  "gc_content": 0, "dna_size": 0,
                  "reference_annotation": 0}
        prov = self.getContext().provenance()
        gfu = GenomeFileUtil(os.environ["SDK_CALLBACK_URL"])
        info = gfu.save_one_genome(
            {"workspace": self.getWsName(), "name": genome_name,
             "data": genome, "provenance": prov})["info"]
        genome_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])
        result = self.getImpl().annotate(self.getContext(),
                                         {"object_ref": "{};{}".format(genome_ref, assembly_ref),
                                          "output_workspace": self.getWsName(),
                                          "output_genome_name": genome_name,
                                          "evalue": None,
                                          "fast": 0,
                                          "gcode": 0,
                                          "genus": "genus",
                                          "kingdom": "Bacteria",
                                          "metagenome": 0,
                                          "mincontiglen": 1,
                                          "norrna": 0,
                                          "notrna": 0,
                                          "rawproduct": 0,
                                          "rfam": 1,
                                          "scientific_name": "Super : diper - name;"
                                          })[0]
        rep = self.getWsClient().get_objects([{"ref": result["report_ref"]}])[0]["data"]
        self.assertTrue("text_message" in rep)
        print("Report:\n" + str(rep["text_message"]))
        genome_ref = self.getWsName() + "/" + genome_name
        genome = self.getWsClient().get_objects([{"ref": genome_ref}])[0]["data"]
        features_to_work = {}
        for feature in genome["features"]:
            features_to_work[feature["id"]] = feature["location"]
        aseq = AssemblySequenceAPI(os.environ["SDK_CALLBACK_URL"], token=self.getContext()["token"])
        dna_sequences = aseq.get_dna_sequences({"requested_features": features_to_work,
                                                "assembly_ref": genome["assembly_ref"]})[
            "dna_sequences"]
        bad_dnas = 0
        for feature in genome["features"]:
            if feature["dna_sequence"] != dna_sequences[feature["id"]]:
                bad_dnas += 1
        self.assertEqual(bad_dnas, 0)

    def Xtest_annotate_contigs_too_big(self):
        """
        simulate a metagenome contig file
        """
        assembly_ref = self.getBogusAssembly()
        genome_name = "Genome.1"
        # This should fail with an error
        with self.assertRaises(ValueError):
            result = self.getImpl().annotate(self.getContext(),
                                             {"object_ref": assembly_ref,
                                              "output_workspace": self.getWsName(),
                                              "output_genome_name": genome_name,
                                              "evalue": None,
                                              "fast": 0,
                                              "gcode": 0,
                                              "genus": "genus",
                                              "kingdom": "Bacteria",
                                              "metagenome": 0,
                                              "mincontiglen": 1,
                                              "norrna": 0,
                                              "notrna": 0,
                                              "rawproduct": 0,
                                              "rfam": 1,
                                              "scientific_name": "Super : diper - name;"
                                              })
