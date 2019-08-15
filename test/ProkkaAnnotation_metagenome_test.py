# -*- coding: utf-8 -*-
import json  # noqa: F401
import os  # noqa: F401
import shutil
import time
import unittest
from configparser import ConfigParser  # py3
from os import environ
from shutil import copyfile

from ProkkaAnnotation.ProkkaAnnotationImpl import ProkkaAnnotation
from ProkkaAnnotation.ProkkaAnnotationServer import MethodContext
from ProkkaAnnotation.authclient import KBaseAuth as _KBaseAuth
from installed_clients.AssemblySequenceAPIClient import AssemblySequenceAPI
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService


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
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.gfu = GenomeFileUtil(cls.callback_url)
        cls.au = AssemblyUtil(cls.callback_url)
        cls.scratch = cls.cfg['scratch']

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

    def create_metagenome(self, gff, fasta, output_name):
        """"""
        fasta_path = self.scratch + "/metagenome.fa"
        gff_path = self.scratch + "/metagenome.gff"
        shutil.copyfile(gff, gff_path)
        shutil.copyfile(fasta, fasta_path)
        metagenome_ref = self.gfu.fasta_gff_to_metagenome({
            "gff_file": {'path': gff_path},
            "fasta_file": {'path': fasta_path},
            "genome_name": output_name,
            "workspace_name": self.getWsName(),
            "generate_missing_genes": True
        })['genome_ref']
        return metagenome_ref

    def create_assembly(self, fasta, output_name):
        """copies to scratch then saves as assembly"""
        fasta_path = self.scratch + "/temp_assembly.fa"
        shutil.copyfile(fasta, fasta_path)
        assembly_ref = self.au.save_assembly_from_fasta({
            "file": {"path": fasta_path},
            "workspace_name": self.getWsName(),
            "assembly_name": output_name,
            "type": "metagenome"
        })
        return assembly_ref

    def _filter_gff(self, gff):
        fasta = []
        save = []
        with open(gff) as f:
            for l in f:
                if '##FASTA' in l:
                    for line in f:
                        fasta.append(line)
                    break
                if '##' in l:
                    continue
                save.append(l)
        gff_path = gff + "_edited.gff"
        with open(gff_path, 'w') as f:
            for l in save:
                f.write(l.strip() + '\n')
        fasta_path = gff  + "edited.fa"
        with open(fasta_path, 'w') as f:
            for l in fasta:
                f.write(l.strip() + '\n')
        return gff_path, fasta_path

    def verify_output(self, ret):
        self.assertTrue('output_metagenome_ref' in ret)
        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)

    @unittest.skip('x')
    def test_annotate_metagenome_from_ama_ref(self):
        ret = self.serviceImpl.annotate_metagenome(self.getContext(), {
            "object_ref": "",
            "output_workspace": self.getWsName(),
            "output_metagenome_name": "test_metagenome_from_metagenome"
        })[0]
        self.verify_output(ret)

    @unittest.skip('x')
    def test_annotate_metagenome_from_annotated_metagenome_assembly_short(self):
        metagenome_gff = "data/short_one.gff"
        metagenome_fasta = "data/short_one.fa"
        ref = self.create_metagenome(metagenome_gff, metagenome_fasta, 'metagenome_metagenome_short')
        ret = self.serviceImpl.annotate_metagenome(self.getContext(), {
            "object_ref": ref,
            "output_workspace": self.getWsName(),
            "output_metagenome_name": "test_metagenome_from_metagenome_short"
        })[0]
        self.verify_output(ret)

    def test_annotate_metagenome_from_annotated_metagenome_assembly(self):
        metagenome_gff = "data/59111.assembled.gff"
        metagenome_fasta = "data/59111.assembled.fna"

        ref = self.create_metagenome(metagenome_gff, metagenome_fasta, 'metagenome_metagenome')
        ret = self.serviceImpl.annotate_metagenome(self.getContext(), {
            "object_ref": ref,
            "output_workspace": self.getWsName(),
            "output_metagenome_name": "test_metagenome_from_metagenome"
        })[0]
        self.verify_output(ret)

    def test_annotate_metagenome_from_assembly(self):
        """"""
        metagenome_fasta = "data/59111.assembled.fna"
        ref = self.create_assembly(metagenome_fasta, 'metagenome_assembly')
        ret = self.serviceImpl.annotate_metagenome(self.getContext(), {
            "object_ref": ref,
            "output_workspace": self.getWsName(),
            "output_metagenome_name": "test_assembly_metagenome",
        })[0]
        self.verify_output(ret)

    def test_validation_integration(self):
        """
        This test does some basic validation tests of the required parameters.
        :return:
        """
        with self.assertRaises(ValueError):
            self.getImpl().annotate_metagenome(self.getContext(), {})
        with self.assertRaises(ValueError):
            self.getImpl().annotate_metagenome(self.getContext(), {"object_workspace": "0"})
        with self.assertRaises(ValueError):
            self.getImpl().annotate_metagenome(self.getContext(), {"object_workspace": "0", "object_ref": 0})
        with self.assertRaises(ValueError):
            self.getImpl().annotate_metagenome(self.getContext(), {"object_workspace": "0", "object_ref": 0,
                                                        "output_metagenome_name": 0})
