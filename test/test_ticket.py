# -*- coding: utf-8 -*-
import json  # noqa: F401
import os  # noqa: F401
import time
import unittest
from uuid import uuid4
from configparser import ConfigParser  # py3
from os import environ

from ProkkaAnnotation.ProkkaAnnotationImpl import ProkkaAnnotation
from ProkkaAnnotation.ProkkaAnnotationServer import MethodContext
from ProkkaAnnotation.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace as workspaceService


class TestTicket(unittest.TestCase):

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

    def test_thing(self):
        self.getImpl().annotate(self.getContext(), {
            "object_ref": "57447/5/2",
            "output_workspace": self.getWsName(),
            "output_genome_name": str(uuid4()),
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
            "scientific_name": "Intrasporangium calvum C5"
        })[0]
