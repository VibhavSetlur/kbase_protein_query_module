"""KBase API utilities for interacting with KBase services and data."""

import os
import uuid
from os.path import exists
from pathlib import Path
from typing import Any, Optional, Union

import requests

from .shared_env_utils import SharedEnvUtils

# TODO: Need to write the callback service and run it on poplar, then write and run tests for this module


class KBCallbackUtils(SharedEnvUtils):
    """Utilities enabling execution of KBase callbacks"""

    def __init__(
        self,
        callback_directory: Optional[Union[str, os.PathLike]] = "/tmp/scratch",
        callback_url: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        """Initialize KBase callback utilities."""
        super().__init__(**kwargs)
        # This should be a unique username for your own callback server
        self._callback_url = callback_url
        self._callback_directory = callback_directory
        self._callback_clients = {}

    def initialize_callback(self):
        """Initialize the callback URL for KBase services.

        This method should be called to set up the callback URL before making any API calls.
        """
        callback_service_url = self.get_config("DevEnv", "callback_service_url")
        if not callback_service_url:
            raise ValueError(
                "Callback service URL not configured. Please set 'callback_service_url' in your configuration."
            )
        response = requests.post(
            f"{callback_service_url}/start",
            json={"token": self.get_token(namespace="kbase")},
        )
        if response.status_code != 200:
            raise RuntimeError(f"Failed to start callback service: {response.text}")
        output = response.json()
        self._callback_url = output["callback_url"]
        os.makedirs(str(Path(self._callback_directory).parent()), exist_ok=True)
        if exists(self._callback_directory):
            os.remove(self._callback_directory)
        os.symlink(
            self._callback_directory + output["directory"],
            "/tmp/" + output["directory"] + "/scratch",
        )
        self.set_config("DevEnv", "call_back_folder", output["callback_folder"])

    def stop_callback(self):
        """Stop the callback service.

        This method should be called to stop the callback service when it is no longer needed.
        """
        callback_service_url = self.get_config("DevEnv", "callback_service_url")
        if not callback_service_url:
            raise ValueError(
                "Callback service URL not configured. Please set 'callback_service_url' in your configuration."
            )
        response = requests.post(
            f"{callback_service_url}/stop",
            json={"token": self.get_token(namespace="kbase")},
        )
        if response.status_code != 200:
            raise RuntimeError(f"Failed to stop callback service: {response.text}")
        self._callback_url = None

    def report_client(self):
        if self._callback_url is None:
            raise ValueError(
                "Either set callback URL if you're using an SDK module, or call the initialize callback function."
            )
        if "KBaseReport" not in self._callback_clients:
            from installed_clients.KBaseReportClient import KBaseReport

            self._callback_clients["KBaseReport"] = KBaseReport(
                self._callback_url, token=self.get_token(namespace="kbase")
            )
        return self._callback_clients["KBaseReport"]

    def dfu_client(self):
        if self._callback_url is None:
            raise ValueError(
                "Either set callback URL if you're using an SDK module, or call the initialize callback function."
            )
        if "DataFileUtil" not in self._callback_clients:
            from installed_clients.DataFileUtilClient import DataFileUtil

            self._callback_clients["DataFileUtil"] = DataFileUtil(
                self._callback_url, token=self.get_token(namespace="kbase")
            )
        return self._callback_clients["DataFileUtil"]

    def gfu_client(self):
        if self._callback_url is None:
            raise ValueError(
                "Either set callback URL if you're using an SDK module, or call the initialize callback function."
            )
        if "GenomeFileUtil" not in self._callback_clients:
            from installed_clients.GenomeFileUtilClient import GenomeFileUtil

            self._callback_clients["GenomeFileUtil"] = GenomeFileUtil(
                self._callback_url, token=self.get_token(namespace="kbase")
            )
        return self._callback_clients["GenomeFileUtil"]

    def afu_client(self):
        if self._callback_url is None:
            raise ValueError(
                "Either set callback URL if you're using an SDK module, or call the initialize callback function."
            )
        if "AssemblyUtil" not in self._callback_clients:
            from installed_clients.AssemblyUtilClient import AssemblyUtil

            self._callback_clients["AssemblyUtil"] = AssemblyUtil(
                self._callback_url, token=self.get_token(namespace="kbase")
            )
        return self._callback_clients["AssemblyUtil"]

    def rast_client(self):
        if "RAST_SDK" not in self._callback_clients:
            from installed_clients.RAST_SDKClient import RAST_SDK

            self._callback_clients["RAST_SDK"] = RAST_SDK(
                self._callback_url, token=self.get_token(namespace="kbase")
            )
        return self._callback_clients["RAST_SDK"]

    def anno_client(self):
        if self._callback_url is None:
            raise ValueError(
                "Either set callback URL if you're using an SDK module, or call the initialize callback function."
            )
        if "cb_annotation_ontology_api" not in self._callback_clients:
            from installed_clients.cb_annotation_ontology_apiClient import (
                cb_annotation_ontology_api,
            )

            self._callback_clients["cb_annotation_ontology_api"] = (
                cb_annotation_ontology_api(
                    self._callback_url, token=self.get_token(namespace="kbase")
                )
            )
        return self._callback_clients["cb_annotation_ontology_api"]

    def devutil_client(self):
        if self._callback_url is None:
            raise ValueError(
                "Either set callback URL if you're using an SDK module, or call the initialize callback function."
            )
        if "KBDevUtils" not in self._callback_clients:
            from installed_clients.chenry_utility_moduleClient import (
                chenry_utility_module,
            )

            self._callback_clients["KBDevUtils"] = chenry_utility_module(
                self._callback_url, token=self.get_token(namespace="kbase")
            )
        return self._callback_clients["KBDevUtils"]

    def annotate_genome_with_rast(self, genome_id, ws=None, output_ws=None):
        if not output_ws:
            output_ws = ws
        rast_client = self.rast_client()
        output = rast_client.annotate_genome(
            {
                "workspace": output_ws,
                "input_genome": genome_id,
                "output_genome": genome_id + ".RAST",
            }
        )
        return output["workspace"] + "/" + output["id"]

    def save_genome_or_metagenome(self, objid, workspace, obj_json):
        self.set_ws(workspace)
        save_output = self.gfu_client().save_one_genome(
            {
                "name": objid,
                "data": obj_json,
                "upgrade": 1,
                "provenance": self.provenance(),
                "hidden": 0,
                "workspace": self.ws_name,
            }
        )
        self.obj_created.append(
            {"ref": self.create_ref(objid, self.ws_name), "description": ""}
        )
        return save_output["info"]

    def save_report_to_kbase(
        self, height=700, message="", warnings=[], file_links=[], summary_height=None
    ):
        """Save a report to KBase with HTML links and file links."""
        if not self.working_dir:
            raise ValueError("Working directory is not set")

        # Prepare HTML files
        os.makedirs(self.working_dir + "/html", exist_ok=True)
        rootDir = self.working_dir + "/html/"
        files = [
            {
                "path": "/kb/module/work/tmp/html/",
                "name": "index.html",
                "description": "HTML report",
            }
        ]
        for dirName, subdirList, fileList in os.walk(rootDir):
            for fname in fileList:
                if fname != "index.html":
                    files.append(
                        {
                            "path": dirName.replace(
                                rootDir, "/kb/module/work/tmp/html/"
                            ),
                            "name": fname,
                            "description": "Files related to HTML report",
                        }
                    )
        report_name = self.method + "-" + str(uuid.uuid4())
        output = self.report_client().create_extended_report(
            {
                "message": message,
                "warnings": warnings,
                "html_links": files,
                "file_links": file_links,
                "direct_html_link_index": 0,
                "html_window_height": height,
                "objects_created": self.obj_created,
                "workspace_name": self.ws_name,
                "report_object_name": report_name,
                "summary_window_height": summary_height,
            }
        )
        return {
            "report_name": report_name,
            "report_ref": output["ref"],
            "workspace_name": self.ws_name,
        }
