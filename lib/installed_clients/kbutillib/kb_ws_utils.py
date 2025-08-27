"""KBase API utilities for interacting with KBase services and data."""

import json
import os
import re
import sys
import time
from typing import Any, Dict, List, Optional, Union

import requests

from .installed_clients.AbstractHandleClient import AbstractHandle as HandleService
from .installed_clients.WorkspaceClient import Workspace
from .shared_env_utils import SharedEnvUtils


class KBWSUtils(SharedEnvUtils):
    """Utilities for interacting with KBase (Department of Energy Systems Biology
    Knowledgebase) APIs and services.

    Provides methods for authentication, data retrieval, workspace operations,
    and other KBase-specific functionality.
    """

    def __init__(
        self, kb_version: Optional[str] = "prod", max_retry: int = 3, **kwargs: Any
    ) -> None:
        """Initialize KBase Workspace utilities."""
        super().__init__(**kwargs)
        self.kb_version = kb_version
        base_url = self.get_base_url_from_version(kb_version)
        self.workspace_url = f"{base_url}/ws"
        self.shock_url = f"{base_url}/shock-api"
        self.hs_url = f"{base_url}/handle_service"
        self.cached_to_obj_path = {}
        self._ws_client = Workspace(
            self.workspace_url, token=self.get_token(namespace="kbase")
        )
        self.max_retry = max_retry
        self.hs_client = HandleService(
            self.hs_url, token=self.get_token(namespace="kbase")
        )
        self.ws_id = None
        self.ws_name = None

    def reset_attributes(self):
        """Resetting workspace elements related to a new method call."""
        super().reset_attributes()
        self.ws_id = None
        self.ws_name = None

    def initialize_call(
        self,
        method: str,
        params: Dict[str, Any],
        print_params: bool = False,
        no_print: List[str] = None,
        no_prov_params: List[str] = None,
    ) -> None:
        """Initialize workspace elements related to a new method call."""
        if not self.initialized:
            # Set workspace if provided
            if "workspace" in params:
                self.set_ws(params["workspace"])
            elif "output_workspace" in params:
                self.set_ws(params["output_workspace"])
        super().initialize_call(method, params, print_params, no_print, no_prov_params)

    def ws_client(self) -> Workspace:
        """Get the Workspace client.

        Returns:
            Workspace client instance
        """
        return self._ws_client

    def set_ws(self, workspace):
        if self.ws_id == workspace or self.ws_name == workspace:
            return
        if not isinstance(workspace, str) or re.search("^\\d+$", workspace) != None:
            if isinstance(workspace, str):
                workspace = int(workspace)
            self.ws_id = workspace
            info = self.ws_client().get_workspace_info({"id": workspace})
            self.ws_name = info[1]
        else:
            self.ws_name = workspace
            info = self.ws_client().get_workspace_info({"workspace": workspace})
            self.ws_id = info[0]

    def get_base_url_from_version(self, version):
        if version == "prod":
            return "https://kbase.us/services"
        elif version == "appdev":
            return "https://appdev.kbase.us/services"
        elif version == "ci":
            return "https://ci.kbase.us/services"
        else:
            self.log_critical("Unknown workspace version: " + version)
            return "https://kbase.us/services"

    def download_blob_file(self, handle_id, file_path):
        headers = {"Authorization": "OAuth " + self.get_token(namespace="kbase")}
        hs = self.hs_client
        handles = hs.hids_to_handles([handle_id])
        shock_id = handles[0]["id"]
        node_url = self.shock_url + "/node/" + shock_id
        r = requests.get(node_url, headers=headers, allow_redirects=True)
        errtxt = ("Error downloading file from shock " + "node {}: ").format(shock_id)
        if not r.ok:
            print(json.loads(r.content)["error"][0])
            return None
        resp_obj = r.json()
        size = resp_obj["data"]["file"]["size"]
        if not size:
            print(f"Node {shock_id} has no file")
            return None
        node_file_name = resp_obj["data"]["file"]["name"]
        attributes = resp_obj["data"]["attributes"]
        # Making the directory if it doesn't exist
        dir = os.path.dirname(file_path)
        os.makedirs(dir, exist_ok=True)
        # Adding filename to the end of the directory
        if os.path.isdir(file_path):
            file_path = os.path.join(file_path, node_file_name)
        with open(file_path, "wb") as fhandle:
            with requests.get(
                node_url + "?download_raw",
                stream=True,
                headers=headers,
                allow_redirects=True,
            ) as r:
                if not r.ok:
                    print(json.loads(r.content)["error"][0])
                    return None
                for chunk in r.iter_content(1024):
                    if not chunk:
                        break
                    fhandle.write(chunk)
        return file_path

    def run_narrative_method(
        self,
        method_name: str,
        parameters: Dict[str, Any],
        workspace_id: Union[int, str],
    ) -> Dict[str, Any]:
        """Run a KBase narrative method (app).

        Args:
            method_name: Name of the method to run
            parameters: Method parameters
            workspace_id: Workspace ID where method should run

        Returns:
            Method execution result
        """
        url = f"{self.kbase_url}/narrative_method_store/rpc"

        data = {
            "version": "1.1",
            "method": "NarrativeMethodStore.run_method",
            "params": [
                {
                    "method": method_name,
                    "params": parameters,
                    "workspace_id": workspace_id,
                }
            ],
        }

        response = self._make_request("POST", url, data=data)
        result = response.json()

        if "error" in result:
            raise Exception(f"Narrative Method error: {result['error']}")

        self.log_info(f"Method '{method_name}' executed successfully")
        return result["result"][0]

    def save_ws_object(self, objid, workspace, obj_json, obj_type):
        self.set_ws(workspace)
        params = {
            "id": self.ws_id,
            "objects": [
                {
                    "data": obj_json,
                    "name": objid,
                    "type": obj_type,
                    "meta": {},
                    "provenance": self.get_provenance(),
                }
            ],
        }
        self.obj_created.append(
            {"ref": self.create_ref(objid, self.ws_name), "description": ""}
        )
        return self.ws_client().save_objects(params)

    def get_provenance(self) -> List[Dict[str, Any]]:
        """Generate provenance information for workspace operations.

        Returns:
            Provenance data list
        """
        return [
            {
                "description": self.method,
                "input_ws_objects": self.input_objects,
                "method": self.method,
                "script_command_line": "",
                "method_params": [self.params],
                "service": self.name,
                "service_ver": self.version,
            }
        ]

    def list_ws_objects(self, wsid_or_ref, type=None, include_metadata=True):
        """List objects in a workspace"""
        ws_client = self.ws_client()
        done = False
        skip = 0
        full_output = {}
        start_after = None
        while not done:
            input = {}
            if type:
                input["type"] = type
            if include_metadata:
                input["includeMetadata"] = 1
            else:
                input["includeMetadata"] = 0
            if isinstance(wsid_or_ref, int):
                input["ids"] = [wsid_or_ref]
                wsid_or_ref = str(wsid_or_ref)
            else:
                input["workspaces"] = [wsid_or_ref]

            if start_after:
                input["startafter"] = start_after
            output = ws_client.list_objects(input)
            start_after = wsid_or_ref + "/" + str(output[-1][0])
            for item in output:
                full_output[item[1]] = item
            if len(output) < 5000:
                done = True
        return full_output

    def process_ws_ids(self, id_or_ref, workspace=None, no_ref=False):
        """IDs should always be processed through this function so we can interchangeably use
        refs, IDs, and names for workspaces and objects
        """
        objspec = {}
        if len(id_or_ref.split(";")) > 1:
            objspec["to_obj_ref_path"] = id_or_ref.split(";")[0:-1]
            id_or_ref = id_or_ref.split(";")[-1]

        if len(id_or_ref.split("/")) > 1:
            if no_ref:
                array = id_or_ref.split("/")
                workspace = array[0]
                id_or_ref = array[1]
            else:
                objspec["ref"] = id_or_ref

        if "ref" not in objspec:
            if isinstance(workspace, int):
                objspec["wsid"] = workspace
            else:
                objspec["workspace"] = workspace
            if isinstance(id_or_ref, int):
                objspec["objid"] = id_or_ref
            else:
                objspec["name"] = id_or_ref
        return objspec

    def get_object_info(self, id_or_ref, ws=None):
        ws_identities = [self.process_ws_ids(id_or_ref, ws)]
        return self.ws_client().get_object_info(ws_identities, 1)[0]

    def get_object(self, id_or_ref, ws=None):
        res = self.ws_get_objects({"objects": [self.process_ws_ids(id_or_ref, ws)]})
        if res is None:
            return None
        return res["data"][0]

    def ws_get_objects(self, args):
        """All functions calling get_objects2 should call this function to ensure they get the retry
        code because workspace periodically times out
        :param args:
        :return:
        """
        tries = 0
        while tries < self.max_retry:
            try:
                return self.ws_client().get_objects2(args)
            except Exception as e:
                self.log_warning(
                    "Workspace get_objects2 call failed [%s:%s - %s]. Trying again! Error: %s"
                    % (sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2], str(e))
                )
                tries += 1
                time.sleep(10)  # Give half second
        self.log_warning(
            "get_objects2 failed after multiple tries: %s" % sys.exc_info()[0]
        )
        raise

    def wsinfo_to_ref(self, info):
        return str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])

    def create_ref(self, id_or_ref, ws=None):
        if isinstance(id_or_ref, int):
            id_or_ref = str(id_or_ref)
        if len(id_or_ref.split("/")) > 1:
            return id_or_ref
        if isinstance(ws, int):
            ws = str(ws)
        return ws + "/" + id_or_ref
