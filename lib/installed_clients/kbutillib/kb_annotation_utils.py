"""Implementation of KBase Annotation Ontology API within these utilities, which enables its use outside the SDK"""

import hashlib
import json
import os
import re
from typing import Any

import pandas as pd

from .kb_ws_utils import KBWSUtils

source_hash = {"MetaCyc": "META", "KEGG": "RO", "BiGG": "BIGG", "Rhea": "RHEA"}

ontology_translation = {
    "KEGGKO": "KO",
    "KEGGRO": "RO",
    "METACYC": "META",
    "SEED": "SSO",
    "TCDB": "TC",
    "MODELSEED": "MSRXN",
    "TIGRFAM": "TIGR",
}

ontology_hash = {
    "antiSMASH-CV": 1,
    "KO": 1,
    "EC": 1,
    "SSO": 1,
    "RO": 1,
    "META": 1,
    "MSRXN": 1,
    "MSCPD": 1,
    "MSCPX": 1,
    "BIGG": 1,
    "BIGGCPD": 1,
    "GO": 1,
    "TC": 1,
    "RHEA": 1,
    "TIGR": 1,
    "PF": 1,
    "PTHR": 1,
}


class KBAnnotationUtils(KBWSUtils):
    """Utilities for managing gene annotation information within KBase genome objects in the workspace"""

    def __init__(self, **kwargs: Any) -> None:
        """Initialize KBase annotation utilities."""
        super().__init__(**kwargs)

        # Ensure required dependencies are available
        self._ensure_dependencies()
        self.annoontology_dir = self.dependencies_dir / "cb_annotation_ontology_api"

        self.object = None
        self.objectinfo = None
        self.type = None
        self.ref = None
        self.msrxn_filter = True
        self.eventarray = []
        self.ftrhash = {}
        self.ftrtypes = {}
        self.alias_hash = {}
        self.object_alias_hash = {}
        self.term_names = {}
        self.ontologies_present = {}
        # Loading filtered reactions
        self.filtered_rxn = {}
        filename = self.annoontology_dir / "data/FilteredReactions.csv"
        filtered_reaction_df = pd.read_csv(filename, sep="\t")
        for index, row in filtered_reaction_df.iterrows():
            self.filtered_rxn[row["id"]] = row["reason"]

    def _ensure_dependencies(self) -> None:
        """Ensure all required dependencies are available."""
        dependencies_ok = True

        if not self.ensure_annotation_ontology_api():
            self.log_error("Failed to obtain Annotation Ontology API dependency")
            dependencies_ok = False

        if not dependencies_ok:
            raise ImportError(
                "Required dependencies for KBAnnotationUtils are not available"
            )

    def reset_attributes(self):
        """Resetting annotation elements related to a new method call."""
        super().reset_attributes()
        self.object = None
        self.objectinfo = None
        self.type = None
        self.ref = None
        self.msrxn_filter = True
        self.eventarray = []
        self.ftrhash = {}
        self.ftrtypes = {}
        self.alias_hash = {}
        self.object_alias_hash = {}
        self.term_names = {}
        self.ontologies_present = {}

    def initialize_call(
        self, method, params, print_params=False, no_print=[], no_prov_params=[]
    ):
        self.reset_attributes()
        super().initialize_call(method, params, print_params, no_print, no_prov_params)

    def get_alias_hash(self, namespace):
        if "MSRXN" not in self.alias_hash:
            filename = self.annoontology_dir / "data/msrxn_hash.json"
            with open(filename) as json_file:
                self.alias_hash["MSRXN"] = json.load(json_file)
        if namespace not in self.alias_hash:
            self.alias_hash[namespace] = {}
            if namespace == "EC":
                filename = self.annoontology_dir / "data/EC_translation.tsv"
                data = ""
                with open(filename) as file:
                    data = file.read()
                lines = data.split("\n")
                lines.pop(0)
                for line in lines:
                    items = line.split("\t")
                    if len(items) >= 2:
                        modelseed = "MSRXN:" + items[0]
                        if modelseed in self.alias_hash["MSRXN"]:
                            modelseed = self.alias_hash["MSRXN"][modelseed][0]
                        if "EC:" + items[1] not in self.alias_hash["EC"]:
                            self.alias_hash["EC"]["EC:" + items[1]] = []
                        self.alias_hash["EC"]["EC:" + items[1]].append(modelseed)
            elif (
                namespace == "META"
                or namespace == "RO"
                or namespace == "BIGG"
                or namespace == "RHEA"
            ):
                filename = self.annoontology_dir / "data/ModelSEED_Reaction_Aliases.txt"
                data = ""
                with open(filename) as file:
                    data = file.read()
                lines = data.split("\n")
                for line in lines:
                    items = line.split("\t")
                    if len(items) >= 3:
                        modelseed = "MSRXN:" + items[0]
                        if modelseed in self.alias_hash["MSRXN"]:
                            modelseed = self.alias_hash["MSRXN"][modelseed][0]
                        source = None
                        if items[2] in source_hash:
                            source = source_hash[items[2]]
                        if source != None:
                            if source not in self.alias_hash:
                                self.alias_hash[source] = {}
                            if items[1] not in self.alias_hash[source]:
                                self.alias_hash[source][source + ":" + items[1]] = []
                            self.alias_hash[source][source + ":" + items[1]].append(
                                modelseed
                            )
            elif namespace == "KO":
                filename = self.annoontology_dir / "data/kegg_95_0_ko_seed.tsv"
                data = ""
                with open(filename) as file:
                    data = file.read()
                lines = data.split("\n")
                lines.pop(0)
                for line in lines:
                    items = line.split("\t")
                    if len(items) >= 2:
                        if items[0] not in self.alias_hash["KO"]:
                            self.alias_hash["KO"]["KO:" + items[0]] = []
                        modelseed_ids = items[1].split(";")
                        id_hash = {}
                        for modelseed in modelseed_ids:
                            modelseed = "MSRXN:" + modelseed
                            if modelseed in self.alias_hash["MSRXN"]:
                                modelseed = self.alias_hash["MSRXN"][modelseed][0]
                            if modelseed not in id_hash:
                                self.alias_hash["KO"]["KO:" + items[0]].append(
                                    modelseed
                                )
                            id_hash[modelseed] = 1
            elif namespace == "SSO":
                sso_template = dict()
                filename = self.annoontology_dir / "data/SSO_reactions.json"
                with open(filename) as json_file:
                    sso_template = json.load(json_file)
                for sso in sso_template:
                    id_hash = {}
                    for modelseed in sso_template[sso]:
                        modelseed = "MSRXN:" + modelseed
                        if modelseed in self.alias_hash["MSRXN"]:
                            modelseed = self.alias_hash["MSRXN"][modelseed][0]
                        if modelseed not in id_hash:
                            if sso not in self.alias_hash["SSO"]:
                                self.alias_hash["SSO"][sso] = []
                            self.alias_hash["SSO"][sso].append(modelseed)
                        id_hash[modelseed] = 1
            elif namespace == "GO":
                go_translation = dict()
                filename = self.annoontology_dir / "data/GO_ontology_translation.json"
                with open(filename) as json_file:
                    go_translation = json.load(json_file)
                for term in go_translation["translation"]:
                    adjusted_term = "GO:" + term
                    if "equiv_terms" in go_translation["translation"][term]:
                        id_hash = {}
                        for rxn_data in go_translation["translation"][term][
                            "equiv_terms"
                        ]:
                            if rxn_data["equiv_term"] != None:
                                modelseed = "MSRXN:" + rxn_data["equiv_term"]
                                if modelseed in self.alias_hash["MSRXN"]:
                                    modelseed = self.alias_hash["MSRXN"][modelseed][0]
                                if adjusted_term not in self.alias_hash["GO"]:
                                    self.alias_hash["GO"][adjusted_term] = []
                                if modelseed not in id_hash:
                                    self.alias_hash["GO"][adjusted_term].append(
                                        modelseed
                                    )
                                id_hash[modelseed] = 1
        return self.alias_hash[namespace]

    def translate_term_to_modelseed(self, term):
        namespace = term.split(":").pop(0)
        output = []
        if namespace == "MSRXN":
            if term not in self.get_alias_hash(namespace):
                output = [term]
            else:
                output = self.get_alias_hash(namespace)[term]
        elif term not in self.get_alias_hash(namespace):
            output = []
        else:
            output = self.get_alias_hash(namespace)[term]
        if self.msrxn_filter:
            new_output = []
            for item in output:
                search = item
                if search[0:6] == "MSRXN:":
                    search = search[6:]
                if search not in self.filtered_rxn:
                    new_output.append(item)
            return new_output
        return output

    def get_annotation_ontology_events(self, params):
        self.initialize_call("get_annotation_ontology_events", params, True)
        params = self.validate_args(
            params,
            [],
            {
                "propagate": True,
                "query_events": None,
                "query_genes": None,
                "object": None,
                "type": None,
                "input_ref": None,
                "filter_msrxn_for_modeling": True,
                "input_workspace": None,
            },
        )
        self.msrxn_filter = params["filter_msrxn_for_modeling"]
        self.process_object(params)
        event_query = None
        if params["query_events"]:
            event_query = {}
            for event in params["query_events"]:
                event_query[event] = 1
        gene_query = None
        if params["query_genes"]:
            gene_query = {}
            for gene in params["query_genes"]:
                gene_query[gene] = 1
        output = {"events": self.eventarray, "feature_types": {}}
        for id in self.ftrhash:
            feature = self.ftrhash[id]
            if gene_query == None or id in gene_query:
                if "ontology_terms" in feature:
                    output["feature_types"][id] = self.ftrtypes[id]
                    self.integrate_terms_from_ftr(id, feature)
                if "cdss" in feature:
                    for cds in feature["cdss"]:
                        subfeature = self.ftrhash[cds]
                        if "ontology_terms" in subfeature:
                            self.integrate_terms_from_ftr(id, subfeature)
                if "parent_gene" in feature and feature["parent_gene"] in self.ftrhash:
                    subfeature = self.ftrhash[feature["parent_gene"]]
                    if "ontology_terms" in subfeature:
                        self.integrate_terms_from_ftr(id, subfeature)
        return output

    def add_annotation_ontology_events(self, params):
        self.initialize_call(
            "add_annotation_ontology_events",
            params,
            True,
            no_print=["events", "object"],
            no_prov_params=["object", "events", "feature_object"],
        )
        params = self.validate_args(
            params,
            ["output_workspace", "events"],
            {
                "provenance": [],
                "overwrite_matching": True,
                "object": None,
                "type": None,
                "input_ref": None,
                "input_workspace": None,
                "output_name": None,
                "save": 1,
            },
        )
        self.process_object(params)
        if not params["output_name"]:
            if self.objectinfo:
                params["output_name"] = self.objectinfo[1]
            else:
                params["output_name"] = self.object["id"]
        if "clear_existing" in params and params["clear_existing"] == 1:
            self.eventarray = []
        output = {"ftrs_not_found": [], "ftrs_found": 0, "terms_not_found": []}
        # Scrolling through new events, stadardizing, and checking for matches
        new_events = []
        for event in params["events"]:
            new_event = self.standardize_event(event)
            new_events.append(new_event)
            match = False
            for i, existing_event in enumerate(self.eventarray):
                # If an existing event has a matching event ID, we overwrite it
                if existing_event["event_id"] == new_event["event_id"]:
                    match = True
                    if params["overwrite_matching"]:
                        self.eventarray[i] = new_event
            if not match:
                self.eventarray.append(new_event)
        # Adding events
        feature_found_hash = {}
        terms_not_found = {}
        for event in new_events:
            for currgene in event["ontology_terms"]:
                genes = []
                if currgene in self.ftrhash:
                    feature_found_hash[currgene] = 1
                    genes = [self.ftrhash[currgene]["id"]]
                elif currgene in self.object_alias_hash:
                    feature_found_hash[currgene] = 1
                    genes = self.object_alias_hash[currgene]
                else:
                    output["ftrs_not_found"].append(currgene)
                for gene in genes:
                    if gene in self.ftrhash:
                        feature = self.ftrhash[gene]
                        add_ftr_output = self.add_feature_ontology_terms(
                            feature, event, currgene
                        )
                        for term in add_ftr_output["terms_not_found"]:
                            terms_not_found[term] = 1
        output["ftrs_found"] = len(feature_found_hash)
        for term in terms_not_found:
            output["terms_not_found"].append(term)
        # Saving object if requested but not if it's an AMA
        if params["save"] == 1:
            save_output = self.save_object(params)
            for key in save_output:
                output[key] = save_output[key]
        else:
            # Returning object if save not requested
            output["object"] = self.object
            output["type"] = self.type
            if "feature_object" in params:
                output["feature_object"] = params["feature_object"]
        return output

    def process_object(self, params):
        if params.get("object"):
            self.object = params["object"]
            self.type = params["type"]
        else:
            res = None
            if "input_workspace" not in params:
                res = self.ws_client().get_objects2(
                    {"objects": [self.process_ws_ids(params["input_ref"], None)]}
                )
            else:
                res = self.ws_client().get_objects2(
                    {
                        "objects": [
                            self.process_ws_ids(
                                params["input_ref"], params["input_workspace"]
                            )
                        ]
                    }
                )
            self.object = res["data"][0]["data"]
            self.objectinfo = res["data"][0]["info"]
            self.type = res["data"][0]["info"][2]
            self.ref = (
                str(res["data"][0]["info"][6])
                + "/"
                + str(res["data"][0]["info"][0])
                + "/"
                + str(res["data"][0]["info"][4])
            )
        self.eventarray = []
        self.ontologies_present = {}
        if "ontology_events" in self.object:
            for event in self.object["ontology_events"]:
                newevent = self.standardize_event(event)
                self.eventarray.append(newevent)
        self.ftrhash = {}
        self.ftrtypes = {}
        self.object_alias_hash = {}
        if "features" in self.object:
            to_remove = []
            for ftr in self.object["features"]:
                if "protein_translation" not in ftr:
                    if "non_coding_features" not in self.object:
                        self.object["non_coding_features"] = []
                    self.object["non_coding_features"].append(ftr)
                    to_remove.append(ftr)
                else:
                    self.ftrhash[ftr["id"]] = ftr
                    self.ftrtypes[ftr["id"]] = "gene"
            for item in to_remove:
                self.object["features"].remove(item)
        if "sequences" in self.object:
            for ftr in self.object["sequences"]:
                self.ftrhash[ftr["id"]] = ftr
                self.ftrtypes[ftr["id"]] = "sequence"
        if "cdss" in self.object:
            for ftr in self.object["cdss"]:
                self.ftrhash[ftr["id"]] = ftr
                self.ftrtypes[ftr["id"]] = "cds"
        # if "mrnas" in self.object:
        #    for ftr in self.object["mrnas"]:
        #        self.ftrhash[ftr["id"]] = ftr
        #        self.ftrtypes[ftr["id"]] = "mrna"
        if "non_coding_features" in self.object:
            for ftr in self.object["non_coding_features"]:
                self.ftrhash[ftr["id"]] = ftr
                self.ftrtypes[ftr["id"]] = "noncoding"
        if "features_handle_ref" in self.object:
            if "feature_object" not in params:
                shock_output = self.dfu_client().shock_to_file(
                    {
                        "handle_id": self.object["features_handle_ref"],
                        "file_path": self.get_config("AnnoOntAPI", "scratch"),
                    }
                )
                os.system("gunzip --force " + shock_output["file_path"])
                shock_output["file_path"] = shock_output["file_path"][0:-3]
                with open(shock_output["file_path"]) as json_file:
                    features = json.load(json_file)
            else:
                features = params["feature_object"]
            for ftr in features:
                self.ftrhash[ftr["id"]] = ftr
                self.ftrtypes[ftr["id"]] = "gene"
        for ftrid in self.ftrhash:
            ftr = self.ftrhash[ftrid]
            if "ontology_terms" not in ftr:
                ftr["ontology_terms"] = {}
            self.upgrade_feature(ftr)
            self.process_feature_aliases(ftr)

    def save_object(self, params):
        # All objects undergo these steps
        # Setting ontology and eventarray in object
        self.object["ontologies_present"] = self.ontologies_present
        self.object["ontology_events"] = self.eventarray
        for event in self.object["ontology_events"]:
            if "ontology_terms" in event:
                event.pop("ontology_terms")

        # If a metagenome, saving features
        if self.type.startswith("KBaseMetagenomes.AnnotatedMetagenomeAssembly"):
            if "feature_object" not in params:
                self.log_critical(
                    "feature_object must exist in order to save the altered metagenome object!"
                )

            json_file_path = (
                self.get_config("AnnoOntAPI", "scratch")
                + self.object["name"]
                + "_features.json"
            )
            with open(json_file_path, "w") as fid:
                json.dump(params["feature_object"], fid)
            json_to_shock = self.dfu_client().file_to_shock(
                {"file_path": json_file_path, "make_handle": 1, "pack": "gzip"}
            )
            # Resetting feature file handle o new value
            self.object["features_handle_ref"] = json_to_shock["handle"]["hid"]
            # Remove json file to avoid disk overload
            os.remove(json_file_path)
            # Removing genbank handle ref because this breaks saving
            self.object.pop("genbank_handle_ref", None)
        elif self.type.startswith("KBaseGenomes.Genome"):
            # Removing genbank handle ref because this breaks saving
            self.object.pop("genbank_handle_ref", None)
            self.check_genome(self.object, self.ref)
        elif self.type.startswith(
            "KBaseSequences.ProteinSequenceSet"
        ) or self.type.startswith("KBaseSequences.DNASequenceSet"):
            pass  # No specific instructions for handling these types yet
        if self.type.startswith("KBaseGenomes.Genome") or self.type.startswith(
            "KBaseMetagenomes.AnnotatedMetagenomeAssembly"
        ):
            save_output = self.save_genome_or_metagenome(
                params["output_name"], params["output_workspace"], self.object
            )
        elif self.type.startswith(
            "KBaseSequences.ProteinSequenceSet"
        ) or self.type.startswith("KBaseSequences.DNASequenceSet"):
            save_output = self.save_ws_object(
                params["output_name"],
                params["output_workspace"],
                self.object,
                self.type,
            )[0]
        output = {}
        output["output_ref"] = self.wsinfo_to_ref(save_output)
        output["output_name"] = str(save_output[1])
        return output

    # Function to standardize ontology tags
    def clean_tag(self, original_tag):
        tag = original_tag.upper()
        if tag not in ontology_hash and tag in ontology_translation:
            tag = ontology_translation[tag]
        return tag

    def clean_term(self, original_term, original_tag, tag):
        term = original_term
        array = term.split(":")
        if len(array) == 1:
            term = tag + ":" + array[0]
        else:
            if array[0].upper() == original_tag.upper() or array[0].upper() == tag:
                array[0] = tag
                term = ":".join(array)
            else:
                term = tag + ":" + ":".join(array)
        return term

    def standardize_event(self, event):
        if "id" not in event and "ontology_id" in event:
            event["id"] = event["ontology_id"]
        if "event_id" not in event:
            event["event_id"] = (
                event["method"]
                + ":"
                + event["method_version"]
                + ":"
                + event["id"]
                + ":"
                + event["timestamp"]
            )
        old_description = None
        if "description" in event:
            old_description = event["description"]
            if (
                event["description"][-1 * len(event["timestamp"]) :]
                != event["timestamp"]
            ):
                event["description"] = event["description"] + ":" + event["timestamp"]
        else:
            event["description"] = (
                event["method"]
                + ":"
                + event["method_version"]
                + ":"
                + event["id"]
                + ":"
                + event["timestamp"]
            )
        if "ontology_id" not in event:
            event["ontology_id"] = event["id"]
        standard_event = {
            "id": self.clean_tag(event["ontology_id"]),
            "event_id": event["event_id"],
            "original_description": old_description,
            "description": event["description"],
            "ontology_id": self.clean_tag(event["ontology_id"]),
            "method": event["method"],
            "method_version": event["method_version"],
            "timestamp": event["timestamp"],
        }
        if "ontology_terms" in event:
            standard_event["ontology_terms"] = event["ontology_terms"]
        return standard_event

    def process_feature_aliases(self, ftr):
        if "aliases" in ftr:
            for alias in ftr["aliases"]:
                if not isinstance(alias, str):
                    alias = alias[1]
                if alias not in self.object_alias_hash:
                    self.object_alias_hash[alias] = []
                self.object_alias_hash[alias].append(ftr["id"])
                if alias.lower() not in self.object_alias_hash:
                    self.object_alias_hash[alias.lower()] = []
                if ftr["id"] not in self.object_alias_hash[alias.lower()]:
                    self.object_alias_hash[alias.lower()].append(ftr["id"])
        if "db_xrefs" in ftr:
            for alias in ftr["db_xrefs"]:
                if alias[1] not in self.object_alias_hash:
                    self.object_alias_hash[alias[1]] = []
                self.object_alias_hash[alias[1]].append(ftr["id"])
                if alias[1].lower() not in self.object_alias_hash:
                    self.object_alias_hash[alias[1].lower()] = []
                if ftr["id"] not in self.object_alias_hash[alias[1].lower()]:
                    self.object_alias_hash[alias[1].lower()].append(ftr["id"])

    def upgrade_feature(self, ftr):
        type = self.ftrtypes[ftr["id"]]
        if "function" in ftr:
            ftr["functions"] = re.split(r"\s*;\s+|\s+[\@\/]\s+", ftr["function"])
            del ftr["function"]
            # Clearing old ontology terms rather than attempting to translate them
            if "ontology_terms" not in ftr:
                ftr["ontology_terms"] = {}
        if type == "gene" and "cdss" not in ftr:
            ftr["cdss"] = []
        if "dna_sequence_length" not in ftr:
            if "location" in ftr:
                ftr["dna_sequence_length"] = ftr["location"][0][3]
            elif "dna_sequence" in ftr:
                ftr["dna_sequence_length"] = len(ftr["dna_sequence"])
            elif "protein_translation" in ftr:
                ftr["dna_sequence_length"] = len(ftr["protein_translation"]) * 3
        if "aliases" in ftr:
            for i in range(0, len(ftr["aliases"])):
                if isinstance(ftr["aliases"][i], str):
                    array = ftr["aliases"][i].split(":")
                    if len(array) > 1:
                        ftr["aliases"][i] = array
                    else:
                        ftr["aliases"][i] = ["Unknown", ftr["aliases"][i]]
        if type == "cdss" and "protein_md5" not in ftr:
            ftr["protein_md5"] = hashlib.md5(
                ftr["protein_translation"].encode()
            ).hexdigest()
        if "md5" not in ftr:
            if "dna_sequence" in ftr:
                ftr["md5"] = hashlib.md5(ftr["dna_sequence"].encode()).hexdigest()
            elif "protein_translation" in ftr:
                ftr["md5"] = hashlib.md5(
                    ftr["protein_translation"].encode()
                ).hexdigest()
            else:
                ftr["md5"] = ""

    def integrate_terms_from_ftr(self, id, feature):
        annotation = False
        if "ontology_terms" in feature:
            for original_tag in feature["ontology_terms"]:
                for original_term in feature["ontology_terms"][original_tag]:
                    tag = self.clean_tag(original_tag)
                    term = self.clean_term(original_term, original_tag, tag)
                    modelseed_ids = self.translate_term_to_modelseed(term)
                    for j, event_index in enumerate(
                        feature["ontology_terms"][original_tag][original_term]
                    ):
                        if event_index == None:
                            for i, event in enumerate(self.eventarray):
                                if event["ontology_id"] == tag:
                                    event_index = i
                                    feature["ontology_terms"][original_tag][
                                        original_term
                                    ][j] = i
                                    break
                        if event_index < len(self.eventarray):
                            if "ontology_terms" not in self.eventarray[event_index]:
                                self.eventarray[event_index]["ontology_terms"] = {}
                            if id not in self.eventarray[event_index]["ontology_terms"]:
                                self.eventarray[event_index]["ontology_terms"][id] = []
                            termdata = {"term": term}
                            if id != feature["id"]:
                                termdata["indirect"] = True
                            if len(modelseed_ids) > 0:
                                termdata["modelseed_ids"] = modelseed_ids
                            if "ontology_evidence" in feature:
                                if term in feature["ontology_evidence"]:
                                    if (
                                        str(event_index)
                                        in feature["ontology_evidence"][term]
                                    ):
                                        termdata["evidence"] = feature[
                                            "ontology_evidence"
                                        ][term][str(event_index)]
                            found = False
                            for outterm in self.eventarray[event_index][
                                "ontology_terms"
                            ][id]:
                                if outterm["term"] == term:
                                    found = True
                            if not found:
                                self.eventarray[event_index]["ontology_terms"][
                                    id
                                ].append(termdata)
                                annotation = True
        return annotation

    def add_feature_ontology_terms(self, feature, event, ftrid):
        term_data = event["ontology_terms"][ftrid]
        event_index = None
        for i, item in enumerate(self.eventarray):
            if item["event_id"] == event["event_id"]:
                event_index = i
        if event_index == None:
            print("Event not found!", ftrid, str(event), str(feature))
        output = {"terms_not_found": []}
        if "ontology_terms" not in feature:
            feature["ontology_terms"] = {}
        if event["ontology_id"] not in feature["ontology_terms"]:
            feature["ontology_terms"][event["ontology_id"]] = {}
        for term in term_data:
            if (
                "indirect" not in term or not term["indirect"]
            ):  # Don't add terms for proteins that might be present because of CDS annotations
                if term["term"].split(":")[0] != event["ontology_id"]:
                    term["term"] = event["ontology_id"] + ":" + term["term"]
                # If this is a SEED role, translate to an SSO
                if (
                    event["ontology_id"] == "SSO"
                    and re.search(r"^SSO:\d+$", term["term"]) == None
                ):
                    term["term"] = re.sub("^SSO:", "", term["term"])
                    terms = re.split(r"\s*;\s+|\s+[\@\/]\s+", term["term"])
                    first = 1
                    for subterm in terms:
                        if first == 1:
                            # Only the first term completes the rest of this code
                            term["term"] = self.translate_rast_function_to_sso(subterm)
                            first = 0
                        else:
                            # Subterms need to be added independently
                            subterm = self.translate_rast_function_to_sso(subterm)
                            if subterm != None:
                                if (
                                    subterm
                                    not in feature["ontology_terms"][
                                        event["ontology_id"]
                                    ]
                                ):
                                    feature["ontology_terms"][event["ontology_id"]][
                                        subterm
                                    ] = []
                                feature["ontology_terms"][event["ontology_id"]][
                                    subterm
                                ].append(event_index)
                                if event["ontology_id"] not in self.ontologies_present:
                                    self.ontologies_present[event["ontology_id"]] = {}
                                self.ontologies_present[event["ontology_id"]][
                                    subterm
                                ] = self.get_term_name(event["ontology_id"], subterm)
                                if (
                                    self.ontologies_present[event["ontology_id"]][
                                        subterm
                                    ]
                                    == "Unknown"
                                ):
                                    output["terms_not_found"].append(subterm)
                if term["term"] == None:
                    continue
                if event["ontology_id"] not in self.ontologies_present:
                    self.ontologies_present[event["ontology_id"]] = {}
                # Dealing with custom names when specified for a term
                if "name" in term or "suffix" in term:
                    if "suffix" in term:
                        term["name"] = (
                            self.get_term_name(event["ontology_id"], term["term"])
                            + term["suffix"]
                        )
                    # If a custom name is specified, for now, we need to make sure the term is unique or the name will be overwritten
                    if (
                        term["term"] in self.ontologies_present[event["ontology_id"]]
                        and self.ontologies_present[event["ontology_id"]][term["term"]]
                        != term["name"]
                    ):
                        index = 1
                        while (
                            term["term"] + ";" + str(index)
                            in self.ontologies_present[event["ontology_id"]]
                            and self.ontologies_present[event["ontology_id"]][
                                term["term"] + ";" + str(index)
                            ]
                            != term["name"]
                        ):
                            index += 1
                        term["term"] = term["term"] + ";" + str(index)
                    self.ontologies_present[event["ontology_id"]][term["term"]] = term[
                        "name"
                    ]
                else:
                    self.ontologies_present[event["ontology_id"]][term["term"]] = (
                        self.get_term_name(event["ontology_id"], term["term"])
                    )
                    if (
                        self.ontologies_present[event["ontology_id"]][term["term"]]
                        == "Unknown"
                    ):
                        output["terms_not_found"].append(term["term"])
                if term["term"] not in feature["ontology_terms"][event["ontology_id"]]:
                    feature["ontology_terms"][event["ontology_id"]][term["term"]] = []
                feature["ontology_terms"][event["ontology_id"]][term["term"]].append(
                    event_index
                )
                if "evidence" in term:
                    if "ontology_evidence" not in feature:
                        feature["ontology_evidence"] = {}
                    if term["term"] not in feature["ontology_evidence"]:
                        feature["ontology_evidence"][term["term"]] = {}
                    feature["ontology_evidence"][term["term"]][str(event_index)] = term[
                        "evidence"
                    ]
        return output

    def check_genome(self, genome, ref=None):
        if "gc_content" in genome and isinstance(genome["gc_content"], str):
            genome["gc_content"] = float(genome["gc_content"])
        if "md5" not in genome:
            genome["md5"] = ""
        if "molecule_type" not in genome:
            genome["molecule_type"] = "DNA"
        if "genome_tiers" not in genome:
            genome["genome_tiers"] = ["ExternalDB", "User"]
        if "feature_counts" not in genome:
            genome["feature_counts"] = {
                "CDS": 0,
                "gene": 0,
                "non-protein_encoding_gene": 0,
            }
            if "cdss" in genome:
                genome["feature_counts"]["CDS"] = len(genome["cdss"])
            if "features" in genome:
                genome["feature_counts"]["gene"] = len(genome["features"])
            if "non_coding_features" in genome:
                genome["feature_counts"]["non-protein_encoding_gene"] = len(
                    genome["non_coding_features"]
                )
        if "cdss" not in genome:
            genome["cdss"] = []
        if "assembly_ref" in genome and not ref == None:
            genome["assembly_ref"] = ref + ";" + genome["assembly_ref"]

    def convert_role_to_searchrole(self, term):
        term = term.lower()
        term = re.sub(r"\s", "", term)
        term = re.sub(r"[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]*", "", term)
        term = re.sub(r"\#.*$", "", term)
        term = re.sub(r"\(ec:*\)", "", term)
        term = re.sub(r"[\(\)\[\],-]", "", term)
        return term

    def translate_rast_function_to_sso(self, input_term):
        # Stripping out SSO prefix if it's present
        input_term = re.sub("^SSO:", "", input_term)
        input_term = self.convert_role_to_searchrole(input_term)
        # Checking for SSO translation file
        if "SEED_ROLE" not in self.alias_hash:
            self.alias_hash["SEED_ROLE"] = {}
            sso_ontology = dict()
            with open(self.annoontology_dir / "data/SSO_dictionary.json") as json_file:
                sso_ontology = json.load(json_file)
            for term in sso_ontology["term_hash"]:
                name = self.convert_role_to_searchrole(
                    sso_ontology["term_hash"][term]["name"]
                )
                self.alias_hash["SEED_ROLE"][name] = term

        # Translating
        if input_term in self.alias_hash["SEED_ROLE"]:
            return self.alias_hash["SEED_ROLE"][input_term]
        else:
            return None

    def get_term_name(self, type, term):
        if type not in self.term_names:
            self.term_names[type] = {}
            if type in [
                "SSO",
                "AntiSmash",
                "EC",
                "TC",
                "META",
                "RO",
                "KO",
                "GO",
                "TIGR",
                "PF",
                "PTHR",
            ]:
                with open(
                    self.annoontology_dir / "data" / f"{type}_dictionary.json"
                ) as json_file:
                    ontology = json.load(json_file)
                    for term in ontology["term_hash"]:
                        self.term_names[type][term] = ontology["term_hash"][term][
                            "name"
                        ]
        if term not in self.term_names[type]:
            return "Unknown"
        return self.term_names[type][term]
