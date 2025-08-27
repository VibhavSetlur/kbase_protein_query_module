"""KBase genome utilities for working with genomic data and annotations."""

from typing import Any, Dict, Optional

from .base_utils import BaseUtils


class KBGenomeUtils(BaseUtils):
    """Utilities for working with KBase genome objects and genomic data.

    Provides methods for genome manipulation, feature extraction, sequence
    analysis, and other genome-specific operations in the KBase environment.
    """

    def __init__(self, **kwargs: Any) -> None:
        """Initialize KBase genome utilities.

        Args:
            **kwargs: Additional keyword arguments passed to BaseUtils
        """
        super().__init__(**kwargs)
        self.genetic_code = self._load_genetic_code()

    def object_to_proteins(self, ref):
        output = self.get_object(ref, self.ws_id)
        self.object_info_hash[ref] = output["info"]
        sequence_list = []
        # TODO: add support for other object types
        if "features" in output["data"]:
            for ftr in output["data"]["features"]:
                if "protein_translation" in ftr:
                    if len(ftr["protein_translation"]) > 0:
                        sequence_list.append([ftr["id"], ftr["protein_translation"]])
        return sequence_list

    def add_annotations_to_object(self, reference, suffix, annotations):
        """Loads specified gene annotation into KBase genome object

        Parameters
        ----------
        string - genome_ref
            KBase workspace reference to genome where annotations should be saved
        string - suffix
            Suffix to be used when saving modified genome back to KBase
        mapping<string gene_id,mapping<string ontology,mapping<string term,{"type":string,"score":float}>>> - annotations
            Annotations to be saved to genome

        Returns:
        -------
        dict

        Raises:
        ------
        """
        ontology_inputs = {}
        for geneid in annotations:
            for ontology in annotations[geneid]:
                if ontology not in ontology_inputs:
                    ontology_inputs[ontology] = {}
                if geneid not in ontology_inputs[ontology]:
                    ontology_inputs[ontology][geneid] = []
                for term in annotations[geneid][ontology]:
                    anno_data = {"term": term}
                    if "scores" in annotations[geneid][ontology][term]:
                        anno_data["evidence"] = {
                            "scores": annotations[geneid][ontology][term]["scores"]
                        }
                    if "name" in annotations[geneid][ontology][term]:
                        anno_data["name"] = (
                            annotations[geneid][ontology][term]["name"] + suffix
                        )
                    ontology_inputs[ontology][geneid].append(anno_data)

        anno_api_input = {
            "input_ref": reference,
            "output_name": self.object_info_hash[reference][1] + suffix,
            "output_workspace": self.ws_id,
            "overwrite_matching": 1,
            "save": 1,
            "provenance": self.provenance(),
            "events": [],
        }
        for ontology in ontology_inputs.keys():
            anno_api_input["events"].append(
                {
                    "ontology_id": ontology,
                    "method": self.name + "." + self.method,
                    "method_version": self.version,
                    "timestamp": self.timestamp,
                    "ontology_terms": ontology_inputs[ontology],
                }
            )
        anno_api_output = self.anno_client().add_annotation_ontology_events(
            anno_api_input
        )
        self.obj_created.append(
            {
                "ref": anno_api_output["output_ref"],
                "description": "Saving annotation for "
                + self.object_info_hash[reference][1],
            }
        )
        return anno_api_output

    def _load_genetic_code(self) -> Dict[str, str]:
        """Load the standard genetic code for translation.

        Returns:
            Dictionary mapping codons to amino acids
        """
        return {
            "TTT": "F",
            "TTC": "F",
            "TTA": "L",
            "TTG": "L",
            "TCT": "S",
            "TCC": "S",
            "TCA": "S",
            "TCG": "S",
            "TAT": "Y",
            "TAC": "Y",
            "TAA": "*",
            "TAG": "*",
            "TGT": "C",
            "TGC": "C",
            "TGA": "*",
            "TGG": "W",
            "CTT": "L",
            "CTC": "L",
            "CTA": "L",
            "CTG": "L",
            "CCT": "P",
            "CCC": "P",
            "CCA": "P",
            "CCG": "P",
            "CAT": "H",
            "CAC": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGT": "R",
            "CGC": "R",
            "CGA": "R",
            "CGG": "R",
            "ATT": "I",
            "ATC": "I",
            "ATA": "I",
            "ATG": "M",
            "ACT": "T",
            "ACC": "T",
            "ACA": "T",
            "ACG": "T",
            "AAT": "N",
            "AAC": "N",
            "AAA": "K",
            "AAG": "K",
            "AGT": "S",
            "AGC": "S",
            "AGA": "R",
            "AGG": "R",
            "GTT": "V",
            "GTC": "V",
            "GTA": "V",
            "GTG": "V",
            "GCT": "A",
            "GCC": "A",
            "GCA": "A",
            "GCG": "A",
            "GAT": "D",
            "GAC": "D",
            "GAA": "E",
            "GAG": "E",
            "GGT": "G",
            "GGC": "G",
            "GGA": "G",
            "GGG": "G",
        }

    def reverse_complement(self, dna_sequence: str) -> str:
        """Get the reverse complement of a DNA sequence.

        Args:
            dna_sequence: DNA sequence string

        Returns:
            Reverse complement sequence
        """
        complement_map = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

        # Handle lowercase and other characters
        complement_map.update({k.lower(): v.lower() for k, v in complement_map.items()})

        complement = "".join(complement_map.get(base, base) for base in dna_sequence)
        return complement[::-1]  # Reverse the string

    def translate_sequence(
        self, dna_sequence: str, genetic_code: Optional[int] = None
    ) -> str:
        """Translate a DNA sequence to amino acids.

        Args:
            dna_sequence: DNA sequence to translate
            genetic_code: Genetic code table number (defaults to standard code)

        Returns:
            Translated amino acid sequence
        """
        # Use standard genetic code for now (could be extended for other codes)
        if genetic_code and genetic_code != 11:
            self.log_warning(
                f"Genetic code {genetic_code} not implemented, using standard code"
            )

        # Ensure sequence length is multiple of 3
        seq_len = len(dna_sequence)
        if seq_len % 3 != 0:
            self.log_warning(
                f"Sequence length {seq_len} not divisible by 3, truncating"
            )
            dna_sequence = dna_sequence[: seq_len - (seq_len % 3)]

        # Translate codons
        amino_acids = []
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i : i + 3].upper()
            amino_acid = self.genetic_code.get(codon, "X")  # 'X' for unknown
            amino_acids.append(amino_acid)

        return "".join(amino_acids)

    def calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content of a DNA sequence.

        Args:
            sequence: DNA sequence string

        Returns:
            GC content as a fraction (0.0 to 1.0)
        """
        if not sequence:
            return 0.0

        sequence = sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        total_bases = len([base for base in sequence if base in "ATGC"])

        return gc_count / total_bases if total_bases > 0 else 0.0
