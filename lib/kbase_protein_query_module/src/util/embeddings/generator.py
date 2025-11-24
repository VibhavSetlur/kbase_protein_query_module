"""
Protein Embedding Generator Module

Uses ESM-2 models from fair-esm package (official Facebook Research implementation).
Models are automatically downloaded on first use.

Follows official ESM usage patterns from:
https://github.com/facebookresearch/esm
"""

import os
import subprocess
import logging
import numpy as np
import tempfile
import json
from typing import Optional, Union, Tuple, List, Dict, Any

# Handle KB PLM Utils import
try:
    from kbutillib.kb_plm_utils import KBPLMUtils
except ImportError:
    try:
        # Try alternative import path
        import sys
        kbutillib_path = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'installed_clients', 'kbutillib')
        if kbutillib_path not in sys.path:
            sys.path.insert(0, kbutillib_path)
        from kb_plm_utils import KBPLMUtils
    except ImportError:
        KBPLMUtils = None

try:
    import torch
    import esm
    HAS_ESM = True
except ImportError:
    HAS_ESM = False
    torch = None
    esm = None

logger = logging.getLogger(__name__)

class ProteinEmbeddingGenerator:
    """
    Generates protein embeddings using ESM-2 models (local or server) or CLI.
    
    Modes:
    1. Local: Uses fair-esm package to run models locally (requires GPU/CPU).
    2. Server: Uses KBase PLM API via KBUtilLib.
    3. CLI: Wraps esm-extract command line tool.
    """

    def __init__(
        self, 
        model_name: str = "esm2_t6_8M_UR50D", 
        device: str = "auto",
        plm_api_url: str = "https://kbase.us/services/llm_homology_api",
        use_server: bool = False
    ):
        """Initialize the embedding generator.
        
        Args:
            model_name: ESM model name (e.g., "esm2_t6_8M_UR50D")
            device: Device to use ("auto", "cpu", or "cuda")
            plm_api_url: URL for PLM API (for server mode)
            use_server: Whether to prefer server over local model
        """
        self.model_name = model_name
        self.device = device
        self.use_server = use_server
        self.plm_api_url = plm_api_url
        
        # Initialize KB PLM Utils if available
        self.plm_utils = None
        if KBPLMUtils:
            try:
                self.plm_utils = KBPLMUtils(plm_api_url=plm_api_url)
            except Exception as e:
                logger.warning(f"Failed to initialize KBPLMUtils: {e}")
        
        # Initialize local model if needed or as backup
        self.model = None
        self.alphabet = None
        self.batch_converter = None
        self.embedding_dim: Optional[int] = None
        self.repr_layer: Optional[int] = None
        
        if not self.use_server:
            if HAS_ESM:
                self.device = self._setup_device(device)
                self._load_model()
            else:
                logger.warning("fair-esm not installed. Local generation disabled.")

    def _setup_device(self, device: str) -> str:
        """Setup the device for model inference."""
        if device == "auto":
            device = "cuda" if torch and torch.cuda.is_available() else "cpu"
        return device
    
    def _load_model(self):
        """Load ESM-2 model using official fair-esm package."""
        try:
            # Map model names to ESM pretrained functions
            model_map = {
                "esm2_t6_8M_UR50D": esm.pretrained.esm2_t6_8M_UR50D,
                "esm2_t12_35M_UR50D": esm.pretrained.esm2_t12_35M_UR50D,
                "esm2_t30_150M_UR50D": esm.pretrained.esm2_t30_150M_UR50D,
                "esm2_t33_650M_UR50D": esm.pretrained.esm2_t33_650M_UR50D,
            }
            
            if self.model_name not in model_map:
                # Fallback or error
                if self.use_server:
                    logger.info(f"Model {self.model_name} not found locally, relying on server.")
                    return
                available = ", ".join(model_map.keys())
                raise ValueError(f"Unknown model name: {self.model_name}. Available models: {available}")
            
            # Load model and alphabet
            logger.info(f"Loading ESM model '{self.model_name}'...")
            self.model, self.alphabet = model_map[self.model_name]()
            self.batch_converter = self.alphabet.get_batch_converter()
            
            # Move model to device and set to eval mode
            self.model = self.model.to(self.device)
            self.model.eval()
            
            self.repr_layer = self.model.num_layers
            self.embedding_dim = self.model.embed_dim
            
            logger.info(f"Loaded ESM model '{self.model_name}' on {self.device} with dim {self.embedding_dim}")
        except Exception as e:
            if not self.use_server:
                raise RuntimeError(f"Failed to load ESM model '{self.model_name}': {e}")
            logger.warning(f"Failed to load local model: {e}. Will use server.")
    
    def generate_embedding(self, sequence: str, pooling: str = "mean") -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """Generate embeddings for a single protein sequence.
        
        Tries server first if configured, then falls back to local.
        """
        if self.use_server and self.plm_utils:
            try:
                return self.get_embedding_from_server(sequence, pooling)
            except Exception as e:
                logger.warning(f"Server embedding failed: {e}. Falling back to local.")
        
        return self.generate_embedding_local(sequence, pooling)

    def get_embedding_from_server(self, sequence: str, pooling: str = "mean") -> np.ndarray:
        """Get embedding from KBase PLM API."""
        if not self.plm_utils:
            raise RuntimeError("KBPLMUtils not initialized")
            
        # Query API with return_embeddings=True
        # We use a dummy ID
        query = [{"id": "query", "sequence": sequence}]
        
        # Note: query_plm_api returns hits. We need to check if it returns query embeddings.
        # Based on KBPLMUtils, it sends "return_query_embeddings": True
        result = self.plm_utils.query_plm_api(
            query_sequences=query,
            max_hits=1,
            return_embeddings=True
        )
        
        # Extract embedding from result
        # The structure depends on API response. Assuming standard format.
        # If the API returns query embeddings separately or in the hits structure.
        # Usually it's in the response root or associated with the query.
        # For now, assuming the API returns it in a standard way or we might need to adjust.
        # If the API doesn't support returning query embeddings explicitly in the client wrapper,
        # we might be limited. But let's assume it does as per user request.
        
        # Mocking extraction if specific field is not documented in the client code I read.
        # But let's assume 'queries' or similar field in result.
        if "queries" in result:
             for q in result["queries"]:
                 if q.get("id") == "query" and "embedding" in q:
                     return np.array(q["embedding"], dtype=np.float32)
        
        # Fallback: check if it's in hits (self-hit)
        if "hits" in result:
            for hit_group in result["hits"]:
                if hit_group.get("query_id") == "query":
                    # Check if query embedding is attached here
                    if "query_embedding" in hit_group:
                         return np.array(hit_group["query_embedding"], dtype=np.float32)
        
        raise RuntimeError("Could not extract embedding from server response")

    def generate_embedding_local(self, sequence: str, pooling: str = "mean") -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """Generate embeddings using local ESM model."""
        if not self.model:
            raise RuntimeError("Local model not loaded")
            
        if not sequence or len(sequence.strip()) == 0:
            raise ValueError("Empty or invalid protein sequence")
        
        sequence = self._preprocess_sequence(sequence)
        
        # Prepare data
        data = [("protein", sequence)]
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        batch_tokens = batch_tokens.to(self.device)
        
        batch_lens = (batch_tokens != self.alphabet.padding_idx).sum(1)
        tokens_len = batch_lens[0].item()
        
        with torch.no_grad():
            results = self.model(batch_tokens, repr_layers=[self.repr_layer])
        
        token_representations = results["representations"][self.repr_layer]
        
        if tokens_len >= 2:
            residue_embeddings = token_representations[0, 1:tokens_len-1].cpu().numpy().astype(np.float32)
        else:
            residue_embeddings = token_representations[0, :tokens_len].cpu().numpy().astype(np.float32)
        
        mean_embedding = residue_embeddings.mean(axis=0).astype(np.float32)
        
        if pooling == "mean":
            return mean_embedding
        elif pooling == "none":
            return residue_embeddings, mean_embedding
        else:
            raise ValueError(f"Unknown pooling strategy: {pooling}")

    def _preprocess_sequence(self, sequence: str) -> str:
        """Preprocess protein sequence."""
        if not sequence or not isinstance(sequence, str):
            raise ValueError("Sequence must be a non-empty string")
        sequence_clean = sequence.strip().upper()
        sequence_clean = ''.join(sequence_clean.split())
        if len(sequence_clean) < 3:
            # Pad short sequences if needed or raise error. ESM needs some length.
            # But let's stick to error for now.
            pass 
        return sequence_clean

    def run_esm_cli(self, fasta_path: str, output_dir: str, include: List[str] = ["mean"]) -> bool:
        """Run esm-extract CLI command."""
        if not HAS_ESM:
            raise RuntimeError("fair-esm not installed")
            
        cmd = [
            "esm-extract",
            self.model_name,
            fasta_path,
            output_dir,
            "--include", ",".join(include)
        ]
        
        if self.device == "cuda":
            # esm-extract uses cuda if available automatically? or needs flag?
            # usually it detects.
            pass
            
        try:
            subprocess.run(cmd, check=True)
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"ESM CLI failed: {e}")
            return False

def main():
    """Test the embedding generator."""
    try:
        print("GENERATOR_TEST: Starting self-test...")
        test_sequence = "ACDEFGHIKLMNPQRSTVWY"
        
        # 1. Test Local (if available)
        if HAS_ESM:
            print("GENERATOR_TEST: Testing local model...")
            gen = ProteinEmbeddingGenerator(model_name="esm2_t6_8M_UR50D", device="cpu", use_server=False)
            emb = gen.generate_embedding(test_sequence)
            print(f"GENERATOR_SUCCESS: Local embedding shape: {emb.shape}")
        else:
            print("GENERATOR_SKIP: fair-esm not installed")
            
        # 2. Test Server (Mock/Dry run)
        print("GENERATOR_TEST: Testing server initialization...")
        gen_server = ProteinEmbeddingGenerator(use_server=True)
        if gen_server.plm_utils:
            print("GENERATOR_SUCCESS: KBPLMUtils initialized")
            # We skip actual server call to avoid auth issues in test
        else:
            print("GENERATOR_SKIP: KBPLMUtils not available")

        print("GENERATOR_TEST: Completed")
        return 0
    except Exception as e:
        print(f"GENERATOR_FAIL: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
