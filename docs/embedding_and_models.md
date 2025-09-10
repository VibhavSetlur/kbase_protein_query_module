# Embedding and Models (Local-only)

Model
- Uses local ESM model at `data/esm2_t6_8M_UR50D_local`
- No model selection in UI/spec; keep a single path

Code
- `processing/embeddings/generator.py` loads tokenizer/model from local path only
- TEST_MODE replaces model with mocks for unit tests

Recommendations
- Keep tensors on CPU unless you verify GPU availability
- Limit sequence length (e.g., 1024) and batch sizes conservatively
- Do not add HF fallbacks; Narrative deploy pulls from Docker image
