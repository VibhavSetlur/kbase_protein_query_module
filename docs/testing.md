# Testing

Unit tests (local)
- Activate env: `conda activate kbase_test`
- Set: `TEST_MODE=true`, optional `KB_AUTH_TOKEN=<token>`
- Run: `python -m pytest` (ignore integration tests if needed)

Integration tests (sdk)
- Run: `kb-sdk test` (uses KBase runtime)
- If validator is down, try `kb-sdk test -s` and retry later

Notes
- Do not add HF model fallbacks; tests use mocks in TEST_MODE
- Keep PYTHONPATH pointing to module root when running locally
