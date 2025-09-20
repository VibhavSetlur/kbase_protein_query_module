# KBase Protein Query Module

A modular protein analysis system for KBase with plug-and-play architecture.

## Architecture

Modular backend with thin KBase facade:

```
KBase Facade (Impl) → WorkflowOrchestrator → Modular Backend
                                           ├── Input Processing
                                           ├── Analysis System  
                                           ├── Output Generation
                                           └── Utilities
```

## Documentation

- [Developer Guide](docs/developer_guide.md) - Architecture and development
- [Analysis Extensions](docs/analysis_extensions.md) - Adding new analyses
- [Testing](docs/testing.md) - Running and writing tests

## Adding New Analyses

1. Create analysis module in `src/analysis/your_analysis/`
2. Register in `src/analysis/config.py`
3. Add output handler in `src/outputs/analysis/your_analysis/`
4. Test with `pytest test/unit_tests/`

## License

MIT License - see [LICENSE](LICENSE) file.