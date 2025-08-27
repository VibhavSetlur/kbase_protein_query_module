# KBase Protein Query Module - UI Documentation

This directory contains the KBase Narrative interface specifications and UI configuration for the Protein Query Module apps.

## Directory Structure

```
ui/
├── narrative/                           # KBase Narrative interface
│   └── methods/                        # App method specifications
│       ├── AssignProteinFamilyWithEmbedding/  # Family assignment app
│       │   ├── spec.json              # Method specification
│       │   ├── display.yaml           # UI display configuration
│       │   └── img/                   # App images (optional)
│       ├── CheckProteinExistence/     # Protein existence checking app
│       │   ├── spec.json
│       │   ├── display.yaml
│       │   └── img/
│       ├── GenerateProteinEmbeddings/ # Protein embedding generation app
│       │   ├── spec.json
│       │   ├── display.yaml
│       │   └── img/
│       ├── FindTopMatchesFromEmbedding/ # Similarity search app
│       │   ├── spec.json
│       │   ├── display.yaml
│       │   └── img/
│       └── SummarizeAndVisualizeResults/ # Results visualization app
│           ├── spec.json
│           ├── display.yaml
│           └── img/
└── README.md                           # This documentation
```

## KBase Apps Overview

The module provides five main KBase apps that form a complete protein analysis workflow:

### 1. Check Protein Existence (`CheckProteinExistence`)
**Purpose**: Quickly check if a UniProt protein exists in the KBase Protein Network

**Input Parameters**:
- `protein_id` (text): UniProt protein ID (e.g., P12345)

**Output**:
- Existence status (boolean)
- Family ID (if found)
- Metadata (if available)
- Report with detailed information

**Usage**: Use this as the first step in your protein analysis workflow to verify protein availability.

### 2. Generate Protein Embeddings (`GenerateProteinEmbeddings`)
**Purpose**: Generate high-dimensional protein embeddings using ESM-2 models

**Input Parameters**:
- `sequence` (textarea): Protein amino acid sequence

**Output**:
- Protein embedding (list of floats)
- Embedding summary
- Sequence length and normalization info
- Report with embedding details

**Usage**: Generate embeddings for downstream analysis like similarity search and family assignment.

### 3. Assign Protein Family (`AssignProteinFamilyWithEmbedding`)
**Purpose**: Assign protein embeddings to families using centroid similarity

**Input Parameters**:
- `embedding_ref` (workspace object): Reference to a previously generated embedding

**Output**:
- Family ID
- Confidence score
- Representative protein ID
- Assignment details

**Usage**: Fast family classification for downstream analysis and network construction.

### 4. Find Top Matches (`FindTopMatchesFromEmbedding`)
**Purpose**: Find similar proteins using protein embeddings

**Input Parameters**:
- `embedding_ref` (workspace object): Reference to a protein embedding
- `family_id` (text): Family to search within
- `top_n` (text, optional): Number of matches to return (default: 10)

**Output**:
- Ranked list of similar proteins
- Similarity scores
- Protein metadata
- Search statistics

**Usage**: Discover similar proteins for network analysis and functional inference.

### 5. Summarize and Visualize Results (`SummarizeAndVisualizeResults`)
**Purpose**: Create interactive visualizations and summary reports

**Input Parameters**:
- `top_matches_result_ref` (workspace object): Reference to search results
- `output_dir` (text, optional): Output directory name

**Output**:
- Interactive network visualization
- Summary statistics
- Downloadable reports
- HTML report with plots

**Usage**: Final step to visualize and analyze protein similarity networks.

## UI Specification Files

### spec.json Structure

Each app has a `spec.json` file that defines:
- **Version**: App version number
- **Authors**: App developers
- **Contact**: Support contact information
- **Categories**: App categorization
- **Parameters**: Input parameter definitions
- **Behavior**: Service mapping and input/output configuration

**Example**:
```json
{
    "ver": "0.0.1",
    "authors": ["yourname"],
    "contact": "https://kbase.us/contact-us/",
    "categories": ["active", "protein", "search"],
    "parameters": [
        {
            "id": "protein_id",
            "field_type": "text",
            "text_options": {
                "validate_as": "regex",
                "regex_constraint": "^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$"
            }
        }
    ],
    "behavior": {
        "service-mapping": {
            "name": "kbase_protein_query_module",
            "method": "check_protein_existence",
            "input_mapping": [...],
            "output_mapping": [...]
        }
    }
}
```

### display.yaml Structure

Each app has a `display.yaml` file that defines:
- **App Name**: Human-readable app name
- **Tooltip**: Brief description
- **Description**: Detailed app description
- **Parameters**: UI parameter labels and hints
- **Suggestions**: Related apps and workflow suggestions
- **Publications**: Relevant publications and links

**Example**:
```yaml
name: Check Protein Existence
tooltip: |
    Quickly check if a UniProt protein exists in the KBase Protein Network.

parameters:
    protein_id:
        ui-name: Protein ID (UniProt)
        short-hint: Enter a valid UniProt accession (e.g., P12345).

description: |
    <p>Enter a UniProt ID to check if it exists in the KBase Protein Network.</p>
    <ul>
    <li><b>Input:</b> UniProt protein ID</li>
    <li><b>Output:</b> Existence status and metadata</li>
    </ul>
```

## Parameter Types

### Text Parameters
```json
{
    "field_type": "text",
    "text_options": {
        "validate_as": "regex",
        "regex_constraint": "^pattern$",
        "regex_error_message": "Invalid format"
    }
}
```

### Textarea Parameters
```json
{
    "field_type": "textarea",
    "text_options": {
        "example": ">ExampleProtein\nMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
    }
}
```

### Workspace Object References
```json
{
    "field_type": "workspaceObjectRef",
    "text_options": {
        "valid_ws_types": ["kbase_protein_query_module.ProteinEmbeddingResult"]
    }
}
```

### Integer Parameters
```json
{
    "field_type": "text",
    "text_options": {
        "validate_as": "int",
        "min_int": 1,
        "max_int": 1000
    }
}
```

## Workflow Integration

### Suggested Workflows

1. **Basic Protein Analysis**:
   ```
   Check Protein Existence → Generate Protein Embeddings → Assign Protein Family
   ```

2. **Similarity Search Workflow**:
   ```
   Generate Protein Embeddings → Find Top Matches → Summarize and Visualize Results
   ```

3. **Complete Analysis**:
   ```
   Check Protein Existence → Generate Protein Embeddings → Assign Protein Family → Find Top Matches → Summarize and Visualize Results
   ```

### App Suggestions

Each app includes suggestions for:
- **Related Apps**: Apps that work well together
- **Next Steps**: Logical next apps in the workflow
- **Categories**: App categorization for discovery

## UI Best Practices

### Parameter Design
- **Clear Labels**: Use descriptive parameter names
- **Helpful Hints**: Provide short and long hints
- **Validation**: Use regex validation for format checking
- **Examples**: Include example values where helpful
- **Defaults**: Provide sensible default values

### Error Handling
- **Input Validation**: Validate parameters before processing
- **Clear Messages**: Provide helpful error messages
- **Graceful Degradation**: Handle missing or invalid data
- **User Feedback**: Show progress and status updates

### Accessibility
- **Keyboard Navigation**: Ensure keyboard accessibility
- **Screen Readers**: Use semantic HTML and ARIA labels
- **Color Contrast**: Ensure sufficient color contrast
- **Responsive Design**: Work on different screen sizes

## Customization

### Adding Images
1. Create an `img/` directory in your app folder
2. Add images (PNG, JPG, GIF supported)
3. Reference images in `display.yaml`:
   ```yaml
   screenshots:
       - img/screenshot1.png
       - img/screenshot2.png
   ```

### Styling
- Use HTML in descriptions for formatting
- Include links to external resources
- Add icons for visual appeal
- Use consistent styling across apps

### Localization
- Support multiple languages if needed
- Use clear, non-technical language
- Provide translations for key terms
- Consider cultural differences

## Testing UI

### Local Testing
```bash
# Test UI specifications
python -c "
import json
import yaml
with open('ui/narrative/methods/CheckProteinExistence/spec.json') as f:
    spec = json.load(f)
print('Spec valid:', spec['ver'])
"
```

### Validation
- **JSON Schema**: Validate spec.json against KBase schema
- **YAML Syntax**: Check display.yaml syntax
- **Parameter Mapping**: Verify input/output mappings
- **Workspace Types**: Confirm workspace object types

## Deployment

### Registration Process
1. **Push to GitHub**: Ensure all UI changes are committed
2. **Register on AppDev**: Use KBase AppDev environment
3. **Test in Development**: Switch to 'D' mode in Apps Panel
4. **Move to Beta**: Migrate to beta for user testing
5. **Request Release**: Submit for production release

### Version Management
- **Update Version**: Increment version in spec.json
- **Update Authors**: Add new contributors
- **Update Contact**: Ensure contact information is current
- **Update Documentation**: Keep descriptions current

## Troubleshooting

### Common Issues

1. **Parameter Validation Errors**
   ```bash
   # Check parameter definitions
   cat ui/narrative/methods/AppName/spec.json | jq '.parameters'
   ```

2. **Output Mapping Errors**
   ```bash
   # Verify output mappings match KIDL spec
   grep -r "output_mapping" ui/narrative/methods/
   ```

3. **Workspace Type Errors**
   ```bash
   # Check workspace object types
   grep -r "valid_ws_types" ui/narrative/methods/
   ```

### Debug Mode
```bash
# Enable debug logging
export KBASE_DEBUG=1
kb-sdk test
```

## Resources

- **[KBase UI Specification](https://kbase.github.io/kb_sdk_docs/references/narrative_app_ui_specification.html)**
- **[Adding UI Elements](https://kbase.github.io/kb_sdk_docs/howtos/add_ui_elements.html)**
- **[KBase Narrative Interface](https://docs.kbase.us/kbase.us)**
- **[KBase Developer Guidelines](https://kbase.github.io/kb_sdk_docs/references/developer_guidelines.html)**

## Contributing

When adding new apps or modifying existing ones:
1. **Follow Naming Conventions**: Use descriptive app names
2. **Update Documentation**: Keep README files current
3. **Test Thoroughly**: Validate all UI specifications
4. **Get Feedback**: Test with potential users
5. **Version Control**: Use descriptive commit messages

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.
