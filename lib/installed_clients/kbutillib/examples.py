"""Example composite utility classes demonstrating the modular design."""

from .kb_callback_utils import KBCallbackUtils
from .kb_genome_utils import KBGenomeUtils
from .kb_model_utils import KBModelUtils
from .kb_sdk_utils import KBSDKUtils
from .ms_biochem_utils import MSBiochemUtils
from .notebook_utils import NotebookUtils
from .shared_env_utils import SharedEnvUtils


class KBaseWorkbench(
    KBCallbackUtils, KBGenomeUtils, KBModelUtils, KBSDKUtils, SharedEnvUtils
):
    """Composite utility class for comprehensive KBase workflows.

    Combines KBase API access, genome utilities, model utilities, SDK utilities,
    and shared environment management for complete KBase data analysis workflows.

    Initialize the KBase workbench with all utilities.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.log_info(
            "KBase Workbench initialized with API, genome, model, and SDK utilities"
        )


class KBaseSDKWorkflow(KBSDKUtils, SharedEnvUtils):
    """Composite utility class for KBase SDK development workflows.

    Combines SDK utilities and shared environment management
    for developing and testing KBase SDK applications.
    """

    def __init__(self, **kwargs):
        """Initialize the SDK workflow environment."""
        super().__init__(**kwargs)
        self.log_info("KBase SDK Workflow initialized")


class NotebookAnalysis(NotebookUtils, MSBiochemUtils, SharedEnvUtils):
    """Composite utility class for interactive notebook-based analysis.

    Combines notebook display utilities, biochemistry database analysis tools,
    and shared environment management for interactive data analysis.
    """

    def __init__(self, **kwargs):
        """Initialize the notebook analysis environment."""
        super().__init__(**kwargs)
        self.log_info("Notebook Analysis environment initialized")


class GenomeMetabolismWorkflow(KBGenomeUtils, KBModelUtils, KBSDKUtils, SharedEnvUtils):
    """Composite utility class for genome-to-model workflows.

    Combines genome analysis, metabolic modeling, SDK utilities,
    and shared environment for complete genome-to-metabolism workflows.
    """

    def __init__(self, **kwargs):
        """Initialize the genome-metabolism workflow utilities."""
        super().__init__(**kwargs)
        self.log_info("Genome-Metabolism Workflow initialized")

    def genome_to_model_pipeline(self, genome_ref: str, workspace_id: str):
        """Example pipeline method combining multiple utilities.

        Args:
            genome_ref: Reference to KBase genome object
            workspace_id: Target workspace for model creation
        """
        self.log_info(f"Starting genome-to-model pipeline for {genome_ref}")

        # Initialize SDK method call
        self.initialize_call(
            "genome_to_model_pipeline",
            {"genome_ref": genome_ref, "workspace_id": workspace_id},
        )

        # Get genome data using KBaseAPI
        genome_data = self.get_object(genome_ref, workspace_id)

        # Analyze genome using KBGenomeUtils
        genome_info = self.parse_genome_object(genome_data)
        features = self.extract_features_by_type(genome_data, "CDS")

        self.log_info(f"Analyzed genome with {len(features)} coding sequences")

        # This would continue with model reconstruction...
        # (Actual model building would require additional KBase services)

        return {
            "genome_info": genome_info,
            "features_analyzed": len(features),
            "status": "pipeline_completed",
        }


class FullUtilityStack(
    KBGenomeUtils,
    MSBiochemUtils,
    KBModelUtils,
    KBSDKUtils,
    NotebookUtils,
    SharedEnvUtils,
):
    """Composite utility class with all available utilities.

    The complete "Swiss Army knife" with access to all utility modules.
    Use this when you need maximum flexibility across all domains.
    """

    def __init__(self, **kwargs):
        """Initialize the full utility stack."""
        super().__init__(**kwargs)
        self.log_info("Full Utility Stack initialized with all available modules")


# Example of how users can create their own custom combinations
class MyCustomUtils(MSBiochemUtils, KBSDKUtils, SharedEnvUtils):
    """Example of a user-defined custom utility combination.

    Users can inherit from any combination of utility classes to create
    their own specialized utility collections.
    """

    def __init__(self, **kwargs):
        """Initialize custom utilities."""
        super().__init__(**kwargs)
        self.log_info("Custom utility combination initialized")

    def my_specialized_workflow(self):
        """Example of adding custom methods to the composite class."""
        self.log_info("Running specialized workflow...")

        # Initialize method call for SDK tracking
        self.initialize_call("my_specialized_workflow", {})

        # Custom workflow combining KBase, MS, and SDK utilities
        env = self.kb_environment()
        self.log_info(f"Running in KBase environment: {env}")

        # Example: save some analysis results
        results = {"environment": env, "timestamp": self.timestamp}
        self.save_json("workflow_results", results)

        return results
