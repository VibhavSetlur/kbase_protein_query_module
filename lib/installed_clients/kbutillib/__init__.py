"""KBUtilLib - Modular utility framework for scientific and development projects."""

# Core utilities - these should always be available
from .base_utils import BaseUtils
from .shared_env_utils import SharedEnvUtils

# Optional modules - wrapped in try-except to handle missing dependencies
try:
    from .notebook_utils import NotebookUtils
except ImportError:
    NotebookUtils = None

try:
    from .kb_ws_utils import KBWSUtils
except ImportError:
    KBWSUtils = None

try:
    from .kb_genome_utils import KBGenomeUtils
except ImportError:
    KBGenomeUtils = None

try:
    from .ms_biochem_utils import MSBiochemUtils
except ImportError:
    MSBiochemUtils = None

try:
    from .kb_model_utils import KBModelUtils
except ImportError:
    KBModelUtils = None

try:
    from .kb_sdk_utils import KBSDKUtils
except ImportError:
    KBSDKUtils = None

try:
    from .argo_utils import ArgoUtils
except ImportError:
    ArgoUtils = None

try:
    from .kb_annotation_utils import KBAnnotationUtils
except ImportError:
    KBAnnotationUtils = None

# Import example composite classes
# Temporarily disabled for testing core functionality
# try:
#     from . import examples
# except ImportError:
#     examples = None
examples = None

__all__ = [
    "ArgoUtilsMSBiochemUtils",
    "BaseUtils",
    "KBAnnotationUtils",
    "KBGenomeUtils",
    "KBModelUtils",
    "KBSDKUtils",
    "KBWSUtils",
    "NotebookUtils",
    "SharedEnvUtils",
    "examples",
]

__version__ = "0.1.0"
