"""Base utility class providing core shared logic for all utility modules."""

import json
import logging
import subprocess
import sys
import os
import time
from genericpath import exists
from pathlib import Path
from typing import Any, Dict, List

import requests

requests.packages.urllib3.disable_warnings()

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

class BaseUtils:
    """Base class for all utility modules in the KBUtilLib framework.

    Provides core shared functionality including logging, error handling,
    and common utility methods that are inherited by all specialized
    utility modules.
    """

    def __init__(self, name="Unknown", log_level: str = "INFO", **kwargs: Any) -> None:
        """Initialize the base utility class."""
        # Add dependencies directory to Python path first thing
        self._setup_dependencies_path()
        self.logger = self._setup_logger(log_level)

        # Allow subclasses to pass additional initialization parameters
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.version = "0.0.0"
        self.name = name
        self.util_directory = script_dir+"/../../"
        self.data_directory = self.util_directory+"/data/"

        # Initialize attributes for tracking provenance on primary method calls
        # self.obj_created = []
        # self.input_objects = []
        # self.method = None
        # self.params = {}
        # self.initialized = False
        # self.timestamp = None
        self.reset_attributes()

    def reset_attributes(self):
        # Initializing stores tracking objects created and input objects
        self.obj_created = []
        self.input_objects = []
        # Initializing attributes tracking method data to support provencance and context
        self.method = None
        self.params = {}
        self.initialized = False
        self.timestamp = None

    def initialize_call(
        self,
        method: str,
        params: Dict[str, Any],
        print_params: bool = False,
        no_print: List[str] = None,
        no_prov_params: List[str] = None,
    ) -> None:
        """This function reiniitializes a module method call for provenance tracking."""
        if no_print is None:
            no_print = []
        if no_prov_params is None:
            no_prov_params = []

        if not self.initialized:
            # Computing timestamp
            ts = time.gmtime()
            self.timestamp = time.strftime("%Y-%m-%d %H:%M:%S", ts)

            self.obj_created = []
            self.input_objects = []
            self.method = method

            # Filter parameters for provenance
            filtered_params = {}
            for key in params:
                if key not in no_prov_params:
                    filtered_params[key] = params[key]

            self.params = filtered_params.copy()
            self.initialized = True

            # Print parameters if requested
            if print_params:
                log_params = filtered_params.copy()
                for item in no_print:
                    if item in log_params:
                        del log_params[item]
                self.log_info(f"{method}: {json.dumps(log_params, indent=2)}")

    def _setup_logger(self, log_level: str) -> logging.Logger:
        """Set up logging for the utility module."""
        logger = logging.getLogger(
            f"{self.__class__.__module__}.{self.__class__.__name__}"
        )
        logger.setLevel(getattr(logging, log_level.upper()))
        logger.propagate = False  # <- stop bubbling to root which causes log messages to show up twice in jupyter notebooks

        # Only add handler if none exists to prevent duplicate logs
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)

        return logger

    def _setup_dependencies_path(self) -> None:
        """Add the dependencies directory to Python path if it exists.
        Also adds specific repository paths for symbolic links.
        """
        # Get the directory containing this module
        current_dir = Path(__file__).parent
        self.dependencies_dir = current_dir / "dependencies"

        if self.dependencies_dir.exists():
            dependencies_str = str(self.dependencies_dir)
            if dependencies_str not in sys.path:
                sys.path.insert(0, dependencies_str)

            # Also add specific paths for repositories that might be symbolic links
            # Add ModelSEEDpy if it exists (the repo root contains the package)
            modelseedpy_path = self.dependencies_dir / "modelseedpy"
            if modelseedpy_path.is_symlink():
                target_repo = str(
                    modelseedpy_path.resolve()
                )  # The ModelSEEDpy repo root
                if target_repo not in sys.path:
                    sys.path.insert(0, target_repo)

            # Add cobrakbase if it exists (need the parent containing the cobrakbase repo)
            cobrakbase_path = self.dependencies_dir / "cobrakbase"
            if cobrakbase_path.is_symlink():
                # cobrakbase symlink points to /Users/chenry/code/cobrakbase/cobrakbase
                # We need /Users/chenry/code/cobrakbase in the path
                target_pkg = cobrakbase_path.resolve()  # Points to the package dir
                target_repo = str(target_pkg.parent)  # The repo root
                if target_repo not in sys.path:
                    sys.path.insert(0, target_repo)

    def _ensure_git_dependency(
        self, module_name: str, git_url: str, branch: str = "master"
    ) -> bool:
        """Ensure a git-based dependency is available in the dependencies directory.

        Args:
            module_name: Name of the module to check for
            git_url: Git URL to clone from if module is missing
            branch: Git branch to checkout (default: master)

        Returns:
            bool: True if module is available, False if failed to obtain
        """
        module_dir = self.dependencies_dir / module_name

        # Check if module is already available and importable
        try:
            __import__(module_name)
            self.log_debug(f"Module {module_name} already available and importable")
            return True
        except ImportError as e:
            self.log_debug(f"Module {module_name} not importable: {e}")
            # Check if it's due to missing dependencies vs missing module
            if "No module named" in str(e) and module_name not in str(e):
                # Module exists but has missing dependencies
                self.log_warning(f"Module {module_name} has missing dependencies: {e}")
                # For now, return False but could be enhanced to install dependencies
                return False

        # Check if module exists in dependencies directory
        if module_dir.exists():
            # Module directory exists, try importing again
            try:
                __import__(module_name)
                return True
            except ImportError:
                # Directory exists but module not working, try updating
                self.log_warning(
                    f"Module {module_name} directory exists but import failed, attempting git pull"
                )
                try:
                    result = subprocess.run(
                        ["git", "pull"],
                        cwd=module_dir,
                        capture_output=True,
                        text=True,
                        timeout=30,
                    )
                    if result.returncode == 0:
                        try:
                            __import__(module_name)
                            return True
                        except ImportError:
                            pass
                except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
                    self.log_warning(f"Git pull failed for {module_name}: {e}")

        # Module not available, try to clone it
        self.log_info(f"Cloning {module_name} from {git_url}")
        try:
            # Ensure dependencies directory exists
            self.dependencies_dir.mkdir(exist_ok=True)

            # Clone the repository
            result = subprocess.run(
                ["git", "clone", "-b", branch, git_url, str(module_dir)],
                capture_output=True,
                text=True,
                timeout=120,  # Increased timeout for large repositories
            )

            if result.returncode == 0:
                self.log_info(f"Successfully cloned {module_name}")
                # Try importing the module
                try:
                    __import__(module_name)
                    return True
                except ImportError as e:
                    self.log_error(f"Cloned {module_name} but import still failed: {e}")
                    return False
            else:
                self.log_error(f"Failed to clone {module_name}: {result.stderr}")
                return False

        except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
            self.log_error(f"Error cloning {module_name}: {e}")
            return False

    def ensure_modelseed_database(self) -> bool:
        """Ensure ModelSEEDDatabase is available."""
        return self._ensure_git_dependency(
            "ModelSEEDDatabase",
            "https://github.com/ModelSEED/ModelSEEDDatabase.git",
            "master",
        )

    def ensure_modelseed_py(self) -> bool:
        """Ensure ModelSEEDpy is available."""
        return self._ensure_git_dependency(
            "modelseedpy", "https://github.com/ModelSEED/ModelSEEDpy.git", "main"
        )

    def ensure_cobra_kbase(self) -> bool:
        """Ensure CobraKBase is available."""
        return self._ensure_git_dependency(
            "cobrakbase", "https://github.com/Fxe/cobrakbase.git", "master"
        )

    def ensure_annotation_ontology_api(self) -> bool:
        """Ensure Annotation Ontology API is available."""
        return self._ensure_git_dependency(
            "cb_annotation_ontology_api",
            "https://github.com/kbaseapps/cb_annotation_ontology_api.git",
            "main",
        )

    def log_info(self, message: str) -> None:
        """Log an info message."""
        self.logger.info(message)

    def log_warning(self, message: str) -> None:
        """Log a warning message."""
        self.logger.warning(message)

    def log_error(self, message: str) -> None:
        """Log an error message."""
        self.logger.error(message)

    def log_debug(self, message: str) -> None:
        """Log a debug message."""
        self.logger.debug(message)

    def log_critical(self, message: str) -> None:
        """Log a critical message."""
        self.logger.critical(message)

    def print_attributes(self, obj=None, properties=True, functions=True):
        """Print attributes and functions of this object (or another object), useful with all the inheritance we're using"""
        if obj is None:
            obj = self
        attributes = dir(obj)
        if properties:
            print("Properties:")
            properties = [
                attr for attr in attributes if not callable(getattr(obj, attr))
            ]
            for property in properties:
                print(f"{property}")
        if functions:
            print("Functions:")
            functions = [attr for attr in attributes if callable(getattr(obj, attr))]
            for func in functions:
                print(f"{func}")

    def validate_args(
        self, params: Dict[str, Any], required: List[str], defaults: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Validate method arguments and apply defaults."""
        for item in required:
            if item not in params:
                raise ValueError(f"Required argument {item} is missing!")

        for key, value in defaults.items():
            if key not in params:
                params[key] = value

        return params

    def transfer_outputs(
        self, output: Dict[str, Any], api_output: Dict[str, Any], key_list: List[str]
    ) -> None:
        """Transfer specified keys from API output to the output dictionary."""
        for key in key_list:
            if key in api_output:
                output[key] = api_output[key]

    def save_util_data(self, name: str, data: Any) -> None:
        """Save data to a JSON file in the notebook data directory."""
        filename = self.data_directory + "/" + name + ".json"
        dir = os.path.dirname(filename)
        os.makedirs(dir, exist_ok=True)
        with open(filename, "w") as f:
            json.dump(data, f, indent=4, skipkeys=True)

    def load_util_data(
        self, name: str, default: Any = None
    ) -> Any:
        """Load data from a JSON file in the notebook data directory."""
        filename = self.data_directory + "/" + name + ".json"
        if not exists(filename):
            if default is None:
                self.log_error(
                    "Requested data " + name + " doesn't exist at " + filename
                )
                raise (
                    ValueError(
                        "Requested data " + name + " doesn't exist at " + filename
                    )
                )
            return default
        with open(filename) as f:
            data = json.load(f)
        return data