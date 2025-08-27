"""Shared environment management for configuration and secrets."""

import os
from configparser import ConfigParser
from pathlib import Path
from typing import Any, Dict, Optional, Union

from .base_utils import BaseUtils


class SharedEnvUtils(BaseUtils):
    """Manages shared environment configuration, secrets, and runtime settings.

    Provides centralized access to configuration files, environment variables,
    authentication tokens, and other shared resources across utility modules.
    """

    def __init__(
        self,
        config_file: Optional[Union[str, Path]] = None,
        token_file: Optional[Union[str, Path]] = Path.home() / ".tokens",
        kbase_token_file: Optional[Union[str, Path]] = Path.home() / ".kbase" / "token",
        token: Optional[Union[str, Dict[str, str]]] = None,
        **kwargs: Any,
    ) -> None:
        """Initialize the shared environment which manages configurations and authentication tokens."""
        super().__init__(**kwargs)
        # Reading config file is specified
        self._config_hash = {}
        self._config_file = None
        if config_file and Path(config_file).exists():
            self._config_file = config_file
            self._config_hash = self.read_config()

        # Reading token file
        self._token_hash = {}
        self._env_vars = {}
        self._token_file = token_file
        self._kbase_token_file = kbase_token_file
        if (token_file and Path(token_file).exists()) or (
            kbase_token_file and Path(kbase_token_file).exists()
        ):
            self._token_hash = self.read_token_file()
        # Loading environment variables
        self.load_environment_variables()

        # Setting KBase token from input argument if provided
        if token is not None:
            if isinstance(token, str):
                # If token is a string, set it directly
                self.set_token(token, namespace="kbase")
            elif isinstance(token, dict):
                # If token is a dictionary, set each key-value pair
                for key, value in token.items():
                    self.set_token(value, key)

    # Reading config file in standard config format
    def read_config(
        self, config_file: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """Read configuration from a file."""
        if config_file is None:
            config_file = self._config_file
        config_path = Path(config_file)
        confighash = {}
        try:
            if config_path.exists():
                config = ConfigParser()
                config.read(config_path)
                for section in config.sections():
                    confighash[section] = {}
                    for nameval in config.items(section):
                        confighash[section][nameval[0]] = nameval[1]
                return confighash
            else:
                self.log_warning(f"Config file not found: {config_path}")
        except Exception as e:
            self.log_error(f"Error parsing config file: {e}")
            raise
        return confighash

    def get_config(self, section, key, default=None):
        """Get a configuration value."""
        if section not in self._config_hash:
            self.log_warning(f"Section '{section}' not found in config")
            return default
        if key not in self._config_hash[section]:
            self.log_warning(f"Key '{key}' not found in section '{section}'")
            return default
        return self._config_hash.get(section).get(key)

    def read_token_file(
        self, token_file: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """Load authentication tokens from a file."""
        # Load a standard text file with token name on one line and token value on the next
        token_hash = {}
        if token_file is None:
            if self._token_file:
                token_file = self._token_file

        token_path = Path(token_file)
        try:
            if token_path.exists():
                self.log_info(f"Loaded {len(token_hash)} tokens from {token_file}")
                with open(token_path) as fh:
                    for line in fh:
                        if line.strip():  # Ignore empty lines
                            key, value = line.strip().split("=", 1)
                            token_hash[key] = value
            if self._kbase_token_file is not None and self._kbase_token_file.exists():
                self.log_info(f"Loaded kbase tokens from {self._kbase_token_file}")
                with open(self._kbase_token_file) as fh:
                    token_hash["kbase"] = fh.read().strip()
            return token_hash
        except Exception as e:
            self.log_error(f"Error loading tokens from {token_file}: {e}")
            raise

    def save_token_file(self, token_file: Optional[Union[str, Path]] = None):
        """Save authentication tokens to a file."""
        if token_file is None:
            token_file = self._token_file
        if token_file is None:
            self.log_warning("No token file specified for saving")
            return

        token_path = Path(token_file)
        try:
            with open(token_path, "w") as fh:
                for key, value in self._token_hash.items():
                    fh.write(f"{key}={value}\n")
            # Also saving kbase tokens in the .kbase directory to ensure compatibility with CobraKBase
            if "kbase" in self._token_hash:
                if self._kbase_token_file is not None:
                    directory = Path(self._kbase_token_file).parent
                    if not directory.exists():
                        directory.mkdir(parents=True, exist_ok=True)
                with open(self._kbase_token_file, "w") as fh:
                    fh.write(self._token_hash["kbase"])
            self.log_info(f"Saved {len(self._token_hash)} tokens to {token_file}")
        except Exception as e:
            self.log_error(f"Error saving tokens to {token_file}: {e}")
            raise

    def set_token(self, token, namespace="kbase", save_file=True):
        """Set a token and save it to the token file."""
        self._token_hash[namespace] = token
        if save_file:
            self.save_token_file()

    # Retrieves specified token from cached data
    def get_token(self, namespace="kbase") -> Any:
        """Retrieve a stored secret."""
        return self._token_hash.get(namespace, None)

    def load_environment_variables(self) -> None:
        """Load relevant environment variables into the shared environment."""
        # Load common environment variables that might be used across utilities
        env_prefixes = ["KB_", "KBASE_", "MS_", "NOTEBOOK_"]

        for key, value in os.environ.items():
            if any(key.startswith(prefix) for prefix in env_prefixes):
                self._env_vars[key] = value

        self.log_debug(f"Loaded {len(self._env_vars)} environment variables")

    def get_env_var(self, key: str, default: Optional[str] = None) -> Optional[str]:
        """Get an environment variable value.

        Args:
            key: Environment variable name
            default: Default value if not found

        Returns:
            Environment variable value or default
        """
        return self._env_vars.get(key, default)

    # Printing environment state for debugging or inspection
    def export_environment(self) -> Dict[str, Any]:
        """Export the current environment state for debugging or inspection."""
        return {
            "config": self._config_hash if self._config_hash else {},
            "config_file": str(self._config_file) if self._config_file else None,
            "kbase_token_file": str(self._kbase_token_file)
            if self._kbase_token_file
            else None,
            "token_file": str(self._token_file) if self._token_file else None,
            "env_vars": self._env_vars if self._env_vars else {},
            "token_keys": list(self._token_hash.keys())
            if self._token_hash
            else [],  # Don't expose secret values
        }
