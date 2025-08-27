"""Notebook utility functions for Jupyter environments and interactive development."""

import json
import os
from genericpath import exists
from typing import Any, Optional, Dict
import pandas as pd

from .base_utils import BaseUtils


class NotebookUtils(BaseUtils):
    """Utility functions for working with Jupyter notebooks and interactive environments.

    Provides tools for notebook manipulation, output formatting, interactive displays,
    and integration with notebook-specific features like widgets and displays.
    """

    def __init__(self, notebook_folder: str, **kwargs: Any) -> None:
        """Initialize notebook utilities.

        Args:
            **kwargs: Additional keyword arguments passed to BaseUtil
        """
        super().__init__(**kwargs)
        self.notebook_folder = notebook_folder
        self.data_dir = os.path.join(notebook_folder, "data")
        os.makedirs(self.data_dir, exist_ok=True)
        self.datacache_dir = os.path.join(notebook_folder, "datacache")
        os.makedirs(self.datacache_dir, exist_ok=True)
        self.output_dir = os.path.join(notebook_folder, "nboutput")
        os.makedirs(self.output_dir, exist_ok=True)

        self.in_notebook = self._detect_notebook_environment()
        if self.in_notebook:
            self.log_info("Notebook environment detected")
        else:
            self.log_debug("Not running in notebook environment")

    def _detect_notebook_environment(self) -> bool:
        """Detect if code is running in a Jupyter notebook environment.

        Returns:
            True if running in a notebook, False otherwise
        """
        try:
            # Check for IPython/Jupyter kernel
            from IPython import get_ipython

            ipython = get_ipython()

            if ipython is None:
                return False

            # Check if we're in a notebook (not just IPython shell)
            return hasattr(ipython, "kernel")
        except ImportError:
            return False

    def display_dataframe(
        self,
        df: Any,
        max_rows: Optional[int] = None,
        max_cols: Optional[int] = None,
        use_interactive: bool = True,
        page_size: int = 25,
        show_search: bool = True,
        show_info: bool = True,
        scrollable: bool = True,
        height: Optional[str] = None,
    ) -> None:
        """Display a DataFrame with enhanced interactive features including pagination and search.

        Args:
            df: DataFrame to display (pandas, polars, etc.)
            max_rows: Maximum number of rows to display (None for all)
            max_cols: Maximum number of columns to display (None for all)
            use_interactive: Use interactive table with pagination and search
            page_size: Number of rows per page when using interactive display
            show_search: Show search functionality
            show_info: Show table information (row count, etc.)
            scrollable: Make table scrollable for large content
            height: Fixed height for the table (e.g., "400px")
        """
        if not self.in_notebook:
            print(df)
            return

        # Convert to pandas if not already
        if not isinstance(df, pd.DataFrame):
            try:
                # Try to convert other dataframe types to pandas
                if hasattr(df, "to_pandas"):
                    df = df.to_pandas()
                elif hasattr(df, "toPandas"):
                    df = df.toPandas()
                else:
                    df = pd.DataFrame(df)
            except Exception:
                # Fallback to simple display
                try:
                    from IPython.display import display

                    display(df)
                except ImportError:
                    print(df)
                return

        # Apply row/column limits if specified
        if max_rows and len(df) > max_rows:
            df = df.head(max_rows)
            self.log_info(f"Displaying first {max_rows} rows of {len(df)} total")

        if max_cols and len(df.columns) > max_cols:
            df = df.iloc[:, :max_cols]
            self.log_info(f"Displaying first {max_cols} columns of {len(df.columns)} total")

        if use_interactive and len(df) > page_size:
            self._display_interactive_dataframe(
                df, page_size, show_search, show_info, scrollable, height
            )
        else:
            self._display_simple_dataframe(df)

    def _display_interactive_dataframe(
        self,
        df: Any,
        page_size: int,
        show_search: bool,
        show_info: bool,
        scrollable: bool,
        height: Optional[str],
    ) -> None:
        """Display DataFrame with interactive features using itables."""
        try:
            # Use itables for interactive display
            import itables
            from IPython.display import display

            # Configure itables options
            itables.options.lengthMenu = [10, 25, 50, 100, -1]
            itables.options.pageLength = page_size
            itables.options.searching = show_search
            itables.options.info = show_info
            itables.options.scrollCollapse = True
            itables.options.scrollX = scrollable
            
            # Configure table width and text wrapping
            itables.options.style = "width:100%; table-layout: fixed;"
            itables.options.columnDefs = [
                {
                    "targets": "_all",
                    "className": "text-wrap"
                }
            ]
            
            if height:
                itables.options.scrollY = height

            # Display with itables (remove conflicting nowrap class)
            display(itables.show(df, classes="display"))
            
            # Add custom CSS for better text wrapping
            from IPython.display import HTML
            display(HTML("""
            <style>
            .dataTable td {
                word-wrap: break-word;
                word-break: break-word;
                white-space: normal !important;
                max-width: 200px;
            }
            .dataTable th {
                word-wrap: break-word;
                word-break: break-word;
                white-space: normal !important;
            }
            .text-wrap {
                white-space: normal !important;
            }
            </style>
            """))
            return

        except ImportError:
            # Clear error message about missing itables
            self.log_error(
                "Interactive table features require 'itables'. "
                "Install with: pip install itables"
            )
            print("\n" + "="*60)
            print("INTERACTIVE FEATURES UNAVAILABLE")
            print("="*60)
            print("To enable interactive tables with pagination and search:")
            print("  pip install itables")
            print("\nAlternatively, use: pip install \"KBUtilLib[notebook]\"")
            print("="*60 + "\n")
            
            # Fall back to simple display
            self._display_simple_dataframe(df)

    def _display_simple_dataframe(self, df: Any) -> None:
        """Display DataFrame using basic IPython display."""
        try:
            from IPython.display import display
            display(df)
        except ImportError:
            print(df)

    def display_json(self, data: Dict[str, Any], indent: int = 2) -> None:
        """Display JSON data with proper formatting in notebook.

        Args:
            data: Dictionary or JSON-serializable data
            indent: Indentation level for pretty printing
        """
        if not self.in_notebook:
            print(json.dumps(data, indent=indent))
            return

        try:
            from IPython.display import JSON, display

            display(JSON(data))
        except ImportError:
            print(json.dumps(data, indent=indent))

    def display_markdown(self, text: str) -> None:
        """Display markdown text in notebook.

        Args:
            text: Markdown formatted text
        """
        if not self.in_notebook:
            print(text)
            return

        try:
            from IPython.display import Markdown, display

            display(Markdown(text))
        except ImportError:
            print(text)

    def display_html(self, html: str) -> None:
        """Display HTML content in notebook.

        Args:
            html: HTML content string
        """
        if not self.in_notebook:
            print(html)
            return

        try:
            from IPython.display import HTML, display

            display(HTML(html))
        except ImportError:
            print(html)

    def create_progress_bar(self, total: int, description: str = "") -> Any:
        """Create a progress bar for notebook environments.

        Args:
            total: Total number of iterations
            description: Description text for the progress bar

        Returns:
            Progress bar object (tqdm if available, None otherwise)
        """
        if not self.in_notebook:
            return None

        try:
            from tqdm.notebook import tqdm

            return tqdm(total=total, desc=description)
        except ImportError:
            self.log_warning("tqdm not available for progress bars")
            return None

    def clear_output(self, wait: bool = False) -> None:
        """Clear notebook cell output.

        Args:
            wait: Whether to wait for new output before clearing
        """
        if not self.in_notebook:
            return

        try:
            from IPython.display import clear_output

            clear_output(wait=wait)
        except ImportError:
            pass

    def create_interactive_widget(self, widget_type: str, **kwargs: Any) -> Any:
        """Create an interactive widget for notebook use.

        Args:
            widget_type: Type of widget ('text', 'slider', 'dropdown', etc.)
            **kwargs: Widget-specific parameters

        Returns:
            Widget object if ipywidgets available, None otherwise
        """
        if not self.in_notebook:
            self.log_warning("Widgets only supported in notebook environments")
            return None

        try:
            import ipywidgets as widgets

            widget_map = {
                "text": widgets.Text,
                "textarea": widgets.Textarea,
                "slider": widgets.IntSlider,
                "float_slider": widgets.FloatSlider,
                "dropdown": widgets.Dropdown,
                "checkbox": widgets.Checkbox,
                "button": widgets.Button,
                "progress": widgets.IntProgress,
            }

            if widget_type in widget_map:
                return widget_map[widget_type](**kwargs)
            else:
                self.log_error(f"Unknown widget type: {widget_type}")
                return None

        except ImportError:
            self.log_warning("ipywidgets not available for interactive widgets")
            return None

    def save(self, name: str, data: Any) -> None:
        """Save data to a JSON file in the notebook data directory."""
        filename = self.datacache_dir + "/" + name + ".json"
        dir = os.path.dirname(filename)
        os.makedirs(dir, exist_ok=True)
        with open(filename, "w") as f:
            json.dump(data, f, indent=4, skipkeys=True)

    def load(
        self, name: str, default: Any = None, kb_type: Optional[str] = None
    ) -> Any:
        """Load data from a JSON file in the notebook data directory."""
        filename = self.datacache_dir + "/" + name + ".json"
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
        if kb_type is not None:
            data = self.kb_object_factory._build_object(kb_type, data, None, None)
        return data

    def list(self):
        """List all JSON files in the notebook data directory."""
        if not exists(self.datacache_dir):
            return []
        files = os.listdir(self.datacache_dir)
        return [x.split(".")[0] for x in files if x.endswith(".json")]

    def exists(self, name: str) -> bool:
        """Check if a JSON file exists in the notebook data directory."""
        return exists(self.datacache_dir + "/" + name + ".json")
