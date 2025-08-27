"""KBase SDK utilities for working with KBase SDK environments and services."""

from typing import Any

from .kb_ws_utils import KBWSUtils


class KBSDKUtils(KBWSUtils):
    """Utilities for working with KBase SDK environments and common SDK operations.

    Provides methods for workspace operations, client management, file handling,
    report generation, and other SDK-specific functionality commonly used in
    KBase applications and services.
    """

    def __init__(self, **kwargs: Any) -> None:
        """Initialize KBase SDK utilities.

        Args:
            module_name: Name of the SDK module
            working_dir: Working directory (defaults to scratch from config)
            module_dir: Module directory path
            callback_url: Callback URL for SDK services
            **kwargs: Additional keyword arguments passed to SharedEnvironment
        """
        super().__init__(**kwargs)

    def build_dataframe_report(self, table, column_list):
        # Convert columns to this format:
        columns = []
        for item in column_list:
            columns.append({"data": item})
        # for index, row in table.iterrows():
        #    pass
        json_str = table.to_json(orient="records")
        # columns=column_list
        html_data = (
            """
<html>
<header>
    <link href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css" rel="stylesheet">
</header>
<body>
<script src="https://code.jquery.com/jquery-3.6.0.slim.min.js" integrity="sha256-u7e5khyithlIdTpu22PHhENmPcRdFiHRjhAuHcs05RI=" crossorigin="anonymous"></script>
<script type="text/javascript" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
<script>
    $(document).ready(function() {
        $('#example').DataTable( {
            "ajax": {
                "url": "data.json"
            },
            "columns": """
            + json.dumps(columns, indent=4)
            + """
        } );
    } );
</script>
</body>
</html>
"""
        )
        os.makedirs(self.working_dir + "/html", exist_ok=True)
        with open(self.working_dir + "/html/index.html", "w") as f:
            f.write(html_data)
        with open(self.working_dir + "/html/data.json", "w") as f:
            f.write(json_str)
