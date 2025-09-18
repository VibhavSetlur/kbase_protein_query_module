"""
Safe deletion utilities.

Prefer moving files to the OS trash (Send2Trash). If unavailable, fallback to
moving into the project-level trash directory.
"""

import os
import shutil
import logging
from typing import Optional

logger = logging.getLogger(__name__)

try:
    from send2trash import send2trash as _send2trash
except Exception:  # pragma: no cover
    _send2trash = None


def move_to_trash(path: str, project_root: Optional[str] = None) -> bool:
    """
    Move file/directory to OS trash if possible; otherwise into project trash/.

    Args:
        path: Absolute or relative path to file or directory
        project_root: Root of project for fallback trash directory

    Returns:
        True if moved successfully, False otherwise
    """
    if not os.path.exists(path):
        logger.warning(f"Path not found for trash: {path}")
        return False

    # Try OS trash first
    if _send2trash is not None:
        try:
            _send2trash(path)
            logger.info(f"Sent to OS trash: {path}")
            return True
        except Exception as e:  # pragma: no cover
            logger.warning(f"Send2Trash failed for {path}: {e}")

    # Fallback to project trash directory
    try:
        root = project_root or os.getcwd()
        trash_dir = os.path.join(root, 'trash')
        os.makedirs(trash_dir, exist_ok=True)

        name = os.path.basename(path.rstrip(os.sep))
        dest = os.path.join(trash_dir, name)

        # Avoid collisions
        suffix = 1
        base_dest = dest
        while os.path.exists(dest):
            dest = f"{base_dest}_{suffix}"
            suffix += 1

        shutil.move(path, dest)
        logger.info(f"Moved to project trash: {path} -> {dest}")
        return True
    except Exception as e:
        logger.error(f"Failed to move to trash: {path}: {e}")
        return False


