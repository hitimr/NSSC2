import os
from pathlib import Path




path = Path(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE_DIR = str(path.parent.absolute()) + "/"
