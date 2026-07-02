from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
# repo root is two levels up from training/tests/
REPO_ROOT = HERE.parents[1]
REAL_PARAM_DIR = REPO_ROOT / "param"


@pytest.fixture
def mini_param_path() -> Path:
    return DATA / "mini.param"


def real_param_files() -> list[Path]:
    if not REAL_PARAM_DIR.is_dir():
        return []
    return sorted(REAL_PARAM_DIR.glob("*.param"))
