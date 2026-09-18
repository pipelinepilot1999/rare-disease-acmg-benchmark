"""Path resolution for the pipeline scripts.

The repository root is found from this file, so a clone runs from wherever it is
checked out rather than from the author's home directory. Everything that is not
in the repository - downloaded resources, Exomiser working directories, the
per-case VCFs - lives outside it and is overridable by environment variable,
because those are tens of gigabytes and are deliberately not tracked.

    RDACMG_HOME       default ~/rdacmg     resources, work, cases, Exomiser output
    RDACMG_REPO       default <this repo>  code, data, results
"""
from __future__ import annotations
import os
from pathlib import Path

REPO = Path(os.environ.get("RDACMG_REPO") or Path(__file__).resolve().parent.parent)
HOME = Path(os.environ.get("RDACMG_HOME") or Path("~/rdacmg").expanduser())

RES = HOME / "resources"
WORK = HOME / "work"
CACHE = WORK / "cache"
CASES = HOME / "cases"
ARMA = HOME / "arma"
JOBS = HOME / "jobs"
BATCHES = HOME / "batches"
LOGS = HOME / "logs"

# Scripts were written against string paths; keep them strings at the boundary.
def s(p) -> str:
    return str(p)


def require(path, what: str) -> str:
    """Fail with an actionable message rather than a bare FileNotFoundError."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(
            f"missing {what}: {p}\n"
            f"  This lives outside the repository. Set RDACMG_HOME, or re-run the\n"
            f"  earlier pipeline step that produces it (see README 'Reproducing')."
        )
    return str(p)
