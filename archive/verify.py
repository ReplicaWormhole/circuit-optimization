"""Check the byte-preserved legacy namespace and its compatibility links."""

import hashlib
import json
import os
from pathlib import Path


def verify(base: Path) -> dict:
    manifest = json.loads((base / "manifest.json").read_text())
    if manifest.get("schema") != 1:
        raise ValueError("unsupported archive manifest schema")
    legacy = base / "legacy_root"
    problems = []
    checked = 0
    for group in ("moved_root_files", "frozen_helper_copies"):
        for name, recorded in manifest[group].items():
            if Path(name).name != name:
                problems.append(f"invalid archive name: {name}")
                continue
            path = legacy / name
            if not path.is_file() or path.is_symlink():
                problems.append(f"missing regular file: {name}")
                continue
            raw = path.read_bytes()
            checked += 1
            if len(raw) != recorded["size"] or hashlib.sha256(raw).hexdigest() != recorded["sha256"]:
                problems.append(f"changed content: {name}")
    for name, target in manifest["compatibility_links"].items():
        path = legacy / name
        if not path.is_symlink() or os.readlink(path) != target or not path.exists():
            problems.append(f"changed compatibility link: {name}")
    return {"checked_files": checked,
            "checked_links": len(manifest["compatibility_links"]),
            "problems": problems}


if __name__ == "__main__":
    result = verify(Path(__file__).resolve().parent)
    print(json.dumps(result, indent=2))
    if result["problems"]:
        raise SystemExit(1)
