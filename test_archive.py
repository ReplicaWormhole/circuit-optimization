"""Focused checks for the legacy archive integrity audit."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from archive.verify import verify


class ArchiveVerificationTests(unittest.TestCase):
    def test_detects_changed_content_and_compatibility_link(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            legacy = base / "legacy_root"
            legacy.mkdir()
            (base / "data").mkdir()
            source = legacy / "saved.json"
            source.write_bytes(b"{}\n")
            (legacy / "data").symlink_to("../data", target_is_directory=True)
            manifest = {
                "schema": 1,
                "moved_root_files": {"saved.json": {
                    "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                    "size": source.stat().st_size}},
                "frozen_helper_copies": {},
                "compatibility_links": {"data": "../data"},
            }
            (base / "manifest.json").write_text(json.dumps(manifest))
            self.assertEqual(verify(base)["problems"], [])
            source.write_bytes(b"changed\n")
            (legacy / "data").unlink()
            (legacy / "data").symlink_to("../other", target_is_directory=True)
            self.assertEqual(verify(base)["problems"], [
                "changed content: saved.json", "changed compatibility link: data"])


if __name__ == "__main__":
    unittest.main()
