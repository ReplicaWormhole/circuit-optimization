import argparse
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import experiment_log as log


class ExperimentLogTests(unittest.TestCase):
    def test_reserve_deduplicate_finish_and_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = (log.ROOT, log.DATABASE, log.CANDIDATES,
                        log.CODE_SNAPSHOTS, log.RUN_WORKSPACES)
            log.ROOT = root
            log.DATABASE = root / "experiments.sqlite3"
            log.CANDIDATES = root / "candidates"
            log.CODE_SNAPSHOTS = root / "code_snapshots"
            log.RUN_WORKSPACES = root / "experiments" / "runs"
            try:
                config = {"seed": 42, "max_iterations": 10}
                code = root / "search.py"
                code.write_text("print('search')\n")
                args = argparse.Namespace(target="V4", method="test-search", agent="test",
                                          description="A test run", config=json.dumps(config),
                                          code=code, parent=None)
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    log.reserve(args)
                run_id = json.loads(output.getvalue())["id"]
                with contextlib.redirect_stdout(output := io.StringIO()):
                    log.workspace(argparse.Namespace(id=run_id, code=code))
                self.assertTrue(json.loads(output.getvalue())["created"])
                plan = root / "experiments" / "runs" / str(run_id) / "PLAN.md"
                self.assertIn(f"experiment_log.py show {run_id}", plan.read_text())
                self.assertEqual((plan.parent / "search.py").read_bytes(), code.read_bytes())
                plan.write_text("edited plan\n")
                with contextlib.redirect_stdout(output := io.StringIO()):
                    log.workspace(argparse.Namespace(id=run_id, code=code))
                self.assertFalse(json.loads(output.getvalue())["created"])
                self.assertEqual(plan.read_text(), "edited plan\n")
                with self.assertRaisesRegex(ValueError, "unknown run id"):
                    log.workspace(argparse.Namespace(id=999, code=None))
                code.write_text("modified search\n")
                with self.assertRaisesRegex(ValueError, "does not match"):
                    log.workspace(argparse.Namespace(id=run_id, code=code))
                code.write_text("print('search')\n")
                with self.assertRaisesRegex(ValueError, "duplicate run"):
                    log.reserve(args)

                failed_args = argparse.Namespace(**{**vars(args), "config": json.dumps({
                    "seed": 43, "max_iterations": 10})})
                failed_output = io.StringIO()
                with contextlib.redirect_stdout(failed_output):
                    log.reserve(failed_args)
                failed_id = json.loads(failed_output.getvalue())["id"]
                with contextlib.redirect_stdout(io.StringIO()):
                    log.finish(argparse.Namespace(id=failed_id, status="failed",
                                                  candidate=None, evidence_level="numerical",
                                                  notes="optimizer crashed"))

                candidate = root / "candidate.json"
                candidate.write_text(json.dumps({"n": 4, "gates": []}))
                end = argparse.Namespace(id=run_id, status="inconclusive",
                                         candidate=candidate, evidence_level="numerical",
                                         notes="identity does not diagonalize V4")
                with contextlib.redirect_stdout(io.StringIO()):
                    log.finish(end)
                with self.assertRaisesRegex(ValueError, "already inconclusive"):
                    log.finish(end)
                with log.connect() as db:
                    row = db.execute("SELECT * FROM attempts WHERE id = ?", (run_id,)).fetchone()
                self.assertEqual(row["cnot_count"], 0)
                self.assertEqual(row["valid_cycle"], 0)
                self.assertTrue((root / "code_snapshots" / row["code_sha256"]).is_file())
                archived = root / row["candidate_path"]
                self.assertTrue(archived.is_file())
                self.assertEqual(json.loads(archived.read_text()), {"n": 4, "gates": []})
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    log.audit(None)
                self.assertEqual(json.loads(output.getvalue())["problems"], [])
                self.assertEqual(json.loads(output.getvalue())["audited_runs"], 2)
                archived.write_text("changed")
                with contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        log.audit(None)
            finally:
                (log.ROOT, log.DATABASE, log.CANDIDATES,
                 log.CODE_SNAPSHOTS, log.RUN_WORKSPACES) = original


if __name__ == "__main__":
    unittest.main()
