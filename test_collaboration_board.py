import concurrent.futures
import tempfile
import unittest
import json
import sqlite3
from pathlib import Path

from collaboration_board import audit, claim, close, connect, propose, submit


class BoardTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / 'board.sqlite3'
        self.db = connect(self.path)
        self.hid = propose(self.db, 'test', 'cx', 'bounded test', 'one run', 'root')

    def tearDown(self):
        self.db.close()
        self.temp.cleanup()

    def test_concurrent_claim_has_one_winner(self):
        def attempt(agent):
            db = connect(self.path)
            try:
                claim(db, self.hid, agent)
                return True
            except ValueError:
                return False
            finally:
                db.close()
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            self.assertEqual(sum(pool.map(attempt, ['a', 'b'])), 1)
        self.assertEqual(self.db.execute("SELECT count(*) FROM events WHERE kind='claim'").fetchone()[0], 1)

    def test_ownership_and_failed_transition_roll_back(self):
        claim(self.db, self.hid, 'a')
        with self.assertRaises(ValueError):
            close(self.db, self.hid, 'b', 'completed', 'unauthorized')
        self.assertEqual(self.db.execute('SELECT status FROM hypotheses').fetchone()[0], 'running')
        close(self.db, self.hid, 'a', 'inconclusive', 'no valid circuit')
        with self.assertRaises(ValueError):
            claim(self.db, self.hid, 'b')

    def test_dependencies_require_completed_work(self):
        child = propose(self.db, 'child', 'cx', 'uses parent', 'one run', 'root', self.hid)
        with self.assertRaises(ValueError):
            claim(self.db, child, 'b')
        claim(self.db, self.hid, 'a')
        close(self.db, self.hid, 'a', 'completed', 'construction available')
        claim(self.db, child, 'b')
        with self.assertRaises(sqlite3.IntegrityError):
            propose(self.db, 'invalid', 'cx', 'missing dependency', 'one run', 'root', 999)

    def test_submission_preserves_original_and_cannot_promote_itself(self):
        claim(self.db, self.hid, 'a')
        source = Path(self.temp.name) / 'candidate.json'
        source.write_text('{"n":4,"gates":[]}')
        with self.assertRaises(ValueError):
            submit(self.db, self.hid, 'b', source, '{}', 'numerical', 'wrong owner')
        self.assertFalse((Path(self.temp.name) / 'submissions').exists())
        sid = submit(self.db, self.hid, 'a', source, '{"compiled_cnot":0}', 'numerical', 'test')
        row = self.db.execute('SELECT * FROM submissions WHERE id=?', (sid,)).fetchone()
        source.write_text('{"changed":true}')
        self.assertEqual(json.loads(Path(row['path']).read_text())['n'], 4)
        with self.assertRaises(ValueError):
            submit(self.db, self.hid, 'a', source, '{}', 'verified', 'not a certificate')
        self.assertEqual(audit(self.db)['audited_submissions'], 1)
        Path(row['path']).write_text('{}')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            audit(self.db)


if __name__ == '__main__':
    unittest.main()
