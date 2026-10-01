"""Transactional hypothesis board; scientific runs remain in experiment_log.py."""
import argparse
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_DB = ROOT / 'collaboration' / 'board.sqlite3'


def connect(path=DEFAULT_DB):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=30)
    db.row_factory = sqlite3.Row
    db.execute('PRAGMA foreign_keys=ON')
    db.executescript('''
        CREATE TABLE IF NOT EXISTS hypotheses (
          id INTEGER PRIMARY KEY, title TEXT NOT NULL, gate_set TEXT NOT NULL,
          proposal TEXT NOT NULL, budget TEXT NOT NULL, dependency INTEGER
          REFERENCES hypotheses(id), status TEXT NOT NULL DEFAULT 'proposed',
          owner TEXT, run_id INTEGER, outcome TEXT);
        CREATE TABLE IF NOT EXISTS events (
          id INTEGER PRIMARY KEY, hypothesis INTEGER REFERENCES hypotheses(id),
          created_at TEXT NOT NULL, author TEXT NOT NULL, kind TEXT NOT NULL,
          recipient TEXT, message TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS submissions (
          id INTEGER PRIMARY KEY, hypothesis INTEGER NOT NULL REFERENCES hypotheses(id),
          gate_set TEXT NOT NULL, author TEXT NOT NULL, path TEXT NOT NULL,
          sha256 TEXT NOT NULL, evidence TEXT NOT NULL, costs TEXT NOT NULL,
          notes TEXT NOT NULL);
    ''')
    return db


def event(db, hid, author, kind, message, recipient=None):
    db.execute('INSERT INTO events(hypothesis,created_at,author,kind,recipient,message) '
               'VALUES(?,?,?,?,?,?)',
               (hid, datetime.now(timezone.utc).isoformat(), author, kind, recipient, message))


def require_owner(db, hid, author):
    row = db.execute('SELECT * FROM hypotheses WHERE id=?', (hid,)).fetchone()
    if row is None or row['status'] != 'running' or row['owner'] != author:
        raise ValueError('hypothesis must be running and owned by this agent')
    return row


def propose(db, title, gate_set, proposal, budget, author, dependency=None):
    with db:
        cursor = db.execute('INSERT INTO hypotheses(title,gate_set,proposal,budget,dependency) '
                            'VALUES(?,?,?,?,?)', (title, gate_set, proposal, budget, dependency))
        event(db, cursor.lastrowid, author, 'proposal', proposal)
    return cursor.lastrowid


def claim(db, hid, author):
    # Lock before reading dependencies; concurrent claims have exactly one winner.
    db.execute('BEGIN IMMEDIATE')
    try:
        row = db.execute('SELECT * FROM hypotheses WHERE id=?', (hid,)).fetchone()
        if row is None or row['status'] != 'proposed':
            raise ValueError('hypothesis unavailable')
        if row['dependency'] is not None:
            dep = db.execute('SELECT status FROM hypotheses WHERE id=?', (row['dependency'],)).fetchone()
            if dep['status'] != 'completed':
                raise ValueError('dependency has not completed')
        db.execute("UPDATE hypotheses SET status='running',owner=? WHERE id=?", (author, hid))
        event(db, hid, author, 'claim', 'claimed')
        db.commit()
    except Exception:
        db.rollback()
        raise


def close(db, hid, author, status, outcome):
    if status not in ('completed', 'inconclusive', 'blocked') or not outcome.strip():
        raise ValueError('close needs a terminal status and a nonempty outcome')
    with db:
        db.execute('BEGIN IMMEDIATE')
        require_owner(db, hid, author)
        db.execute('UPDATE hypotheses SET status=?,outcome=? WHERE id=?', (status, outcome, hid))
        event(db, hid, author, status, outcome)


def submit(db, hid, author, path, costs, evidence, notes):
    if evidence not in ('numerical', 'exact_claim'):
        raise ValueError('submissions are proposals, never automatically verified incumbents')
    data = Path(path).resolve().read_bytes()
    json.loads(data)  # Ensure the submission is a JSON artifact.
    costs = json.dumps(json.loads(costs), sort_keys=True, allow_nan=False)
    sha = hashlib.sha256(data).hexdigest()
    with db:
        db.execute('BEGIN IMMEDIATE')
        row = require_owner(db, hid, author)
        directory = Path(db.execute('PRAGMA database_list').fetchone()['file']).parent / 'submissions'
        directory.mkdir(parents=True, exist_ok=True)
        dest = directory / (sha + '.json')
        # Exclusive creation prevents simultaneous writers overwriting a snapshot.
        try:
            with dest.open('xb') as stream:
                stream.write(data)
        except FileExistsError:
            if dest.read_bytes() != data:
                raise ValueError('content-addressed artifact differs from its hash')
        try:
            stored_path = str(dest.relative_to(ROOT))
        except ValueError:
            stored_path = str(dest)
        cursor = db.execute('INSERT INTO submissions(hypothesis,gate_set,author,path,sha256,'
                            'evidence,costs,notes) VALUES(?,?,?,?,?,?,?,?)',
                            (hid, row['gate_set'], author, stored_path, sha,
                             evidence, costs, notes))
        event(db, hid, author, 'submission', str(cursor.lastrowid))
    return cursor.lastrowid


def audit(db):
    rows = db.execute('SELECT * FROM submissions').fetchall()
    for row in rows:
        path = Path(row['path'])
        if not path.is_absolute():
            path = ROOT / path
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ValueError(f"submission {row['id']}: content hash mismatch")
        json.loads(data)
    return {'audited_submissions': len(rows), 'scope': 'content hashes and JSON only; no proof validation'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, default=DEFAULT_DB)
    subs = parser.add_subparsers(dest='command', required=True)
    p = subs.add_parser('propose')
    for field in ('title', 'gate-set', 'proposal', 'budget', 'agent'):
        p.add_argument('--' + field, required=True)
    p.add_argument('--depends-on', type=int)
    for name in ('claim', 'close', 'handoff', 'link-run', 'submit'):
        p = subs.add_parser(name)
        p.add_argument('id', type=int)
        p.add_argument('--agent', required=True)
        if name == 'close':
            p.add_argument('--status', choices=['completed', 'inconclusive', 'blocked'], required=True)
            p.add_argument('--outcome', required=True)
        elif name == 'handoff':
            p.add_argument('--to', required=True)
            p.add_argument('--message', required=True)
        elif name == 'link-run':
            p.add_argument('--run', type=int, required=True)
        elif name == 'submit':
            p.add_argument('--candidate', type=Path, required=True)
            p.add_argument('--costs', required=True, help='JSON native/depth/CNOT cost record')
            p.add_argument('--evidence', choices=['numerical', 'exact_claim'], required=True)
            p.add_argument('--notes', required=True)
    subs.add_parser('list')
    subs.add_parser('events')
    subs.add_parser('submissions')
    subs.add_parser('audit')
    args = parser.parse_args()
    db = connect(args.db)
    try:
        if args.command == 'propose':
            print(propose(db, args.title, args.gate_set, args.proposal, args.budget,
                          args.agent, args.depends_on))
        elif args.command == 'claim':
            claim(db, args.id, args.agent)
        elif args.command == 'close':
            close(db, args.id, args.agent, args.status, args.outcome)
        elif args.command == 'submit':
            print(submit(db, args.id, args.agent, args.candidate, args.costs,
                         args.evidence, args.notes))
        elif args.command == 'link-run':
            with sqlite3.connect('file:' + str(ROOT / 'experiments.sqlite3') + '?mode=ro', uri=True) as ledger:
                run = ledger.execute('SELECT agent FROM attempts WHERE id=?', (args.run,)).fetchone()
            if run is None or run[0] != args.agent:
                raise ValueError('run must exist and belong to this agent')
            with db:
                db.execute('BEGIN IMMEDIATE')
                row = require_owner(db, args.id, args.agent)
                if row['run_id'] is not None:
                    raise ValueError('run already linked; create a new bounded hypothesis')
                db.execute('UPDATE hypotheses SET run_id=? WHERE id=?', (args.run, args.id))
                event(db, args.id, args.agent, 'run', str(args.run))
        elif args.command == 'handoff':
            with db:
                event(db, args.id, args.agent, 'handoff', args.message, args.to)
        elif args.command == 'audit':
            print(json.dumps(audit(db)))
        else:
            table = {'list': 'hypotheses', 'events': 'events', 'submissions': 'submissions'}[args.command]
            for row in db.execute('SELECT * FROM ' + table + ' ORDER BY id'):
                print(json.dumps(dict(row)))
    except (ValueError, OSError, sqlite3.Error) as error:
        parser.error(str(error))
    finally:
        db.close()


if __name__ == '__main__':
    main()
