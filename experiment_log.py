"""Concurrent local ledger for four-qubit circuit-search attempts.

Reserve a run before computing; finish it once a candidate or failure is known.
The SQLite database is the authoritative attempt index. Candidate JSON is
archived under candidates/ by content hash, so later edits cannot alter evidence.
"""

import argparse
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from check_circuit import evaluate


ROOT = Path(__file__).resolve().parent
DATABASE = ROOT / "experiments.sqlite3"
CANDIDATES = ROOT / "candidates"
CODE_SNAPSHOTS = ROOT / "code_snapshots"
RUN_WORKSPACES = ROOT / "experiments" / "runs"


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def connect():
    db = sqlite3.connect(DATABASE, timeout=30)
    db.row_factory = sqlite3.Row
    db.execute("""CREATE TABLE IF NOT EXISTS attempts (
        id INTEGER PRIMARY KEY,
        run_key TEXT NOT NULL UNIQUE,
        created_at TEXT NOT NULL,
        finished_at TEXT,
        target TEXT NOT NULL,
        method TEXT NOT NULL,
        agent TEXT NOT NULL,
        description TEXT NOT NULL,
        config_json TEXT NOT NULL,
        code_sha256 TEXT,
        parent_id INTEGER REFERENCES attempts(id),
        status TEXT NOT NULL,
        candidate_sha256 TEXT,
        candidate_path TEXT,
        topology_sha256 TEXT,
        cnot_count INTEGER,
        off_diagonal_error REAL,
        total_spin_error REAL,
        valid_cycle INTEGER,
        valid_spin INTEGER,
        evidence_level TEXT,
        notes TEXT
    )""")
    return db


def reserve(args):
    config = json.loads(args.config)
    if not isinstance(config, dict):
        raise ValueError("--config must be a JSON object; include seed and search limits")
    code_data = args.code.read_bytes() if args.code else None
    code_hash = digest(code_data) if code_data is not None else None
    if code_data is not None:
        CODE_SNAPSHOTS.mkdir(exist_ok=True)
        snapshot = CODE_SNAPSHOTS / code_hash
        try:
            with snapshot.open("xb") as stream:
                stream.write(code_data)
        except FileExistsError:
            if snapshot.read_bytes() != code_data:
                raise ValueError("archived code hash collision or modification")
    key_data = {"target": args.target, "method": args.method,
                "config": config, "code_sha256": code_hash, "parent_id": args.parent}
    run_key = digest(canonical(key_data).encode())
    with connect() as db:
        try:
            cursor = db.execute("""INSERT INTO attempts
                (run_key, created_at, target, method, agent, description,
                 config_json, code_sha256, parent_id, status)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 'running')""",
                (run_key, now(), args.target, args.method, args.agent,
                 args.description, canonical(config), code_hash, args.parent))
        except sqlite3.IntegrityError as error:
            existing = db.execute("SELECT id, status FROM attempts WHERE run_key = ?",
                                  (run_key,)).fetchone()
            if existing:
                raise ValueError(f"duplicate run: id={existing['id']} status={existing['status']}") from error
            raise
    print(json.dumps({"id": cursor.lastrowid, "run_key": run_key, "status": "running"}))


def finish(args):
    with connect() as db:
        row = db.execute("SELECT * FROM attempts WHERE id = ?", (args.id,)).fetchone()
        if row is None:
            raise ValueError(f"unknown run id {args.id}")
        if row["status"] != "running":
            raise ValueError(f"run {args.id} is already {row['status']}")
        fields = {"finished_at": now(), "status": args.status,
                  "evidence_level": args.evidence_level, "notes": args.notes}
        if args.candidate:
            candidate = json.loads(args.candidate.read_text())
            if candidate.get("n") != 4:
                raise ValueError("candidate must use four qubits")
            candidate_data = (canonical(candidate) + "\n").encode()
            candidate_hash = digest(candidate_data)
            CANDIDATES.mkdir(exist_ok=True)
            archive = CANDIDATES / f"{candidate_hash}.json"
            try:
                with archive.open("xb") as stream:
                    stream.write(candidate_data)
            except FileExistsError:
                if archive.read_bytes() != candidate_data:
                    raise ValueError("archived candidate hash collision or modification")
            result = evaluate(candidate)
            topology = [(gate["control"], gate["target"])
                        for gate in candidate["gates"] if gate["gate"].lower() == "cx"]
            fields.update(candidate_sha256=candidate_hash,
                          candidate_path=str(archive.relative_to(ROOT)),
                          topology_sha256=digest(canonical(topology).encode()),
                          cnot_count=result["cnot_count"],
                          off_diagonal_error=result["off_diagonal_error"],
                          total_spin_error=result["total_spin_off_diagonal_error"],
                          valid_cycle=int(result["valid_diagonalizer"]),
                          valid_spin=int(result["diagonalizes_total_spin"]))
        columns = ", ".join(f"{name} = ?" for name in fields)
        updated = db.execute(f"UPDATE attempts SET {columns} WHERE id = ? AND status = 'running'",
                             (*fields.values(), args.id))
        if updated.rowcount != 1:
            raise ValueError(f"run {args.id} was finished concurrently")
    print(json.dumps({"id": args.id, "status": args.status, **fields}))


def show(args):
    with connect() as db:
        row = db.execute("SELECT * FROM attempts WHERE id = ?", (args.id,)).fetchone()
    if row is None:
        raise ValueError(f"unknown run id {args.id}")
    result = dict(row)
    result["config"] = json.loads(result.pop("config_json"))
    print(json.dumps(result, indent=2))


def list_attempts(args):
    with connect() as db:
        rows = db.execute("""SELECT id, status, target, method, agent, cnot_count,
            off_diagonal_error, valid_cycle, valid_spin, topology_sha256
            FROM attempts ORDER BY id DESC""").fetchall()
    for row in rows:
        print(json.dumps(dict(row)))


def export_attempts(_args):
    with connect() as db:
        rows = db.execute("SELECT * FROM attempts ORDER BY id").fetchall()
    for row in rows:
        result = dict(row)
        result["config"] = json.loads(result.pop("config_json"))
        print(json.dumps(result, sort_keys=True))


def summary(_args):
    with connect() as db:
        total = db.execute("SELECT COUNT(*) FROM attempts").fetchone()[0]
        running = db.execute("SELECT COUNT(*) FROM attempts WHERE status = 'running'").fetchone()[0]
        topologies = db.execute("SELECT COUNT(DISTINCT topology_sha256) FROM attempts").fetchone()[0]
        best = db.execute("""SELECT id, cnot_count, candidate_path
            FROM attempts WHERE valid_cycle = 1 ORDER BY cnot_count, id LIMIT 1""").fetchone()
    print(json.dumps({"attempts": total, "running": running,
                      "candidate_topologies": topologies,
                      "best_valid_cycle": dict(best) if best else None}, indent=2))


def workspace(args):
    """Create an editable home for a reserved run without changing the ledger."""
    with connect() as db:
        row = db.execute("SELECT id, description, code_sha256 FROM attempts WHERE id = ?",
                         (args.id,)).fetchone()
    if row is None:
        raise ValueError(f"unknown run id {args.id}; reserve the run first")

    code_data = None
    if args.code is not None:
        code_data = args.code.read_bytes()
        if row["code_sha256"] is None or digest(code_data) != row["code_sha256"]:
            raise ValueError("--code does not match the reserved code snapshot")

    directory = RUN_WORKSPACES / str(args.id)
    try:
        directory.mkdir(parents=True)
        created = True
    except FileExistsError:
        if not directory.is_dir():
            raise ValueError(f"workspace path is not a directory: {directory}")
        created = False
    plan = directory / "PLAN.md"
    template = (
        f"# Run {args.id}: {row['description']}\n\n"
        f"Ledger record: `python3 experiment_log.py show {args.id}`\n\n"
        "## Hypothesis and scope\n\n"
        "## Frozen inputs and limits\n\n"
        "## Reproduction command\n\n"
        "## Outputs and independent checks\n\n"
        "## Conclusion and limitations\n"
    )
    try:
        with plan.open("x") as stream:
            stream.write(template)
    except FileExistsError:
        pass
    if code_data is not None:
        destination = directory / args.code.name
        try:
            with destination.open("xb") as stream:
                stream.write(code_data)
        except FileExistsError:
            if destination.read_bytes() != code_data:
                raise ValueError(f"workspace code differs from reserved snapshot: {destination}")
    print(json.dumps({"id": args.id, "workspace": str(directory.relative_to(ROOT)),
                      "created": created}))


def audit(_args):
    with connect() as db:
        rows = db.execute("SELECT * FROM attempts").fetchall()
    problems = []
    checked_candidates = 0
    for row in rows:
        key_data = {"target": row["target"], "method": row["method"],
                    "config": json.loads(row["config_json"]),
                    "code_sha256": row["code_sha256"], "parent_id": row["parent_id"]}
        if digest(canonical(key_data).encode()) != row["run_key"]:
            problems.append(f"run {row['id']}: run key does not match configuration")
        if row["code_sha256"]:
            snapshot = CODE_SNAPSHOTS / row["code_sha256"]
            if not snapshot.is_file() or digest(snapshot.read_bytes()) != row["code_sha256"]:
                problems.append(f"run {row['id']}: code snapshot missing or changed")
        if row["candidate_path"] is None:
            continue
        checked_candidates += 1
        path = ROOT / row["candidate_path"]
        if not path.is_file():
            problems.append(f"run {row['id']}: missing candidate {path}")
            continue
        raw = path.read_bytes()
        if digest(raw) != row["candidate_sha256"]:
            problems.append(f"run {row['id']}: candidate content hash changed")
            continue
        result = evaluate(json.loads(raw))
        candidate = json.loads(raw)
        topology = [(gate["control"], gate["target"])
                    for gate in candidate["gates"] if gate["gate"].lower() == "cx"]
        if digest(canonical(topology).encode()) != row["topology_sha256"]:
            problems.append(f"run {row['id']}: topology hash changed")
        if result["cnot_count"] != row["cnot_count"]:
            problems.append(f"run {row['id']}: CNOT count changed")
        if bool(result["valid_diagonalizer"]) != bool(row["valid_cycle"]):
            problems.append(f"run {row['id']}: cycle validity changed")
        if bool(result["diagonalizes_total_spin"]) != bool(row["valid_spin"]):
            problems.append(f"run {row['id']}: spin validity changed")
        if abs(result["off_diagonal_error"] - row["off_diagonal_error"]) > 1e-12:
            problems.append(f"run {row['id']}: off-diagonal residual changed")
    print(json.dumps({"audited_runs": len(rows), "audited_candidates": checked_candidates,
                      "problems": problems}, indent=2))
    if problems:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    start = subparsers.add_parser("reserve")
    start.add_argument("--target", choices=["V4", "V4+S2"], default="V4")
    start.add_argument("--method", required=True)
    start.add_argument("--agent", required=True)
    start.add_argument("--description", required=True)
    start.add_argument("--config", required=True, help="JSON object with seed and limits")
    start.add_argument("--code", type=Path)
    start.add_argument("--parent", type=int)
    start.set_defaults(action=reserve)
    end = subparsers.add_parser("finish")
    end.add_argument("id", type=int)
    end.add_argument("--status", choices=["complete", "failed", "inconclusive"], required=True)
    end.add_argument("--candidate", type=Path)
    end.add_argument("--evidence-level", choices=["numerical", "symbolic", "analytic"],
                     default="numerical")
    end.add_argument("--notes", default="")
    end.set_defaults(action=finish)
    detail = subparsers.add_parser("show")
    detail.add_argument("id", type=int)
    detail.set_defaults(action=show)
    subparsers.add_parser("list").set_defaults(action=list_attempts)
    subparsers.add_parser("export").set_defaults(action=export_attempts)
    subparsers.add_parser("summary").set_defaults(action=summary)
    run_space = subparsers.add_parser("workspace", help="create a directory for a reserved run")
    run_space.add_argument("id", type=int)
    run_space.add_argument("--code", type=Path, help="copy and verify the reserved search script")
    run_space.set_defaults(action=workspace)
    subparsers.add_parser("audit").set_defaults(action=audit)
    args = parser.parse_args()
    try:
        args.action(args)
    except (ValueError, OSError, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
