"""Render persisted operational team messages, not reconstructed dialogue."""
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
db = sqlite3.connect('file:' + str(ROOT/'collaboration/board.sqlite3') + '?mode=ro', uri=True)
db.row_factory = sqlite3.Row
rows = db.execute("SELECT * FROM events WHERE hypothesis>=20 AND kind IN ('proposal','claim','handoff','completed','inconclusive','failed') ORDER BY id").fetchall()
db.close()
lines = ['# Researcher chats — current coordinated round', '',
         'These are verbatim messages persisted on the shared board. Event IDs and',
         'timestamps provide provenance. Assignments and outcomes are labelled separately',
         'from researcher handoffs. Private internal reasoning is not included.', '',
         'The individual researcher CHATS.md files also retain their live sent/received',
         'messages. This export does not invent dialogue or claim every tool call is a chat.', '']
for r in rows:
    target = ' → ' + r['recipient'] if r['recipient'] else ''
    lines += [f"## Event {r['id']} — {r['author']}{target}", '',
              f"{r['created_at']} · hypothesis {r['hypothesis']} · {r['kind']}", '']
    lines += ['> ' + s for s in r['message'].splitlines()]
    lines += ['']
(ROOT/'collaboration/CHATS.md').write_text('\n'.join(lines).rstrip()+'\n')
print(f'Exported {len(rows)} recorded events to collaboration/CHATS.md')
