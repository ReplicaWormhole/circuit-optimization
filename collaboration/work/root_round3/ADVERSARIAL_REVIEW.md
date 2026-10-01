# Coordinator adversarial closeout

Scope: code, math, computation, explicit-circuit transcription and workflow.
The current work changes evidence and next actions; it does not complete the
persistent optimum objective. Requirements are retained in GOAL_AUDIT.md.

## Confirmed issues corrected

- Initial architecture motivation claimed one mutation could form an aba
  precursor. Two are needed: [01,32,12] to [01,12,01]. Reports and actual chat
  corrections retain this distinction. Frozen direct12 proposal is not changed
  retroactively; run374 explicitly uses the corrected bridge.
- Readable circuit initially omitted the meaning of numeric D labels. Exporter
  and artifact now state0,1,2,3 mean1,i,-1,-i. No gate/metadata change.

## Checks and scope

- All22 focused base/native tests pass. Aggregate run emitted SQLite
  unclosed-connection ResourceWarnings from the existing test suite. No test
  or database helper behavior was changed in this work; warnings are recorded
  separately from successful assertions.
- Independent parsing compares all73 rendered rows and exact rotation metadata,
  controls/targets, source hash and output-label root mapping against the
  accepted prescription. Pass. This is transcription evidence, not a second
  exact circuit implementation.
- Decorated identity reviewed in chronological matrix order, with both interior
  commutations required. Product-local centralizer proof was independently
  confirmed; failed optimizer runs do not exclude that class.
- Joint halfshift proof uses distinct input wires, nonzero off-diagonal axis
  components, consistent inverse-diagonalizer conjugation and multiplicities
  (10,6). Valid necessary condition; no surviving topology excluded.
- Two runs reserved/snapshotted before compute, finite seeds/topologies and
  one-thread limits recorded. Run373 elapsed5.864 seconds;374 elapsed10.649,
  both within their declared240-second budgets. All saved native gate lists
  independently rejected. Bridge complete matrices agree up to phase to
  5.36e-16 and3.56e-16; this verifies conversion, not diagonalization.
- Ledger audit374 runs,189 archived candidates, no problems; no running rows.
  Board auditfour submission hashes/JSONs valid. Database audits do not prove
  mathematical claims or liveness.
- Actual operational chats are preserved with sender/recipient/provenance;
  erroneous initial messages remain alongside correction messages. No invented
  dialogue or private internal reasoning is supplied.

No remaining confirmed defect. Current exact upper13 and global lower6 leave
optimality unproved. Old exact13 implementation caveats and gate-set distinctions
remain. Keep the goal active, not complete or blocked.

All new research JSON, scripts, reports, chats, code snapshots and SQLite files
are intentional retained artifacts. Caches/journals ignored. No dependencies,
destructive cleanup, remote modifications or pushes. Focused local commit after
review; repository state checked afterward.
