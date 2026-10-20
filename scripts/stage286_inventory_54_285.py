#!/usr/bin/env python3
"""Stage 286: read-only inventory of results/ JSON records for stages >= 54.

Extracts status / claim-like text / blocker / dependency metadata from every
result JSON produced by stages 54-285 and writes a single machine-readable
inventory plus a compact aggregate summary.  Performs no computation on the
scientific content and never writes to any pre-existing file.
"""
import json, os, re, sys, hashlib
from collections import defaultdict

ROOT = os.path.expanduser("~/projects/recomb2027_novelty")
RES = os.path.join(ROOT, "results")

STAGE_RE = re.compile(r"^stage(\d+)_(.+)\.json$")
TEXT_KEYS = ("title", "summary", "claim", "statement", "theorem", "lemma",
             "finding", "conclusion", "verdict", "decision", "objective",
             "question", "goal", "description", "result", "headline",
             "scope", "purpose", "topic", "name")

def first_text(obj, depth=0):
    """Collect claim-like text fields from the top two levels of a record."""
    out = {}
    if not isinstance(obj, dict) or depth > 1:
        return out
    for k, v in obj.items():
        lk = k.lower()
        if any(t == lk or lk.endswith("_" + t) or lk.startswith(t + "_") for t in TEXT_KEYS):
            if isinstance(v, str) and v.strip():
                out[k] = v.strip()[:400]
            elif isinstance(v, list) and v and all(isinstance(x, str) for x in v):
                out[k] = " | ".join(v)[:400]
    return out

records = []
bad = []
for fn in sorted(os.listdir(RES)):
    m = STAGE_RE.match(fn)
    if not m:
        continue
    stage = int(m.group(1))
    if stage < 54:
        continue
    stem = m.group(2)
    path = os.path.join(RES, fn)
    size = os.path.getsize(path)
    rec = {"stage": stage, "stem": stem, "file": fn, "size": size}
    try:
        with open(path) as fh:
            obj = json.load(fh)
    except Exception as e:
        bad.append((fn, repr(e)[:120]))
        rec["parse_error"] = repr(e)[:120]
        records.append(rec)
        continue
    if isinstance(obj, dict):
        rec["keys"] = sorted(obj.keys())[:40]
        st = obj.get("status")
        rec["status"] = st if isinstance(st, str) else (json.dumps(st)[:120] if st is not None else None)
        bl = obj.get("blockers")
        if isinstance(bl, list):
            rec["n_blockers"] = len(bl)
            rec["blockers"] = [json.dumps(b)[:300] for b in bl[:5]]
        elif bl is not None:
            rec["n_blockers"] = 1
            rec["blockers"] = [json.dumps(bl)[:300]]
        else:
            rec["n_blockers"] = 0
        ch = obj.get("checks")
        if isinstance(ch, list):
            rec["n_checks"] = len(ch)
            sts = [x.get("status") for x in ch if isinstance(x, dict)]
            rec["check_statuses"] = dict((s, sts.count(s)) for s in set(sts) if s)
        deps = obj.get("dependencies_sha256") or obj.get("dependencies")
        if isinstance(deps, dict):
            rec["deps"] = sorted(deps.keys())[:40]
        rec["text"] = first_text(obj)
        # nested one level for text
        nested = {}
        for k, v in obj.items():
            if isinstance(v, dict):
                t = first_text(v, 1)
                if t:
                    nested[k] = t
        if nested:
            rec["text_nested"] = nested
    elif isinstance(obj, list):
        rec["keys"] = ["<list:%d>" % len(obj)]
        if obj and isinstance(obj[0], dict):
            rec["text"] = first_text(obj[0])
    records.append(rec)

out = {"root": ROOT, "n_records": len(records), "parse_failures": bad, "records": records}
dest = os.path.join(RES, "stage286_inventory_54_285.json")
assert not os.path.exists(dest), "refusing to overwrite " + dest
with open(dest, "w") as fh:
    json.dump(out, fh, indent=1, sort_keys=True)

# ---- aggregates ----
by_stage = defaultdict(list)
for r in records:
    by_stage[r["stage"]].append(r["stem"])
fam_re = re.compile(r"(_case[_\-]?\d+|_\d+|_part\d+|_shard\d+|_chunk\d+)$")
fams = defaultdict(lambda: {"n": 0, "stages": set()})
for r in records:
    f = fam_re.sub("", r["stem"])
    fams[f]["n"] += 1
    fams[f]["stages"].add(r["stage"])
statuses = defaultdict(int)
for r in records:
    statuses[(r.get("status") or "<none>")[:80]] += 1

print("N_RECORDS", len(records))
print("N_STAGES", len(by_stage), "range", min(by_stage), max(by_stage))
print("PARSE_FAILURES", len(bad))
print("TOTAL_BYTES", sum(r["size"] for r in records))
print("N_FAMILIES", len(fams))
print("N_WITH_BLOCKERS", sum(1 for r in records if r.get("n_blockers")))
print("---TOP_STATUSES---")
for s, n in sorted(statuses.items(), key=lambda x: -x[1])[:40]:
    print(n, "\t", s)
print("---DEST---", dest, os.path.getsize(dest))

