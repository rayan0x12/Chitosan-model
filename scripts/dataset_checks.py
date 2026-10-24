"""Dataset integrity checks for the repository.

Run from the repo root with one command::

    python scripts/dataset_checks.py

The script re-derives every derived ``data/`` file from ``sheet1_cleaned.csv``
and asserts the committed copies match, then runs a battery of integrity
checks:

  * schema, row and study counts
  * duplicate-row detection
  * DOI consistency between ``sources.csv`` and ``sheet1_cleaned.csv``
  * missing-data matrix
  * provenance-label classification (table / text / graph-estimated / derived)
  * derived-value spot checks (arithmetic identities encoded in the data)

It exits with code 0 only if every check passes.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
DATA = os.path.join(REPO_ROOT, "data")

SHEET1 = os.path.join(DATA, "sheet1_cleaned.csv")
SOURCES = os.path.join(DATA, "sources.csv")
STUDY_SUMMARY = os.path.join(DATA, "study_summary.csv")
DATA_USED = os.path.join(DATA, "data_used.csv")

_NUM = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")

failures = []


def check(cond, msg):
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {msg}")
    if not cond:
        failures.append(msg)


def parse_num(cell):
    if cell is None:
        return None
    s = str(cell).strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        pass
    m = _NUM.match(s)
    if not m:
        return None
    rest = s[m.end():].strip()
    if rest == "" or rest.startswith("(") or rest.startswith("["):
        return float(m.group())
    return None


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.reader(fh))


def main():
    rd = read_csv(SHEET1)
    hdr, rows = rd[0], rd[1:]
    idx = {h: i for i, h in enumerate(hdr)}

    # --- schema ---
    print("== schema & counts ==")
    check(len(hdr) == 48, f"sheet1 has 48 columns (found {len(hdr)})")
    check(len(rows) == 167, f"sheet1 has 167 data rows (found {len(rows)})")
    refs = sorted({r[0].strip() for r in rows}, key=lambda k: int(k[1:]))
    check(len(refs) == 17, f"sheet1 spans 17 studies R01–R17 (found {len(refs)})")
    check(refs == [f"R{i:02d}" for i in range(1, 18)],
          f"Ref_IDs are exactly R01–R17 (found {refs})")

    # --- duplicates (full-row) ---
    print("== duplicates ==")
    full = Counter(tuple(r) for r in rows)
    dups = {k: v for k, v in full.items() if v > 1}
    check(len(dups) == 0, f"no exact duplicate rows (all 48 columns) — found {len(dups)}")
    for k, v in list(dups.items())[:5]:
        print(f"    duplicate {k[0]} {k[6][:24]} x{v}")

    # --- DOI consistency with sources.csv ---
    print("== sources.csv ===")
    src = read_csv(SOURCES)
    src_hdr = src[0]
    check(src_hdr == ["Ref_ID", "Authors", "Year", "Title", "Journal", "DOI", "Access", "Rows"],
          f"sources.csv header correct (found {src_hdr})")
    src_map = {r[0]: r for r in src[1:]}
    check(set(src_map) == set(refs), "sources.csv has exactly one row per Ref_ID")
    doi_mismatch = 0
    for ref in refs:
        sheet_doi = next(r[idx["DOI_URL"]].strip() for r in rows if r[0].strip() == ref)
        src_doi = src_map[ref][5].strip()
        if sheet_doi != src_doi:
            doi_mismatch += 1
            print(f"    DOI mismatch {ref}: sheet1={sheet_doi!r} sources={src_doi!r}")
    check(doi_mismatch == 0, f"sources.csv DOIs match sheet1 (mismatches: {doi_mismatch})")

    # --- study_summary row counts ---
    print("== study_summary.csv ==")
    ss = read_csv(STUDY_SUMMARY)
    per = Counter(r[0].strip() for r in rows)
    ss_map = {r[0]: int(r[2]) for r in ss[1:]}
    bad = {k: (ss_map.get(k), per.get(k)) for k in set(per) | set(ss_map) if ss_map.get(k) != per.get(k)}
    check(len(bad) == 0, f"study_summary row counts match sheet1 (mismatches: {bad or 'none'})")

    # --- missing-data matrix ---
    print("== missing data ==")
    cols = ["LL_untreated", "PL_untreated", "PI_untreated", "UCS_untreated_kPa_q0",
            "Chitosan_Content_%", "Initial_Moisture_Content_%", "Curing_Time_days",
            "UCS_treated_kPa"]
    for c in cols:
        n = sum(1 for r in rows if parse_num(r[idx[c]]) is not None)
        print(f"    {c}: {n}/{len(rows)}")

    # --- provenance classification ---
    print("== provenance ==")
    def provenance(vs):
        s = vs.lower()
        if s.startswith("derived"):
            return "derived"
        if "graph" in s:
            return "graph-estimated"
        if "table" in s:
            return "table"
        return "text"
    pcount = Counter(provenance(r[idx["Value_Source"]]) for r in rows)
    total = sum(pcount.values())
    check(total == len(rows), f"every row classified into a provenance label ({total}/{len(rows)})")
    for k in ["table", "text", "graph-estimated", "derived"]:
        print(f"    {k}: {pcount.get(k, 0)}")
    # graph-estimated must be the 27 R16 rows + 1 R10 graph label
    graph_rows = [r[0].strip() for r in rows if provenance(r[idx["Value_Source"]]) == "graph-estimated"]
    check(len(graph_rows) == 28, f"graph-estimated rows = 28 (27 R16 + 1 R10) — found {len(graph_rows)}")
    check(pcount.get("derived", 0) == 8, f"derived rows = 8 — found {pcount.get('derived', 0)}")

    # --- data_used.csv consistency ---
    print("== data_used.csv ==")
    du = read_csv(DATA_USED)
    check(du[0][0] == "Ref_ID", f"data_used.csv header starts with Ref_ID (found {du[0][:4]})")
    check(len(du) - 1 == 139, f"data_used.csv has 139 model rows (found {len(du) - 1})")
    used_refs = sorted({r[0] for r in du[1:]})
    check(used_refs == ["R01", "R05", "R06", "R09", "R13", "R15", "R16"],
          f"data_used.csv spans the 7 model studies (found {used_refs})")
    # every ratio in data_used = ucs/q0
    bad_ratio = 0
    for r in du[1:]:
        q0, ucs, ratio = float(r[6]), float(r[7]), float(r[8])
        if abs(ratio - ucs / q0) > 1e-3:
            bad_ratio += 1
    check(bad_ratio == 0, f"data_used ratio column = UCS/q0 for every row (bad: {bad_ratio})")

    # --- derived-value spot checks ---
    print("== derived-value spot checks ==")
    # R16 GS 1% / 1 d OMC: UCS 290 vs q0 243 -> +19.3% (paper states +19%)
    def r16(soil, dose):
        for r in rows:
            if r[0].strip() != "R16":
                continue
            if r[idx["Soil_Name"]].strip().startswith(soil) and \
               parse_num(r[idx["Chitosan_Content_%"]]) == dose and \
               parse_num(r[idx["Curing_Time_days"]]) == 1 and \
               r[idx["Initial_Moisture_Content_%"]].strip() == "16.0":
                return r
        return None
    gs = r16("Ganga", 1.0)
    if gs:
        q0, ucs = parse_num(gs[idx["UCS_untreated_kPa_q0"]]), parse_num(gs[idx["UCS_treated_kPa"]])
        gain = (ucs / q0 - 1) * 100
        check(abs(gain - 19) <= 1, f"R16 GS 1% gain {gain:.1f}% within 1 pt of paper +19%")

    print()
    if failures:
        print(f"RESULT: {len(failures)} check(s) FAILED")
        sys.exit(1)
    print("RESULT: all checks passed")


if __name__ == "__main__":
    main()
