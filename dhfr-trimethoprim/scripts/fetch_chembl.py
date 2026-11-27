"""Pull DHFR inhibition data from ChEMBL for the selectivity benchmark.

Only Ki is kept. IC50 is deliberately excluded: by Cheng-Prusoff,

    Ki = IC50 / (1 + [S]/Km)

so an IC50 depends on the substrate concentration the assay used. Pooling IC50
values across labs that used different [S] compares numbers that are not the same
quantity. ChEMBL rarely records [S], so the conversion cannot be done after the fact.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

BASE = "https://www.ebi.ac.uk/chembl/api/data"

TARGETS = {
    "ecoli": "CHEMBL1809",   # Escherichia coli DHFR
    "human": "CHEMBL202",    # Homo sapiens DHFR
}

FIELDS = [
    "molecule_chembl_id", "canonical_smiles", "standard_type", "standard_relation",
    "standard_value", "standard_units", "pchembl_value", "assay_chembl_id",
    "target_chembl_id", "document_chembl_id", "assay_description",
]


def fetch_all(target_id: str, standard_type: str = "Ki", page: int = 1000) -> list[dict]:
    rows, offset = [], 0
    while True:
        q = urllib.parse.urlencode({
            "target_chembl_id": target_id,
            "standard_type": standard_type,
            "format": "json",
            "limit": page,
            "offset": offset,
        })
        url = f"{BASE}/activity?{q}"
        with urllib.request.urlopen(url, timeout=120) as r:
            data = json.load(r)
        batch = data.get("activities", [])
        rows.extend({k: a.get(k) for k in FIELDS} for a in batch)
        total = data["page_meta"]["total_count"]
        offset += len(batch)
        print(f"  {target_id} {standard_type}: {offset}/{total}", file=sys.stderr)
        if offset >= total or not batch:
            return rows
        time.sleep(0.2)


def main() -> None:
    out = {}
    for name, tid in TARGETS.items():
        print(f"fetching {name} ({tid})", file=sys.stderr)
        out[name] = fetch_all(tid)

    dest = RAW / "chembl_dhfr_ki.json"
    dest.write_text(json.dumps(out, indent=1))
    for name, rows in out.items():
        print(f"{name}: {len(rows)} Ki records")
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
