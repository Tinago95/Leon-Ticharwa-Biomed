# Decisions

Every judgement call and the reason for it. Newest last.

---

## 2026-11-27 — Human target is P00374, not P00382

**Decision.** Use `P00374` as the human DHFR comparator. Keep `P00382` in `data/raw/` but
do not use it as the human target.

**Why.** The project plan's page 6 table gives `P00382` as human DHFR, 187 aa. It is not:
UniProt returns *Dihydrofolate reductase type 1*, *Escherichia coli*, 157 aa, gene `dhfrI` —
a plasmid-borne trimethoprim-resistance enzyme. `P00374` is *Homo sapiens* DHFR at the
187 aa the plan expects.

Left unfixed this would have compared *E. coli* against drug-resistant *E. coli* and
reported it as bacterial-vs-human selectivity. The run would have completed and the
numbers would have looked reasonable, which is what makes it dangerous.

**Kept anyway.** `P00382` is a real resistance enzyme, so it is useful for the page 10
resistance stretch task — just not here.

See [`identifiers.md`](identifiers.md) and `notebooks/00_verify_sources.ipynb`.

---

## 2026-11-27 — Redocking gate uses 2W9G only

**Decision.** Run the Phase 4b validation gate against `2W9G`. Do not substitute `2W9H`.

**Why.** The plan names "2W9G / 2W9H" interchangeably as *S. aureus* DHFR with
trimethoprim and NADPH. Only `2W9G` has the cofactor — `2W9H` contains trimethoprim and
1,2-ethanediol and no NADPH. The Phase 4a recipe retains NADPH in the receptor, so 2W9H
cannot satisfy it.

---

## 2026-11-27 — Environment is a venv, not conda

**Decision.** `.venv` at the repo root, dependencies in
[`requirements.txt`](../requirements.txt), rather than the `environment.yml` the plan
specifies.

**Why.** No conda, mamba, pyenv, poetry or uv on this machine — only Homebrew Python
3.13.2. A venv reproduces the pinned-version discipline the plan is actually asking for.

**Open.** AutoDock Vina, Meeko and ProLIF are not yet installed and are not pip-only on
macOS ARM. Resolve before Phase 3.

---

## Still to decide

- Trimethoprim protonation state for the main runs (plan recommends protonated; run
  neutral as comparison and record the score and pose change).
- What to do with `MN` and `BME` in `1RX2` — the plan's receptor recipe does not mention
  them. Default is to strip both as crystallisation additives.
- Docking box centre: compute from the crystal ligand centroid, do not guess.
- Random seeds and exhaustiveness for the sensitivity repeats.
