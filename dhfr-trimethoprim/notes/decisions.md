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

Installed and verified: RDKit 2026.09.2, BioPython 1.88, Meeko 0.8.0, ProLIF 2.2.3,
AutoDock Vina 1.2.7.

---

## 2026-11-27 — Vina 1.2.7 built from source, not the 1.1.2 binary

**Decision.** Build AutoDock Vina **1.2.7** from the `ccsb-scripps/AutoDock-Vina` source
release and install to `~/.local/bin`. Do not use the downloaded 1.1.2 macOS binary.

**Why.** Three independent reasons, any one sufficient:

1. **macOS deleted the 1.1.2 binary.** The Catalina x86_64 build from 2019 is unsigned and
   Chrome-quarantined. It vanished from `~/Downloads` between inspection and execution —
   not moved to Trash, just gone. `vina_split` survived; `vina` did not. Re-extracting
   would repeat it.
2. **It is the wrong version.** The plan cites Eberhardt et al. 2021, which is Vina **1.2**.
   1.1.2 is Trott & Olson 2010.
3. **The downloaded tarball contains binaries only, no source**, so "build 1.1.2 from
   source" was not possible without fetching source separately — and 1.1.2's source
   targets Boost 1.41 (2009), which will not compile cleanly against the installed
   Boost 1.87.

**How.** `brew install boost` (1.87.0), then `make` in `build/mac/release` of the 1.2.7
source — its Makefile already resolves `brew --prefix`. Produces a native **arm64**
binary, no Rosetta. Build steps are in [`../scripts/build_vina.sh`](../scripts/build_vina.sh).

**Note.** `pip install vina` does **not** work here: its `setup.py` searches only
`conda env`, `/usr/local/include` and `/usr/include`, and ignores `BOOST_ROOT`, `CPATH`
and `CPPFLAGS`. Homebrew on Apple Silicon is `/opt/homebrew`, so it never finds Boost.
The command-line executable is what the plan actually invokes, so the Python bindings
are not needed.

---

## 2026-11-27 — Phase 4b command needs rewriting for Vina 1.2

**Decision.** Drop `--log` from the docking commands; redirect stdout instead.

**Why.** The plan's Phase 4b command ends `--out results/redock_2w9g.pdbqt --log
results/redock_2w9g.log`. **`--log` was removed in Vina 1.2** — zero occurrences in
`vina --help`. The command as written fails. Use:

```
vina --receptor ... --ligand ... --out results/redock_2w9g.pdbqt > results/redock_2w9g.log
```

This is the plan citing Vina 1.2 while using 1.1.2 syntax.

---

## Open question — rotatable-bond count

The plan's Phase 3 says to "check the output reports the expected rotatable bonds (the
methylene bridge and methoxy groups - roughly five to six)".

A smoke test through the real path — RDKit embed from SMILES, then Meeko 0.8.0 — gives
**7 rotatable branches**, not five or six. Three methoxy torsions plus the methylene
bridge plus the ring-bridge bond accounts for it. Not an error in the chemistry, but the
plan's stated range will look like a failure when it is not. Expect 7.

---

## Still to decide

- Trimethoprim protonation state for the main runs (plan recommends protonated; run
  neutral as comparison and record the score and pose change).
- What to do with `MN` and `BME` in `1RX2` — the plan's receptor recipe does not mention
  them. Default is to strip both as crystallisation additives.
- Docking box centre: compute from the crystal ligand centroid, do not guess.
- Random seeds and exhaustiveness for the sensitivity repeats.
