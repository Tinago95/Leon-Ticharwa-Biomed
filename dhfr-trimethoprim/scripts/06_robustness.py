"""Stress-test the selectivity benchmark.

Reads results/benchmark_compounds.csv (written by notebook 05) and asks whether any
of its conclusions survive scrutiny. Runs in seconds: no docking is repeated.

    python scripts/06_robustness.py

Writes results/robustness.csv.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from scipy import stats

RDLogger.DisableLog("rdApp.*")

ROOT = pathlib.Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
REF = "CHEMBL22"          # trimethoprim
SIM_CUT = 0.4             # "similar to the reference"
rng = np.random.default_rng(0)

findings: list[dict] = []


def record(question: str, answer: str, number: str = "") -> None:
    findings.append({"question": question, "number": number, "answer": answer})
    print(f"  {answer}")


def partial_rho(x, y, z) -> float:
    """Spearman between x and y with z held fixed."""
    rxy = stats.spearmanr(x, y).statistic
    rxz = stats.spearmanr(x, z).statistic
    ryz = stats.spearmanr(y, z).statistic
    return (rxy - rxz * ryz) / np.sqrt((1 - rxz**2) * (1 - ryz**2))


def boot_ci(x, y, n=10000):
    x, y = np.asarray(x), np.asarray(y)
    idx = np.arange(len(x))
    vals = [stats.spearmanr(x[s], y[s]).statistic
            for s in (rng.choice(idx, len(idx), replace=True) for _ in range(n))]
    vals = [v for v in vals if not np.isnan(v)]
    return np.percentile(vals, [2.5, 97.5])


def main() -> None:
    df = pd.read_csv(RESULTS / "benchmark_compounds.csv", index_col=0)

    fps = {}
    for cid, row in df.iterrows():
        m = Chem.MolFromSmiles(row.smiles)
        if m is not None:
            fps[cid] = AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)

    df["sim"] = pd.Series({c: DataStructs.TanimotoSimilarity(fps[REF], f)
                           for c, f in fps.items()})
    df["abs_err"] = (df.ddG_calc - df.ddG_exp).abs()
    df["abs_sel"] = df.ddG_exp.abs()
    others = df.drop(REF)

    print("\n1. HOW WIDE ARE THE CONFIDENCE INTERVALS?")
    for label, a, b in (("selectivity", df.ddG_calc, df.ddG_exp),
                        ("potency, E. coli", df.score_ecoli, df.dG_exp_ecoli),
                        ("potency, human", df.score_human, df.dG_exp_human)):
        rho = stats.spearmanr(a, b).statistic
        lo, hi = boot_ci(a, b)
        record("confidence interval",
               f"{label:18s} rho {rho:+.2f}  95% CI [{lo:+.2f}, {hi:+.2f}]"
               f"{'  spans zero' if lo < 0 < hi else ''}",
               f"{rho:+.2f}")

    print("\n2. IS THE POTENCY SIGNAL JUST MOLECULAR SIZE?")
    for e in ("ecoli", "human"):
        raw = stats.spearmanr(df[f"score_{e}"], df[f"dG_exp_{e}"]).statistic
        par = partial_rho(df[f"score_{e}"], df[f"dG_exp_{e}"], df.heavy_atoms)
        size = stats.spearmanr(df[f"score_{e}"], df.heavy_atoms).statistic
        record("size confound",
               f"{e:7s} score-vs-size {size:+.2f} | potency {raw:+.2f} "
               f"-> {par:+.2f} with size fixed",
               f"{par:+.2f}")

    print("\n3. DOES LIGAND EFFICIENCY RESCUE IT? (no: shared divisor)")
    le = stats.spearmanr(-df.score_ecoli / df.heavy_atoms,
                         -df.dG_exp_ecoli / df.heavy_atoms).statistic
    null = [stats.spearmanr(
                -rng.normal(df.score_ecoli.mean(), df.score_ecoli.std(), len(df)) / df.heavy_atoms,
                -rng.normal(df.dG_exp_ecoli.mean(), df.dG_exp_ecoli.std(), len(df)) / df.heavy_atoms,
            ).statistic for _ in range(5000)]
    record("ligand efficiency artefact",
           f"LE rho {le:+.2f}; same statistic on independent noise: "
           f"median {np.median(null):+.2f} (our value at the "
           f"{(np.array(null) < le).mean():.0%} percentile)",
           f"{le:+.2f}")

    print("\n4. HOW MANY INDEPENDENT COMPOUNDS ARE THERE REALLY?")
    ids, unassigned, clusters = list(fps), set(fps), []
    while unassigned:
        seed = unassigned.pop()
        group = {seed}
        for other in list(unassigned):
            if DataStructs.TanimotoSimilarity(fps[seed], fps[other]) > 0.7:
                group.add(other)
                unassigned.discard(other)
        clusters.append(group)
    identical = sum(1 for i, a in enumerate(ids) for b in ids[i + 1:]
                    if DataStructs.TanimotoSimilarity(fps[a], fps[b]) == 1.0)
    record("effective sample size",
           f"{len(ids)} compounds form {len(clusters)} clusters at Tanimoto 0.7; "
           f"{identical} pairs share an identical fingerprint",
           str(len(clusters)))

    print("\n5. SIMILARITY TO THE REFERENCE COMPOUND")
    o = others.sim
    record("similarity spread",
           f"Tanimoto to trimethoprim: median {o.median():.2f}, "
           f"range {o.min():.2f}-{o.max():.2f}, none above 0.7",
           f"{o.median():.2f}")

    print("\n6. DOES DOCKING DO BETTER ON COMPOUNDS LIKE THE REFERENCE?")
    near, far = others[others.sim >= SIM_CUT], others[others.sim < SIM_CUT]
    p_naive = stats.mannwhitneyu(near.abs_err, far.abs_err).pvalue
    record("applicability domain, naive",
           f"near {near.abs_err.mean():.2f} vs far {far.abs_err.mean():.2f} kcal/mol, "
           f"p = {p_naive:.3f}  (looks real)",
           f"p={p_naive:.3f}")

    # The confound: near-reference compounds are also the non-selective ones,
    # and a method that predicts near-zero selectivity gets those right for free.
    low = others[others.abs_sel < 1.0]
    ln, lf = low[low.sim >= SIM_CUT], low[low.sim < SIM_CUT]
    p_matched = stats.mannwhitneyu(ln.abs_err, lf.abs_err).pvalue
    record("applicability domain, matched",
           f"near |selectivity| {near.abs_sel.mean():.2f} vs far {far.abs_sel.mean():.2f} "
           f"kcal/mol, so the groups differ in difficulty, not just chemistry",
           "")
    record("applicability domain, matched",
           f"restricted to |selectivity| < 1: near {ln.abs_err.mean():.2f} vs "
           f"far {lf.abs_err.mean():.2f}, p = {p_matched:.3f}  (effect gone)",
           f"p={p_matched:.3f}")
    record("applicability domain, matched",
           f"similarity vs error: rho {stats.spearmanr(others.sim, others.abs_err).statistic:+.2f} "
           f"-> {partial_rho(others.sim, others.abs_err, others.abs_sel):+.2f} "
           f"with difficulty held fixed",
           "")

    print("\n7. DOES DOCKING BEAT A TRIVIAL NULL MODEL?")
    r_diff, _ = stats.spearmanr(others.abs_sel, others.abs_err)
    record("error is driven by difficulty",
           f"|true selectivity| vs |error|: rho = {r_diff:+.2f} "
           f"(docking predicts a narrow band, so error tracks distance from it)",
           f"{r_diff:+.2f}")
    mae_dock = others.abs_err.mean()
    mae_zero = others.ddG_exp.abs().mean()
    mae_mean = (others.ddG_exp - others.ddG_calc.mean()).abs().mean()
    record("null model comparison",
           f"docking {mae_dock:.2f} | always predict 0.00: {mae_zero:.2f} | "
           f"always predict {others.ddG_calc.mean():+.2f}: {mae_mean:.2f} kcal/mol",
           f"{mae_dock:.2f} vs {mae_zero:.2f}")
    verdict = ("predicting zero selectivity for everything is MORE accurate "
               "than running the docking") if mae_zero < mae_dock else \
              "docking beats the trivial null"
    record("null model verdict", verdict, "")

    out = pd.DataFrame(findings)
    out.to_csv(RESULTS / "robustness.csv", index=False)
    print(f"\nwrote {RESULTS / 'robustness.csv'}")


if __name__ == "__main__":
    main()
