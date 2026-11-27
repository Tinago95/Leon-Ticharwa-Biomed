# Verified identifiers

Canonical reference for this project. **Use this table, not the one on page 6 of the
project-plan PDF**: that table contains four errors, listed at the bottom.

Every row was re-derived from the files in `data/raw/` or from a live RCSB API query by
[`notebooks/00_verify_sources.ipynb`](../notebooks/00_verify_sources.ipynb)
(20 checks, 5 failures). Raw output: [`results/source_verification.csv`](../results/source_verification.csv).

Last verified: 2026-11-27.

---

## Proteins

| System | UniProt | Length | Organism | Gene | Notes |
|---|---|---|---|---|---|
| *E. coli* DHFR, **main target** | `P0ABQ4` | 159 aa | *Escherichia coli* K12 | `folA` | 18.0 kDa (computed). Asp27 and Met20 both confirmed present at those positions. |
| Human DHFR, **selectivity comparison** | `P00374` | 187 aa | *Homo sapiens* | `DHFR` | **Corrected**: the PDF gives `P00382`, which is a different protein entirely. |
| *S. aureus* DHFR, docking validation | not specified | 157 aa | *Staphylococcus aureus* |, | Work from the 2W9G coordinates; no UniProt accession verified for this project. |

### Do not use `P00382`

| Accession | What the PDF says | What it actually is |
|---|---|---|
| `P00382` | Human DHFR, 187 aa | **Dihydrofolate reductase type 1**, *E. coli*, 157 aa, gene `dhfrI` |

`P00382` is a plasmid-borne **trimethoprim-resistance** enzyme. Using it as the "human"
comparator would silently compare *E. coli* against drug-resistant *E. coli* and produce a
selectivity result that looks plausible and means nothing. This is the one error in the PDF
that would invalidate a headline result rather than just a detail.

`P00382` is not stored in `data/raw/`: it is genuinely interesting for the resistance
stretch task on page 10, just not as the human target.

---

## Structures

Ligand codes are the actual chemical components present in each entry, read from RCSB.

| PDB | Organism | Ligands present | Cofactor state | Use |
|---|---|---|---|---|
| `1RX2` | *E. coli* | FOL, NAP, **MN, BME** | NADP+ (oxidised) | Phase 2 superposition target |
| `7DFR` | *E. coli* | FOL, NAP | NADP+ (oxidised) | Alternative *E. coli* holoenzyme |
| `1U72` | *Homo sapiens* | MTX, **NDP** | **NADPH (reduced)** | Human receptor, methotrexate control |
| `2W3A` | *Homo sapiens* | TOP, NDP, GOL, SO4 | NADPH (reduced) | Human DHFR + trimethoprim |
| `2W9G` | *S. aureus* | TOP, NDP | NADPH (reduced) | **Redocking gate (Phase 4b)**: has the cofactor |
| `2W9H` | *S. aureus* | TOP, **EDO** | **none** | Apo-cofactor; do **not** use for the gate |

**Ligand codes:** `FOL` folic acid · `NAP` NADP+ · `NDP` NADPH · `TOP` trimethoprim ·
`MTX` methotrexate · `MN` manganese(II) · `BME` β-mercaptoethanol ·
`EDO` 1,2-ethanediol · `GOL` glycerol · `SO4` sulfate

### Receptor-prep consequences

- **`1RX2` carries Mn²⁺ and β-mercaptoethanol** beyond the folate and NADP+ the PDF
  mentions. Phase 4a says "remove waters and the bound ligand, except keep NADPH/NADP+".
  That instruction does not cover MN and BME. Decide explicitly and record it; both are
  crystallisation additives and are normally stripped.
- **`2W9H` has no cofactor.** Phase 4b's gate assumes trimethoprim *and* NADPH. Only
  `2W9G` satisfies that. Running the gate on 2W9H would prepare a receptor missing the
  cofactor the whole recipe depends on.
- **`1U72` is NADPH, not NADP+.** Oxidised and reduced nicotinamide differ in charge and
  geometry at the nicotinamide ring. If the receptor recipe is meant to be identical across
  systems, note that the *E. coli* entries supply NADP+ while the human entry supplies NADPH.

### Numbering

`1RX2` chain A is numbered **1–159 and matches UniProt `P0ABQ4` numbering directly**:
residue 27 is Asp, residue 20 is Met. No offset correction is needed for the Phase 5
Asp27 contact test.

One sequence difference: **position 37 is Asp in `1RX2`, Asn in `P0ABQ4`.** Expect one
mismatch when aligning the AlphaFold model to the crystal; it is not a parsing bug.

---

## Ligand

| Property | Value |
|---|---|
| Compound | Trimethoprim |
| PubChem CID | `5578` |
| Formula | C14H18N4O3 |
| Atoms | 39 total, **21 heavy** + 18 hydrogens |
| SMILES | `COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC` |

The SMILES string hardcoded in the PDF (`COc1cc(cc(c1OC)OC)Cc2cnc(nc2N)N`) is
**correct**: RDKit canonicalises it and the PubChem string to the identical form. Safe
to use as written.

---

## AlphaFold

Current release is **v6**; v4 URLs return 404. File URL pattern:

```
https://alphafold.ebi.ac.uk/files/AF-{ACCESSION}-F1-model_v6.pdb
https://alphafold.ebi.ac.uk/files/AF-{ACCESSION}-F1-model_v6.cif
https://alphafold.ebi.ac.uk/files/AF-{ACCESSION}-F1-confidence_v6.json
https://alphafold.ebi.ac.uk/files/AF-{ACCESSION}-F1-predicted_aligned_error_v6.json
```

The EBI CDN returns **403 to a default curl User-Agent**. Send a browser User-Agent or the
download silently writes a 134-byte HTML error page with a `.pdb` extension. Check file
sizes after fetching.

PDB and mmCIF residue counts agree for all three models (159 / 187 / 157).

---

## Sequence identity

The PDF states ~26% identity between human and *E. coli* DHFR, and "roughly 74% of
sequence differs".

A BLOSUM62 global alignment (gap open −11, extend −1) over `P0ABQ4` and `P00374` gives
**32.3%** (51 identities across 158 aligned columns).

This figure is sensitive to alignment method and parameters, and published values span
roughly 26–35%, so 26% is not wrong so much as imprecise. **Quote your own computed number
and state the method** rather than citing 26% as if it were measured.

---

## Summary of PDF errors

| # | Page | Claim | Reality | Severity |
|---|---|---|---|---|
| 1 | 6 | `P00382` is human DHFR, 187 aa | *E. coli* DHFR type 1, 157 aa | **Invalidates the selectivity result** |
| 2 | 6 | `1U72` = methotrexate + NADP+ | Contains NADPH (NDP) | Affects receptor prep |
| 3 | 6 | `2W9G / 2W9H` = trimethoprim + NADPH | 2W9H has no cofactor | Breaks the redocking gate if 2W9H chosen |
| 4 | 7 | `1RX2` = folate + NADP+ | Also Mn²⁺ and β-mercaptoethanol | Undefined prep step |
| 5 | 9 | Phase 4b command uses `vina --log` | `--log` was removed in Vina 1.2 | **Command fails as written** |

Error 5: the plan cites Vina 1.2 (Eberhardt 2021) but writes 1.1.2 syntax. Redirect
stdout instead. See [`decisions.md`](decisions.md).

### Verified correct

Gene `folA`; 159 residues and 18.0 kDa (computed) for *E. coli* DHFR; Asp27 as the
anchoring acidic residue; Met20 at position 20; the trimethoprim SMILES; all six PDB
organisms; `2W3A` as human DHFR with trimethoprim and NADPH; the DHF + NADPH + H⁺ → THF +
NADP+ reaction; the co-trimoxazole sequential-blockade pharmacology; and the methotrexate
contrast.

The **"central eight-stranded beta sheet flanked by four alpha helices"** description is
confirmed from `1RX2`'s own header annotation: one sheet (`A`) declaring 8 strands, and
4 HELIX records.

### Not checked

The journal citations on page 13 (Matthews 1985, Sawaya & Kraut 1997, Heaslet 2009,
Eberhardt 2021, Jumper 2021) were **not verified**: authors, years and titles still need
confirming against the literature before they go in a write-up.