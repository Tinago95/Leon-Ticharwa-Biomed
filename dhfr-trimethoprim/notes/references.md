# Verified references

Every citation checked against Crossref, PubMed and the RCSB primary-citation records
on 2026-11-27. Volumes, pages and DOIs below were read from those records, not copied
from the project plan.

All five citations in the original plan are **real and substantially correct**. Two
details needed fixing, and one mapping was missing.

---

## Methods

**AutoDock Vina 1.2**, cited correctly by the plan.

> Eberhardt J, Santos-Martins D, Tillack AF, Forli S. *AutoDock Vina 1.2.0: New Docking
> Methods, Expanded Force Field, and Python Bindings.* Journal of Chemical Information
> and Modeling **61**(8):3891–3898 (2021). doi:10.1021/acs.jcim.1c00203

**AlphaFold2**, cited correctly by the plan.

> Jumper J, Evans R, Pritzel A, *et al.* *Highly accurate protein structure prediction
> with AlphaFold.* Nature **596**:583–589 (2021). doi:10.1038/s41586-021-03819-2

---

## Structures

**Matthews 1985**, real, but the plan's description is incomplete.

> Matthews DA, Bolin JT, Burridge JM, *et al.* *Refined crystal structures of
> Escherichia coli and chicken liver dihydrofolate reductase containing bound
> trimethoprim.* Journal of Biological Chemistry **260**(1):381–391 (1985).
> doi:10.1016/s0021-9258(18)89743-5

The plan describes this as covering *E. coli* DHFR. The paper covers **both *E. coli*
and chicken liver** DHFR, and the cross-species comparison is much of its point, which
makes it more relevant to this project's selectivity question than the plan suggests.

**Sawaya & Kraut 1997**, real; the plan gives no journal, so here it is.

> Sawaya MR, Kraut J. *Loop and subdomain movements in the mechanism of Escherichia coli
> dihydrofolate reductase: crystallographic evidence.* Biochemistry **36**(3):586–603
> (1997). PMID 9012674

This is the primary citation for **1RX2**, as the plan implies. It is also the paper
behind the Met20 loop story, so it is the right reference for any claim about loop
movement.

**Heaslet 2009**, cited correctly by the plan.

> Heaslet H, Harris M, Fahnoe K, *et al.* *Structural comparison of chromosomal and
> exogenous dihydrofolate reductase from Staphylococcus aureus in complex with the
> potent inhibitor trimethoprim.* Proteins **76**(3):706–717 (2009). PMID 19280600

Primary citation for **2W9G**, the redocking-gate structure.

---

## Missing from the plan

**7DFR is not the Matthews paper.** The plan lists `7DFR` among its structures and
Matthews 1985 among its references, which invites the assumption that one cites the
other. It does not:

> Bystroff C, Oatley SJ, Kraut J. *Crystal structures of Escherichia coli dihydrofolate
> reductase: the NADP+ holoenzyme and the folate·NADP+ ternary complex.* Biochemistry
> **29**(13):3263–3277 (1990). PMID 2185835

**2W3A is unpublished.** RCSB lists its citation as "To be Published", by Leung AKW,
Ross LJ, Zywno-Van Ginkel S *et al.* If you cite 2W3A, cite the **PDB entry**, not a
paper.

**1U72 has its own citation**, which the plan does not give:

> Cody V, Luft JR, Pangborn W. *Understanding the role of Leu22 variants in methotrexate
> resistance: comparison of wild-type and Leu22Arg variant mouse and human dihydrofolate
> reductase ternary crystal complexes with methotrexate and NADPH.* Acta
> Crystallographica Section D **61**(2):147–155 (2005). PMID 15681865

---

## Citation map

Which paper goes with which structure, so nothing is misattributed in the write-up:

| PDB | Primary citation | Used here for |
|---|---|---|
| `1RX2` | Sawaya & Kraut 1997, Biochemistry | AlphaFold comparison, *E. coli* docking |
| `7DFR` | Bystroff, Oatley & Kraut 1990, Biochemistry | mentioned, not used |
| `2W9G` | Heaslet *et al.* 2009, Proteins | the redocking gate |
| `2W9H` | Heaslet *et al.* 2009, Proteins | not used, no cofactor |
| `1U72` | Cody, Luft & Pangborn 2005, Acta Cryst D | human receptor |
| `2W3A` | unpublished, cite the PDB entry | not used |

---

## Still worth adding

The write-up makes claims about method accuracy that currently rest on general knowledge
rather than a citation. Before publishing, find a reference for:

- Typical RMSE of relative FEP on congeneric series (~1–1.5 kcal/mol)
- Typical accuracy of docking for absolute affinity (~2–3 kcal/mol, weak rank order)
- Trimethoprim's selectivity for bacterial over human DHFR (~10³–10⁵ fold)

The third is partly covered by this project's own ChEMBL data, which gives 904-fold for
trimethoprim from measured Kᵢ values, so that number can be reported as derived here
rather than cited.
