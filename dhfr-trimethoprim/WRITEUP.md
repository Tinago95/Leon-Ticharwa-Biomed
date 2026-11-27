# I tried to predict why an antibiotic is safe. It didn't work, and that's the result.

**A beginner's computational drug discovery project on DHFR and trimethoprim.**

*Draft. Numbers are live from the repository; prose is mine to rewrite in my own voice
before posting.*

---

## Why I did this

I'm preparing for a Bachelor of Health and Medical Science and, beyond it, graduate-entry
medicine, with the aim of working somewhere between medicine and drug discovery. I have a
programming background and limited formal chemistry, so I wanted to start building the
structural-biology half before the degree does.

I picked a problem where the answer is already known and published, so that at every stage
I could check my work against reality instead of against my own expectations.

---

## What I built

An end-to-end docking pipeline for one bacterial enzyme and one real antibiotic:
dihydrofolate reductase (DHFR) from *Escherichia coli*, and trimethoprim.

Sequence from UniProt, predicted structure from AlphaFold, comparison against a crystal
structure, ligand preparation from SMILES, molecular docking with AutoDock Vina 1.2.7,
interaction analysis, and a benchmark against experimental binding data from ChEMBL.

Six notebooks, each ending in a pass check that has to hold before the next one runs.
Everything reproduces from a clean checkout.

**Repository:** `github.com/Tinago95/Leon-Ticharwa-Biomed`

DHFR recycles dihydrofolate back into tetrahydrofolate, the carrier cells use to build
DNA. Block it and a dividing bacterium stops dividing. Trimethoprim blocks the bacterial
enzyme roughly a thousand times harder than the human one, which is the whole reason it
is a usable antibiotic rather than a poison.

That thousand-fold difference is what I tried to reproduce computationally.

---

## How I know it works

Before trusting any prediction, I validated the protocol against a known experimental
answer.

I took `2W9G`, a crystal structure of *S. aureus* DHFR with trimethoprim bound, deleted
the drug, rebuilt it from scratch from its SMILES string, and docked it back in. If the
protocol works, the molecule should land where crystallography found it.

**It landed 0.44 Å from the crystal pose, reproduced across three random seeds.**

The conventional threshold is 2 Å. Rebuilding from SMILES rather than reusing the crystal
coordinates matters: otherwise you hand the docking program part of the answer and the
test proves nothing.

I also checked the prediction against experiment upstream of that. The AlphaFold model of
*E. coli* DHFR agreed with crystal structure `1RX2` to **0.361 Å** backbone RMSD. The
prediction was never the weak link.

---

## What I found

With a validated protocol, I docked trimethoprim into bacterial and human DHFR:

| System | Vina score (kcal/mol) |
|---|---|
| *E. coli* crystal (1RX2) | −7.38 |
| *E. coli* AlphaFold model | −7.70 |
| Human (1U72) | **−7.99** |
| Methotrexate → human *(control)* | −10.11 |

Human DHFR scored **better**. The difference, 0.61 kcal/mol, points the wrong way.

Since every factor of ten in binding corresponds to about 1.42 kcal/mol, trimethoprim's
real ~900-fold selectivity needs roughly **4.2 kcal/mol**. I measured 0.61, backwards,
against a seed-to-seed spread of only 0.06. The search was precise and the answer was
wrong, which turns out to be the most useful thing I learned all project: **precision is
not accuracy**.

### One compound is an anecdote

Rather than stop there, I pulled every compound ChEMBL has measured against both enzymes:
25 of them, spanning seven orders of magnitude of selectivity, with trimethoprim itself
among them at a measured 904-fold.

I used only Kᵢ, never IC₅₀. By Cheng–Prusoff, an IC₅₀ depends on the substrate
concentration the assay used, and ChEMBL rarely records it, so pooling IC₅₀ across labs
would compare quantities that aren't the same thing.

| Test | Spearman ρ |
|---|---|
| Ranking potency within *E. coli* | **+0.54** |
| Ranking potency within human | +0.11 |
| **Predicting selectivity between them** | **−0.21** (p = 0.32) |

Docking got the direction of selectivity right for **1 compound out of 14**. A coin flip
would manage about 7. Its predictions spanned 23% of the experimental range and never
crossed zero, meaning the method cannot express a bacterial-selective compound at all.

The trimethoprim result wasn't bad luck.

### Where exactly it breaks

A failed prediction is only interesting if you can say why. So I measured what the docked
poses actually touched:

| Target | Anchor residue | Distance |
|---|---|---|
| *E. coli* | Asp27 | 3.09 Å |
| Human | Glu30 | 2.86 Å |

Both poses contact the correct anchoring residue. **The geometry is right; only the
ranking is wrong.** The failure is in the scoring function, not in conformational
sampling.

That fits a known design decision: AutoDock Vina deliberately ignores partial charges and
handles electrostatics only implicitly. Trimethoprim's selectivity is largely
electrostatic, resting on how well each pocket stabilises a protonated ring against a
buried carboxylate. The one piece of physics that matters most here is the piece the
method doesn't model.

---

## What went wrong along the way

**The project plan I started from had six factual errors.** I only found them because I
wrote a notebook that re-derives every identifier from source instead of trusting the
document. The worst: it gave `P00382` as human DHFR. That accession is actually an
*E. coli* trimethoprim-**resistance** enzyme. Left uncaught, my "bacterial versus human"
comparison would have been *E. coli* versus drug-resistant *E. coli*. It would have run
fine and produced plausible numbers. The correct accession is `P00374`.

**I put a docking box in empty space.** AlphaFold models arrive in their own coordinate
frame. Mine sat 51 Å from the crystal structure, so reusing the crystal's box centre
searched solvent. I caught it only because I checked how many protein atoms fell inside
the box. There is now an assertion that fails loudly if it ever recurs.

**macOS silently deleted my docking software.** The prebuilt AutoDock Vina 1.1.2 binary is
unsigned and from 2019; it vanished from `~/Downloads` between my inspecting it and
running it, not even to the Trash. I built 1.2.7 from source instead, which was the right
version anyway.

**I made an error in my own documentation** and had to correct it: I recorded trimethoprim
as having 39 heavy atoms. It has 21 heavy atoms and 39 total including hydrogens. The
notebook now asserts both so it can't drift again.

---

## What this does not show

This is not drug discovery. No new compound, no patient-relevant claim.

A docking score is not a measured binding affinity, and I have not reported one as if it
were. The benchmark covers 25 compounds from heterogeneous assays against one target pair:
small enough that every correlation here has a wide confidence interval, and not a general
statement about docking as a technique. It worked for pose prediction and partially for
potency ranking; it failed at one specific task that is known to be hard.

I also haven't tested the obvious next steps. MM-GBSA rescoring would put explicit charges
and desolvation back in and might close the gap. Relative free-energy perturbation through
a thermodynamic cycle is the rigorous answer and is well beyond what I can currently run or
validate.

---

## What I'm learning next

Acetylcholinesterase and donepezil, the Alzheimer's target. Everything transfers: the repo
layout, the preparation recipe, the redocking-first discipline. What gets harder is the
point of doing it second. The target is human, the binding site is a deep narrow gorge
rather than an open cleft, donepezil is larger and more flexible, and the limits of
docking as a model of a brain drug need even more careful wording.

Before that, the gaps I want to close: the physical chemistry behind scoring functions,
enough statistics to put proper confidence intervals on a correlation of 25 points, and
the thermodynamics to understand what a free-energy perturbation is actually computing.

---

## The honest summary

I built a pipeline, validated it against a known answer to 0.44 Å, and then used it to
ask a question it turned out to be incapable of answering. I know it's incapable because I
checked it against 25 experimental measurements rather than stopping at one convenient
result, and I know *where* it fails, because the poses were right while the ranking was
wrong.

That's a negative result. It is also the only kind of result I'd have trusted, because the
alternative was a number I had no way to check.

---

*Code, data and full methods: `github.com/Tinago95/Leon-Ticharwa-Biomed`. Every figure
regenerates from a script; every number above traces to a CSV in `results/`. A guide to
the background concepts, written for readers new to the field, is in `GUIDE.md`.*
