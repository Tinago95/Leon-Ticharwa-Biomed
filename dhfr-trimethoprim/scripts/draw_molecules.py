"""Render the molecules in this project as PNGs for embedding in GUIDE.md.

GitHub renders PNG reliably in markdown, so PNG rather than SVG.
Run:  python scripts/draw_molecules.py
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "molecules"
OUT.mkdir(parents=True, exist_ok=True)

MOLECULES = {
    "trimethoprim": (
        "COc1cc(cc(c1OC)OC)Cc2cnc(nc2N)N",
        "Trimethoprim - antibacterial, targets bacterial DHFR",
    ),
    "trimethoprim_protonated": (
        "COc1cc(cc(c1OC)OC)Cc2c[nH+]c(N)nc2N",
        "Trimethoprim, protonated at N1 - the form Asp27 binds (34% at pH 7.4)",
    ),
    "methotrexate": (
        "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1",
        "Methotrexate - anticancer, targets human DHFR",
    ),
    "dihydrofolate": (
        "Nc1nc2NCC(CNc3ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc3)=Nc2c(=O)[nH]1",
        "Dihydrofolate (DHF) - the substrate DHFR acts on",
    ),
    "tetrahydrofolate": (
        "C1C(NC2=C(N1)N=C(NC2=O)N)CNC3=CC=C(C=C3)C(=O)N[C@@H](CCC(=O)O)C(=O)O",
        "Tetrahydrofolate (THF) - the product, the active one-carbon carrier",
    ),
    "folic_acid": (
        "Nc1nc2ncc(CNc3ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc3)nc2c(=O)[nH]1",
        "Folic acid - the vitamin, seen in the 1RX2 crystal as FOL",
    ),
    "nadph": (
        "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]2O[C@@H]"
        "(n3cnc4c(N)ncnc43)[C@H](OP(=O)(O)O)[C@@H]2O)[C@@H](O)[C@H]1O",
        "NADPH - the cofactor, reduced form (PDB code NDP). Carries the hydride.",
    ),
}


def draw(name: str, smiles: str, size=(560, 340)) -> Path | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  SKIP {name}: SMILES would not parse")
        return None
    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    opts.bondLineWidth = 2
    opts.additionalAtomLabelPadding = 0.08
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()

    path = OUT / f"{name}.png"
    path.write_bytes(drawer.GetDrawingText())
    return path



def main() -> None:
    print(f"writing to {OUT}")
    for name, (smiles, caption) in MOLECULES.items():
        path = draw(name, smiles)
        if path:
            mol = Chem.MolFromSmiles(smiles)
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            print(f"  {path.name:34s} {formula:16s} {mol.GetNumAtoms():3d} heavy")

    # Side-by-side comparison of the two drugs.
    pair = [Chem.MolFromSmiles(MOLECULES[k][0]) for k in ("trimethoprim", "methotrexate")]
    for m in pair:
        AllChem.Compute2DCoords(m)
    img = Draw.MolsToGridImage(
        pair, molsPerRow=2, subImgSize=(420, 300),
        legends=["Trimethoprim -> bacterial DHFR", "Methotrexate -> human DHFR"],
        returnPNG=False,
    )
    img.save(OUT / "drug_comparison.png")
    print(f"  {'drug_comparison.png':34s} two drugs, one enzyme family")


if __name__ == "__main__":
    main()
