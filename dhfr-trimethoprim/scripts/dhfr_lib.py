"""Shared helpers for the DHFR + trimethoprim pipeline.

Imported by the notebooks. Keeping the logic here rather than in cells means the
receptor recipe is literally identical across crystal, AlphaFold and human targets,
which is what makes the cross-system comparison meaningful.
"""

from __future__ import annotations

import json
import subprocess
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- paths

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PREP = ROOT / "data" / "prepared"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

for _d in (RAW, PREP, RESULTS, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)

# Browser UA: the EBI CDN 403s curl/urllib defaults and silently returns an
# HTML error page with whatever extension you asked for.
_UA = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
)

# Trimethoprim. Verified against PubChem CID 5578 - canonicalises identically.
TMP_SMILES = "COc1cc(cc(c1OC)OC)Cc2cnc(nc2N)N"
# Protonated at N1 of the diaminopyrimidine: the form Asp27 actually binds.
TMP_SMILES_PROTONATED = "COc1cc(cc(c1OC)OC)Cc2c[nH+]c(N)nc2N"
MTX_SMILES = (
    "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1"
)

# Docking box: 22 A per side, per the project plan.
BOX_SIZE = (22.0, 22.0, 22.0)


# ---------------------------------------------------------------- fetching


def fetch(url: str, dest: Path, *, min_bytes: int = 1000) -> Path:
    """Download unless present. Guards against the silent-403 trap."""
    if dest.exists() and dest.stat().st_size >= min_bytes:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read()
    if len(data) < min_bytes:
        raise RuntimeError(f"{url} returned only {len(data)} bytes - likely an error page")
    dest.write_bytes(data)
    return dest


def fetch_pdb(pdb_id: str) -> Path:
    pdb_id = pdb_id.upper()
    return fetch(f"https://files.rcsb.org/download/{pdb_id}.pdb", RAW / f"{pdb_id}.pdb")


def fetch_alphafold(acc: str, version: int = 6) -> dict[str, Path]:
    base = "https://alphafold.ebi.ac.uk/files"
    want = {
        "pdb": f"AF-{acc}-F1-model_v{version}.pdb",
        "cif": f"AF-{acc}-F1-model_v{version}.cif",
        "confidence": f"AF-{acc}-F1-confidence_v{version}.json",
        "pae": f"AF-{acc}-F1-predicted_aligned_error_v{version}.json",
    }
    return {k: fetch(f"{base}/{fn}", RAW / fn) for k, fn in want.items()}


# ---------------------------------------------------------------- PDB surgery

_WATER = {"HOH", "WAT", "DOD"}
# Crystallisation additives that should not be in a docking receptor.
ADDITIVES = {"BME", "EDO", "GOL", "SO4", "PO4", "ACT", "MPD", "PEG", "DMS", "CL", "NA", "MN"}
COFACTORS = {"NDP", "NAP", "NAD", "NAI"}


def split_pdb(
    pdb_path: Path,
    chain: str = "A",
    *,
    keep_cofactor: bool = True,
    drop_ligands: tuple[str, ...] = (),
) -> tuple[list[str], dict[str, list[str]]]:
    """Return (receptor_lines, {resname: ligand_lines}) for one chain.

    Waters and crystallisation additives are always dropped. The cofactor is kept
    in the receptor when keep_cofactor is True - trimethoprim binds alongside it,
    not instead of it, and removing it lets the ligand occupy the cofactor's space.
    """
    receptor: list[str] = []
    ligands: dict[str, list[str]] = {}

    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[21] != chain:
            continue
        resname = line[17:20].strip()

        if resname in _WATER or resname in ADDITIVES:
            continue
        if resname in drop_ligands:
            ligands.setdefault(resname, []).append(line)
            continue
        if resname in COFACTORS and not keep_cofactor:
            continue
        receptor.append(line)

    return receptor, ligands


def incomplete_residues(lines: list[str], min_atoms: int = 4) -> list[str]:
    """Residues with too few atoms to match a template (disordered termini)."""
    counts: dict[tuple[str, int], int] = {}
    names: dict[tuple[str, int], str] = {}
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        key = (line[21], int(line[22:26]))
        counts[key] = counts.get(key, 0) + 1
        names[key] = line[17:20].strip()
    return [f"{c}:{n}" for (c, n), k in counts.items() if k < min_atoms and names[(c, n)] != "GLY"]


def altloc_residues(lines: list[str]) -> list[str]:
    return sorted({f"{l[21]}:{int(l[22:26])}" for l in lines if l[16] != " "})


def centroid(lines: list[str]) -> np.ndarray:
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
    return xyz.mean(axis=0)


# ---------------------------------------------------------------- preparation


@dataclass
class Receptor:
    name: str
    pdbqt: Path
    source: Path
    cofactor: bool
    deleted: list[str]


def prepare_receptor(
    pdb_path: Path,
    name: str,
    *,
    chain: str = "A",
    keep_cofactor: bool = True,
    drop_ligands: tuple[str, ...] = (),
    default_altloc: str = "A",
    python_bin: str | None = None,
) -> tuple[Receptor, dict[str, list[str]]]:
    """One receptor recipe, used for every target."""
    rec_lines, ligands = split_pdb(
        pdb_path, chain, keep_cofactor=keep_cofactor, drop_ligands=drop_ligands
    )
    raw = PREP / f"{name}_receptor_raw.pdb"
    raw.write_text("\n".join(rec_lines) + "\nEND\n")

    bad = incomplete_residues(rec_lines)
    out_base = PREP / f"{name}_receptor"

    mk = Path(python_bin).parent / "mk_prepare_receptor.py" if python_bin else "mk_prepare_receptor.py"
    cmd = [str(mk), "--read_pdb", str(raw), "-o", str(out_base), "-p",
           "--default_altloc", default_altloc]
    for res in bad:
        cmd += ["-d", res]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    pdbqt = out_base.with_suffix(".pdbqt")
    if not pdbqt.exists():
        raise RuntimeError(f"receptor prep failed for {name}:\n{proc.stdout}\n{proc.stderr}")

    return Receptor(name, pdbqt, pdb_path, keep_cofactor, bad), ligands


def superpose_onto(
    mobile_pdb: Path,
    target_pdb: Path,
    out_pdb: Path,
    *,
    target_chain: str = "A",
) -> tuple[Path, float, int]:
    """Fit `mobile` onto `target` on matched CA atoms and write the moved coordinates.

    Required before an apo predicted model can share a docking box with a crystal
    structure: AlphaFold models come in their own coordinate frame, tens of angstroms
    away from the deposited crystal. Without this the box lands in empty solvent.
    """
    import warnings

    from Bio.PDB import PDBIO, PDBParser, Superimposer
    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    warnings.simplefilter("ignore", PDBConstructionWarning)
    parser = PDBParser(QUIET=True)

    def cas(path, chain_id=None):
        model = parser.get_structure("s", path)[0]
        chain = model[chain_id] if chain_id else next(iter(model))
        return {r.id[1]: r["CA"] for r in chain if r.id[0] == " " and "CA" in r}

    mobile_struct = parser.get_structure("mobile", mobile_pdb)
    mob = cas(mobile_pdb)
    tgt = cas(target_pdb, target_chain)
    common = sorted(set(mob) & set(tgt))
    if len(common) < 20:
        raise RuntimeError(f"only {len(common)} shared residues - numbering mismatch?")

    sup = Superimposer()
    sup.set_atoms([tgt[i] for i in common], [mob[i] for i in common])
    sup.apply(list(mobile_struct.get_atoms()))

    io = PDBIO()
    io.set_structure(mobile_struct)
    io.save(str(out_pdb))
    return out_pdb, float(sup.rms), len(common)


def prepare_ligand(smiles: str, name: str, *, seed: int = 42) -> Path:
    """SMILES -> 3D -> PDBQT. Never from crystal coordinates: building from SMILES
    is what keeps the redocking test honest."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad SMILES for {name}: {smiles}")

    # ChEMBL often stores salts (parent + counterion). Dock the parent only:
    # the counterion is not part of the binding species.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        raise RuntimeError(f"embedding failed for {name}")
    AllChem.MMFFOptimizeMolecule(mol)

    setup = MoleculePreparation().prepare(mol)[0]
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise RuntimeError(f"meeko failed for {name}: {err}")

    out = PREP / f"{name}.pdbqt"
    out.write_text(pdbqt)
    return out


def n_rotatable(pdbqt: Path) -> int:
    return sum(1 for l in pdbqt.read_text().splitlines() if l.startswith("BRANCH"))


# ---------------------------------------------------------------- docking


@dataclass
class DockResult:
    run: str
    receptor: str
    ligand: str
    seed: int
    exhaustiveness: int
    top_score: float
    scores: list[float]
    out_pdbqt: Path
    log: Path


def run_vina(
    receptor: Path,
    ligand: Path,
    center: np.ndarray,
    run_name: str,
    *,
    seed: int = 42,
    exhaustiveness: int = 32,
    num_modes: int = 9,
    size: tuple[float, float, float] = BOX_SIZE,
    vina_bin: str = "vina",
) -> DockResult:
    out = RESULTS / f"{run_name}.pdbqt"
    log = RESULTS / f"{run_name}.log"

    cmd = [
        vina_bin,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--out", str(out),
    ]
    # Vina 1.2 removed --log; capture stdout instead.
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log.write_text(proc.stdout + proc.stderr)
    if not out.exists():
        raise RuntimeError(f"vina failed for {run_name}:\n{proc.stdout}\n{proc.stderr}")

    scores = parse_scores(log)
    return DockResult(
        run_name, receptor.stem, ligand.stem, seed, exhaustiveness,
        scores[0] if scores else float("nan"), scores, out, log,
    )


def parse_scores(log: Path) -> list[float]:
    scores, in_table = [], False
    for line in log.read_text().splitlines():
        if line.strip().startswith("-----+"):
            in_table = True
            continue
        if in_table:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                scores.append(float(parts[1]))
            elif scores:
                break
    return scores


# ---------------------------------------------------------------- RMSD


def pose_rmsds(dock_pdbqt: Path, reference_lines: list[str], smiles: str) -> list[float]:
    """Symmetry-aware heavy-atom RMSD of each docked pose against a crystal pose.

    Uses RDKit CalcRMS against a template-corrected reference, NOT Vina's own
    rmsd columns, which only compare poses with each other.
    """
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign

    ref_path = PREP / "_ref_tmp.pdb"
    ref_path.write_text("\n".join(reference_lines) + "\nEND\n")

    template = neutralise(Chem.MolFromSmiles(smiles))
    ref = Chem.MolFromPDBFile(str(ref_path), removeHs=True, sanitize=False)
    ref = AllChem.AssignBondOrdersFromTemplate(template, ref)

    pm = PDBQTMolecule.from_file(str(dock_pdbqt), skip_typing=True)
    probe = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
    # A protonated ligand and a neutral crystal reference share an identical
    # heavy-atom skeleton but will not substructure-match while the charge differs.
    # RMSD here is heavy-atom only, so neutralise before matching.
    probe = neutralise(probe)

    return [rdMolAlign.CalcRMS(probe, ref, prbId=i) for i in range(probe.GetNumConformers())]


def neutralise(mol):
    """Strip formal charges and explicit H counts, keeping the heavy-atom skeleton."""
    from rdkit import Chem

    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetFormalCharge() != 0:
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


# ---------------------------------------------------------------- misc


def toolchain() -> dict[str, str]:
    import Bio
    import meeko
    import rdkit

    out = {"rdkit": rdkit.__version__, "biopython": Bio.__version__, "meeko": meeko.__version__}
    try:
        v = subprocess.run(["vina", "--version"], capture_output=True, text=True).stdout
        out["vina"] = v.strip().splitlines()[0] if v.strip() else "?"
    except FileNotFoundError:
        out["vina"] = "NOT FOUND"
    return out


def save_json(obj, path: Path) -> Path:
    path.write_text(json.dumps(obj, indent=2, default=str))
    return path
