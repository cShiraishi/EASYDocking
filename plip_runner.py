"""
plip_runner.py — Helper script to run PLIP analysis and output JSON results.
This script is called from app.py via subprocess using the vina-rdkit conda environment
which has PLIP + openbabel installed (compatible with Python 3.11).

Usage:
    /path/to/vina-rdkit/python plip_runner.py <complex.pdb>

Output: JSON to stdout with interaction data.
"""
import sys
import json

def run_plip(complex_pdb_path):
    try:
        from plip.structure.preparation import PDBComplex
    except ImportError:
        print(json.dumps({"error": "PLIP not installed in this environment"}))
        sys.exit(1)

    result = {
        "bsid": None,
        "interactions": {
            "hbonds":        [],  # {residue, dist, angle, donor}
            "hydrophobic":   [],  # {residue, dist}
            "pistacking":    [],  # {residue, dist, type}
            "pication":      [],  # {residue, dist}
            "saltbridges":   [],  # {residue, dist}
            "halogenbonds":  [],  # {residue, dist}
            "metal":         [],  # {residue, dist}
            "waterbridges":  [],  # {residue, dist}
        },
        "error": None
    }

    try:
        mol = PDBComplex()
        mol.load_pdb(complex_pdb_path)
        mol.analyze()

        # Find LIG binding site
        bsid = None
        for key in mol.interaction_sets.keys():
            if 'LIG' in key:
                bsid = key
                break

        if bsid is None and mol.interaction_sets:
            bsid = list(mol.interaction_sets.keys())[0]

        if bsid is None:
            result["error"] = "No binding site found"
            print(json.dumps(result))
            return

        result["bsid"] = bsid
        interactions = mol.interaction_sets[bsid]

        # Hydrogen bonds
        for hb in interactions.hbonds_pdon + interactions.hbonds_ldon:
            result["interactions"]["hbonds"].append({
                "residue": f"{hb.restype}{hb.resnr}{hb.reschain}",
                "dist": round(float(hb.distance_ad), 2),
                "angle": round(float(hb.angle), 1),
                "donor": hb.protisdon
            })

        # Hydrophobic contacts
        for hc in interactions.hydrophobic_contacts:
            result["interactions"]["hydrophobic"].append({
                "residue": f"{hc.restype}{hc.resnr}{hc.reschain}",
                "dist": round(float(hc.distance), 2)
            })

        # Pi-stacking
        for ps in interactions.pistacking:
            result["interactions"]["pistacking"].append({
                "residue": f"{ps.restype}{ps.resnr}{ps.reschain}",
                "dist": round(float(ps.distance), 2),
                "type": str(ps.type)
            })

        # Pi-cation
        for pc in interactions.pication_laro + interactions.pication_paro:
            result["interactions"]["pication"].append({
                "residue": f"{pc.restype}{pc.resnr}{pc.reschain}",
                "dist": round(float(pc.distance), 2)
            })

        # Salt bridges
        for sb in interactions.saltbridges_lneg + interactions.saltbridges_pneg:
            result["interactions"]["saltbridges"].append({
                "residue": f"{sb.restype}{sb.resnr}{sb.reschain}",
                "dist": round(float(sb.distance), 2)
            })

        # Halogen bonds
        for xb in interactions.halogen_bonds:
            result["interactions"]["halogenbonds"].append({
                "residue": f"{xb.restype}{xb.resnr}{xb.reschain}",
                "dist": round(float(xb.distance), 2)
            })

        # Metal complexes
        for mb in interactions.metal_complexes:
            result["interactions"]["metal"].append({
                "residue": f"{mb.restype}{mb.resnr}{mb.reschain}",
                "dist": round(float(mb.distance), 2)
            })

        # Water bridges
        for wb in interactions.water_bridges:
            result["interactions"]["waterbridges"].append({
                "residue": f"{wb.restype}{wb.resnr}{wb.reschain}",
                "dist": round(float(wb.distance_aw), 2)
            })

    except Exception as e:
        result["error"] = str(e)

    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: plip_runner.py <complex.pdb>"}))
        sys.exit(1)
    run_plip(sys.argv[1])
