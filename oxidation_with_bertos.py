import os
import subprocess
from pymatgen.core.structure import Structure
from pymatgen.analysis.bond_valence import BVAnalyzer

# === CONFIG ===
#need to modify :D
BERTOS_DIR = r"C:\Users\95224\Documents\GitHub\BERTOS"
BERTOS_MODEL = os.path.join(BERTOS_DIR, "trained_models", "ICSD_CN")
GET_OS_SCRIPT = os.path.join(BERTOS_DIR, "getOS.py")
CIF_PATH = r"C:\Users\95224\Neutron Project\Datasets\Cif\mp-186.cif"


def bertos_predict(cif_path):
    try:
        structure = Structure.from_file(cif_path)
        formula = structure.composition.reduced_formula

        result = subprocess.run(
            ["python", GET_OS_SCRIPT, "--i", formula],
            capture_output=True,
            text=True,
            check=True
        )

        print("\n[BERTOS Output]")
        print(result.stdout)

        # Parse predicted oxidation states
        os_dict = {}
        lines = result.stdout.strip().splitlines()
        parsing = False
        for line in lines:
            if "Predicted Oxidation States" in line:
                parsing = True
                continue
            if parsing and ":" in line:
                elem, ox = line.strip().split(":")
                os_dict[elem.strip()] = int(ox.strip())

        ox_list = [os_dict.get(site.specie.symbol, None) for site in structure]
        return ox_list  # this becomes oxidation_set

    except subprocess.CalledProcessError as e:
        print("BERTOS failed:", e.stderr)
        return [None] * len(Structure.from_file(cif_path))  # fallback


# === Wrapper that tries BVAnalyzer first ===
def get_valences_with_fallback(structure, cif_path=None):
    try:
        bv = BVAnalyzer()
        decorated = bv.get_oxi_state_decorated_structure(structure)
        oxidation_states = [site.oxi_state for site in decorated]
        return oxidation_states
    except Exception as e:
        print("BVAnalyzer failed:", e)
        if cif_path is None:
            raise ValueError("Must provide cif_path for BERTOS fallback")
        return bertos_predict(cif_path)  # returns oxidation states list




if __name__ == "__main__":
    import sys
    cif_path = sys.argv[1] if len(sys.argv) > 1 else CIF_PATH

    structure = Structure.from_file(cif_path)
    print(f"[Structure loaded] Composition: {structure.composition}")

    valences = get_valences_with_fallback(structure, cif_path)

    print("\n[Final Valence Result]")
    print(valences)