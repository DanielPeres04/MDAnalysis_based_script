import MDAnalysis as mda
from MDAnalysis.analysis import rms
import sys
import pandas as pd
from matplotlib import pyplot as plt

START_TIME = 0
STOP_TIME = 50

# the protein key is used to select residues in all monomers
# the protein is parsed into different monomers automatically
PROTEIN_KEY = "resid 1-50"

# if more ligands are needed, add them using the reffered notation
# specific atom IDs can selected for the ligands 
LIGANDS_KEYS = ["resname SAM", "resname NAS"]

# choose functions from function selection pool. Functions are applied from left to right
# function options: 
FUNCTIONS = [compute_RMSD_fit_protein, ]


def analyse_rgyr(u, atom_list, step: int = 100):
    """
    Computes the gyration radius for each corresponding frame and simulation time. The outputted values are dependent on the inputed step
    """
    u.trajectory[0].frame
    reduced_step = range(0,len(u.trajectory)+1,100)
    rgyr_data = []
    for st in reduced_step: 
        rgyr_data.append((u.trajectory.time, atom_list.radius_of_gyration()))
        #print(f"Frame = {frame:3d}; Time = {time:4.0f} ps; Rgyr = {rgyr:.4f} A")
    u.trajectory[0]
    return rgyr_data


def parse_protein(atom_list):
    """
    Creates a dict divided by complexes and adds a different monomer to each. For instance, a dimer dict \
    will have 2 complexes and a tetramer dict will have 4 complexes.
    
    :param atom_list: An AtomGroup object containing atoms of a protein.
    """
    complexes = {}
    current_complex = 1
    current_resid_number = -1
    current_index = -1
    for i in range(len(atom_list)):
        atom_info = atom_list[i]
        current_index = int(atom_info.index)
        resid_number = atom_info.resid
        #update the chain number if a new chain has been found
        if resid_number < current_resid_number:
            current_complex += 1
        current_resid_number = resid_number
        #add a chain to the dict if it does not exist
        if f"{current_complex}" not in complexes:
            complexes[f"{current_complex}"] = {}
        if "protein_chain" not in complexes[f"{current_complex}"]:
            complexes[f"{current_complex}"] ["protein_chain"] = [current_index]
        else:
            complexes[f"{current_complex}"] ["protein_chain"].append(current_index)

    return complexes

def add_ligands_to_complexes(u, complexes_dict, ligand_key_list):
    """
    Adds the ligands of each chain complex in the system to a pre-existing complexes-dict
    
    :param u: universe
    :param complexes_dict: a complexes dictionaty that contains protein chains for each chain complex
    :param first_lai: the first atom id of the first appearing ligand in the topology file
    """
    final_dict = complexes_dict.copy()

    for ligand_key in ligand_key_list:
        ligand_atom_group = u.select_atoms(ligand_key)
        current_resid = ligand_atom_group[0].resid
        ligand_name = ligand_atom_group[0].resname
        current_complex = 1

        for atom_info in ligand_atom_group:
            current_index = int(atom_info.index)

            if atom_info.resid > current_resid:
                current_complex += 1
                current_resid = atom_info.resid

            if "ligands" not in final_dict[f"{current_complex}"]:
                final_dict[f"{current_complex}"]["ligands"] = {}
                final_dict[f"{current_complex}"]["ligands"][ligand_name] = [current_index]

            elif ligand_name not in final_dict[f"{current_complex}"]["ligands"]:
                final_dict[f"{current_complex}"]["ligands"][ligand_name] = [current_index]

            else:
                final_dict[f"{current_complex}"]["ligands"][ligand_name].append(current_index)

        current_complex = 1
        
    return final_dict

def compute_RMSD_fit_protein(u, full_dict, complex_id, ligand_names):
    """
    Computes the RMSD of a inputed complex with fitting on the complex protein_chain, from a universe and a dictionary of all atoms in the system
    
    :param u: Universe
    :param full_dict: disctionary containing atom ids, separated into complexes. Each complex contains \
    the specific monomer and ligands atoms.
    :param complex_n: the identifier of a complex
    """

    complex_atom_group, complex_keys_list = get_complex_data(u, full_dict, complex_id)

    r = rms.RMSD(complex_atom_group, complex_atom_group, select = "name CA", groupselections = complex_keys_list)
    r.run()

    run_df = pd.DataFrame(r.results.rmsd, columns = ["Frame", "Time(ns)", "Backbone", "Alpha Carbons"] +  ligand_names)

    run_df.plot(x = "Frame", y = ["Alpha Carbons"] + ligand_names); plt.legend(); plt.show()

    return run_df


def compute_RMSD_fit_ligand(u, full_dict, complex_id, ligand_names ):
    """
    Computes the RMSD of a inputed complex with fitting on the ligand, from a universe and a dictionary of all atoms in the system
    
    :param u: Universe
    :param full_dict: disctionary containing atom ids, separated into complexes. Each complex contains \
    the specific monomer and ligands atoms.
    :param complex_n: the identifier of a complex
    """

    complex_atom_group, complex_keys_list = get_complex_data(u, full_dict, complex_id)

    r = rms.RMSD(complex_atom_group, complex_atom_group, select = "name CA", groupselections = complex_keys_list)
    r.run()

    run_df = pd.DataFrame(r.results.rmsd, columns = ["Frame", "Time(ns)", "Backbone", "Alpha Carbons"] +  ligand_names)

    run_df.plot(x = "Frame", y = ["Alpha Carbons"] + ligand_names); plt.legend(); plt.show()

    return run_df

def ids_to_key(list):
    out = ""
    for item in list:
        out += f"{item} "
    return f"index {out.strip()}"

def complex_to_ids(full_dict, complex_id):
    """
    Docstring for complex_to_keys
    
    :param full_dict: Description
    :param complex_id: Description
    """
    ids_list= []
    chosen_complex = full_dict[f"{complex_id}"]
    for k in chosen_complex.keys():
        if k == "protein_chain":
            ids_list.append(chosen_complex[k])
        else:
            for kk in chosen_complex[k].keys():
                ids_list.append(chosen_complex[k][kk])
    return ids_list

def get_complex_data(u, full_dict, complex_id):
    """
    From a universe, a dictionary containing all atoms in the system and a specific complex id, \
    outputs the atomids of the complex, the atom_group and the separate keys for all the groups in \
    the complex.
    
    :param full_dict: Description
    :param complex_id: Description
    """
    complex_groups_atomids = complex_to_ids(full_dict, complex_id)

    all_complex_atomids = []
    for group_list in complex_groups_atomids:
        for atom in group_list:
            all_complex_atomids.append(atom)
    complex_atom_group = u.atoms[all_complex_atomids]

    complex_keys_list = []
    for group_atomids in complex_groups_atomids:
        group_keys = ids_to_key(group_atomids)
        complex_keys_list.append(group_keys)

    return complex_atom_group, complex_keys_list



def main():

    #create universe
    U = mda.Universe(sys.argv[1],sys.argv[2])
       
    #select the mobile domain of the protein specifically
    protein_full = U.select_atoms(PROTEIN_KEY)
    #analyse the radius of gyration
    protein_dict = parse_protein(protein_full)
    #select the ligands specifically
    full_dict = add_ligands_to_complexes(U, protein_dict, LIGANDS_KEYS)


    complex_id = "1"
    #compute the RMSD for each complex
    ligand_names = list(full_dict[complex_id]["ligands"].keys())
    
    chain_A_df = compute_RMSD_fit_protein(U, full_dict, complex_id, ligand_names)


if __name__ == "__main__":
    main()
